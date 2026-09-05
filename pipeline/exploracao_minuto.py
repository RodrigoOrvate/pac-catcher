"""PASSO 0 - exploracao_minuto.py"""
import os, sys, argparse
sys.path.insert(0, r"C:\acoplamento_theta-gamma\SCRIPT\pipeline")
import numpy as np
from scipy.signal import welch, hilbert, butter, filtfilt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ns2_utils import carrega_dados, fatia_janela

def filtra_sinal(sinal, lowcut, highcut, fs, order=3):
    nyq = 0.5 * fs
    b, a = butter(order, [lowcut/nyq, highcut/nyq], btype="band")
    return filtfilt(b, a, sinal)

def detecta_arquivos(pasta):
    return sorted([f for f in os.listdir(pasta) if f.endswith(".ns2")])

def calcula_psd_minuto(pasta, arquivo, canal, t_ini_s, t_fim_s):
    path = os.path.join(pasta, arquivo)
    dados, fs, canal_ids = carrega_dados(path)
    idx = canal_ids.index(canal)
    sinal = fatia_janela(dados, fs, t_ini_s, t_fim_s)[:, idx]
    sinal = filtra_sinal(sinal, 4, 100, fs, order=3)
    f, p = welch(sinal, fs=fs, nperseg=int(1.2*fs), noverlap=int(1.2*fs)//2, nfft=4000)
    return f, p

def calcula_timeline_mi(pasta, arquivo, canal, t_ini_s, t_fim_s, passo_s=5):
    path = os.path.join(pasta, arquivo)
    dados, fs, canal_ids = carrega_dados(path)
    idx = canal_ids.index(canal)
    ts, zs = [], []
    for ini in np.arange(t_ini_s, t_fim_s - 10, passo_s):
        fim = ini + 10
        sinal = fatia_janela(dados, fs, ini, fim)[:, idx]
        s_th = filtra_sinal(sinal, 4, 8, fs, order=3)
        s_g = filtra_sinal(sinal, 30, 80, fs, order=3)
        fase = np.angle(hilbert(s_th))
        env = np.abs(hilbert(s_g))
        # MI simplificado (amostra 1/10 para velocidade)
        fs_down = max(1, fs // 10)
        fase_s = fase[::10]
        env_s = env[::10]
        bins = np.linspace(-np.pi, np.pi, 19)
        bin_idx = np.clip(np.digitize(fase_s, bins) - 1, 0, 17)
        # MI observado
        p_obs = np.bincount(bin_idx, weights=env_s, minlength=18) / np.bincount(bin_idx, minlength=18).clip(1, None)
        p_obs = np.clip(p_obs, 1e-12, 1)
        mi_obs = np.sum(p_obs * np.log(p_obs * 18))
        # Surrogate simples
        shift = max(10, len(env_s)//10)
        env_shift = np.roll(env_s, shift)
        p_s = np.bincount(bin_idx, weights=env_shift, minlength=18) / np.bincount(bin_idx, minlength=18).clip(1, None)
        p_s = np.clip(p_s, 1e-12, 1)
        mi_s = np.sum(p_s * np.log(p_s * 18))
        z = (mi_obs - mi_s) / max(mi_s * 0.1, 0.01) if mi_s > 0 else 0.0
        ts.append(ini + 5)
        zs.append(z)
    return np.array(ts), np.array(zs)

def calcula_timeline_theta_power(pasta, arquivo, canal, t_ini_s, t_fim_s, passo_s=5):
    path = os.path.join(pasta, arquivo)
    dados, fs, canal_ids = carrega_dados(path)
    idx = canal_ids.index(canal)
    ts, powers = [], []
    for ini in np.arange(t_ini_s, t_fim_s - 10, passo_s):
        fim = ini + 10
        sinal = fatia_janela(dados, fs, ini, fim)[:, idx]
        s_th = filtra_sinal(sinal, 4, 8, fs, order=3)
        f, p = welch(s_th, fs=fs, nperseg=int(1.2*fs), noverlap=int(1.2*fs)//2, nfft=4000)
        mask = (f >= 4) & (f <= 8)
        powers.append(np.log10(np.mean(p[mask])) if mask.sum()>0 else -3)
        ts.append(ini + 5)
    return np.array(ts), np.array(powers)

def plota_painel(minuto_idx, pasta, arquivo, canal, t_ini_s, t_fim_s, saida_path):
    os.makedirs(os.path.dirname(saida_path) or ".", exist_ok=True)
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))
    f, p = calcula_psd_minuto(pasta, arquivo, canal, t_ini_s, t_fim_s)
    ax = axes[0]
    ax.semilogy(f, p, color="#1f77b4", lw=1.5, label="PSD")
    ax.axvspan(4, 8, color="#2ca02c", alpha=0.15, label="Theta 4-8 Hz")
    ax.axvspan(30, 80, color="#ff7f0e", alpha=0.15, label="Gamma 30-80 Hz")
    for fl in [60, 120, 180]:
        ax.axvline(fl, color="#555555", ls="--", lw=0.8, alpha=0.7)
    ax.set_xlim(0, 200)
    ax.set_xlabel("Frequencia (Hz)", fontsize=11)
    ax.set_ylabel("PSD (uV2/Hz)", fontsize=11)
    ax.set_title(f"Minuto {minuto_idx} â€” Canal {canal} â€” PSD media (60s)", fontsize=13, weight="bold")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.3)
    ts, zs = calcula_timeline_mi(pasta, arquivo, canal, t_ini_s, t_fim_s, passo_s=5)
    ax = axes[1]
    ax.plot(ts - t_ini_s, zs, marker="o", markersize=3, color="#9467bd", lw=1.5, label="MI z-score")
    ax.axhline(3, color="#d62728", ls="--", lw=1, label="z=3 (limiar)")
    ax.fill_between(ts - t_ini_s, 0, zs, alpha=0.15, color="#9467bd")
    ax.set_ylabel("MI z-score", fontsize=11)
    ax.set_title("Timeline MI (janelas 10s / passo 5s)", fontsize=11)
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(True, alpha=0.3)
    ts, th = calcula_timeline_theta_power(pasta, arquivo, canal, t_ini_s, t_fim_s, passo_s=5)
    ax = axes[2]
    ax.plot(ts - t_ini_s, th, marker="s", markersize=3, color="#2ca02c", lw=1.5, label="Theta 4-8 Hz (log10)")
    ax.set_ylabel("Theta Power (log10 uV2)", fontsize=11)
    ax.set_xlabel("Tempo dentro do minuto (s)", fontsize=11)
    ax.set_title("Timeline Theta Power", fontsize=11)
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(saida_path, dpi=110, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {saida_path}")


# Ranking rÃ¡pido por potÃªncia (sem MI completo) â€” top N canais por minuto
# Cada minuto pode ter canais diferentes; evita varrer 32 canais com MI completo.
def ranking_canais_potencia(pasta, arquivo, t_ini_s, t_fim_s, n_top=5):
    path = __import__("os").path.join(pasta, arquivo)
    dados, fs, canal_ids = carrega_dados(path)
    potencias = {}
    for ch in canal_ids:
        idx = canal_ids.index(ch)
        sinal = fatia_janela(dados, fs, t_ini_s, t_fim_s)[:, idx]
        s_th = filtra_sinal(sinal, 4, 8, fs, order=3)
        s_g = filtra_sinal(sinal, 30, 80, fs, order=3)
        f_th, p_th = welch(s_th, fs=fs, nperseg=int(1.2*fs), noverlap=int(1.2*fs)//2, nfft=4000)
        f_g, p_g = welch(s_g, fs=fs, nperseg=int(1.2*fs), noverlap=int(1.2*fs)//2, nfft=4000)
        p_th_mean = np.mean(p_th[(f_th>=4)&(f_th<=8)])
        p_g_mean = np.mean(p_g[(f_g>=30)&(f_g<=80)])
        potencias[ch] = p_th_mean + p_g_mean
    top = sorted(potencias.items(), key=lambda x: x[1], reverse=True)[:n_top]
    return [ch for ch,_ in top]


# Proteção de duração real — evita overrun quando arquivo é ~298s (não 300)
def duracao_arquivo(pasta, arquivo):
    path = os.path.join(pasta, arquivo)
    dados, fs, _ = carrega_dados(path)
    return dados.shape[0] / fs

def main():
    ap = argparse.ArgumentParser(description="PASSO 0: timeline por minuto")
    ap.add_argument("--pasta_ns2", required=True)
    ap.add_argument("--canal", default="chan20")
    ap.add_argument("--saida_dir", default="test_minutos")
    ap.add_argument("--top_n", type=int, default=5, help="Top N canais/minuto (0=fixed --canal)")
    args = ap.parse_args()
    os.makedirs(args.saida_dir, exist_ok=True)
    arquivos = detecta_arquivos(args.pasta_ns2)
    if not arquivos:
        print("[ERRO] Nenhum .ns2"); return
    total_min = len(arquivos) * 5
    print(f"= PASSO 0 â€” Timeline por minuto =")
    print(f"Pasta: {args.pasta_ns2}  Arquivos: {arquivos}  Canal fixo: {args.canal}  Minutos: {total_min}  top_n={args.top_n}")
    for minuto in range(1, total_min + 1):
        arquivo_idx = (minuto - 1) // 5
        arquivo = arquivos[min(arquivo_idx, len(arquivos) - 1)]
        t_ini_local = ((minuto - 1) % 5) * 60
        t_fim_local = t_ini_local + 60
        # Limita à duração real do arquivo (variação ~2 s por .ns2)
        dur_real = duracao_arquivo(args.pasta_ns2, arquivo)
        t_fim_local = min(t_fim_local, max(dur_real - 0.5, 0))
        if args.top_n > 0:
            base = os.path.join(args.saida_dir, f"minuto_{minuto:03d}")
            os.makedirs(base, exist_ok=True)
            top = ranking_canais_potencia(args.pasta_ns2, arquivo, t_ini_local, t_fim_local, n_top=args.top_n)
            print(f"\n-> Minuto {minuto:02d} (arquivo {arquivo}, local {t_ini_local}-{t_fim_local}s) â€” top: {top}")
            for ch in top:
                nome_png = f"{ch}.png"
                saida_path = os.path.join(base, nome_png)
                plota_painel(minuto, args.pasta_ns2, arquivo, ch, t_ini_local, t_fim_local, saida_path)
        else:
            nome_png = f"minuto_{minuto:03d}.png"
            saida_path = os.path.join(args.saida_dir, nome_png)
            print(f"\n-> Minuto {minuto:02d} (arquivo {arquivo}, local {t_ini_local}-{t_fim_local}s)")
            plota_painel(minuto, args.pasta_ns2, arquivo, args.canal, t_ini_local, t_fim_local, saida_path)
    print(f"\nOK: {total_min} minutos processados em {args.saida_dir}/")
if __name__ == "__main__":
    main()
