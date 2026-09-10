"""PASSO 0 - exploracao_minuto.py"""
import os, sys, argparse
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
from scipy.signal import welch, hilbert, butter, filtfilt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pac_core.io import carrega_dados, fatia_janela
from pipeline.etapa1_triagem.triagem_pac import mi_com_surrogates  # z-score com 200 surrogates (correto)

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

def calcula_timeline_mi_triplo(pasta, arquivo, canal, t_ini_s, t_fim_s,
                                passo_s=5, n_surr=100):
    """
    Timeline de MI z-score para os TRES pares de banda:
      theta_gamma (4-8 x 30-80 Hz)   — exploração locomotora
      theta_hg    (4-8 x 80-150 Hz)  — estado misto
      theta_hfo   (4-8 x 150-250 Hz) — repouso / SWR-associado

    Retorna: ts, z_tg, z_thg, z_thfo (arrays com mesmo comprimento)

    Usa mi_com_surrogates() com n_surr=100 surrogates (mais rápido que 200
    mas estátisticamente válido para painel visual interativo).
    """
    path = os.path.join(pasta, arquivo)
    dados, fs, canal_ids = carrega_dados(path)
    idx = canal_ids.index(canal)
    rng = np.random.default_rng(42)
    nyq = fs * 0.5
    ts, z_tg, z_thg, z_thfo, ratio_hfo_g = [], [], [], [], []

    for ini in np.arange(t_ini_s, t_fim_s - 10, passo_s):
        fim = ini + 10
        sinal = fatia_janela(dados, fs, ini, fim)[:, idx].astype(float)

        # Fase (theta) — calculada UMA vez, reutilizada para os 3 pares
        s_th = filtra_sinal(sinal, 4, 8, fs, order=3)
        fase = np.angle(hilbert(s_th))

        # Par 1: Theta x Gamma
        env = np.abs(hilbert(filtra_sinal(sinal, 30, 80, fs)))
        _, z1, _, _, _ = mi_com_surrogates(fase, env, fs, n_surr=n_surr, rng=rng)

        # Par 2: Theta x HG
        env = np.abs(hilbert(filtra_sinal(sinal, 80, 150, fs)))
        _, z2, _, _, _ = mi_com_surrogates(fase, env, fs, n_surr=n_surr, rng=rng)

        # Par 3: Theta x HFO (respeita Nyquist)
        hfo_hi = min(250, nyq * 0.95)
        if 150 < hfo_hi:
            env = np.abs(hilbert(filtra_sinal(sinal, 150, hfo_hi, fs)))
            _, z3, _, _, _ = mi_com_surrogates(fase, env, fs, n_surr=n_surr, rng=rng)
        else:
            z3 = float('nan')

        # Ratio HFO/Gamma (proxy de harmônico de Gamma — baixo = HFO genuino)
        p_g   = np.mean(filtra_sinal(sinal, 30, 80, fs) ** 2) + 1e-12
        p_hfo = np.mean(filtra_sinal(sinal, 150, hfo_hi, fs) ** 2) + 1e-12 if 150 < hfo_hi else 0.0
        ratio = p_hfo / p_g

        ts.append(ini + 5)
        z_tg.append(z1)
        z_thg.append(z2)
        z_thfo.append(z3)
        ratio_hfo_g.append(ratio)

    return (np.array(ts), np.array(z_tg), np.array(z_thg),
            np.array(z_thfo), np.array(ratio_hfo_g))


# Mantido para compatibilidade com outros scripts que chamam calcula_timeline_mi
def calcula_timeline_mi(pasta, arquivo, canal, t_ini_s, t_fim_s, passo_s=5,
                         n_surr=100):
    ts, z_tg, _, _, _ = calcula_timeline_mi_triplo(
        pasta, arquivo, canal, t_ini_s, t_fim_s, passo_s, n_surr)
    return ts, z_tg

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

def plota_painel(minuto_idx, pasta, arquivo, canal, t_ini_s, t_fim_s, saida_path,
                  n_surr=100):
    """
    Painel de 4 subplots por minuto:
      1. PSD (4 bandas sombreadas: Theta, Gamma, HG, HFO)
      2. Timeline MI z-score — 3 pares sobrepostos (Theta-Gamma, Theta-HG, Theta-HFO)
         Interpretacao: TG alto = exploracao; THFO alto = repouso/SWR
      3. Timeline Theta Power (log10) — proxy de oscilação ativa
      4. ratio_hfo_gamma — proxy de harmônico Gamma->HFO (alto = suspeito)
    """
    os.makedirs(os.path.dirname(saida_path) or ".", exist_ok=True)
    fig, axes = plt.subplots(4, 1, figsize=(13, 14))
    fig.suptitle(f"Minuto {minuto_idx} — Canal {canal} — {t_ini_s:.0f}-{t_fim_s:.0f}s",
                  fontsize=14, weight="bold", y=0.98)

    # --- Painel 1: PSD ---
    f, p = calcula_psd_minuto(pasta, arquivo, canal, t_ini_s, t_fim_s)
    ax = axes[0]
    ax.semilogy(f, p, color="#1f77b4", lw=1.5, label="PSD")
    ax.axvspan(4,   8,   color="#2ca02c", alpha=0.18, label="Theta 4-8 Hz")
    ax.axvspan(30,  80,  color="#ff7f0e", alpha=0.15, label="Gamma 30-80 Hz")
    ax.axvspan(80,  150, color="#9467bd", alpha=0.13, label="HG 80-150 Hz")
    ax.axvspan(150, 250, color="#e377c2", alpha=0.12, label="HFO 150-250 Hz")
    for fl in [60, 120, 180]:
        ax.axvline(fl, color="#555555", ls="--", lw=0.8, alpha=0.6)
    ax.set_xlim(0, 260)
    ax.set_xlabel("Frequencia (Hz)", fontsize=10)
    ax.set_ylabel("PSD (uV2/Hz)", fontsize=10)
    ax.set_title("PSD media da janela (4 bandas PAC marcadas)", fontsize=11)
    ax.legend(loc="upper right", fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3)

    # --- Painel 2: Timeline MI z-score (3 pares) ---
    ts, z_tg, z_thg, z_thfo, ratio = calcula_timeline_mi_triplo(
        pasta, arquivo, canal, t_ini_s, t_fim_s, passo_s=5, n_surr=n_surr)
    t_rel = ts - t_ini_s
    ax = axes[1]
    ax.plot(t_rel, z_tg,   marker="o", ms=3, color="#ff7f0e", lw=1.8,
            label="z Theta x Gamma (exploracao)")
    ax.plot(t_rel, z_thg,  marker="s", ms=3, color="#9467bd", lw=1.8,
            label="z Theta x HG")
    ax.plot(t_rel, z_thfo, marker="^", ms=3, color="#e377c2", lw=1.8,
            label="z Theta x HFO (repouso/SWR)")
    ax.axhline(3, color="#d62728", ls="--", lw=1, label="z=3 (limiar)")
    ax.axhline(0, color="#aaaaaa", ls="-", lw=0.5)
    ax.set_ylabel("MI z-score", fontsize=10)
    ax.set_title("Timeline MI — 3 pares (janelas 10s / passo 5s)", fontsize=11)
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(True, alpha=0.3)

    # --- Painel 3: Theta Power ---
    ts_th, th = calcula_timeline_theta_power(pasta, arquivo, canal, t_ini_s, t_fim_s)
    ax = axes[2]
    ax.plot(ts_th - t_ini_s, th, marker="s", ms=3, color="#2ca02c", lw=1.5,
            label="Theta 4-8 Hz (log10 uV2)")
    ax.set_ylabel("Theta Power\n(log10 uV2)", fontsize=10)
    ax.set_title("Timeline Theta Power — proxy de estado ativo", fontsize=11)
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(True, alpha=0.3)

    # --- Painel 4: ratio HFO/Gamma ---
    ax = axes[3]
    ax.plot(t_rel, ratio, marker="D", ms=3, color="#8c564b", lw=1.3,
            label="Potencia HFO / Gamma")
    ax.axhline(0.3, color="#d62728", ls="--", lw=1,
               label="0.3 (alerta harmonico Gamma->HFO)")
    ax.set_ylabel("ratio HFO/Gamma", fontsize=10)
    ax.set_xlabel("Tempo dentro da janela (s)", fontsize=10)
    ax.set_title("ratio_hfo_gamma — baixo = HFO genuino; alto = provavel harmonico",
                  fontsize=11)
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(True, alpha=0.3)

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig(saida_path, dpi=110, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {saida_path}")


# Ranking rápido por potência (sem MI completo) — top N canais por minuto
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
        s_hg = filtra_sinal(sinal, 80, 150, fs, order=3)
        s_hfo = filtra_sinal(sinal, 150, 250, fs, order=3)
        f_th, p_th = welch(s_th, fs=fs, nperseg=int(1.2*fs), noverlap=int(1.2*fs)//2, nfft=4000)
        f_g, p_g = welch(s_g, fs=fs, nperseg=int(1.2*fs), noverlap=int(1.2*fs)//2, nfft=4000)
        f_hg, p_hg = welch(s_hg, fs=fs, nperseg=int(1.2*fs), noverlap=int(1.2*fs)//2, nfft=4000)
        f_hfo, p_hfo = welch(s_hfo, fs=fs, nperseg=int(1.2*fs), noverlap=int(1.2*fs)//2, nfft=4000)
        p_th_mean = np.mean(p_th[(f_th>=4)&(f_th<=8)])
        p_g_mean = np.mean(p_g[(f_g>=30)&(f_g<=80)])
        p_hg_mean = np.mean(p_hg[(f_hg>=80)&(f_hg<=150)])
        p_hfo_mean = np.mean(p_hfo[(f_hfo>=150)&(f_hfo<=250)])
        potencias[ch] = p_th_mean + p_g_mean + p_hg_mean + p_hfo_mean
    top = sorted(potencias.items(), key=lambda x: x[1], reverse=True)[:n_top]
    return [ch for ch,_ in top]


# Prote��o de dura��o real � evita overrun quando arquivo � ~298s (n�o 300)
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
    print(f"= PASSO 0 — Timeline por minuto =")
    print(f"Pasta: {args.pasta_ns2}  Arquivos: {arquivos}  Canal fixo: {args.canal}  Minutos: {total_min}  top_n={args.top_n}")
    for minuto in range(1, total_min + 1):
        arquivo_idx = (minuto - 1) // 5
        arquivo = arquivos[min(arquivo_idx, len(arquivos) - 1)]
        t_ini_local = ((minuto - 1) % 5) * 60
        t_fim_local = t_ini_local + 60
        # Limita � dura��o real do arquivo (varia��o ~2 s por .ns2)
        dur_real = duracao_arquivo(args.pasta_ns2, arquivo)
        t_fim_local = min(t_fim_local, max(dur_real - 0.5, 0))
        if args.top_n > 0:
            base = os.path.join(args.saida_dir, f"minuto_{minuto:03d}")
            os.makedirs(base, exist_ok=True)
            top = ranking_canais_potencia(args.pasta_ns2, arquivo, t_ini_local, t_fim_local, n_top=args.top_n)
            print(f"\n-> Minuto {minuto:02d} (arquivo {arquivo}, local {t_ini_local}-{t_fim_local}s) — top: {top}")
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
