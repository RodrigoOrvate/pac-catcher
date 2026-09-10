"""
audita_transientes.py (23/08/2026)
==========================================
Auditoria da hipótese de ARTEFATO POR TRANSIENTES broadband nos vencedores
de PAC theta-gamma.

Hipótese em teste (levantada na revisão do vencedor 2, chan22 @ 003 20-30s):
"espigões" de ruído mecânico/EMG passam pelos dois filtros passa-banda e
geram pico simultâneo de fase (theta) e amplitude (gamma), inflando o MI.
A nula de deslocamento circular NÃO protege contra isso -- o deslocamento
quebra justamente a coincidência dos transientes, então transientes
coincidentes produzem z alto. Esta auditoria discrimina com:

  1. CONTEXTO de amplitude: o RMS/pico da janela é excepcional para o
     canal (vs. janelas de 10s do próprio arquivo) e para a sessão
     (rank vs. os 32 canais na mesma janela)? Kurtose do bruto e da banda
     gamma (transiente -> kurtose alta).
  2. DESPIKE: interpolar linearmente as amostras fora de k*sigma_robusta
     (mediana/MAD da própria janela, guarda de +-50 ms para o ringing dos
     filtros) e recomputar o z do par de pico com a MESMA nula (mesmos
     deslocamentos, rng semeado). Se o z colapsa, o acoplamento morava nos
     transientes; se sobrevive, nao era dirigido por espigoes.
  3. SUB-JANELAS de 2s: acoplamento genuíno aparece na maioria das sub-
     janelas; dirigido por 1-2 episódios concentra todo o z nelas.
  4. CONTROLES para validar o proprio teste:
     - chan20 @ 001 85-95 (rearing) e chan10 @ 001 30-40 (grooming):
       controles "bons" (validados em comportamento e respiracao);
     - chan32 @ 002 95-105 (imovel): controle ARTEFATO conhecido
       (janela rejeitada #3 -- transientes ritmados, canal spiky de borda).
     O veredito sobre chan22 so vale se o teste separar esses dois grupos.
  5. VIZINHOS do caso disputado (chan24/26/28 na mesma janela): fonte
     cortical local tem gradiente; cabo batendo/EMG difuso nao respeita
     gradiente entre canais vizinhos.

Uso (casos da sessão 08/07 como default — nada precisa ser passado):

    python audita_transientes.py
    python audita_transientes.py --pasta_ns2 "../Basal antes da infusao" --saida_dir auditoria

Nova sessão, por CLI (casos específicos NUNCA entram no código). Formato de
--casos: "rotulo,arquivo,canal,ini,fim,fp,fa" separados por ';'. Ex.:

    python audita_transientes.py \
        --pasta_ns2 "<SESSAO>/<BASAL>" --saida_dir "<SESSAO>/RESULTADOS/auditoria" \
        --casos "ev02_disputado,20240711-121046-001.ns2,chan26,65,75,5,70;ev15_limpo,20240711-121046-003.ns2,chan24,275,285,6,40" \
        --vizinhos ""   # pula o bloco de vizinhos (ele é específico de cada sessão)

Saídas: <saida_dir>/auditoria_transientes.csv (z por condição),
        <saida_dir>/contexto_amplitude.csv (estatísticas de contexto)
        e um PNG por caso (bruto com transientes marcados, gamma+envelope,
        histograma de fase, z por condição/sub-janela).

Requer: neo, numpy, scipy, pandas, matplotlib (mesmos do pipeline).
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
import pandas as pd
import scipy.signal as signal
import scipy.stats as stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ns2_utils import le_ns2, fatia_janela
from comodulogram import _mi_de_bin_idx
from pac_core.filtering import filtra_sinal, aplica_notch

# Mesmos parâmetros do comodulogram.py / triagem_pac.py
N_BINS = 18
N_SURR = 200
NOTCH_HZ = 60.0

# Casos: o disputado, dois controles bons e um controle de artefato conhecido.
# fp/fa = par de pico (Hz) de cada um, conforme robustez_parametros.csv / CLAUDE.md.
CASOS = [
    dict(rotulo="chan22@20-30_disputado", arquivo="20240708-123605-003.ns2",
         canal="chan22", ini=20.0, fim=30.0, fp=5, fa=35),
    dict(rotulo="chan20@85-95_controle_bom", arquivo="20240708-123605-001.ns2",
         canal="chan20", ini=85.0, fim=95.0, fp=7, fa=75),
    dict(rotulo="chan10@30-40_controle_bom", arquivo="20240708-123605-001.ns2",
         canal="chan10", ini=30.0, fim=40.0, fp=8, fa=60),
    dict(rotulo="chan32@95-105_controle_artefato", arquivo="20240708-123605-002.ns2",
         canal="chan32", ini=95.0, fim=105.0, fp=8, fa=40),
]

# Vizinhos imediatos do caso disputado (mesma janela, mesmo par de pico)
VIZINHOS = [("chan24", 6.8), ("chan26", 5.8), ("chan28", 5.9)]  # (canal, z publicado)


# ==========================================
# NÚCLEO (mesma matemática do comodulogram.py)
# ==========================================

def mi_z_par(lfp, fs, f_fase, f_amp, n_bins=N_BINS, n_surr=N_SURR, rng=None):
    """
    MI z-scoredo na célula (f_fase, f_amp) -- exatamente a mesma matemática
    de calcula_comodulograma_z(), mas só para UM par de frequências.
    rng semeado -> os MESMOS deslocamentos de surrogate em todas as
    condições (baseline vs despike), então a diferença de z é só do sinal.
    """
    if rng is None:
        rng = np.random.default_rng(42)
    bins = np.linspace(-np.pi, np.pi, n_bins + 1)

    lfp_fase = filtra_sinal(lfp, f_fase - 1.0, f_fase + 1.0, fs)
    fase = np.angle(signal.hilbert(lfp_fase))
    bin_idx = np.clip(np.digitize(fase, bins) - 1, 0, n_bins - 1)

    lfp_amp = filtra_sinal(lfp, f_amp - 5.0, f_amp + 5.0, fs)
    env = np.abs(signal.hilbert(lfp_amp))

    mi_obs = _mi_de_bin_idx(bin_idx, env, n_bins)

    n = lfp.size
    shift_min = int(1.0 * fs)  # >= 1s de deslocamento, igual à nula do triagem
    if n <= 2 * shift_min:
        shift_min = max(1, n // 10)
    deslocamentos = rng.integers(shift_min, n - shift_min, size=n_surr)
    mi_surr = np.array([_mi_de_bin_idx(bin_idx, np.roll(env, d), n_bins)
                        for d in deslocamentos])

    dp = mi_surr.std()
    z = float((mi_obs - mi_surr.mean()) / dp) if dp > 0 else 0.0
    return z, float(mi_obs), env, bin_idx


def remove_transientes(lfp, fs, k=6.0, guarda_ms=50.0):
    """
    Marca amostras fora de med +- k*sigma_robusta (1.4826*MAD da PRÓPRIA
    janela), dilata a máscara +-guarda_ms (ringing dos filtros) e substitui
    por interpolação linear a partir das amostras limpas vizinhas.
    Retorna (sinal_limpo, mascara_ruins, limiar, frac_removida).
    """
    med = np.median(lfp)
    mad = np.median(np.abs(lfp - med)) * 1.4826
    limiar = k * mad
    ruins = np.abs(lfp - med) > limiar

    raio = int(guarda_ms / 1000.0 * fs)
    if raio > 0 and ruins.any():
        kernel = np.ones(2 * raio + 1, dtype=int)
        ruins = np.convolve(ruins.astype(int), kernel, mode="same") > 0

    y = lfp.copy()
    if ruins.any():
        idx_ok = np.flatnonzero(~ruins)
        if len(idx_ok) > 1:
            y[ruins] = np.interp(np.flatnonzero(ruins), idx_ok, lfp[idx_ok])
    frac = float(ruins.mean())
    return y, ruins, float(limiar), frac


# ==========================================
# CONTEXTO DE AMPLITUDE
# ==========================================

def contexto_amplitude(lfp_bruto, coluna_canal_inteira, fs, dados_janela_todos,
                       idx_canal, fa):
    """
    Estatísticas que dizem se a janela é excepcional:
    - vs o próprio canal: percentil do RMS da janela entre as janelas de 10s
      não-sobrepostas do arquivo inteiro; nº de janelas do arquivo com pico
      >= ao pico desta janela;
    - vs a sessão: rank do RMS entre os 32 canais na MESMA janela;
    - forma: kurtose do bruto (excesso) e kurtose do sinal filtrado em
      fa+-5 Hz (transiente -> kurtose alta nas duas).
    """
    rms_janela = float(np.sqrt(np.mean(lfp_bruto ** 2)))
    pico_janela = float(np.max(np.abs(lfp_bruto)))

    n_win = int(len(coluna_canal_inteira) // (10 * fs))
    rms_arquivo = np.array([
        np.sqrt(np.mean((coluna_canal_inteira[i * int(10 * fs):(i + 1) * int(10 * fs)]) ** 2))
        for i in range(n_win)
    ])
    pico_arquivo = np.array([
        np.max(np.abs(coluna_canal_inteira[i * int(10 * fs):(i + 1) * int(10 * fs)]))
        for i in range(n_win)
    ])
    pct_rms = float(stats.percentileofscore(rms_arquivo, rms_janela))
    n_janelas_pico_maior = int(np.sum(pico_arquivo >= pico_janela))

    rms_todos = np.sqrt(np.mean(dados_janela_todos.astype(float) ** 2, axis=0))
    ordem = np.argsort(rms_todos)[::-1]
    rank = int(np.flatnonzero(ordem == idx_canal)[0]) + 1

    lfp_n = aplica_notch(lfp_bruto, fs, linha_hz=NOTCH_HZ)
    kurt_bruto = float(stats.kurtosis(lfp_bruto))
    gama = filtra_sinal(lfp_n, fa - 5.0, fa + 5.0, fs)
    kurt_gama = float(stats.kurtosis(gama))

    return {
        "rms_janela": rms_janela,
        "pico_abs_janela": pico_janela,
        "min_janela": float(np.min(lfp_bruto)),
        "max_janela": float(np.max(lfp_bruto)),
        "pct_rms_vs_arquivo": pct_rms,
        "n_janelas_10s_com_pico_maior": n_janelas_pico_maior,
        "rank_rms_entre_32_canais": rank,
        "kurtose_bruto": kurt_bruto,
        "kurtose_gama": kurt_gama,
    }


# ==========================================
# FIGURA POR CASO
# ==========================================

def figura_caso(lfp_n, fs, ruins, env, hist_fase, rotulo_z, titulo, caminho_png):
    t = np.arange(lfp_n.size) / fs
    gama = filtra_sinal(lfp_n, 30.0, 150.0, fs)  # só para desenho

    fig, axs = plt.subplots(2, 2, figsize=(13, 8))

    ax = axs[0, 0]
    ax.plot(t, lfp_n, linewidth=0.4, color="tab:blue")
    if ruins.any():
        ax.fill_between(t, lfp_n.min(), lfp_n.max(), where=ruins,
                        color="red", alpha=0.25, linewidth=0)
    ax.set_title("LFP (notch 60) — vermelho = amostras removidas no despike")
    ax.set_xlabel("Tempo (s)")

    ax = axs[0, 1]
    ax.plot(t, gama, linewidth=0.4, color="tab:orange", alpha=0.6)
    ax.plot(t, env, linewidth=1.2, color="darkred")
    if ruins.any():
        ax.fill_between(t, gama.min(), gama.max(), where=ruins,
                        color="red", alpha=0.25, linewidth=0)
    ax.set_title("Gamma 30–150 Hz + envelope")
    ax.set_xlabel("Tempo (s)")

    ax = axs[1, 0]
    ax.bar(np.arange(len(hist_fase)), hist_fase, color="tab:purple")
    ax.axhline(1.0, color="k", linestyle="--", linewidth=0.8)
    ax.set_title("Amplitude de γ normalizada por bin de fase θ (média = 1)")
    ax.set_xlabel("Bin de fase (0–17, −π → π)")
    ax.set_ylabel("Amplitude relativa")

    ax = axs[1, 1]
    nomes = list(rotulo_z.keys())
    vals = [rotulo_z[n] for n in nomes]
    cores = ["tab:green" if "baseline" in n else
             ("tab:red" if "despike" in n else "tab:gray") for n in nomes]
    ax.bar(range(len(vals)), vals, color=cores)
    ax.axhline(0.0, color="k", linewidth=0.8)
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(nomes, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("z do MI (mesma nula)")
    ax.set_title("z do par de pico por condição")

    fig.suptitle(titulo, fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(caminho_png, dpi=150)
    plt.close(fig)


# ==========================================
# MAIN
# ==========================================

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pasta_ns2", default="../Basal antes da infusao",
                    help="Pasta com os .ns2 originais")
    ap.add_argument("--saida_dir", default="auditoria",
                    help="Diretório de saída (CSVs e PNGs)")
    ap.add_argument("--casos", default=None,
                    help='Casos no formato "rotulo,arquivo,canal,ini,fim,fp,fa" '
                         "separados por ';'. Omitir = casos default da sessão 08/07.")
    ap.add_argument("--vizinhos", default="chan24:6.8,chan26:5.8,chan28:5.9",
                    help='"canal:z_publicado,...". Passar "" para pular o bloco.')
    ap.add_argument("--vizinhos_ini", type=float, default=20.0)
    ap.add_argument("--vizinhos_fim", type=float, default=30.0)
    ap.add_argument("--vizinhos_fp", type=int, default=5)
    ap.add_argument("--vizinhos_fa", type=int, default=35)
    args = ap.parse_args()

    if args.casos:
        casos = []
        for item in args.casos.split(";"):
            rot, arq, can, ini, fim, fp, fa = item.split(",")
            casos.append(dict(rotulo=rot, arquivo=arq, canal=can,
                              ini=float(ini), fim=float(fim), fp=int(fp), fa=int(fa)))
    else:
        casos = CASOS
    vizinhos = []
    if args.vizinhos.strip():
        for item in args.vizinhos.split(","):
            can, zpub = item.split(":")
            vizinhos.append((can.strip(), float(zpub)))

    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

    os.makedirs(args.saida_dir, exist_ok=True)
    cache_arquivos = {}

    def carrega(arquivo):
        if arquivo not in cache_arquivos:
            dados, fs, nomes = le_ns2(os.path.join(args.pasta_ns2, arquivo))
            cache_arquivos[arquivo] = (dados, fs, {str(n): i for i, n in enumerate(nomes)})
        return cache_arquivos[arquivo]

    linhas_z = []
    linhas_ctx = []

    for caso in casos:
        dados, fs, mapa = carrega(caso["arquivo"])
        idx = mapa[caso["canal"]]
        lfp = fatia_janela(dados[:, idx], fs, caso["ini"], caso["fim"]).astype(float)
        lfp_n = aplica_notch(lfp, fs, linha_hz=NOTCH_HZ)
        coluna = dados[:, idx]

        dados_janela = fatia_janela(dados, fs, caso["ini"], caso["fim"])
        ctx = contexto_amplitude(lfp, coluna, fs, dados_janela, idx, caso["fa"])
        ctx.update({"caso": caso["rotulo"], "arquivo": caso["arquivo"],
                    "canal": caso["canal"]})
        linhas_ctx.append(ctx)

        print(f"\n=== {caso['rotulo']} ({caso['arquivo']}, pico "
              f"{caso['fp']}x{caso['fa']} Hz) ===")
        print(f"  bruto: min={ctx['min_janela']:.0f} max={ctx['max_janela']:.0f} "
              f"RMS={ctx['rms_janela']:.0f} | percentil RMS vs arquivo="
              f"{ctx['pct_rms_vs_arquivo']:.0f}% | "
              f"janelas 10s com pico maior: {ctx['n_janelas_10s_com_pico_maior']}/"
              f"{int(len(coluna) // (10 * fs))} | rank RMS/32 canais="
              f"{ctx['rank_rms_entre_32_canais']} | "
              f"kurtose bruto={ctx['kurtose_bruto']:.1f} "
              f"gama={ctx['kurtose_gama']:.1f}")

        # 1) baseline (idêntico ao publicado)
        z0, mi0, env, bin_idx = mi_z_par(lfp_n, fs, caso["fp"], caso["fa"])
        print(f"  baseline:           z={z0:6.2f}  MI={mi0:.4f}")
        linhas_z.append(dict(caso=caso["rotulo"], canal=caso["canal"],
                             condicao="baseline", z=z0, mi=mi0, frac_removida=0.0))

        # 2) despike k=6 e k=5
        rotulo_z = {"baseline": z0}
        for k in (6.0, 5.0):
            limpo, ruins, limiar, frac = remove_transientes(lfp_n, fs, k=k)
            zk, mik, _, _ = mi_z_par(limpo, fs, caso["fp"], caso["fa"])
            nome = f"despike_k{int(k)}"
            print(f"  {nome} (frac={frac*100:5.2f}%): z={zk:6.2f}  MI={mik:.4f}  "
                  f"limiar=+-{limiar:.0f}")
            rotulo_z[nome] = zk
            linhas_z.append(dict(caso=caso["rotulo"], canal=caso["canal"],
                                 condicao=nome, z=zk, mi=mik, frac_removida=frac))

        # 3) sub-janelas de 2s
        n_sub = 5
        dur = lfp_n.size // n_sub
        zs_sub = []
        for s in range(n_sub):
            pedaco = lfp_n[s * dur:(s + 1) * dur]
            zsub, misub, _, _ = mi_z_par(pedaco, fs, caso["fp"], caso["fa"])
            zs_sub.append(zsub)
            nome = f"sub{s + 1}_{caso['ini'] + s * dur / fs:.0f}-{caso['ini'] + (s + 1) * dur / fs:.0f}s"
            rotulo_z[nome] = zsub
            linhas_z.append(dict(caso=caso["rotulo"], canal=caso["canal"],
                                 condicao=nome, z=zsub, mi=misub, frac_removida=np.nan))
        print(f"  sub-janelas 2s:     z=" + ", ".join(f"{z:5.2f}" for z in zs_sub))

        # histograma de fase (desritivo) na condição baseline
        cont = np.bincount(bin_idx, minlength=N_BINS)
        soma = np.bincount(bin_idx, weights=env, minlength=N_BINS)
        media_bins = np.divide(soma, cont, out=np.zeros(N_BINS), where=cont > 0)
        hist_fase = media_bins / media_bins.sum() * N_BINS
        print(f"  excesso máx de fase: {hist_fase.max():.2f}x a média "
              f"(bin {int(np.argmax(hist_fase))})")

        # figura
        limpo6, ruins6, _, _ = remove_transientes(lfp_n, fs, k=6.0)
        figura_caso(
            lfp_n, fs, ruins6, env, hist_fase, rotulo_z,
            f"{caso['rotulo']} — {caso['arquivo']} — pico {caso['fp']}×{caso['fa']} Hz\n"
            f"baseline z={z0:.1f} | despiked k6 z={rotulo_z['despike_k6']:.1f} | "
            f"kurtose bruto={ctx['kurtose_bruto']:.1f}",
            os.path.join(args.saida_dir, f"auditoria_{caso['rotulo']}.png"),
        )

    # 5) vizinhos do caso disputado (mesma janela/par de pico)
    if vizinhos:
        print(f"\n=== vizinhos ({args.vizinhos_ini:g}-{args.vizinhos_fim:g}s, "
              f"pico {args.vizinhos_fp}x{args.vizinhos_fa} Hz) ===")
        dados, fs, mapa = carrega(casos[0]["arquivo"])
        for canal, z_pub in vizinhos:
            idx = mapa[canal]
            lfp = fatia_janela(dados[:, idx], fs,
                               args.vizinhos_ini, args.vizinhos_fim).astype(float)
            lfp_n = aplica_notch(lfp, fs, linha_hz=NOTCH_HZ)
            z0, mi0, _, _ = mi_z_par(lfp_n, fs, args.vizinhos_fp, args.vizinhos_fa)
            limpo, ruins, _, frac = remove_transientes(lfp_n, fs, k=6.0)
            zk, mik, _, _ = mi_z_par(limpo, fs, args.vizinhos_fp, args.vizinhos_fa)
            print(f"  {canal}: baseline z={z0:5.2f} (publicado ~{z_pub}) -> "
                  f"despike k6 z={zk:5.2f} (frac={frac * 100:.2f}%)")
            linhas_z.append(dict(caso="vizinho_" + canal, canal=canal,
                                 condicao="baseline", z=z0, mi=mi0, frac_removida=0.0))
            linhas_z.append(dict(caso="vizinho_" + canal, canal=canal,
                                 condicao="despike_k6", z=zk, mi=mik, frac_removida=frac))

    pd.DataFrame(linhas_z).to_csv(
        os.path.join(args.saida_dir, "auditoria_transientes.csv"), index=False)
    pd.DataFrame(linhas_ctx).to_csv(
        os.path.join(args.saida_dir, "contexto_amplitude.csv"), index=False)
    print(f"\nCSVs salvos em {args.saida_dir}/ "
          f"(auditoria_transientes.csv, contexto_amplitude.csv)")


if __name__ == "__main__":
    main()
