"""
comodulogram.py
==========================================
Comodulograma fase-amplitude -- suporta os três pares:
  theta_gamma  : Theta 4-8 Hz x Gamma  30-80  Hz
  theta_hg     : Theta 4-8 Hz x HG     80-150 Hz
  theta_hfo    : Theta 4-8 Hz x HFO   150-250 Hz

DOIS MODOS:

  1. MODO LOTE (--csv): lê o resultados_refinados.csv produzido pelo
     refina_candidatos.py e gera um PNG por candidato, com o MI
     Z-SCOREADO contra surrogates de deslocamento circular.
     Usa --par para decidir qual quadrante destacar no mapa.

  2. JANELA ÚNICA: exploração rápida de um trecho específico.

EM AMBOS OS MODOS, --fdr_q 0.05 aplica Benjamini-Hochberg sobre o mapa
e contorna em preto as células significantes.

Frequências: cobre 4-14 Hz (fase) x 30-260 Hz (amplitude), abrangendo
os três pares de uma vez. O argumento --par apenas decide qual quadrante
é usado para extrair o z do pico e classificar o cluster FDR.

Uso:
    # Modo lote para theta_hg:
    python comodulogram.py --csv resultados_refinados.csv \
        --pasta_ns2 "../Basal" --par theta_hg --saida_dir comodulogramas_hg

    # Janela única, theta_hfo, notch 60 Hz:
    python comodulogram.py --arquivo "../Basal/sessao.ns2" \
        --canal 17 --inicio 205 --fim 215 --par theta_hfo --notch 60

Requer: neo, numpy, scipy, pandas, matplotlib
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import scipy.signal as signal
import scipy.stats as stats
from scipy import ndimage
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ns2_utils import le_ns2, carrega_dados, fatia_janela
from triagem_pac import BAND_PAIRS

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

# Grelha de frequencias: cobre os 3 pares em uma unica passagem
# fase: 4-14 Hz (theta e suas bordas)
# amp:  30-260 Hz (gamma + HG + HFO)
FASES_DEFAULT = np.arange(4, 15, 1)            # 11 pontos
AMPS_DEFAULT  = np.concatenate([               # grade mais densa onde importa
    np.arange(30,  80,  5),                    # Gamma: 10 pts
    np.arange(80,  155, 5),                    # HG:    15 pts
    np.arange(155, 265, 10),                   # HFO:   11 pts
])                                             # total 36 pts


# ==========================================
# NÚCLEO DE CÁLCULO (mesmo MI dos outros scripts)
# ==========================================

def filtra_sinal(sinal_in, lowcut, highcut, fs, order=3):
    nyq = 0.5 * fs
    low = max(lowcut / nyq, 1e-6)
    high = min(highcut / nyq, 0.999)
    b, a = signal.butter(order, [low, high], btype="bandpass")
    return signal.filtfilt(b, a, sinal_in)


def _mi_de_bin_idx(bin_idx, envelope, n_bins):
    soma_bins = np.bincount(bin_idx, weights=envelope, minlength=n_bins)
    cont_bins = np.bincount(bin_idx, minlength=n_bins)
    media_bins = np.divide(soma_bins, cont_bins, out=np.zeros(n_bins), where=cont_bins > 0)
    soma = np.sum(media_bins)
    if soma <= 0:
        return 0.0
    P = media_bins / soma
    H = -np.sum(P * np.log(P + 1e-10))
    return (np.log(n_bins) - H) / np.log(n_bins)


def aplica_notch(sinal_in, fs, freqs_notch, q_factor=30.0):
    """
    Rejeita a frequência da rede elétrica e seus harmônicos (iirnotch).
    Aceita uma lista de frequências (ex: [60, 120, 180, 240]).
    Necessário porque picos de linha dentro da banda de amplitude inflam
    o MI simulando acoplamento onde não há.
    """
    if not freqs_notch:
        return sinal_in
        
    nyq = 0.5 * fs
    out = sinal_in
    
    # Suporta tanto um número (comportamento antigo, mas aqui usamos lista) quanto lista
    if not isinstance(freqs_notch, (list, tuple, np.ndarray)):
        freqs_notch = [freqs_notch]
        
    for f in freqs_notch:
        if f >= nyq * 0.98:
            continue
        b, a = signal.iirnotch(f / nyq, q_factor)
        out = signal.filtfilt(b, a, out)
    return out


def calcula_comodulograma_z(lfp_ativo, fs, fases_freq, amps_freq,
                            n_surr=200, n_bins=18, rng=None, notch_hz=None,
                            retorna_mi=False):
    """
    Matriz de z-scores (n_amps x n_fases): para cada par (freq de fase,
    freq de amplitude), compara o KL-MI observado contra n_surr surrogates
    gerados por deslocamento circular do envelope de amplitude.

    notch_hz: se informado (ex.: 60), aplica rejeição na rede elétrica e
    harmônicos ANTES de qualquer filtragem -- rode sempre uma vez com e
    uma vez sem para saber quanto do acoplamento depende do ruído de linha.

    retorna_mi: se True, devolve também o mapa de MI observado e a matriz
    de MIs dos surrogates (n_amps x n_fases x n_surr) -- necessários para o
    p-valor por célula do FDR.

    Eficiência: cada frequência é filtrada UMA vez (fase e amplitude);
    os surrogates só rolam o envelope já pronto, que é barato -- mesmo
    padrão vetorizado do mi_com_surrogates() do triagem_pac.py.
    """
    if rng is None:
        rng = np.random.default_rng()

    if notch_hz:
        lfp_ativo = aplica_notch(lfp_ativo, fs, freqs_notch=notch_hz)

    bins = np.linspace(-np.pi, np.pi, n_bins + 1)

    fases_por_freq = []
    for f_fase in fases_freq:
        lfp_fase = filtra_sinal(lfp_ativo, f_fase - 1.0, f_fase + 1.0, fs)
        fase = np.angle(signal.hilbert(lfp_fase))
        fases_por_freq.append(np.clip(np.digitize(fase, bins) - 1, 0, n_bins - 1))

    envelopes_por_freq = []
    for f_amp in amps_freq:
        lfp_amp = filtra_sinal(lfp_ativo, f_amp - 5.0, f_amp + 5.0, fs)
        envelopes_por_freq.append(np.abs(signal.hilbert(lfp_amp)))

    n = lfp_ativo.size
    shift_min = int(1.0 * fs)  # >= 1s de deslocamento, igual à nula do triagem
    if n <= 2 * shift_min:
        shift_min = max(1, n // 10)
    deslocamentos = rng.integers(shift_min, n - shift_min, size=n_surr)

    mi_obs_mapa = np.zeros((len(amps_freq), len(fases_freq)))
    mi_surr_mapa = np.zeros((len(amps_freq), len(fases_freq), n_surr))
    for i, bin_idx in enumerate(fases_por_freq):
        for j, env in enumerate(envelopes_por_freq):
            mi_obs_mapa[j, i] = _mi_de_bin_idx(bin_idx, env, n_bins)
            for k, desloc in enumerate(deslocamentos):
                mi_surr_mapa[j, i, k] = _mi_de_bin_idx(bin_idx, np.roll(env, desloc), n_bins)

    media = mi_surr_mapa.mean(axis=2)
    dp = mi_surr_mapa.std(axis=2)
    z_mapa = np.divide(mi_obs_mapa - media, dp,
                       out=np.zeros_like(mi_obs_mapa), where=dp > 0)

    if retorna_mi:
        return z_mapa, mi_obs_mapa, mi_surr_mapa
    return z_mapa


def p_valores_por_celula(mi_obs_mapa, mi_surr_mapa):
    """
    p-valor por célula do mapa: cauda da Gama ajustada aos surrogates da
    PRÓPRIA célula -- mesma escolha do refina_candidatos.py (KL-MI é
    não-negativo e assimétrico; assumir normalidade subestimaria a cauda).
    Cai para o p empírico se o ajuste não for possível.
    """
    p = np.empty(mi_obs_mapa.shape)
    for j in range(mi_obs_mapa.shape[0]):
        for i in range(mi_obs_mapa.shape[1]):
            mi_obs = mi_obs_mapa[j, i]
            mi_surr = mi_surr_mapa[j, i]
            mi_surr_pos = mi_surr[mi_surr > 0]
            if len(mi_surr_pos) < 20 or mi_obs <= 0:
                p[j, i] = np.mean(mi_surr >= mi_obs)
                continue
            try:
                forma, _, escala = stats.gamma.fit(mi_surr_pos, floc=0)
                p[j, i] = float(np.clip(
                    stats.gamma.sf(mi_obs, forma, loc=0, scale=escala), 0.0, 1.0))
            except Exception:
                p[j, i] = np.mean(mi_surr >= mi_obs)
    return p


def bh_fdr_mapa(p_mapa, alpha=0.05):
    """
    Benjamini-Hochberg sobre TODAS as células do mapa -- a família aqui é o
    próprio comodulograma ("esta célula se destaca no mapa DESTA janela?").
    Válido sob dependência positiva (PRDS): células vizinhas compartilham os
    mesmos sinais filtrados, que é exatamente o caso aqui.
    """
    p = np.asarray(p_mapa, dtype=float)
    forma = p.shape
    p_flat = p.ravel()
    n = p_flat.size
    ordem = np.argsort(p_flat)
    passou = p_flat[ordem] <= (np.arange(1, n + 1) / n) * alpha
    sig = np.zeros(n, dtype=bool)
    if np.any(passou):
        corte = np.max(np.where(passou))
        sig[ordem[:corte + 1]] = True
    return sig.reshape(forma)


def resume_cluster_fdr(mascara_sig, fases_freq, amps_freq, par="theta_gamma"):
    """
    Quantifica o padrao das celulas significativas dentro do quadrante
    do par ativo (nao sempre Theta-Gamma). O discriminador e' a fracao
    de celulas sig dentro do quadrante correto (frac_sig_par).
    """
    cfg = BAND_PAIRS.get(par, BAND_PAIRS["theta_gamma"])
    fase_band = cfg["fase"]
    amp_band  = cfg["amp"]

    n_sig = int(mascara_sig.sum())
    if n_sig == 0:
        return {"n_sig": 0, "n_sig_par": 0, "frac_sig_par": 0.0,
                "maior_cluster": 0}

    mask_par = ((fases_freq[None, :] >= fase_band[0]) &
                (fases_freq[None, :] <= fase_band[1]) &
                (amps_freq[:, None]  >= amp_band[0])  &
                (amps_freq[:, None]  <= amp_band[1]))
    n_sig_par = int((mascara_sig & mask_par).sum())

    rotulos, _ = ndimage.label(mascara_sig)
    maior_cluster = 0
    if rotulos.max() > 0:
        tamanhos = ndimage.sum(mascara_sig, rotulos,
                               range(1, rotulos.max() + 1))
        maior_cluster = int(np.max(tamanhos))

    return {"n_sig": n_sig, "n_sig_par": n_sig_par,
            "frac_sig_par": n_sig_par / n_sig,
            "maior_cluster": maior_cluster}


def z_pico_par(z_mapa, fases_freq, amps_freq, par="theta_gamma"):
    """Maior z dentro do quadrante do par ativo + onde ocorre."""
    cfg = BAND_PAIRS.get(par, BAND_PAIRS["theta_gamma"])
    mask_fase = (fases_freq >= cfg["fase"][0]) & (fases_freq <= cfg["fase"][1])
    mask_amp  = (amps_freq  >= cfg["amp"][0])  & (amps_freq  <= cfg["amp"][1])
    sub = z_mapa[np.ix_(mask_amp, mask_fase)]
    ii, jj = np.unravel_index(np.argmax(sub), sub.shape)
    return (
        float(sub[ii, jj]),
        float(fases_freq[mask_fase][jj]),
        float(amps_freq[mask_amp][ii]),
    )


# ==========================================
# VISUALIZAÇÃO
# ==========================================

def plota_comodulograma_z(z_mapa, fases_freq, amps_freq, titulo, caminho_png,
                          par="theta_gamma", mascara_fdr=None):
    """
    Heatmap divergente (RdBu_r, centro em z=0). Destaca em tracejado branco
    o quadrante do par ativo (nao sempre Theta-Gamma).
    """
    cfg = BAND_PAIRS.get(par, BAND_PAIRS["theta_gamma"])
    fase_band = cfg["fase"]
    amp_band  = cfg["amp"]

    vmax = max(float(np.max(np.abs(z_mapa))), 1.0)

    fig, ax = plt.subplots(figsize=(10, 7))
    X, Y = np.meshgrid(fases_freq, amps_freq)
    pcm = ax.pcolormesh(X, Y, z_mapa, shading="auto", cmap="RdBu_r",
                        norm=mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax))
    cb = fig.colorbar(pcm, ax=ax)
    cb.set_label("z do MI (vs surrogates)", fontsize=11)

    if mascara_fdr is not None and mascara_fdr.any():
        ax.contour(X, Y, mascara_fdr.astype(float), levels=[0.5],
                   colors="black", linewidths=1.6)

    # Quadrante do par ativo (tracejado branco)
    for x in fase_band:
        ax.axvline(x, color="white", linestyle="--", linewidth=0.9, alpha=0.7)
    for y in amp_band:
        ax.axhline(y, color="white", linestyle="--", linewidth=0.9, alpha=0.7)

    # Linhas horizontais separando os 3 pares (cinza claro)
    for sep in [80, 150]:
        if sep <= amps_freq[-1]:
            ax.axhline(sep, color="#cccccc", linestyle=":", linewidth=0.7, alpha=0.5)
    ax.text(fases_freq[-1] * 0.98, 55,  "Gamma",  color="#ff7f0e", ha="right", fontsize=8)
    ax.text(fases_freq[-1] * 0.98, 110, "HG",     color="#9467bd", ha="right", fontsize=8)
    if amps_freq[-1] > 150:
        ax.text(fases_freq[-1] * 0.98, 200, "HFO", color="#e377c2", ha="right", fontsize=8)

    ax.set_title(titulo, fontsize=12)
    ax.set_xlabel("Frequencia da Fase (Hz)", fontsize=11)
    ax.set_ylabel("Frequencia da Amplitude (Hz)", fontsize=11)
    ax.set_xticks(fases_freq[::2])
    ax.set_yticks(np.arange(30, amps_freq[-1] + 1, 20))

    fig.tight_layout()
    fig.savefig(caminho_png, dpi=150)
    plt.close(fig)


# ==========================================
# MODO LOTE
# ==========================================

def roda_lote(args):
    print(f"Lendo {args.csv} ...")
    df = pd.read_csv(args.csv)

    # Filtra por par se especificado
    if args.par and "par" in df.columns:
        df = df[df["par"] == args.par].copy()
        print(f"  -> {len(df)} candidatos do par '{args.par}'")
    elif args.par and "par" not in df.columns:
        print(f"  [info] Coluna 'par' nao encontrada; usando todos os candidatos.")

    if "veredito" in df.columns and args.veredito_prefixo:
        df = df[df["veredito"].astype(str).str.startswith(args.veredito_prefixo)].copy()
    print(f"  -> {len(df)} candidatos apos filtro de veredito ('{args.veredito_prefixo}')")

    if args.top_n:
        df = df.nlargest(args.top_n, "z_score_refinado")
        print(f"  -> mantendo os {len(df)} melhores por z_score_refinado")

    if len(df) == 0:
        print("Nada a processar.")
        return

    os.makedirs(args.saida_dir, exist_ok=True)

    fases_freq = FASES_DEFAULT
    amps_freq  = AMPS_DEFAULT
    rng = np.random.default_rng(42)
    par_ativo = args.par

    resumo = []
    for arquivo, grupo in df.groupby("arquivo"):
        caminho = os.path.join(args.pasta_ns2, arquivo)
        print(f"\nCarregando {arquivo} ({len(grupo)} candidatos) ...")
        dados, fs, nomes_canais = le_ns2(caminho)
        mapa_canal = {str(nome): i for i, nome in enumerate(nomes_canais)}

        for n_i, row in enumerate(grupo.itertuples(), 1):
            canal_idx = mapa_canal.get(str(row.canal))
            if canal_idx is None:
                print(f"  aviso: canal {row.canal} não encontrado em {arquivo}, pulando")
                continue

            lfp = fatia_janela(
                dados[:, canal_idx], fs, row.janela_ini_s, row.janela_fim_s
            ).astype(float)

            if args.fdr_q:
                z_mapa, mi_obs_mapa, mi_surr_mapa = calcula_comodulograma_z(
                    lfp, fs, fases_freq, amps_freq,
                    n_surr=args.n_surr, rng=rng, notch_hz=args.notch,
                    retorna_mi=True,
                )
            else:
                z_mapa = calcula_comodulograma_z(
                    lfp, fs, fases_freq, amps_freq,
                    n_surr=args.n_surr, rng=rng, notch_hz=args.notch,
                )
            z_pico, f_pico, a_pico = z_pico_par(z_mapa, fases_freq, amps_freq, par=par_ativo)
            
            # MI bruto na celula do pico do par ativo
            i_fp = int(np.argmin(np.abs(fases_freq - f_pico)))
            j_ap = int(np.argmin(np.abs(amps_freq  - a_pico)))
            mi_pico = float(mi_obs_mapa[j_ap, i_fp]) if args.fdr_q else float('nan')

            mascara_fdr = None
            info_fdr = None
            if args.fdr_q:
                p_mapa = p_valores_por_celula(mi_obs_mapa, mi_surr_mapa)
                mascara_fdr = bh_fdr_mapa(p_mapa, alpha=args.fdr_q)
                info_fdr = resume_cluster_fdr(mascara_fdr, fases_freq, amps_freq,
                                               par=par_ativo)
                if info_fdr["n_sig"] == 0:
                    classe_fdr = "nada sobrevive ao FDR"
                elif info_fdr["n_sig"] >= 2 and info_fdr["frac_sig_par"] >= 0.5:
                    classe_fdr = f"concentrado em {par_ativo}"
                else:
                    classe_fdr = f"esparso/fora de {par_ativo}"

            stem = os.path.splitext(arquivo)[0]
            nome_png = (f"{stem}_{row.canal}_{row.janela_ini_s:g}-"
                        f"{row.janela_fim_s:g}s_{par_ativo}_zcomodo.png")
            caminho_png = os.path.join(args.saida_dir, nome_png)

            titulo = (f"{row.canal} @ {row.janela_ini_s:g}-{row.janela_fim_s:g}s | {arquivo}\n"
                      f"par: {par_ativo} | pico z={z_pico:.2f} ({f_pico:g} Hz x {a_pico:g} Hz)")
            if args.fdr_q and info_fdr:
                titulo += (f"\nFDR q={args.fdr_q:g}: {info_fdr['n_sig']} celulas sig "
                           f"({info_fdr['n_sig_par']} no quadrante), "
                           f"maior cluster={info_fdr['maior_cluster']} | {classe_fdr}")
            plota_comodulograma_z(z_mapa, fases_freq, amps_freq, titulo, caminho_png,
                                  par=par_ativo, mascara_fdr=mascara_fdr)

            linha_resumo = {
                "arquivo": arquivo,
                "canal": row.canal,
                "janela_ini_s": row.janela_ini_s,
                "janela_fim_s": row.janela_fim_s,
                "par": par_ativo,
                "z_score_refinado": getattr(row, 'z_score_refinado', float('nan')),
                "z_pico_par": z_pico,
                "fase_pico_hz": f_pico,
                "amp_pico_hz": a_pico,
                "mi_pico": mi_pico,
                "png": nome_png,
            }
            if args.fdr_q:
                linha_resumo.update({
                    "n_sig_fdr": info_fdr["n_sig"],
                    "n_sig_fdr_par": info_fdr["n_sig_par"],
                    "maior_cluster_fdr": info_fdr["maior_cluster"],
                    "classe_fdr": classe_fdr,
                })
            resumo.append(linha_resumo)

            sufixo = f" | FDR: {classe_fdr} ({info_fdr['n_sig']} células)" if args.fdr_q else ""
            print(f"  [{n_i}/{len(grupo)}] {nome_png}  (pico {par_ativo} z={z_pico:.2f}{sufixo})", end="\r")

    resumo_df = pd.DataFrame(resumo).sort_values("z_pico_par", ascending=False)
    caminho_resumo = os.path.join(args.saida_dir, "resumo_comodulogramas.csv")
    resumo_df.to_csv(caminho_resumo, index=False)

    print(f"\n{len(resumo_df)} comodulogramas salvos em {args.saida_dir}/")
    print(f"Resumo salvo: {caminho_resumo}")
    if args.fdr_q and "classe_fdr" in resumo_df.columns:
        print("\nDistribuicao das classes FDR:")
        print(resumo_df["classe_fdr"].value_counts().to_string())
    print(f"\nTop 10 por pico z no quadrante {par_ativo}:")
    print(resumo_df.head(10)[[
        "canal", "janela_ini_s", "janela_fim_s",
        "z_pico_par", "fase_pico_hz", "amp_pico_hz", "png"
    ]].to_string(index=False))


# ==========================================
# MODO JANELA ÚNICA
# ==========================================

def roda_janela_unica(args):
    print(f"Carregando {args.arquivo} ...")
    dados_brutos, fs, canal_ids = carrega_dados(
        args.arquivo, n_canais_bin=args.n_canais, fs_bin=args.fs
    )
    print(f"  -> {dados_brutos.shape[1]} canais, fs={fs} Hz")

    dados_janela = fatia_janela(dados_brutos, fs, args.inicio, args.fim)
    lfp_ativo = dados_janela[:, args.canal].astype(float)
    nome_canal = canal_ids[args.canal] if args.canal < len(canal_ids) else args.canal

    fases_freq = FASES_DEFAULT
    amps_freq  = AMPS_DEFAULT
    par_ativo  = args.par

    print(f"Calculando mapa z-scoredo para canal {nome_canal} "
          f"({args.inicio:.0f}-{args.fim:.0f}s) | par={par_ativo} ...")
    if args.fdr_q:
        z_mapa, mi_obs_mapa, mi_surr_mapa = calcula_comodulograma_z(
            lfp_ativo, fs, fases_freq, amps_freq,
            n_surr=args.n_surr, rng=np.random.default_rng(42), notch_hz=args.notch,
            retorna_mi=True,
        )
    else:
        z_mapa = calcula_comodulograma_z(
            lfp_ativo, fs, fases_freq, amps_freq,
            n_surr=args.n_surr, rng=np.random.default_rng(42), notch_hz=args.notch,
        )

    z_pico, f_pico, a_pico = z_pico_par(z_mapa, fases_freq, amps_freq, par=par_ativo)
    titulo = (f"Comodulograma (z vs surrogates) | par={par_ativo}\n"
              f"Canal {nome_canal}, {os.path.basename(args.arquivo)}, "
              f"{args.inicio:.0f}-{args.fim:.0f}s | pico z={z_pico:.2f} "
              f"({f_pico:g}x{a_pico:g} Hz)")

    mascara_fdr = None
    if args.fdr_q:
        p_mapa = p_valores_por_celula(mi_obs_mapa, mi_surr_mapa)
        mascara_fdr = bh_fdr_mapa(p_mapa, alpha=args.fdr_q)
        info = resume_cluster_fdr(mascara_fdr, fases_freq, amps_freq, par=par_ativo)
        if info["n_sig"] == 0:
            classe = "nada sobrevive ao FDR"
        elif info["n_sig"] >= 2 and info["frac_sig_par"] >= 0.5:
            classe = f"concentrado em {par_ativo}"
        else:
            classe = f"esparso/fora de {par_ativo}"
        titulo += (f"\nFDR q={args.fdr_q:g}: {info['n_sig']} celulas sig "
                   f"({info['n_sig_par']} no quadrante) | {classe}")
        print(f"FDR q={args.fdr_q:g}: {info}")

    caminho_png = args.saida_png or f"comodulograma_{par_ativo}_janela_unica.png"
    plota_comodulograma_z(z_mapa, fases_freq, amps_freq, titulo, caminho_png,
                          par=par_ativo, mascara_fdr=mascara_fdr)
    print(f"Pico no quadrante {par_ativo}: z={z_pico:.2f} ({f_pico:g} Hz x {a_pico:g} Hz)")
    plt.show()


# ==========================================
# MAIN
# ==========================================

def main():
    try:  # console Windows pode estar em cp1252; Θ/Γ quebram o print
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arquivo",
                    help="(janela única) Caminho do .ns2 (direto) ou .bin (legado)")
    ap.add_argument("--csv",
                    help="(modo lote) CSV gerado pelo refina_candidatos.py")
    ap.add_argument("--pasta_ns2", default=".",
                    help="(modo lote) Pasta com os .ns2 originais")
    ap.add_argument("--par",
                    default="theta_gamma",
                    choices=list(BAND_PAIRS.keys()),
                    help="Par a destacar no mapa: theta_gamma (default), "
                         "theta_hg, theta_hfo. Nao altera o calculo — "
                         "apenas o quadrante marcado e o z do pico reportado.")
    ap.add_argument("--veredito_prefixo", default="Candidato robusto",
                    help="(modo lote) Prefixo do veredito a filtrar "
                         "(default: 'Candidato robusto'; vazio = todos)")
    ap.add_argument("--top_n", type=int, default=None,
                    help="(modo lote) Processar só os N melhores por z_score_refinado")
    ap.add_argument("--saida_dir", default="comodulogramas",
                    help="(modo lote) Diretório de saída dos PNGs")
    ap.add_argument("--canal", type=int, default=0,
                    help="(janela única) Índice do canal (default: 0)")
    ap.add_argument("--n_canais", type=int, default=1,
                    help="(janela única .bin legado) Nº de canais no arquivo")
    ap.add_argument("--fs", type=float, default=1000.0,
                    help="(janela única .bin legado) Taxa de amostragem")
    ap.add_argument("--inicio", type=float, default=130.0, help="Início da janela (s)")
    ap.add_argument("--fim", type=float, default=140.0, help="Fim da janela (s)")
    ap.add_argument("--n_surr", type=int, default=200, help="Número de surrogates por par")
    ap.add_argument("--notch", type=float, nargs="+", default=[60.0, 120.0, 180.0, 240.0], metavar="HZ",
                    help="Frequência(s) da rede elétrica a notchar "
                         "(ex.: --notch 60 120 180 240).")
    ap.add_argument("--fdr_q", type=float, default=None, metavar="Q",
                    help="Se informado (ex.: 0.05), aplica Benjamini-Hochberg "
                         "sobre as células do mapa, contorna as significantes "
                         "no PNG e classifica: cluster focal ΘΓ vs esparso.")
    ap.add_argument("--saida_png", default=None,
                    help="(janela única) Onde salvar o PNG antes de mostrar")
    args = ap.parse_args()

    if args.csv:
        if not args.pasta_ns2 or args.pasta_ns2 == ".":
            print("Erro: no modo lote forneça --pasta_ns2 com os .ns2 originais.")
            raise SystemExit(1)
        roda_lote(args)
    elif args.arquivo:
        roda_janela_unica(args)
    else:
        print("Erro: forneça --csv (modo lote) ou --arquivo (janela única).")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
