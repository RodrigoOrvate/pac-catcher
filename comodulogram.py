"""
comodulogram.py (atualizado)
==========================================
Comodulograma Theta-Gamma -- agora lendo .ns2 diretamente, sem
precisar do extrator.exe. Ainda aceita .bin legado.

DOIS MODOS:

  1. MODO LOTE (--csv): lê o resultados_refinados.csv produzido pelo
     refina_candidatos.py e gera um PNG por candidato, com o MI
     Z-SCOREADO contra surrogates de deslocamento circular (a mesma
     nula do triagem_pac.py). Motivo: MI bruto não tem escala
     interpretável -- 0.03 é muito ou pouco? -- e o mapa bruto pinta
     qualquer ruído como "hot". O z-score dá origem significante ao
     heatmap (z=0 == nível de acaso, negativo = abaixo do acaso), o
     que pede colormap DIVERGENTE (RdBu_r centrado em 0) em vez do
     jet arco-íris: azul = abaixo do acaso, vermelho = acima,
     branco/claro = nada.

  2. JANELA ÚNICA: exploração rápida de um trecho específico, agora
     também z-scoredo, pelos mesmos motivos (era o MI bruto criticado
     na cabeça do triagem_pac.py).

EM AMBOS OS MODOS, --fdr_q 0.05 adiciona a correção de múltiplas
comparações SOBRE O MAPA: p-valor por célula (Gama ajustada aos
surrogates da própria célula), Benjamini-Hochberg sobre as 275 células,
contorno preto nas células significantes e classificação do padrão --
"cluster focal ΘΓ" (acoplamento) vs "esparso/fora de ΘΓ" (transientes
ritmados, como na coluna de 8 Hz dos canais ruins).

Uso:
    # Modo lote: um PNG por 'Candidato robusto' do refinamento
    python comodulogram.py --csv resultados_refinados.csv \
        --pasta_ns2 "../Basal antes da infusao" \
        --saida_dir comodulogramas --n_surr 200

    # Só os N melhores por z_score_refinado da triagem
    python comodulogram.py --csv resultados_refinados.csv \
        --pasta_ns2 "../Basal antes da infusao" --top_n 10

    # Outro veredito (ex.: revisar os de ruído comum suspeito)
    python comodulogram.py --csv resultados_refinados.csv \
        --pasta_ns2 "../Basal antes da infusao" \
        --veredito_prefixo "Revisar"

    # Janela única, canal 17 (índice; chan18), janela 205-215s
    python comodulogram.py --arquivo "../Basal antes da infusao/20240708-123605-003.ns2" \
        --canal 17 --inicio 205 --fim 215

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

from ns2_utils import le_ns2, carrega_dados, fatia_janela


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


def aplica_notch(sinal_in, fs, linha_hz=60.0, harmonicos=2, q_factor=30.0):
    """
    Rejeita a frequência da rede elétrica e seus harmônicos (iirnotch).
    Necessário porque um pico de 60 Hz (e o harmônico em 120 Hz) dentro da
    banda de amplitude infla o MI nas células 55/60/65 Hz do comodulograma
    e pode simular acoplamento Theta-Gamma onde não há.
    """
    nyq = 0.5 * fs
    out = sinal_in
    for h in range(1, harmonicos + 1):
        f = linha_hz * h
        if f >= nyq * 0.98:
            break
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
        lfp_ativo = aplica_notch(lfp_ativo, fs, linha_hz=notch_hz)

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


def resume_cluster_fdr(mascara_sig, fases_freq, amps_freq,
                       theta_band=(4, 8), gamma_band=(30, 80)):
    """
    Quantifica o PADRÃO das células significativas: acoplamento genuíno e
    estreito ocupa POUCAS células, mas concentradas no quadrante Theta-Gamma;
    transientes ritmados espalham células significativas por várias amplitudes
    (a "coluna" de 8 Hz). O discriminador é a fração de células significativas
    dentro de ΘΓ (frac_sig_tg), não o tamanho do cluster.
    """
    n_sig = int(mascara_sig.sum())
    if n_sig == 0:
        return {"n_sig": 0, "n_sig_tg": 0, "frac_sig_tg": 0.0,
                "maior_cluster": 0}

    mask_tg = ((fases_freq[None, :] >= theta_band[0]) &
               (fases_freq[None, :] <= theta_band[1]) &
               (amps_freq[:, None] >= gamma_band[0]) &
               (amps_freq[:, None] <= gamma_band[1]))
    n_sig_tg = int((mascara_sig & mask_tg).sum())

    rotulos, _ = ndimage.label(mascara_sig)
    maior_cluster = 0
    if rotulos.max() > 0:
        tamanhos = ndimage.sum(mascara_sig, rotulos,
                               range(1, rotulos.max() + 1))
        maior_cluster = int(np.max(tamanhos))

    return {"n_sig": n_sig, "n_sig_tg": n_sig_tg,
            "frac_sig_tg": n_sig_tg / n_sig,
            "maior_cluster": maior_cluster}


def z_pico_theta_gamma(z_mapa, fases_freq, amps_freq,
                       theta_band=(4, 8), gamma_band=(30, 80)):
    """Maior z dentro do retângulo clássico Theta-Gamma + onde ele ocorre."""
    mask_fase = (fases_freq >= theta_band[0]) & (fases_freq <= theta_band[1])
    mask_amp = (amps_freq >= gamma_band[0]) & (amps_freq <= gamma_band[1])
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
                          theta_band=(4, 8), gamma_band=(30, 80), mascara_fdr=None):
    """
    Heatmap divergente (RdBu_r, centro em z=0): vermelho = acoplamento
    acima do acaso, azul = abaixo, claro = nada. O piso da escala é ±1
    para não estourar quando o mapa todo for ruído puro.
    mascara_fdr: se informado, contorna em preto as células significantes
    após a correção Benjamini-Hochberg sobre o mapa.
    """
    vmax = max(float(np.max(np.abs(z_mapa))), 1.0)

    fig, ax = plt.subplots(figsize=(9, 6))
    X, Y = np.meshgrid(fases_freq, amps_freq)
    pcm = ax.pcolormesh(X, Y, z_mapa, shading="auto", cmap="RdBu_r",
                        norm=mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax))
    cb = fig.colorbar(pcm, ax=ax)
    cb.set_label("z do MI (vs surrogates)", fontsize=11)

    if mascara_fdr is not None and mascara_fdr.any():
        ax.contour(X, Y, mascara_fdr.astype(float), levels=[0.5],
                   colors="black", linewidths=1.6)

    # Caixa tracejada marcando o quadrante Theta-Gamma clássico
    for x in theta_band:
        ax.axvline(x, color="white", linestyle="--", linewidth=0.8, alpha=0.6)
    for y in gamma_band:
        ax.axhline(y, color="white", linestyle="--", linewidth=0.8, alpha=0.6)

    ax.set_title(titulo, fontsize=12)
    ax.set_xlabel("Frequência da Fase (Hz)", fontsize=11)
    ax.set_ylabel("Frequência da Amplitude (Hz)", fontsize=11)
    ax.set_xticks(fases_freq)
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

    if "veredito" in df.columns and args.veredito_prefixo:
        df = df[df["veredito"].astype(str).str.startswith(args.veredito_prefixo)].copy()
    print(f"  -> {len(df)} candidatos após filtro de veredito "
          f"('{args.veredito_prefixo}')")

    if args.top_n:
        df = df.nlargest(args.top_n, "z_score_refinado")
        print(f"  -> mantendo os {len(df)} melhores por z_score_refinado")

    if len(df) == 0:
        print("Nada a processar.")
        return

    os.makedirs(args.saida_dir, exist_ok=True)

    fases_freq = np.arange(4, 15, 1)
    amps_freq = np.arange(30, 155, 5)
    rng = np.random.default_rng(42)

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
                z_mapa, mi_obs_mapa, _ = calcula_comodulograma_z(
                    lfp, fs, fases_freq, amps_freq,
                    n_surr=args.n_surr, rng=rng, notch_hz=args.notch,
                    retorna_mi=True,
                )
            z_pico, f_pico, a_pico = z_pico_theta_gamma(z_mapa, fases_freq, amps_freq)
            # MI bruto (KL-MI, sem normalizar) na célula do pico ΘΓ
            i_fp = int(np.argmin(np.abs(fases_freq - f_pico)))
            j_ap = int(np.argmin(np.abs(amps_freq - a_pico)))
            mi_pico = float(mi_obs_mapa[j_ap, i_fp])

            mascara_fdr = None
            info_fdr = None
            if args.fdr_q:
                p_mapa = p_valores_por_celula(mi_obs_mapa, mi_surr_mapa)
                mascara_fdr = bh_fdr_mapa(p_mapa, alpha=args.fdr_q)
                info_fdr = resume_cluster_fdr(mascara_fdr, fases_freq, amps_freq)
                if info_fdr["n_sig"] == 0:
                    classe_fdr = "nada sobrevive ao FDR"
                elif (info_fdr["n_sig"] >= 2
                      and info_fdr["frac_sig_tg"] >= 0.5):
                    classe_fdr = "concentrado em ΘΓ"
                else:
                    classe_fdr = "esparso/fora de ΘΓ"

            stem = os.path.splitext(arquivo)[0]
            nome_png = f"{stem}_{row.canal}_{row.janela_ini_s:g}-{row.janela_fim_s:g}s_zcomodo.png"
            caminho_png = os.path.join(args.saida_dir, nome_png)

            titulo = (f"{row.canal} @ {row.janela_ini_s:g}–{row.janela_fim_s:g}s — {arquivo}\n"
                      f"pico ΘΓ: z={z_pico:.2f} ({f_pico:g} Hz × {a_pico:g} Hz)")
            if args.fdr_q:
                titulo += (f"\nFDR q={args.fdr_q:g}: {info_fdr['n_sig']} células sig "
                           f"({info_fdr['n_sig_tg']} em ΘΓ), maior cluster="
                           f"{info_fdr['maior_cluster']} — {classe_fdr}")
            plota_comodulograma_z(z_mapa, fases_freq, amps_freq, titulo, caminho_png,
                                  mascara_fdr=mascara_fdr)

            linha_resumo = {
                "arquivo": arquivo,
                "canal": row.canal,
                "janela_ini_s": row.janela_ini_s,
                "janela_fim_s": row.janela_fim_s,
                "z_score_refinado": row.z_score_refinado,
                "z_pico_theta_gamma": z_pico,
                "fase_pico_hz": f_pico,
                "amp_pico_hz": a_pico,
                "mi_pico": mi_pico,
                "png": nome_png,
            }
            if args.fdr_q:
                linha_resumo.update({
                    "n_sig_fdr": info_fdr["n_sig"],
                    "n_sig_fdr_theta_gamma": info_fdr["n_sig_tg"],
                    "maior_cluster_fdr": info_fdr["maior_cluster"],
                    "classe_fdr": classe_fdr,
                })
            resumo.append(linha_resumo)

            sufixo = f" | FDR: {classe_fdr} ({info_fdr['n_sig']} células)" if args.fdr_q else ""
            print(f"  [{n_i}/{len(grupo)}] {nome_png}  (pico ΘΓ z={z_pico:.2f}{sufixo})", end="\r")

    resumo_df = pd.DataFrame(resumo).sort_values("z_pico_theta_gamma", ascending=False)
    caminho_resumo = os.path.join(args.saida_dir, "resumo_comodulogramas.csv")
    resumo_df.to_csv(caminho_resumo, index=False)

    print(f"\n\n{len(resumo_df)} comodulogramas salvos em {args.saida_dir}/")
    print(f"Resumo salvo: {caminho_resumo}")
    if args.fdr_q and "classe_fdr" in resumo_df.columns:
        print("\nDistribuição das classes FDR:")
        print(resumo_df["classe_fdr"].value_counts().to_string())
    print("\nTop 10 por pico z dentro do quadrante Theta-Gamma:")
    print(resumo_df.head(10)[[
        "canal", "janela_ini_s", "janela_fim_s",
        "z_pico_theta_gamma", "fase_pico_hz", "amp_pico_hz", "png"
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

    fases_freq = np.arange(4, 15, 1)
    amps_freq = np.arange(30, 155, 5)

    print(f"Calculando o mapa z-scoredo para o canal {nome_canal} "
          f"({args.inicio:.0f}-{args.fim:.0f}s) ...")
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

    z_pico, f_pico, a_pico = z_pico_theta_gamma(z_mapa, fases_freq, amps_freq)
    titulo = (f"Comodulograma Theta-Gamma (z vs surrogates)\n"
              f"Canal {nome_canal}, {os.path.basename(args.arquivo)}, "
              f"{args.inicio:.0f}-{args.fim:.0f}s — pico ΘΓ: z={z_pico:.2f}")

    mascara_fdr = None
    if args.fdr_q:
        p_mapa = p_valores_por_celula(mi_obs_mapa, mi_surr_mapa)
        mascara_fdr = bh_fdr_mapa(p_mapa, alpha=args.fdr_q)
        info = resume_cluster_fdr(mascara_fdr, fases_freq, amps_freq)
        if info["n_sig"] == 0:
            classe = "nada sobrevive ao FDR"
        elif info["n_sig"] >= 2 and info["frac_sig_tg"] >= 0.5:
            classe = "concentrado em ΘΓ"
        else:
            classe = "esparso/fora de ΘΓ"
        titulo += (f"\nFDR q={args.fdr_q:g}: {info['n_sig']} células sig "
                   f"({info['n_sig_tg']} em ΘΓ), maior cluster="
                   f"{info['maior_cluster']} — {classe}")
        print(f"FDR q={args.fdr_q:g}: {info}")

    plota_comodulograma_z(z_mapa, fases_freq, amps_freq, titulo,
                          args.saida_png or "comodulograma_janela_unica.png",
                          mascara_fdr=mascara_fdr)
    print(f"Pico no quadrante ΘΓ: z={z_pico:.2f} ({f_pico:g} Hz × {a_pico:g} Hz)")
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
    ap.add_argument("--notch", type=float, default=None, metavar="HZ",
                    help="Frequência da rede elétrica a notchar (ex.: 60), "
                         "com harmônicos. Rode com e sem para comparar.")
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
