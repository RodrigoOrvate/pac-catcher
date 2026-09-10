"""
figura_apresentacao.py
==========================================
Etapa 5 (opcional) da validação PAC: figura de apresentação dos VENCEDORES.

Para cada janela/canal vencedor (etapas 1-3), uma figura com:
  ESQUERDA (séries alinhadas no tempo, janela completa):
    1. LFP bruto (com notch 60 Hz)
    2. theta filtrado no par de pico (fase_pico ± 1 Hz)
    3. gamma filtrado no par de pico (amp_pico ± 5 Hz) + envelope
    4. espectrograma (STFT) com as bandas de theta e gamma marcadas
  DIREITA:
    5. comodulograma z-scoredo (mesmo cálculo do comodulogram.py) com o
       par de pico marcado
    6. distribuição POLAR fase×amplitude (o clássico do Tort et al. 2010):
       amplitude média de gamma em 18 bins de fase do theta, normalizada
       pela média -- a "modulação" que o MI quantifica

Todos os números (MI, z) usam a mesma nula de deslocamento circular do
pipeline; o comodulograma é recalculado aqui com notch.

Uso (SCRIPT central do estudo — uma chamada por sessão):
    python figura_apresentacao.py \
        --pasta_ns2 "<sessao>/Basal antes da infusao" \
        --vencedores "<sessao>/vencedores.csv" \
        [--saida_dir "<sessao>/figuras"]

vencedores.csv (uma linha por vencedor; colunas):
    rotulo,arquivo,canal,inicio_s,fim_s,fase_pico_hz,amp_pico_hz[,comportamento]
(comportamento é opcional; entra no título da figura)

Saída: <saida_dir>/<rotulo>_<canal>.png
Requer: neo, numpy, scipy, matplotlib
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import scipy.signal as signal
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from pipeline.etapa3_comodulograma.comodulogram import (calcula_comodulograma_z,
                          _mi_de_bin_idx, FASES_DEFAULT, AMPS_DEFAULT)
from pipeline.etapa7_validacao.robustez_parametros import mi_z_par
from pac_core.io import le_ns2, fatia_janela
from pac_core.filtering import filtra_sinal, aplica_notch

N_SURR = 200
N_BINS = 18


def plota_comodo_no_eixo(ax, z_mapa, fases_freq, amps_freq, f_pico, a_pico):
    """Comodulograma z num eixo dado: RdBu_r centrado em 0, caixa do par de
    pico e marcador no pico."""
    vmax = max(np.abs(z_mapa).max(), 1.0)
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    pc = ax.pcolormesh(fases_freq, amps_freq, z_mapa, cmap="RdBu_r",
                       norm=norm, shading="auto", rasterized=True)
    cb = ax.figure.colorbar(pc, ax=ax, pad=0.02)
    cb.set_label("z do MI (vs surrogates)", fontsize=9)
    cb.ax.tick_params(labelsize=8)
    # caixa do par de pico (bandas usadas nas séries temporais), recortada
    # ao domínio da grade (p.ex. pico em 5 Hz x 35 Hz encosta nas bordas)
    x0 = max(f_pico - 1.0, fases_freq.min())
    y0 = max(a_pico - 5.0, amps_freq.min())
    x1 = min(f_pico + 1.0, fases_freq.max())
    y1 = min(a_pico + 5.0, amps_freq.max())
    ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0,
                               fill=False, ec="black", lw=1.2, ls="--"))
    ax.plot(f_pico, a_pico, marker="+", color="black", ms=14, mew=2)
    ax.set_xlabel("Frequência da fase (Hz)", fontsize=10)
    ax.set_ylabel("Frequência da amplitude (Hz)", fontsize=10)
    ax.tick_params(labelsize=9)
    return pc


def plota_polar(ax, fase_theta, env_gamma, n_bins=N_BINS):
    """Distribuição polar fase×amplitude (Tort et al. 2010): amplitude média
    de gamma por bin de fase do theta, normalizada pela média geral.
    Devolve o MI normalizado (mesma fórmula do pipeline)."""
    bins = np.linspace(-np.pi, np.pi, n_bins + 1)
    bin_idx = np.clip(np.digitize(fase_theta, bins) - 1, 0, n_bins - 1)
    mi = _mi_de_bin_idx(bin_idx, env_gamma, n_bins)

    medias = np.array([env_gamma[bin_idx == k].mean() for k in range(n_bins)])
    razoes = medias / medias.mean()

    centros = np.linspace(-np.pi, np.pi, n_bins, endpoint=False) + np.pi / n_bins
    cores = plt.cm.Purples(0.45 + 0.4 * (razoes - razoes.min()) /
                           (razoes.max() - razoes.min() + 1e-12))
    ax.bar(centros, razoes, width=2 * np.pi / n_bins * 0.92, bottom=0,
           color=cores, edgecolor="white", linewidth=0.5, align="center")

    # referência: razão = 1 (amplitude média; sem modulação = círculo plano)
    ref = np.linspace(0, 2 * np.pi, 200)
    ax.plot(ref, np.ones_like(ref), "--", color="#555555", lw=1.0)
    ax.set_theta_zero_location("E")
    ax.set_theta_direction(1)
    ax.set_xticks(np.linspace(0, 2 * np.pi, 5, endpoint=False))
    ax.set_xticklabels(["0", "π/2", "π", "3π/2", ""], fontsize=9)
    ax.set_ylim(0, max(1.25, razoes.max() * 1.1))
    ax.tick_params(axis="y", labelsize=8)
    ax.set_ylabel("Amplitude de γ normalizada", fontsize=9, labelpad=25)
    return mi


def figura_vencedor(rotulo, arquivo, canal, inicio, fim, f_pico, a_pico,
                    comportamento, pasta_ns2, par_nome="theta_gamma", saida_dir="figuras"):
    rng = np.random.default_rng(42)
    caminho = f"{pasta_ns2}/{arquivo}"
    print(f"Gerando figura: {canal} @ {inicio:g}-{fim:g}s ...")

    dados, fs, nomes = le_ns2(caminho)
    idx = {str(n): i for i, n in enumerate(nomes)}[canal]
    lfp = aplica_notch(fatia_janela(dados[:, idx], fs, inicio, fim).astype(float),
                       fs, 60.0)
    del dados
    t = np.arange(lfp.size) / fs

    theta = filtra_sinal(lfp, f_pico - 1.0, f_pico + 1.0, fs)
    gamma = filtra_sinal(lfp, a_pico - 5.0, a_pico + 5.0, fs)
    env_gamma = np.abs(signal.hilbert(gamma))
    fase_theta = np.angle(signal.hilbert(theta))

    z_mi = mi_z_par(lfp, fs, f_pico, a_pico, rng=rng)

    # comodulograma completo (mesmo cálculo do lote; notch já aplicado)
    fases_freq = FASES_DEFAULT
    amps_freq = AMPS_DEFAULT
    z_mapa = calcula_comodulograma_z(lfp, fs, fases_freq, amps_freq,
                                     n_surr=N_SURR, n_bins=N_BINS, rng=rng,
                                     notch_hz=None)

    # ----------------------------------
    # FIGURA
    # ----------------------------------
    fig = plt.figure(figsize=(13, 9.5))
    gs = fig.add_gridspec(4, 2, width_ratios=[1.75, 1.0],
                          height_ratios=[1, 1, 1, 1.35],
                          hspace=0.55, wspace=0.28)
    cor_txt = "#333333"
    ax_raw = fig.add_subplot(gs[0, 0])
    ax_th = fig.add_subplot(gs[1, 0], sharex=ax_raw)
    ax_ga = fig.add_subplot(gs[2, 0], sharex=ax_raw)
    ax_spec = fig.add_subplot(gs[3, 0], sharex=ax_raw)
    ax_com = fig.add_subplot(gs[0:2, 1])
    ax_pol = fig.add_subplot(gs[2:4, 1], projection="polar")

    ax_raw.plot(t, lfp, lw=0.5, color="#1f77b4")
    ax_raw.set_ylabel("LFP (notch 60 Hz)", fontsize=9, color=cor_txt)
    ax_raw.set_title(f"PAC {par_nome} VALIDADO — {canal}, {arquivo.replace('.ns2', '')}\n"
                     f"{comportamento}\n"
                     f"pico {f_pico:g} Hz × {a_pico:g} Hz, z={z_mi:.1f}",
                     fontsize=11, color=cor_txt)

    ax_th.plot(t, theta, lw=0.9, color="#9467bd")
    ax_th.set_ylabel(f"Fase {f_pico:g}±1 Hz", fontsize=9, color=cor_txt)

    ax_ga.plot(t, gamma, lw=0.5, color="#ff7f0e", alpha=0.65)
    ax_ga.plot(t, env_gamma, lw=1.3, color="#7f2704")
    ax_ga.set_ylabel(f"Amp {a_pico:g}±5 Hz\n+ envelope", fontsize=8.5, color=cor_txt)

    f_st, t_st, sxx = signal.spectrogram(lfp, fs, nperseg=512, noverlap=384)
    ax_spec.pcolormesh(t_st, f_st, sxx, shading="auto", cmap="magma",
                       vmin=0, vmax=np.percentile(sxx, 99), rasterized=True)
    ax_spec.axhline(f_pico, color="white", lw=0.8, ls="--", alpha=0.8)
    ax_spec.axhspan(a_pico - 5, a_pico + 5, color="white", alpha=0.12)
    ax_spec.axhspan(4, 12, color="white", alpha=0.08)
    ax_spec.set_ylim(0, 150)
    ax_spec.set_ylabel("Espectrograma (Hz)", fontsize=9, color=cor_txt)
    ax_spec.set_xlabel("Tempo (s)", fontsize=10, color=cor_txt)

    plota_comodo_no_eixo(ax_com, z_mapa, fases_freq, amps_freq, f_pico, a_pico)
    ax_com.set_title("Comodulograma (z do MI)", fontsize=10, color=cor_txt)

    mi = plota_polar(ax_pol, fase_theta, env_gamma)
    ax_pol.set_title("Fase × Amplitude",
                     fontsize=10, color=cor_txt, pad=8)
    # MI/z ABAIXO do círculo polar (acima colidiria com o xlabel do comodulograma)
    ax_pol.text(0.5, -0.16, f"MI = {mi:.3f}  (z = {z_mi:.1f} vs surrogates)",
                transform=ax_pol.transAxes, ha="center", va="top",
                fontsize=10, color=cor_txt)

    for ax in (ax_raw, ax_th, ax_ga, ax_spec):
        ax.set_xlim(0, t[-1])
        ax.grid(alpha=0.15)
        ax.tick_params(labelsize=9, colors=cor_txt)
    plt.setp(ax_raw.get_xticklabels(), visible=False)
    plt.setp(ax_th.get_xticklabels(), visible=False)
    plt.setp(ax_ga.get_xticklabels(), visible=False)

    saida = os.path.join(saida_dir, f"{rotulo}_{canal}.png")
    os.makedirs(saida_dir, exist_ok=True)
    fig.savefig(saida, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  salva: {saida}  (MI={mi:.3f}, z={z_mi:.2f})")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pasta_ns2", required=True,
                    help="Pasta com os .ns2 da sessão")
    ap.add_argument("--vencedores", required=True,
                    help="CSV da sessão: rotulo,arquivo,canal,inicio_s,fim_s,"
                         "fase_pico_hz,amp_pico_hz[,comportamento]")
    ap.add_argument("--saida_dir", default="figuras")
    args = ap.parse_args()

    venc = pd.read_csv(args.vencedores)
    for _, r in venc.iterrows():
        comp = r.get("comportamento")
        par_nome = r.get("par", "theta_gamma")
        figura_vencedor(str(r["rotulo"]), str(r["arquivo"]), str(r["canal"]),
                        float(r["inicio_s"]), float(r["fim_s"]),
                        float(r["fase_pico_hz"]), float(r["amp_pico_hz"]),
                        "" if pd.isna(comp) else str(comp),
                        args.pasta_ns2, par_nome, args.saida_dir)


if __name__ == "__main__":
    main()
