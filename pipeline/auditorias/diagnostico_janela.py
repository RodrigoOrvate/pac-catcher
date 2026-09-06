"""
diagnostico_janela.py
==========================================
Figura diagnóstica de uma janela: separa acoplamento de verdade de
artefato de transiente/respiração.

POR QUE EXISTE: o comodulograma z-scoredo pode mostrar uma COLUNA de fase
(uma freq de fase acoplando a muitas amplitudes) em vez de um cluster
focal theta-gamma. Isso acontece quando o sinal de "fase" não é senoide
limpa, mas tem transientes afiados (potenciais respiratórios, ondas
agudas, spikes): cada transiente é um burst de banda larga, e se os
transientes ocorrem numa fase preferencial, TODAS as bandas de amplitude
acendem juntas (Aru et al. 2015; Cole & Voytek 2017). O z vs surrogates
confirma que a modulação existe, mas não distingue "oscilação aninhada"
de "transientes ritmados".

O QUE FAZ, para um arquivo/canal/janela:
  Plota 5 painéis alinhados no tempo:
    1. traço bruto
    2. banda respiratória (0.5-3 Hz) -- deflexões lentas grandes =
       artefato/potencial respiratório; rato em repouso respira a ~1-2 Hz,
       sniffing a ~6-10 Hz (que CAI DENTRO do theta -- daí o perigo)
    3. banda theta (4-12 Hz)
    4. banda gamma (30-80 Hz) + envelope
    5. espectrograma (STFT)
  E calcula, com a mesma nula de surrogates do pipeline:
    - MI theta(4-8) x gamma(30-80)   -- o acoplamento original
    - MI resp(0.5-3) x gamma(30-80)  -- se também alto, respiração explica
    - MI resp(0.5-3) x theta(4-8)    -- respiração modula o "theta"?
    - fator de crista do theta filtrado (senoide ~1.4; transientes >>2)

Uso:
    python diagnostico_janela.py --arquivo "../Basal antes da infusao/20240708-123605-002.ns2" \
        --canal chan32 --inicio 95 --fim 105 --notch 60

Requer: neo, numpy, scipy, matplotlib
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt

from ns2_utils import le_ns2, fatia_janela


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


def mi_z(fase, envelope, fs, n_surr=200, n_bins=18, rng=None):
    """MI observado + z contra surrogates de deslocamento circular
    (mesma nula do triagem_pac.py)."""
    if rng is None:
        rng = np.random.default_rng()
    bins = np.linspace(-np.pi, np.pi, n_bins + 1)
    bin_idx = np.clip(np.digitize(fase, bins) - 1, 0, n_bins - 1)
    mi_obs = _mi_de_bin_idx(bin_idx, envelope, n_bins)

    n = len(envelope)
    shift_min = int(1.0 * fs)
    if n <= 2 * shift_min:
        shift_min = max(1, n // 10)
    deslocamentos = rng.integers(shift_min, n - shift_min, size=n_surr)
    mi_surr = np.array([_mi_de_bin_idx(bin_idx, np.roll(envelope, d), n_bins)
                        for d in deslocamentos])
    dp = np.std(mi_surr)
    z = (mi_obs - np.mean(mi_surr)) / dp if dp > 0 else 0.0
    return mi_obs, z


def aplica_notch(sinal_in, fs, freqs_notch, q_factor=30.0):
    if not freqs_notch:
        return sinal_in
    nyq = 0.5 * fs
    out = sinal_in
    if not isinstance(freqs_notch, (list, tuple, np.ndarray)):
        freqs_notch = [freqs_notch]
    for f in freqs_notch:
        if f >= nyq * 0.98:
            continue
        b, a = signal.iirnotch(f / nyq, q_factor)
        out = signal.filtfilt(b, a, out)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arquivo", required=True, help="Caminho do .ns2")
    ap.add_argument("--canal", required=True, help="Nome do canal (ex.: chan32)")
    ap.add_argument("--inicio", type=float, required=True, help="Início da janela (s)")
    ap.add_argument("--fim", type=float, required=True, help="Fim da janela (s)")
    ap.add_argument("--notch", type=float, nargs="+", default=None, metavar="HZ",
                    help="Frequência(s) da rede a notchar (ex.: 60 120 180)")
    ap.add_argument("--n_surr", type=int, default=200)
    ap.add_argument("--saida_png", default=None)
    args = ap.parse_args()

    rng = np.random.default_rng(42)

    print(f"Carregando {args.arquivo} ...")
    dados, fs, nomes = le_ns2(args.arquivo)
    mapa = {str(n): i for i, n in enumerate(nomes)}
    if args.canal not in mapa:
        raise SystemExit(f"Canal '{args.canal}' não encontrado. "
                         f"Disponíveis: chan1..chan{dados.shape[1]}")
    idx = mapa[args.canal]

    lfp = fatia_janela(dados[:, idx], fs, args.inicio, args.fim).astype(float)
    if args.notch:
        lfp = aplica_notch(lfp, fs, freqs_notch=args.notch)
    t = np.arange(lfp.size) / fs

    # ==========================================
    # BANDAS DE INTERESSE
    # ==========================================
    resp = filtra_sinal(lfp, 0.5, 3.0, fs)          # respiração (repouso 1-2 Hz)
    theta = filtra_sinal(lfp, 4.0, 12.0, fs)        # theta largo, p/ visual
    theta_mi = filtra_sinal(lfp, 4.0, 8.0, fs)      # theta do pipeline
    gamma = filtra_sinal(lfp, 30.0, 80.0, fs)
    env_gamma = np.abs(signal.hilbert(gamma))
    env_theta = np.abs(signal.hilbert(theta_mi))

    fase_theta = np.angle(signal.hilbert(theta_mi))
    fase_resp = np.angle(signal.hilbert(resp))

    # ==========================================
    # NÚMEROS: qual "fase" organiza a amplitude de gamma?
    # ==========================================
    mi_tg, z_tg = mi_z(fase_theta, env_gamma, fs, args.n_surr, rng=rng)
    mi_rg, z_rg = mi_z(fase_resp, env_gamma, fs, args.n_surr, rng=rng)
    mi_rt, z_rt = mi_z(fase_resp, env_theta, fs, args.n_surr, rng=rng)

    # Fator de crista do theta filtrado: senoide pura ~1.4; com transientes >>2
    crista_theta = float(np.max(np.abs(theta_mi)) / np.sqrt(np.mean(theta_mi ** 2)))

    print(f"\n=== {args.canal} @ {args.inicio:g}-{args.fim:g}s "
          f"({os.path.basename(args.arquivo)}{' | notch ' + str(args.notch) if args.notch else ''}) ===")
    print(f"MI theta(4-8) x gamma(30-80):  MI={mi_tg:.4f}  z={z_tg:.2f}   <- acoplamento original")
    print(f"MI resp(0.5-3) x gamma(30-80): MI={mi_rg:.4f}  z={z_rg:.2f}   <- respiração explica gamma?")
    print(f"MI resp(0.5-3) x theta(4-8):   MI={mi_rt:.4f}  z={z_rt:.2f}   <- respiração modula o 'theta'?")
    print(f"Fator de crista do theta filtrado: {crista_theta:.2f}  "
          f"(senoide ~1.4; >2 indica transientes/ondas agudas)")
    print(f"Amplitude pico-a-pico da banda respiratória: "
          f"{np.ptp(resp):.0f} unidades (comparar com o traço bruto)")

    # ==========================================
    # FIGURA (5 painéis alinhados no tempo)
    # ==========================================
    fig, eixos = plt.subplots(5, 1, figsize=(12, 11), sharex=True,
                              gridspec_kw={"height_ratios": [2, 1.2, 1.2, 1.2, 1.6]})
    ax1, ax2, ax3, ax4, ax5 = eixos
    cor_txt = "#333333"

    ax1.plot(t, lfp, lw=0.6, color="#1f77b4")
    ax1.set_ylabel("Bruto", fontsize=10, color=cor_txt)
    ax1.set_title(f"{args.canal} @ {args.inicio:g}–{args.fim:g}s — {os.path.basename(args.arquivo)}"
                  f"{' | notch ' + str(args.notch) + ' Hz' if args.notch else ''}", fontsize=12)

    ax2.plot(t, resp, lw=0.8, color="#2ca02c")
    ax2.set_ylabel("Resp 0.5–3 Hz", fontsize=10, color=cor_txt)

    ax3.plot(t, theta, lw=0.8, color="#9467bd")
    ax3.set_ylabel("Theta 4–12 Hz", fontsize=10, color=cor_txt)

    ax4.plot(t, gamma, lw=0.6, color="#ff7f0e", alpha=0.7)
    ax4.plot(t, env_gamma, lw=1.2, color="#7f2704")
    ax4.set_ylabel("Gamma 30–80 Hz\n+ envelope", fontsize=9, color=cor_txt)

    f_st, t_st, sxx = signal.spectrogram(lfp, fs, nperseg=512, noverlap=384)
    ax5.pcolormesh(t_st, f_st, sxx, shading="auto", cmap="magma",
                   vmin=0, vmax=np.percentile(sxx, 99))
    ax5.axhline(8, color="white", lw=0.7, ls="--", alpha=0.7)
    ax5.axhline(30, color="white", lw=0.7, ls="--", alpha=0.7)
    ax5.axhline(80, color="white", lw=0.7, ls="--", alpha=0.7)
    ax5.set_ylabel("Espectrograma\n(Hz)", fontsize=10, color=cor_txt)
    ax5.set_xlabel("Tempo (s)", fontsize=11, color=cor_txt)
    ax5.set_ylim(0, 150)

    for ax in eixos:
        ax.grid(alpha=0.15)
        ax.tick_params(labelsize=9, colors=cor_txt)

    resumo = (f"MI θ×γ z={z_tg:.1f} | MI resp×γ z={z_rg:.1f} | "
              f"MI resp×θ z={z_rt:.1f} | crista θ={crista_theta:.2f}")
    fig.suptitle(resumo, fontsize=11, y=0.995, color=cor_txt)
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    caminho = args.saida_png or (f"diagnostico_{args.canal}_{args.inicio:g}-{args.fim:g}s.png")
    fig.savefig(caminho, dpi=150)
    plt.close(fig)
    print(f"\nFigura salva: {caminho}")


if __name__ == "__main__":
    main()
