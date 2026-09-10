"""
audita_segmentos.py (23/08/2026, CLI em 24/08/2026)
==========================================
Complemento do audita_transientes.py: localiza NO TEMPO onde mora o
acoplamento do caso disputado. Para cada sub-segmento, recomputa o mapa z
completo (mesma matemática e mesma nula do comodulogram.py) e reporta o z
na célula do pico original e o pico do quadrante Theta-Gamma do próprio
segmento.

Uso (o caso da sessão 08/07 é o default — nada precisa ser passado):

    python audita_segmentos.py

Nova sessão, por CLI (valores específicos NUNCA entram no código):

    python audita_segmentos.py \
        --pasta "<SESSAO>/<BASAL>" \
        --arquivo 20240711-121046-001.ns2 \
        --canal chan26 \
        --segmentos 65-70,70-75,65-75 \
        --fp 5 --fa 70
"""
import sys
import argparse
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from ns2_utils import le_ns2, fatia_janela
from comodulogram import aplica_notch, calcula_comodulograma_z, z_pico_par

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pasta", default="../Basal antes da infusao",
                    help="Pasta com os .ns2 (padrão: caso 08/07)")
    ap.add_argument("--arquivo", default="20240708-123605-003.ns2",
                    help=".ns2 dentro de --pasta")
    ap.add_argument("--canal", default="chan22")
    ap.add_argument("--segmentos", default="20-26,26-30,28-30,20-30",
                    help="Lista ini-fim separada por vírgula, ex.: '65-70,70-75,65-75'")
    ap.add_argument("--fp", type=int, default=5, help="Fase de pico original (Hz)")
    ap.add_argument("--fa", type=int, default=35, help="Amplitude de pico original (Hz)")
    ap.add_argument("--n_surr", type=int, default=200)
    args = ap.parse_args()

    segmentos = []
    for parte in args.segmentos.split(","):
        ini, fim = parte.strip().split("-")
        segmentos.append((float(ini), float(fim)))

    fases_freq = np.arange(4, 15, 1)
    amps_freq = np.arange(30, 155, 5)

    dados, fs, nomes = le_ns2(f"{args.pasta}/{args.arquivo}")
    idx = {str(n): i for i, n in enumerate(nomes)}[args.canal]

    print(f"{args.canal} @ {args.arquivo} — z na célula do pico original "
          f"({args.fp}x{args.fa} Hz) e pico Theta-Gamma de cada segmento "
          f"(mesma nula, {args.n_surr} surrogates):\n")
    for ini, fim in segmentos:
        lfp = fatia_janela(dados[:, idx], fs, ini, fim).astype(float)
        lfp_n = aplica_notch(lfp, fs, linha_hz=60.0)
        z_mapa, mi_obs, _ = calcula_comodulograma_z(
            lfp_n, fs, fases_freq, amps_freq, n_surr=args.n_surr,
            rng=np.random.default_rng(42), notch_hz=None, retorna_mi=True)
        i_fp = int(np.argmin(np.abs(fases_freq - args.fp)))
        j_ap = int(np.argmin(np.abs(amps_freq - args.fa)))
        z_celula = float(z_mapa[j_ap, i_fp])
        z_pico, f_pico, a_pico = z_pico_par(z_mapa, fases_freq, amps_freq)
        print(f"  {ini:5g}-{fim:<5g} ({fim - ini:.0f}s): z({args.fp}x{args.fa})={z_celula:6.2f} | "
              f"pico TG do segmento: z={z_pico:5.2f} @ {f_pico:g}x{a_pico:g} Hz | "
              f"MI na celula={mi_obs[j_ap, i_fp]:.4f}")


if __name__ == "__main__":
    main()
