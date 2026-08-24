"""
audita_grooming_robustez.py
==========================================
Robustez da janela de GROOMING puro do vencedor 2 (003.ns2 @ 20-27s, par
5x35 Hz), criada em 24/08/2026 quando o usuario ancorou no video:
15:14-15:21 = grooming -> ns2 20-27s; 15:22-15:29 = walking/sniffing/rearing.

O z=9,4 da janela original 20-30s era inflado pela escolha de janela (mistura
de estados + nula estreita); o nucleo robusto e z~4 durante o grooming.
Este script verifica se o nucleo sobrevive ao mesmo batizado dos outros
vencedores: sweep de n_bins + MVL (sem bins) no chan22, e consistencia nos
vizinhos chan24/26/28. Mesma nula de deslocamento circular; com notch 60 Hz.

Uso:
    python audita_grooming_robustez.py

Requer: neo, numpy, scipy
"""

import sys

import numpy as np

from comodulogram import aplica_notch
from robustez_parametros import mi_z_par, mvl_z_par
from ns2_utils import le_ns2, fatia_janela

PASTA_NS2 = "../Basal antes da infusao"
ARQUIVO = "20240708-123605-003.ns2"
INI, FIM = 20.0, 27.0          # grooming puro segundo o video
F_FASE, F_AMP = 5.0, 35.0      # par canônico do vencedor 2


def main():
    try:  # console Windows pode estar em cp1252; símbolos quebram o print
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

    dados, fs, nomes = le_ns2(f"{PASTA_NS2}/{ARQUIVO}")
    mapa = {str(n): i for i, n in enumerate(nomes)}

    # --- chan22: sweep de n_bins + MVL -------------------------------
    lfp = aplica_notch(
        fatia_janela(dados[:, mapa["chan22"]], fs, INI, FIM).astype(float),
        fs, 60.0)
    print(f"chan22 @ {INI:g}-{FIM:g}s (grooming), par {F_FASE:g}x{F_AMP:g} Hz:")
    for nb in (10, 12, 15, 18, 24, 30):
        z = mi_z_par(lfp, fs, F_FASE, F_AMP, n_bins=nb,
                     rng=np.random.default_rng(42))
        print(f"  n_bins={nb:>2}: z={z:5.2f}")
    z_mvl = mvl_z_par(lfp, fs, F_FASE, F_AMP, rng=np.random.default_rng(42))
    print(f"  MVL (sem bins): z={z_mvl:5.2f}")

    # --- vizinhos do cluster: mesma janela, mesmíssimo par -----------
    print("vizinhos (n_bins=18):")
    for ch in ("chan24", "chan26", "chan28"):
        lfp_c = aplica_notch(
            fatia_janela(dados[:, mapa[ch]], fs, INI, FIM).astype(float),
            fs, 60.0)
        z = mi_z_par(lfp_c, fs, F_FASE, F_AMP, rng=np.random.default_rng(42))
        print(f"  {ch}: z={z:5.2f}")


if __name__ == "__main__":
    main()
