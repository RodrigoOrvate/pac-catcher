"""
audita_segmentos.py (23/08/2026)
==========================================
Complemento do audita_transientes.py: localiza NO TEMPO onde mora o
acoplamento do caso disputado (chan22 @ 003 20-30s, pico 5x35 Hz).
Para cada sub-segmento, recomputa o mapa z completo (mesma matemática e
mesma nula do comodulogram.py) e reporta o z na célula do pico original
e o pico do quadrante Theta-Gamma do próprio segmento.

Uso:
    python audita_segmentos.py
"""
import sys

import numpy as np
import scipy.signal as signal

from ns2_utils import le_ns2, fatia_janela
from comodulogram import aplica_notch, calcula_comodulograma_z, z_pico_theta_gamma

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

PASTA = "../Basal antes da infusao"
ARQUIVO = "20240708-123605-003.ns2"
CANAL = "chan22"
SEGMENTOS = [(20, 26), (26, 30), (28, 30), (20, 30)]
FP, FA = 5, 35

fases_freq = np.arange(4, 15, 1)
amps_freq = np.arange(30, 155, 5)

dados, fs, nomes = le_ns2(f"{PASTA}/{ARQUIVO}")
idx = {str(n): i for i, n in enumerate(nomes)}[CANAL]

print(f"{CANAL} @ {ARQUIVO} — z na célula do pico original ({FP}x{FA} Hz) "
      f"e pico Theta-Gamma de cada segmento (mesma nula, 200 surrogates):\n")
for ini, fim in SEGMENTOS:
    lfp = fatia_janela(dados[:, idx], fs, ini, fim).astype(float)
    lfp_n = aplica_notch(lfp, fs, linha_hz=60.0)
    z_mapa, mi_obs, _ = calcula_comodulograma_z(
        lfp_n, fs, fases_freq, amps_freq, n_surr=200,
        rng=np.random.default_rng(42), notch_hz=None, retorna_mi=True)
    i_fp = int(np.argmin(np.abs(fases_freq - FP)))
    j_ap = int(np.argmin(np.abs(amps_freq - FA)))
    z_celula = float(z_mapa[j_ap, i_fp])
    z_pico, f_pico, a_pico = z_pico_theta_gamma(z_mapa, fases_freq, amps_freq)
    print(f"  {ini:2d}-{fim:2d}s ({fim - ini:.0f}s): z({FP}x{FA})={z_celula:6.2f} | "
          f"pico TG do segmento: z={z_pico:5.2f} @ {f_pico:g}x{a_pico:g} Hz | "
          f"MI na celula={mi_obs[j_ap, i_fp]:.4f}")
