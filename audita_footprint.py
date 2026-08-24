"""
audita_footprint.py (23/08/2026)
==========================================
Complemento do audita_transientes.py: mede a PEGADA ESPACIAL do acoplamento
-- o z na célula do pico em TODOS os 32 canais da mesma janela.

Lógica do discriminador: fonte cortical local produz gradiente suave
(poucos canais vizinhos significantes); artefato difuso (respiração/
sniffing, movimento, volume conduzido, cabo) aparece simultaneamente em
muitos canais distantes. A janela rearing (001 @ 85-95, quarteto
chan18/20/30/32 validado) serve de referência do que é pegada "local".

Uso:
    python audita_footprint.py
"""
import sys

import numpy as np

from ns2_utils import le_ns2, fatia_janela
from comodulogram import aplica_notch, filtra_sinal
from audita_transientes import mi_z_par

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

PASTA = "../Basal antes da infusao"

# (arquivo, janela, par de pico, rótulo)
JANELAS = [
    ("20240708-123605-003.ns2", (20, 30), (5, 35),
     "DISPUTADO: ativo 15:14-15:29"),
    ("20240708-123605-001.ns2", (85, 95), (7, 75),
     "REFERÊNCIA: rearing 6:46-6:56 (validado)"),
    ("20240708-123605-001.ns2", (30, 40), (8, 60),
     "REFERÊNCIA: grooming 5:51-6:06 (validado)"),
]

cache = {}
for arquivo, (ini, fim), (fp, fa), rotulo in JANELAS:
    if arquivo not in cache:
        dados, fs, nomes = le_ns2(f"{PASTA}/{arquivo}")
        cache[arquivo] = (dados, fs, [str(n) for n in nomes])
    dados, fs, nomes = cache[arquivo]

    zs = []
    for i, canal in enumerate(nomes):
        lfp = fatia_janela(dados[:, i], fs, ini, fim).astype(float)
        lfp_n = aplica_notch(lfp, fs, linha_hz=60.0)
        z, _, _, _ = mi_z_par(lfp_n, fs, fp, fa)
        zs.append(z)
    zs = np.array(zs)
    ordem = np.argsort(zs)[::-1]

    print(f"\n=== {rotulo} — {arquivo} @ {ini}-{fim}s — célula {fp}x{fa} Hz ===")
    print("  z por canal (ordem de gravação):")
    for linha in range(0, 32, 8):
        print("   " + "  ".join(f"{nomes[i]:>6s}:{zs[i]:5.1f}"
                                for i in range(linha, linha + 8)))
    top = [f"{nomes[i]}({zs[i]:.1f})" for i in ordem[:6]]
    n3 = int(np.sum(zs >= 3.0))
    print(f"  top-6: {', '.join(top)}")
    print(f"  canais com z>=3: {n3}/32")
