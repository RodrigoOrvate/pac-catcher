"""Fonte unica de verdade para o caminho-base do workspace fora do repo git
(a pasta acoplamento_theta-gamma, que contem SCRIPT/, LAC_NOCI/, etc. como
irmas -- ver .claude/CLAUDE.md). Nunca faca hardcode de C:\\acoplamento_theta-gamma
ou D:\\acoplamento_theta-gamma num script novo -- importe daqui.

Sobrescreva com a variavel de ambiente ACOPLAMENTO_BASE se o workspace
estiver montado em outro lugar (ex.: apos mover a pasta de disco, como em
13/09/2026 quando foi copiada de C: para D: por falta de espaco)."""
import os

BASE_WORKSPACE = os.environ.get("ACOPLAMENTO_BASE", r"D:\acoplamento_theta-gamma")
BASE_LAC_NOCI = os.path.join(BASE_WORKSPACE, "LAC_NOCI")
BASE_RESULTADOS_MESTRADO = os.path.join(BASE_WORKSPACE, "RESULTADOS_MESTRADO")
BASE_SCRIPT = os.path.join(BASE_WORKSPACE, "SCRIPT")
