"""
audita_footprint.py (23/08/2026, CLI em 30/08/2026)
==========================================
Complemento do audita_transientes.py: mede a PEGADA ESPACIAL do acoplamento
-- o z na célula do pico em TODOS os 32 canais da mesma janela.

Lógica do discriminador: fonte cortical local produz gradiente suave
(poucos canais vizinhos significantes); artefato difuso (respiração/
sniffing, movimento, volume conduzido, cabo) aparece simultaneamente em
muitos canais distantes. É a ferramenta certa para decidir se um "cluster"
de canais de mesmo pico é UM evento co-detectado (pegada similiar, mesma
fonte) ou fontes distintas.

Refactor de audita_footprint.py (lógica preservada) → por CLI (regra de ouro:
nada de sessão no código). Formato de --casos:
    "rotulo,arquivo,ini,fim,fp,fa" separados por ';'.

    python audita_footprint.py --pasta_ns2 "<SESSAO>/<BASAL>" \
        --casos "ev1,2024...-001.ns2,135,145,8,50;ev2,2024...-002.ns2,285,295,5,20" \
        [--saida footprint.csv]
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

from pac_core.io import le_ns2, fatia_janela
from pipeline.auditorias.audita_transientes import mi_z_par
from pipeline.etapa3_comodulograma.comodulogram import aplica_notch


def _carrega(path):
    dados, fs, nomes = le_ns2(path)
    return dados, fs, [str(n) for n in nomes]


def footprint_de_janela(dados, fs, nomes, ini, fim, fp, fa, notch_hz=60.0):
    """
    Núcleo: z do MI na célula (fp,fa) para TODOS os canais em `nomes`,
    na mesma janela [ini,fim). Recebe os dados do arquivo já carregados
    em memória (dados, fs, nomes) -- não relê arquivo.

    Devolve dict: {"nomes", "zs" (ndarray), "top6" (lista de (nome,z)),
    "n_z_ge_3" (int)}.
    """
    zs = []
    for i, canal in enumerate(nomes):
        lfp = fatia_janela(dados[:, i], fs, ini, fim).astype(float)
        lfp_n = aplica_notch(lfp, fs, linha_hz=notch_hz)
        z, _, _, _ = mi_z_par(lfp_n, fs, fp, fa)
        zs.append(z)
    zs = np.array(zs)
    ordem = np.argsort(zs)[::-1]
    top6 = [(nomes[i], zs[i]) for i in ordem[:6]]
    n_z_ge_3 = int(np.sum(zs >= 3.0))

    return {"nomes": nomes, "zs": zs, "top6": top6, "n_z_ge_3": n_z_ge_3}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pasta_ns2", required=True)
    ap.add_argument("--csv", required=True,
                    help="'rotulo,arquivo,ini,fim,fp,fa;...'")
    ap.add_argument("--saida", default=None,
                    help="opcional: grava footprint.csv (rotulo, canal, z)")
    args = ap.parse_args()

    casos = []
    for c in args.csv.split(";"):
        if not c.strip():
            continue
        rotulo, arquivo, ini, fim, fp, fa = [x.strip() for x in c.split(",")]
        casos.append((rotulo, arquivo, int(ini), int(fim), float(fp), float(fa)))

    cache = {}
    linhas = []
    for rotulo, arquivo, ini, fim, fp, fa in casos:
        path = os.path.join(args.pasta_ns2, arquivo)
        if path not in cache:
            cache[path] = _carrega(path)
        dados, fs, nomes = cache[path]

        res = footprint_de_janela(dados, fs, nomes, ini, fim, fp, fa, notch_hz=60.0)
        zs = res["zs"]

        n_canais = len(nomes)
        print(f"\n=== {rotulo} — {arquivo} @ {ini}-{fim}s — célula {fp:.0f}x{fa:.0f} Hz "
              f"— {n_canais} canais ===")
        print("  z por canal (ordem de gravação):")
        for linha in range(0, n_canais, 8):
            print("   " + "  ".join(f"{nomes[i]:>6s}:{zs[i]:5.1f}"
                                    for i in range(linha, min(linha + 8, n_canais))))
        top = [f"{nome}({z:.1f})" for nome, z in res["top6"]]
        n3 = res["n_z_ge_3"]
        print(f"  top-6: {', '.join(top)}")
        print(f"  canais com z>=3: {n3}/32")

        for i, canal in enumerate(nomes):
            linhas.append({"rotulo": rotulo, "canal": canal, "z": round(zs[i], 2)})

    if args.saida:
        import pandas as pd
        pd.DataFrame(linhas).to_csv(args.saida, index=False)
        print(f"\nSalvo: {args.saida}")


if __name__ == "__main__":
    main()
