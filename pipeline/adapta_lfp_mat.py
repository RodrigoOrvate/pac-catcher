"""adapta_lfp_mat.py
==================
Adapta arquivos .mat para o pipeline de PAC multi-banda.

Diagnóstico de escala: imprime min/max/std/mean do sinal carregado ANTES
de qualquer processamento. Isso é essencial para detectar problemas de
normalização que causam z-scores absurdos (overflow).

Shim de compatibilidade: a leitura de .mat foi movida para pac_core/io.py.
Este módulo continua existindo para preservar a CLI (`--mat/--saida_csv/
--normaliza/--verbose`) e o nome `carrega_lfp_mat` que outros scripts
importam. Novo código deve importar diretamente de `pac_core.io`.

Funções exportadas para uso em outros scripts:
  carrega_lfp_mat(path, chave)   — carrega uma variável do .mat
  normaliza_sinal(arr)           — z-score robusto (median/MAD)
  info_sinal(arr, nome)          — imprime diagnóstico de escala
  fs_do_mat(m)                   — detecta fs no .mat
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pac_core.io import (  # noqa: F401 - re-export para compatibilidade
    le_mat as carrega_lfp_mat,
    normaliza_sinal,
    info_sinal,
    fs_do_mat,
)

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


# ============================================================
# Modo CLI: diagnóstico + exportação para CSV
# ============================================================

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--mat", required=True, help="Arquivo .mat de entrada")
    ap.add_argument("--canal_mat", default=None,
                    help="Variável do .mat a usar (default: auto-detecta)")
    ap.add_argument("--saida_csv", default=None,
                    help="Caminho de saída CSV (default: tmp_<nome>.csv)")
    ap.add_argument("--normaliza", action="store_true",
                    help="Aplica normalização robusta (median/MAD) antes de exportar")
    ap.add_argument("--verbose", action="store_true",
                    help="Imprime diagnóstico de escala")
    args = ap.parse_args()

    if not os.path.isfile(args.mat):
        print(f"[ERRO] Arquivo não encontrado: {args.mat}")
        sys.exit(1)

    print(f"Carregando: {args.mat}")
    try:
        arr, fs, chave_usada = carrega_lfp_mat(args.mat, args.canal_mat)
    except KeyError as e:
        print(f"[ERRO] {e}")
        sys.exit(1)

    print(f"  Chave usada: '{chave_usada}'  |  fs detectada: {fs:.0f} Hz")
    info_sinal(arr, nome=chave_usada, fs=fs)

    if args.normaliza:
        print("  Aplicando normalização robusta (median/MAD)...")
        arr_out = normaliza_sinal(arr)
        print(f"  Após normalização: min={arr_out.min():.4g}  max={arr_out.max():.4g}  "
              f"std={arr_out.std():.4g}")
    else:
        arr_out = arr
        if arr.std() > 5000 or arr.std() < 0.001:
            print("  [AVISO] Considera usar --normaliza para corrigir escala antes do pipeline.")

    # Saída CSV
    nome_base = os.path.splitext(os.path.basename(args.mat))[0]
    saida = args.saida_csv or f"tmp_{nome_base.lower()}.csv"
    np.savetxt(saida, arr_out.reshape(-1, 1), delimiter=",", fmt="%.8f")
    print(f"\nExportado: {saida}")
    print(f"  {len(arr_out)} amostras, {len(arr_out)/fs:.1f} s @ {fs:.0f} Hz")
    print(f"  1 coluna (canal único)")
    if args.normaliza:
        print("  [normalizado: median=0, MAD=1]")

    # Documenta no README se existir
    readme_path = os.path.join(os.path.dirname(args.mat), "README.md")
    if os.path.isfile(readme_path):
        try:
            with open(readme_path, "a", encoding="utf-8") as f:
                f.write(f"\n\n## Exportação adapta_lfp_mat ({nome_base})\n")
                f.write(f"- Chave: `{chave_usada}` | fs: {fs:.0f} Hz | "
                        f"duração: {len(arr)/fs:.1f} s\n")
                f.write(f"- Normalizado: {'sim (median/MAD)' if args.normaliza else 'não'}\n")
                f.write(f"- Saída: `{saida}`\n")
        except Exception:
            pass


if __name__ == "__main__":
    main()
