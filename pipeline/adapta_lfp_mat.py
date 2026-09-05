"""adapta_lfp_mat.py
==================
Adapta arquivos .mat para o pipeline de PAC multi-banda.

Diagnóstico de escala: imprime min/max/std/mean do sinal carregado ANTES
de qualquer processamento. Isso é essencial para detectar problemas de
normalização que causam z-scores absurdos (overflow).

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
import scipy.io

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


# ============================================================
# Funções utilitárias exportáveis
# ============================================================

def fs_do_mat(m):
    """
    Detecta a taxa de amostragem em um dicionário do scipy.io.loadmat.
    Tenta variáveis comuns: fs, Fs, srate, samplingRate.
    Retorna 1000.0 se não encontrar.
    """
    for chave in ["fs", "Fs", "srate", "samplingRate", "SR"]:
        if chave in m:
            val = m[chave]
            try:
                return float(np.squeeze(val))
            except Exception:
                pass
    return 1000.0


def carrega_lfp_mat(path, chave=None):
    """
    Carrega um LFP de um arquivo .mat.

    Se 'chave' é None, tenta na ordem:
      lfp, LFP, lfpBruto, lfpHG, lfpHFO, signal, data

    Retorna (array_1d_float64, fs_float, chave_usada).
    """
    m = scipy.io.loadmat(path)
    chaves_candidatas = (
        [chave] if chave else
        ["lfp", "LFP", "lfpBruto", "lfpHG", "lfpHFO", "signal", "data"]
    )
    for c in chaves_candidatas:
        if c in m:
            arr = np.squeeze(m[c]).astype(np.float64)
            fs = fs_do_mat(m)
            return arr, fs, c
    raise KeyError(
        f"Nenhuma das chaves {chaves_candidatas} encontrada em {path}.\n"
        f"Chaves disponíveis: {[k for k in m if not k.startswith('_')]}"
    )


def info_sinal(arr, nome="sinal", fs=None):
    """
    Imprime diagnóstico de escala do sinal. Essencial para detectar
    problemas de normalização antes de passar pelo pipeline.
    """
    duracao = len(arr) / fs if fs else None
    print(f"  [{nome}]")
    print(f"    shape: {arr.shape}, dtype: {arr.dtype}")
    if duracao:
        print(f"    duração: {duracao:.1f} s @ {fs:.0f} Hz")
    print(f"    min={arr.min():.4g}  max={arr.max():.4g}")
    print(f"    mean={arr.mean():.4g}  std={arr.std():.4g}")
    med = np.median(arr)
    mad = np.median(np.abs(arr - med))
    print(f"    median={med:.4g}  MAD={mad:.4g}")
    # Alerta de escala
    if arr.std() > 5000:
        print(f"    *** AVISO: std={arr.std():.1f} muito alto — provável escala em nV ou "
              f"contagem ADC. Considere --normaliza ou converta para µV.")
    elif arr.std() < 0.001:
        print(f"    *** AVISO: std={arr.std():.6f} muito baixo — provável escala em V.")
    else:
        print(f"    Escala aparenta ser µV (std razoável para LFP).")


def normaliza_sinal(arr):
    """
    Normalização robusta: (x - median) / MAD.
    Mantém a forma do sinal, remove outliers de offset, escala para ~unidades.
    Preferida ao z-score mean/std quando há outliers (artefatos).
    """
    med = np.median(arr)
    mad = np.median(np.abs(arr - med))
    if mad < 1e-12:
        # Fallback para std se MAD for zero (sinal constante)
        std = arr.std()
        if std < 1e-12:
            return arr - med  # sinal constante
        return (arr - med) / std
    return (arr - med) / mad


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
