"""
audita_skewness.py — Teste de skewness do theta (harmônicos / integer ratio)

Verifica se a banda de fase (theta) de cada candidato é senoide limpa ou tem
forma assimétrica (theta não-senoidal). Quando a razão amp/fase é um inteiro
(ex.: 5×20 = 4.0, 5×35, 5×45), uma theta assimétrica pode gerar MI espúrio
(cada ciclo carrega um harmônico na banda de amplitude). Skewness |g|>limiar
= SUSPEITO; ~0 = senoide simétrica (ok).

Refactor de test_skewness.py (lógica preservada) → agora por CLI, lendo o
vencedores.csv da sessão (regra de ouro: nada de sessão entra no código).

    python audita_skewness.py --csv "<sessao>/RESULTADOS/vencedores.csv" \
        --pasta_ns2 "<sessao>/<BASAL>" --saida "<sessao>/RESULTADOS/skewness.csv"
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import signal

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

from ns2_utils import carrega_dados


def calculate_skewness(data):
    """Coeficiente padronizado de Fisher-Pearson (skewness amostral)."""
    mean = np.mean(data)
    std = np.std(data)
    if std == 0:
        return 0.0
    n = len(data)
    if n < 3:
        return np.nan
    return (n / ((n - 1) * (n - 2))) * np.sum(((data - mean) / std) ** 3)


def theta_skewness_for_window(file_path, channel_name, start_s, end_s, fs_banda=(4, 8), order=4):
    dados, fs, canal_ids = carrega_dados(file_path)
    if channel_name not in canal_ids:
        raise ValueError(f"{channel_name} nao encontrado em {canal_ids}")
    chan_idx = canal_ids.index(channel_name)

    idx_inicio = int(start_s * fs)
    idx_fim = int(end_s * fs)
    chan_data = dados[idx_inicio:idx_fim, chan_idx]
    n_amostras = len(chan_data)

    nyq = 0.5 * fs
    lo, hi = fs_banda
    b, a = signal.butter(order, [lo / nyq, hi / nyq], btype='band')
    theta_filtered = signal.filtfilt(b, a, chan_data)

    skew = calculate_skewness(theta_filtered)
    return skew, n_amostras


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True, help="vencedores.csv da sessão")
    ap.add_argument("--pasta_ns2", required=True, help="pasta com os .ns2")
    ap.add_argument("--saida", default="skewness.csv",
                    help="onde gravar o CSV de resultado")
    ap.add_argument("--limiar", type=float, default=0.5,
                    help="|skewness| acima disto = SUSPEITO (default 0.5)")
    args = ap.parse_args()

    df = pd.read_csv(args.csv)

    linhas = []
    print(f"{'Canal':<10} | {'Janela':<16} | {'N':<6} | {'Skew':<10} | Veredito")
    print("-" * 70)

    for _, r in df.iterrows():
        canal = str(r["canal"])
        arquivo = str(r["arquivo"])
        path = os.path.join(args.pasta_ns2, arquivo)
        ini = float(r.get("janela_ini_s", r.get("inicio_s", 0)))
        fim = float(r.get("janela_fim_s", r.get("fim_s", 0)))

        try:
            skew, n = theta_skewness_for_window(path, canal, ini, fim)
        except Exception as e:
            print(f"{canal:<10} | ERROR: {e}")
            linhas.append({"rotulo": r.get("rotulo"), "arquivo": arquivo, "canal": canal,
                           "janela_ini_s": ini, "janela_fim_s": fim,
                           "janela": f"{ini:.0f}-{fim:.0f}s", "janela_tipo": "full",
                           "n": np.nan, "skewness": np.nan, "veredito_skew": "ERROR",
                           "motivo": str(e)})
            continue

        veredito = "SUSPECT (Asymmetric)" if abs(skew) > args.limiar else "CLEAN (Symmetric)"
        if n < 500:
            veredito += " [n baixo, checar manualmente]"
        print(f"{canal:<10} | {ini:.0f}-{fim:.0f}s (full) | {n:<6} | {skew:<10.4f} | {veredito}")
        linhas.append({"rotulo": r.get("rotulo"), "arquivo": arquivo, "canal": canal,
                       "janela_ini_s": ini, "janela_fim_s": fim,
                       "janela": f"{ini:.0f}-{fim:.0f}s", "janela_tipo": "full",
                       "n": n, "skewness": round(skew, 4), "veredito_skew": veredito,
                       "motivo": ""})

        # ilha fina, se presente no CSV (colunas inicio_ilha_s/fim_ilha_s)
        if "inicio_ilha_s" in df.columns and "fim_ilha_s" in df.columns \
                and not pd.isna(r.get("inicio_ilha_s")):
            ii, fi = float(r["inicio_ilha_s"]), float(r["fim_ilha_s"])
            try:
                skew_i, n_i = theta_skewness_for_window(
                    path, canal, ii, fi, order=2)  # ordem menor p/ poucas amostras
            except Exception:
                continue
            v_i = "SUSPECT (Asymmetric)" if abs(skew_i) > args.limiar else "CLEAN (Symmetric)"
            if n_i < 500:
                v_i += " [n baixo]"
            print(f"{canal:<10} | {ii:.0f}-{fi:.0f}s (ilha) | {n_i:<6} | {skew_i:<10.4f} | {v_i}")
            linhas.append({"rotulo": r.get("rotulo"), "arquivo": arquivo, "canal": canal,
                           "janela_ini_s": ii, "janela_fim_s": fi,
                           "janela": f"{ii:.0f}-{fi:.0f}s", "janela_tipo": "island",
                           "n": n_i, "skewness": round(skew_i, 4),
                           "veredito_skew": v_i, "motivo": ""})

    out = pd.DataFrame(linhas)
    out.to_csv(args.saida, index=False)
    print(f"\nSalvo: {args.saida} ({len(out)} linhas)")


if __name__ == "__main__":
    main()
