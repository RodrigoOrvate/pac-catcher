"""
extrair_picos.py — Re-verificação do pico de cada vencedor com banda restrita + FDR
==========================================
Integração da extração de picos no pipeline, CORRIGIDA contra os dois vícios
do atualizar_picos.py antigo:
  1. Argmax cru infla o z por múltiplas comparações (mesmo vício do
     double-dipping temporal) -> agora o pico só é aceito se a célula
     vencedora for significativa sob FDR de Benjamini-Hochberg sobre TODO o
     mapa (bh_fdr_mapa, mesma nula do comodulogram.py).
  2. Busca ampla (4-14 x 20-150) deixa o pico escapar da banda clássica
     Theta-Gamma -> por padrão procura só em 4-8 x 30-80 Hz; a busca ampla
     é opt-in via --busca_ampla.

Pra cada linha do vencedores.csv da sessão, grava a comparação lado-a-lado
do pico antigo (colunas fase_pico_hz/amp_pico_hz) versus o pico novo dentro
da banda (restrita + FDR), em extracao_picos_v1_vs_v2.csv. Serve para decidir
se o pico reportado se mantém (confirmado), se migra para dentro da banda,
ou se não sobrevive ao FDR.

    python extrair_picos.py --csv "<sessao>/RESULTADOS/vencedores.csv" \
        --pasta_ns2 "<sessao>/<BASAL>" \
        --saida "<sessao>/RESULTADOS/extracao_picos_v1_vs_v2.csv"
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

from ns2_utils import le_ns2, fatia_janela
from triagem_pac import BAND_PAIRS
from comodulogram import (aplica_notch, calcula_comodulograma_z,
                          p_valores_por_celula, bh_fdr_mapa,
                          resume_cluster_fdr, z_pico_par,
                          FASES_DEFAULT, AMPS_DEFAULT)

FASES_PADRAO = FASES_DEFAULT
AMPS_PADRAO = AMPS_DEFAULT


def _pico_in_banda(z_mapa, fases_freq, amps_freq, theta_band, gamma_band):
    mask_fase = (fases_freq >= theta_band[0]) & (fases_freq <= theta_band[1])
    mask_amp = (amps_freq >= gamma_band[0]) & (amps_freq <= gamma_band[1])
    sub = z_mapa[np.ix_(mask_amp, mask_fase)]
    ii, jj = np.unravel_index(np.argmax(sub), sub.shape)
    return (float(sub[ii, jj]),
            float(fases_freq[mask_fase][jj]),
            float(amps_freq[mask_amp][ii]))


def _z_na_celula(z_mapa, fases_freq, amps_freq, fp, fa):
    i = int(np.argmin(np.abs(amps_freq - fa)))
    j = int(np.argmin(np.abs(fases_freq - fp)))
    return float(z_mapa[i, j])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True, help="vencedores.csv da sessão")
    ap.add_argument("--pasta_ns2", required=True, help="pasta com os .ns2")
    ap.add_argument("--saida", default="extracao_picos_v1_vs_v2.csv")
    ap.add_argument("--busca_ampla", action="store_true",
                    help="expande a busca original para 4-14 x 20-250 Hz (legacy)")
    ap.add_argument("--n_surr", type=int, default=200,
                    help="surrogates por célula (default 200)")
    ap.add_argument("--fdr_q", type=float, default=0.05, help="q do FDR (default 0.05)")
    ap.add_argument("--par", default="theta_gamma",
                    choices=list(BAND_PAIRS.keys()),
                    help="Qual par usar para a banda padrão (default: theta_gamma)")
    args = ap.parse_args()

    cfg = BAND_PAIRS.get(args.par, BAND_PAIRS["theta_gamma"])
    if args.busca_ampla:
        theta_band = (4.0, 14.0)
        gamma_band = (20.0, 250.0)
    else:
        theta_band = cfg["fase"]
        gamma_band = cfg["amp"]

    df = pd.read_csv(args.csv)

    linhas = []
    print(f"{'Canal':<8} | {'pico antigo':<14} | {'pico novo (in-banda)':<18} | "
          f"FDR | veredito")
    print("-" * 78)

    for _, r in df.iterrows():
        canal = str(r["canal"])
        arquivo = str(r["arquivo"])
        path = os.path.join(args.pasta_ns2, arquivo)
        ini, fim = float(r["inicio_s"]), float(r["fim_s"])
        fp_old, fa_old = float(r["fase_pico_hz"]), float(r["amp_pico_hz"])

        dados, fs, nomes = le_ns2(path)
        idx = {str(n): i for i, n in enumerate(nomes)}[canal]
        lfp = fatia_janela(dados[:, idx], fs, ini, fim).astype(float)
        lfp_n = aplica_notch(lfp, fs, linha_hz=60.0)

        z_mapa, mi_obs, mi_surr = calcula_comodulograma_z(
            lfp_n, fs, FASES_PADRAO, AMPS_PADRAO, n_surr=args.n_surr,
            rng=np.random.default_rng(42), notch_hz=None, retorna_mi=True)
        p_mapa = p_valores_por_celula(mi_obs, mi_surr)
        mask_fdr = bh_fdr_mapa(p_mapa, alpha=args.fdr_q)
        frac = resume_cluster_fdr(mask_fdr, FASES_PADRAO, AMPS_PADRAO, par=args.par)

        z_old = _z_na_celula(z_mapa, FASES_PADRAO, AMPS_PADRAO, fp_old, fa_old)
        z_new, fp_new, fa_new = _pico_in_banda(z_mapa, FASES_PADRAO, AMPS_PADRAO,
                                               theta_band, gamma_band)
        j_new = int(np.argmin(np.abs(FASES_PADRAO - fp_new)))
        i_new = int(np.argmin(np.abs(AMPS_PADRAO - fa_new)))
        fdr_sig_novo = bool(mask_fdr[i_new, j_new])
        n_sig = frac["n_sig"]
        frac_par = frac.get("frac_sig_par", frac.get("frac_sig_tg", 0.0))

        in_band_old = (theta_band[0] <= fp_old <= theta_band[1] and
                       gamma_band[0] <= fa_old <= gamma_band[1])
        if not in_band_old:
            verdict = "pico_antigo_fora_da_banda"
        elif fdr_sig_novo:
            verdict = "confirmado_in_banda"
        elif n_sig == 0:
            verdict = "nao_sig_fdr"
        else:
            verdict = "pico_migrou_ou_nao_sig"

        print(f"{canal:<8} | {fp_old:.0f}x{fa_old:.0f} (z={z_old:.2f}) | "
              f"{fp_new:.0f}x{fa_new:.0f} (z={z_new:.2f}) | "
              f"{'sim' if fdr_sig_novo else 'nao'} | {verdict}")

        linhas.append({
            "rotulo": r.get("rotulo"), "canal": canal, "arquivo": arquivo,
            "inicio_s": ini, "fim_s": fim, "par": args.par,
            "fase_antiga_hz": fp_old, "amp_antiga_hz": fa_old, "z_antiga": round(z_old, 2),
            "fase_nova_hz": fp_new, "amp_nova_hz": fa_new, "z_nova": round(z_new, 2),
            "fdr_sig_nova": fdr_sig_novo, "n_sig_fdr": n_sig,
            "frac_sig_par": round(frac_par, 3), "veredito": verdict,
        })

    out = pd.DataFrame(linhas)
    out.to_csv(args.saida, index=False)
    print(f"\nSalvo: {args.saida} ({len(out)} linhas)")


if __name__ == "__main__":
    main()
