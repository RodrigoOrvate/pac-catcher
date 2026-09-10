"""
audita_harmonico.py - Teste de razao harmonica teta-gama via FOOOF e Rigidez de Fase

Complementa audita_skewness.py: em vez de julgar teta apenas pela forma de onda
(skewness), usa o FOOOF (modificado, Kuhn et al. 2026 - LFP_FOOOF) para extrair
um cf de teta refinado numa janela de CONTEXTO longa, e testa se o par
(fase_pico_hz, amp_pico_hz) do vencedor e consistente com um harmonico
inteiro desse teta, tanto em FREQUENCIA quanto em FASE (PLV).

A rigidez de fase (PLV alto entre phi_gamma e n*phi_theta) e o discriminador
mais forte entre harmônicos matemáticos e acoplamentos genuínos.

    python audita_harmonico.py --csv "<sessao>/RESULTADOS/vencedores.csv" \
        --pasta_ns2 "<sessao>/<BASAL>" \
        --saida "<sessao>/RESULTADOS/harmonico.csv" \
        --janela_contexto_s 45

==========================================================================
ATENCAO - Codigo experimental em calibracao:
Esta e uma extensao original (Razao Freq + PLV + Qualidade FOOOF) NAO
publicada nem validada formalmente.

RISCOS E LIMITACOES (Ordem de Prioridade):
1. VALIDACAO SINTETICA: O script nao foi testado contra sinais sinteticos
   (Teta assimetrico com harmonico vs Teta+Gamma independentes). a taxa de
   falsos positivos/negativos e desconhecida.
2. MODELO FOOOF: O 'aperiodic_mode=knee' e uma aproximacao. Para canais de
   DG, modelos '2exp' ou '3exp' (LFP_FOOOF original) sao necessarios para
   reduzir erro de estimativa de cf (< 0.1%).
3. ESTACIONARIEDADE: Em imobilidade, o teta pode ocorrer em rajadas. Janelas
   de 45s podem misturar estados, distorcendo o cf.
4. CALIBRACAO: Limiares de erro_ajuste (0.15), skewness (0.5) e PLV (0.8)
   sao chutes iniciais.

Trate saidas como rotulo/auditoria, NAO como filtro.
==========================================================================
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.signal import welch, butter, filtfilt, hilbert

# Dependency: pip install fooof (or specparam)
try:
    from fooof import FOOOF
except ImportError:
    sys.stderr.write(
        "[ERRO] Pacote 'fooof' nao encontrado. Instale com:\n"
        "    pip install fooof\n"
    )
    sys.exit(1)

# Path setup for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # pipeline/auditorias/
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # pipeline/

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

from ns2_utils import carrega_dados, fatia_janela
from audita_skewness import skewness_de_sinal  # REUSO, nao duplicacao
from linha_noise_kuhn import aplica_modo  # REUSO: limpeza de linha Kuhn (60Hz)


from utils_harmonico import (
    extrai_cf_teta_fooof,
    compute_plv_harmonico,
    testa_razao_harmonica,
    calcula_n_max
)


def avalia_harmonico(sinal_ctx, sinal_cand, fs, amp_pico, skew=None,
                     tol_rel=0.1, limiar_plv=0.8, limiar_skew=0.5,
                     modo_preprocesso="hibrido", f_linha=60.0):
    """
    Núcleo: extrai o teta refinado via FOOOF de `sinal_ctx` (janela de
    contexto longa) e testa se (cf_teta, amp_pico) é consistente com um
    harmônico inteiro, em frequência (razão) e fase (PLV, calculado sobre
    `sinal_cand` -- a janela do candidato). `skew` é opcional (se None,
    o teste de assimetria fica sempre falso).

    Devolve dict: {"cf_teta_fooof", "erro_ajuste_fooof", "plv_harmonico",
    "ordem", "veredito_harmonico"}.
    """
    res_fooof = extrai_cf_teta_fooof(sinal_ctx, fs,
                                     modo_preprocesso=modo_preprocesso,
                                     f_linha=f_linha)

    if not res_fooof["qualidade_ok"]:
        veredito = "SEM_REFERENCIA_TETA"
        plv_val = np.nan
        ordem = None
    else:
        n_max = calcula_n_max(res_fooof["cf_teta"], amp_pico, tol_rel * res_fooof["cf_teta"])
        suspeito, ordem, desvio, ambiguo, candidatos = testa_razao_harmonica(
            res_fooof["cf_teta"], amp_pico, tol_rel, n_max=n_max)

        plv_val = np.nan
        if suspeito and ordem is not None:
            plv_val = compute_plv_harmonico(sinal_cand, fs,
                                           res_fooof["cf_teta"], amp_pico,
                                           ordem)

        skew_alto = skew is not None and abs(skew) > limiar_skew
        plv_alto = plv_val is not None and plv_val > limiar_plv

        if ambiguo:
            veredito = f"AMBIGUO_MULTIPLOS_N ({len(candidatos)} candidatos: {[c[0] for c in candidatos]})"
        elif suspeito and skew_alto and plv_alto:
            veredito = f"SUSPEITO_HARMONICO_FORTE ({ordem}x, PLV={plv_val:.2f})"
        elif suspeito and plv_alto:
            veredito = f"REVISAR_FASE_TRAVADA ({ordem}x, PLV={plv_val:.2f})"
        elif suspeito:
            veredito = f"REVISAR_RAZAO_INTEIRA ({ordem}x)"
        elif skew_alto:
            veredito = "REVISAR_TETA_ASSIMETRICO"
        else:
            veredito = "CLEAN"

    return {
        "cf_teta_fooof": res_fooof["cf_teta"],
        "erro_ajuste_fooof": res_fooof["erro_ajuste"],
        "plv_harmonico": plv_val,
        "ordem": ordem,
        "veredito_harmonico": veredito,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True, help="vencedores.csv da sessao")
    ap.add_argument("--pasta_ns2", required=True, help="pasta com os .ns2")
    ap.add_argument("--saida", default="harmonico.csv", help="caminho de saida")
    ap.add_argument("--janela_contexto_s", type=float, default=45.0,
                    help="janela de contexto p/ FOOOF (default 45s)")
    ap.add_argument("--tol_rel", type=float, default=0.1,
                    help="tolerancia relativa da razao harmonica (default 0.1 = 10%%)")
    ap.add_argument("--limiar_plv", type=float, default=0.8,
                    help="limiar de PLV para harmonio forte (default 0.8)")
    ap.add_argument("--limiar_skew", type=float, default=0.5,
                    help="|skewness| acima disto conta como assimetria")
    ap.add_argument("--modo_preprocesso", default="hibrido",
                    choices=["gaussiana", "cirurgica", "hibrido"],
                    help="limpeza de linha Kuhn antes do FOOOF (default hibrido)")
    ap.add_argument("--f_linha", type=float, default=60.0,
                    help="frequencia da rede eletrica (BR=60, EU=50)")
    args = ap.parse_args()

    df = pd.read_csv(args.csv)
    linhas = []

    print(f"{'Canal':<10} | {'Janela':<14} | {'cf_teta':<8} | {'PLV':<6} | {'skew':<7} | Veredito")
    print("-" * 85)

    for _, r in df.iterrows():
        canal, arquivo = str(r["canal"]), str(r["arquivo"])
        path = os.path.join(args.pasta_ns2, arquivo)
        ini = float(r.get("janela_ini_s", r.get("inicio_s", 0)))
        fim = float(r.get("janela_fim_s", r.get("fim_s", 0)))
        fase_pico, amp_pico = float(r["fase_pico_hz"]), float(r["amp_pico_hz"])

        centro = (ini + fim) / 2
        ctx_ini = max(0, centro - args.janela_contexto_s / 2)
        ctx_fim = centro + args.janela_contexto_s / 2

        try:
            dados, fs, canal_ids = carrega_dados(path)
            chan_idx = canal_ids.index(canal)
            # Sinal para FOOOF (Contexto longo)
            sinal_ctx = fatia_janela(dados, fs, ctx_ini, ctx_fim)[:, chan_idx]
            # Sinal para PLV (Janela do Candidato)
            sinal_cand = fatia_janela(dados, fs, ini, fim)[:, chan_idx]
        except Exception as e:
            linhas.append({"rotulo": r.get("rotulo"), "arquivo": arquivo, "canal": canal,
                            "par": r.get("par"),
                            "janela_ini_s": ini, "janela_fim_s": fim,
                            "janela": f"{ini:.0f}-{fim:.0f}s",
                            "cf_teta_fooof": None, "veredito_harmonico": "ERROR", "motivo": str(e)})
            print(f"{canal:<10} | {ini:.0f}-{fim:.0f}s | ERROR: {e}")
            continue

        try:
            # sinal_cand ja e exatamente a mesma fatia que
            # theta_skewness_for_window(path, canal, ini, fim) releria do
            # zero (mesmo carrega_dados + mesmo fatiamento) -- reusar evita
            # I/O redundante.
            skew, _ = skewness_de_sinal(sinal_cand, fs)
        except Exception:
            skew = None

        res_h = avalia_harmonico(sinal_ctx, sinal_cand, fs, amp_pico, skew=skew,
                                 tol_rel=args.tol_rel, limiar_plv=args.limiar_plv,
                                 limiar_skew=args.limiar_skew,
                                 modo_preprocesso=args.modo_preprocesso,
                                 f_linha=args.f_linha)
        res_fooof = {"cf_teta": res_h["cf_teta_fooof"], "erro_ajuste": res_h["erro_ajuste_fooof"]}
        plv_val = res_h["plv_harmonico"]
        veredito = res_h["veredito_harmonico"]

        cf_str = f"{res_fooof['cf_teta']:.2f}" if res_fooof['cf_teta'] is not None else "n/a"
        plv_str = f"{plv_val:.2f}" if not np.isnan(plv_val) else "n/a"
        skw_str = f"{skew:+.2f}" if skew is not None else "n/a"
        print(f"{canal:<10} | {ini:.0f}-{fim:.0f}s | {cf_str:<8} | {plv_str:<6} | {skw_str:<7} | {veredito}")

        linhas.append({
            "rotulo": r.get("rotulo"), "arquivo": arquivo, "canal": canal,
            "par": r.get("par"),
            "janela_ini_s": ini, "janela_fim_s": fim,
            "janela": f"{ini:.0f}-{fim:.0f}s",
            "cf_teta_fooof": res_fooof["cf_teta"],
            "erro_ajuste_fooof": res_fooof["erro_ajuste"],
            "plv_harmonico": plv_val,
            "skewness_janela": skew, "veredito_harmonico": veredito,
        })

    pd.DataFrame(linhas).to_csv(args.saida, index=False)
    print(f"\nSalvo: {args.saida} ({len(linhas)} linhas)")
    print("Lembrete: este CSV e de AUDITORIA. Nao substitui vencedores.csv.")


if __name__ == "__main__":
    main()
