"""
audita_harmonico_hfo.py - Teste de razao harmonica Gamma→HFO via FOOOF e PLV

Complementa audita_harmonico.py: enquanto aquele testa se Gamma é harmônico
de Theta, este testa se HFO (150-250 Hz) é harmônico de Gamma (30-80 Hz).

Mecanismo: Theta assimétrico gera harmônicos no Gamma (bem documentado);
analogamente, Gamma assimétrico pode gerar energia em HFO. O detector de HFO
por limiar de energia simples não distingue esse falso positivo de um ripple
genuíno — este script é a auditoria para essa distinção.

Arquitetura:
  1. FOOOF sobre PSD ampla (4–300 Hz, se fs permitir) para estimar cf_gamma
     (frequência central do pico de Gamma no ajuste aperiodico+periódico).
  2. Teste de razão: f_hfo ≈ n × cf_gamma (tolerância relativa 15%).
  3. PLV(n × phi_gamma, phi_hfo) — harmônico matemático tem fase travada (PLV ≈ 1);
     acoplamento genuíno tem PLV mais baixo e variável.

Veredito:
  CLEAN                   → não é harmônico de Gamma por nenhum critério
  REVISAR_RAZAO_INTEIRA   → razão suspeita, PLV não calculado ou baixo
  REVISAR_FASE_TRAVADA    → razão suspeita + PLV alto (>0.8)
  SUSPEITO_HARMONICO_FORTE → razão + PLV alto + potência HFO/Gamma ratio alto
  SEM_REFERENCIA_GAMMA    → FOOOF não detectou Gamma com qualidade suficiente

CLI:
    python audita_harmonico_hfo.py \
        --csv resultados_triagem_mat.csv \
        --mat DADOS_EXEMPLO/LFP_HG_HFO.mat \
        --saida harmonico_hfo.csv
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.signal import welch, butter, filtfilt, hilbert
import scipy.io

try:
    from fooof import FOOOF
except ImportError:
    sys.stderr.write(
        "[ERRO] Pacote 'fooof' nao encontrado. Instale com:\n"
        "    pip install fooof\n"
    )
    sys.exit(1)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


from utils_harmonico import (
    extrai_cf_gamma_fooof,
    extrai_cf_teta_fooof,
    testa_razao_harmonica,
    calcula_n_max,
    compute_plv_harmonico,
    calcula_ratio_banda
)


# ============================================================
# Carrega .mat ou LFP bruto
# ============================================================

def _carrega_sinal_janela(origem, canal, ini, fim):
    """
    Tenta carregar de .mat (lfpHG/lfpHFO) ou de ns2_utils.
    Retorna (sinal_1d, fs).
    """
    if origem.endswith('.mat'):
        m = scipy.io.loadmat(origem)
        # Procura sinal bruto (lfp, LFP, lfpBruto) ou canal especificado
        for chave in [canal, 'lfp', 'LFP', 'lfpBruto', 'lfpHG', 'lfpHFO']:
            if chave in m:
                arr = np.squeeze(m[chave]).astype(float)
                fs = float(m.get('fs', m.get('Fs', [[1000.0]]))[0][0])
                n_ini = int(ini * fs)
                n_fim = int(fim * fs)
                return arr[n_ini:n_fim], fs
        raise KeyError(f"Nenhuma variável de sinal encontrada em {origem}")
    else:
        # ns2 via ns2_utils
        try:
            from ns2_utils import carrega_dados, fatia_janela
            dados, fs, canal_ids = carrega_dados(origem)
            idx = canal_ids.index(canal)
            return fatia_janela(dados, fs, ini, fim)[:, idx].astype(float), fs
        except Exception as e:
            raise RuntimeError(f"Erro ao carregar {origem}: {e}")


# ============================================================
# MAIN
# ============================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    # Fonte de candidatos
    ap.add_argument("--csv", required=True,
                    help="CSV de triagem com colunas: janela_ini_s, janela_fim_s, "
                         "z_theta_hfo (ou similar), arquivo/canal")
    # Fonte do sinal
    ap.add_argument("--mat", default=None,
                    help="Arquivo .mat com o LFP (bruto de banda larga preferido)")
    ap.add_argument("--pasta_ns2", default=None,
                    help="Pasta com .ns2 (alternativa ao --mat)")
    # Parâmetros
    ap.add_argument("--canal_mat", default="lfpBruto",
                    help="Variável do .mat a usar (default: lfpBruto)")
    ap.add_argument("--z_corte", type=float, default=2.0,
                    help="z mínimo de theta_hfo para entrar na auditoria (default 2.0)")
    ap.add_argument("--janela_contexto_s", type=float, default=30.0,
                    help="Janela de contexto para o FOOOF (default 30s)")
    ap.add_argument("--tol_rel", type=float, default=0.15,
                    help="Tolerância relativa da razão harmônica (default 0.15)")
    ap.add_argument("--limiar_plv", type=float, default=0.8,
                    help="PLV acima disto = fase travada (default 0.8)")
    ap.add_argument("--limiar_ratio", type=float, default=0.3,
                    help="Potência HFO/Gamma acima disto reforça suspeita (default 0.3)")
    ap.add_argument("--f_linha", type=float, default=60.0,
                    help="Freq. rede elétrica para limpeza antes do FOOOF (default 60)")
    ap.add_argument("--saida", default="harmonico_hfo.csv")
    args = ap.parse_args()

    df = pd.read_csv(args.csv)
    linhas = []

    # Detecta coluna de z para HFO
    z_col = next((c for c in ["z_theta_hfo", "z_hfo", "z_score"] if c in df.columns), None)
    if z_col is None:
        print("[AVISO] Nenhuma coluna z_theta_hfo encontrada — auditando todas as janelas.")
        candidatos = df
    else:
        candidatos = df[df[z_col] >= args.z_corte].copy()
        print(f"Candidatos com {z_col} >= {args.z_corte}: {len(candidatos)} de {len(df)} janelas")

    print(f"{'Janela':<14} | {'cf_gamma':<9} | {'ordem':^5} | {'PLV':^6} | "
          f"{'ratio':^6} | Veredito")
    print("-" * 80)

    for _, r in candidatos.iterrows():
        ini = float(r.get("janela_ini_s", r.get("inicio_s", 0)))
        fim = float(r.get("janela_fim_s", r.get("fim_s", ini + 10)))
        arquivo = str(r.get("arquivo", ""))
        canal_str = str(r.get("canal", args.canal_mat))

        # Determina origem do sinal
        if args.mat:
            origem = args.mat
            canal_load = args.canal_mat
        elif args.pasta_ns2 and arquivo:
            origem = os.path.join(args.pasta_ns2, arquivo)
            canal_load = canal_str
        else:
            print(f"[ERRO] Forneça --mat ou --pasta_ns2 para carregar o sinal.")
            return

        # Janela de contexto para FOOOF
        centro = (ini + fim) / 2
        ctx_ini = max(0, centro - args.janela_contexto_s / 2)
        ctx_fim = centro + args.janela_contexto_s / 2

        try:
            sinal_ctx, fs = _carrega_sinal_janela(origem, canal_load, ctx_ini, ctx_fim)
            sinal_cand, _  = _carrega_sinal_janela(origem, canal_load, ini, fim)
        except Exception as e:
            print(f"  {ini:.0f}-{fim:.0f}s: ERRO ao carregar — {e}")
            linhas.append({"janela": f"{ini:.0f}-{fim:.0f}s",
                           "arquivo": arquivo, "canal": canal_str,
                           "par": r.get("par", "theta_hfo"),
                           "janela_ini_s": ini, "janela_fim_s": fim,
                           "veredito_harmonico_hfo": "ERRO_CARGA", "motivo": str(e)})
            continue

        # FOOOF para cf_gamma e cf_teta
        res_gamma = extrai_cf_gamma_fooof(sinal_ctx, fs, f_linha=args.f_linha)
        res_teta = extrai_cf_teta_fooof(sinal_ctx, fs, f_linha=args.f_linha)

        # Potência ratio HFO/Gamma na janela candidata
        p_gamma = calcula_ratio_banda(sinal_cand, fs, (30, 80), (30, 80)) # trick to just get p_gamma
        ratio = calcula_ratio_banda(sinal_cand, fs, (150, min(250, fs * 0.45)), (30, 80))

        if not res_gamma["qualidade_ok"]:
            veredito = "SEM_REFERENCIA_GAMMA"
            ordem_str = "n/a"
            plv_str = "n/a"
            linhas.append({
                "janela": f"{ini:.0f}-{fim:.0f}s", "cf_gamma_fooof": None,
                "erro_fooof": res_gamma["erro_ajuste"], "ordem": None, "plv": None,
                "ratio_hfo_gamma_audit": round(ratio, 4), "veredito_harmonico_hfo": veredito,
                "arquivo": arquivo, "canal": canal_str,
                "par": r.get("par", "theta_hfo"),
                "janela_ini_s": ini, "janela_fim_s": fim,
            })
            print(f"  {ini:.0f}-{fim:.0f}s     | n/a       | n/a   | n/a    | "
                  f"{ratio:.3f}  | {veredito}")
            continue

        cf_gamma = res_gamma["cf_gamma"]
        cf_teta = res_teta["cf_teta"] if res_teta["qualidade_ok"] else None
        
        # Frequência HFO representativa: usa potência média da banda HFO
        # Estima o pico HFO por Welch simples sobre o candidato
        f_welch, psd_cand = welch(sinal_cand, fs=fs,
                                   nperseg=min(len(sinal_cand), int(2 * fs)),
                                   nfft=4 * int(fs))
        mask_hfo = (f_welch >= 100) & (f_welch <= min(300, fs * 0.45))
        if mask_hfo.sum() > 0:
            f_hfo_pico = float(f_welch[mask_hfo][np.argmax(psd_cand[mask_hfo])])
        else:
            f_hfo_pico = 150.0

        tol_abs_gamma = args.tol_rel * cf_gamma
        n_max_gamma = calcula_n_max(cf_gamma, f_hfo_pico, tol_abs_gamma)
        susp_g, ordem_g, desvio_g, ambig_g, cands_g = testa_razao_harmonica(
            cf_gamma, f_hfo_pico, args.tol_rel, n_max=n_max_gamma)
            
        # Teste Teta -> HFO (NOVO)
        susp_t = False
        if cf_teta is not None:
            tol_abs_teta = args.tol_rel * cf_teta
            n_max_teta = calcula_n_max(cf_teta, f_hfo_pico, tol_abs_teta)
            susp_t, ordem_t, desvio_t, ambig_t, cands_t = testa_razao_harmonica(
                cf_teta, f_hfo_pico, args.tol_rel, n_max=n_max_teta)

        plv_g = np.nan
        if susp_g and ordem_g is not None:
            try:
                plv_g = compute_plv_harmonico(sinal_cand, fs, cf_gamma, f_hfo_pico, ordem_g, bw_fase=5.0)
            except Exception:
                plv_g = np.nan
                
        plv_t = np.nan
        if susp_t and ordem_t is not None:
            try:
                plv_t = compute_plv_harmonico(sinal_cand, fs, cf_teta, f_hfo_pico, ordem_t, bw_fase=2.0)
            except Exception:
                plv_t = np.nan

        plv_alto_g = not np.isnan(plv_g) and plv_g > args.limiar_plv
        plv_alto_t = not np.isnan(plv_t) and plv_t > args.limiar_plv
        ratio_alto = ratio > args.limiar_ratio

        if ambig_g or (susp_t and ambig_t):
            len_g = len(cands_g) if susp_g else 0
            len_t = len(cands_t) if susp_t else 0
            veredito = f"AMBIGUO_MULTIPLOS_N ({len_g} em Gama, {len_t} em Teta)"
            ordem_str = f"g:{ordem_g} t:{ordem_t if susp_t else 'n/a'}"
            plv_str2 = f"g:{plv_g:.2f}"
        elif susp_t and ordem_t > 30 and not (plv_alto_t and ratio_alto):
            veredito = f"REVISAR_ORDEM_ALTA_TETA_SEM_CORROBORACAO ({ordem_t}x)"
            ordem_str = str(ordem_t)
            plv_str2 = f"{plv_t:.3f}"
        elif susp_g and ordem_g > 30 and not (plv_alto_g and ratio_alto):
            veredito = f"REVISAR_ORDEM_ALTA_GAMA_SEM_CORROBORACAO ({ordem_g}x)"
            ordem_str = str(ordem_g)
            plv_str2 = f"{plv_g:.3f}"
        elif (susp_g and plv_alto_g and ratio_alto) or (susp_t and plv_alto_t and ratio_alto):
            qual = "Gama" if susp_g else "Teta"
            ordem_f = ordem_g if susp_g else ordem_t
            plv_f = plv_g if susp_g else plv_t
            veredito = f"SUSPEITO_HARMONICO_FORTE_{qual} ({ordem_f}x, PLV={plv_f:.2f})"
            ordem_str = str(ordem_f)
            plv_str2 = f"{plv_f:.3f}"
        elif (susp_g and plv_alto_g) or (susp_t and plv_alto_t):
            qual = "Gama" if susp_g else "Teta"
            ordem_f = ordem_g if susp_g else ordem_t
            plv_f = plv_g if susp_g else plv_t
            veredito = f"REVISAR_FASE_TRAVADA_{qual} ({ordem_f}x, PLV={plv_f:.2f})"
            ordem_str = str(ordem_f)
            plv_str2 = f"{plv_f:.3f}"
        elif susp_g or susp_t:
            qual = "Gama" if susp_g else "Teta"
            ordem_f = ordem_g if susp_g else ordem_t
            plv_f = plv_g if susp_g else plv_t
            veredito = f"REVISAR_RAZAO_INTEIRA_{qual} ({ordem_f}x)"
            ordem_str = str(ordem_f)
            plv_str2 = f"{plv_f:.3f}" if not np.isnan(plv_f) else "n/a"
        else:
            veredito = "CLEAN"
            ordem_str = "n/a"
            plv_str2 = "n/a"

        cf_str = f"{cf_gamma:.1f} (T:{cf_teta:.1f})" if cf_teta else f"{cf_gamma:.1f}"
        print(f"  {ini:.0f}-{fim:.0f}s     | {cf_str:<9} | {ordem_str:^5} | "
              f"{plv_str2:<6} | {ratio:.3f}  | {veredito}")

        linhas.append({
            "janela": f"{ini:.0f}-{fim:.0f}s",
            "arquivo": arquivo, "canal": canal_str,
            "par": r.get("par", "theta_hfo"),
            "janela_ini_s": ini, "janela_fim_s": fim,
            "cf_gamma_fooof": round(cf_gamma, 2) if cf_gamma else None,
            "cf_teta_fooof_hfo": round(cf_teta, 2) if cf_teta else None,
            "erro_fooof": round(res_gamma["erro_ajuste"], 4) if res_gamma["erro_ajuste"] else None,
            "expoente_gamma_fooof": round(res_gamma["expoente_gamma"], 4) if res_gamma.get("expoente_gamma") else None,
            "knee_gamma_fooof": round(res_gamma["knee_gamma"], 4) if res_gamma.get("knee_gamma") else None,
            "f_hfo_pico": round(f_hfo_pico, 1),
            "ordem_harmonico_gama": ordem_g if susp_g else None,
            "ordem_harmonico_teta": ordem_t if susp_t else None,
            "plv_gamma_hfo": round(plv_g, 4) if not np.isnan(plv_g) else None,
            "plv_teta_hfo": round(plv_t, 4) if not np.isnan(plv_t) else None,
            "ratio_hfo_gamma_audit": round(ratio, 4),
            "veredito_harmonico_hfo": veredito,
        })

    out = pd.DataFrame(linhas)
    out.to_csv(args.saida, index=False)
    print(f"\nSalvo: {args.saida} ({len(out)} linhas)")
    print("Nota: SEM_REFERENCIA_GAMMA = FOOOF nao detectou pico de Gamma na janela de contexto.")
    print("CLEAN = nao e' harmonico de Gamma por nenhum dos tres criterios.")


if __name__ == "__main__":
    main()
