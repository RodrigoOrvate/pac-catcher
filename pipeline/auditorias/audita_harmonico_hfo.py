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


# ============================================================
# FOOOF para cf_gamma
# ============================================================

def extrai_cf_gamma_fooof(sinal, fs, fit_range=None, gamma_cf_bounds=(25, 90),
                           gamma_bw_limits=(4, 30), min_peak_height=0.05,
                           nperseg_s=1.2, aperiodic_mode='knee', max_n_peaks=6,
                           preprocessar_linha=True, f_linha=60.0):
    """
    Estima a frequência central (cf) do pico de Gamma via FOOOF.

    Diferente de cf_teta (estreito, ~2-5 Hz BW), Gamma é mais largo (10-30 Hz BW)
    e menos estável — o FOOOF pode encontrar sub-picos dentro da banda. Este
    método escolhe o pico de MAIOR AMPLITUDE dentro de gamma_cf_bounds.

    fit_range: default (4, min(300, fs/2 * 0.95)) — cobrir todo o espectro
               disponível para ancorar a curva aperiodica.

    Retorna dict com:
      cf_gamma, gamma_detectado, erro_ajuste, qualidade_ok, n_picos,
      cf_gamma_alternativo (lista de todos os picos em gamma_cf_bounds)
    """
    if fit_range is None:
        f_max = min(300.0, fs * 0.5 * 0.95)
        fit_range = (4.0, f_max)

    nperseg = int(nperseg_s * fs)
    nfft = max(nperseg, 4 * int(fs))  # resolução >= 0.25 Hz

    freqs, psd = welch(sinal, fs=fs, window='hann',
                        nperseg=nperseg, noverlap=nperseg // 2, nfft=nfft)

    # Limpeza de linha simples: zera bins em 60/120/180 Hz no PSD
    if preprocessar_linha:
        for h in [1, 2, 3]:
            f_linha_h = f_linha * h
            if f_linha_h >= fs / 2:
                break
            mask_linha = np.abs(freqs - f_linha_h) < 2.0
            if mask_linha.any():
                # Interpola linear nos bins vizinhos
                idx = np.where(mask_linha)[0]
                i0, i1 = max(0, idx[0] - 2), min(len(psd) - 1, idx[-1] + 2)
                if i0 < i1:
                    psd[idx] = np.interp(freqs[idx],
                                          [freqs[i0], freqs[i1]],
                                          [psd[i0], psd[i1]])

    fm = FOOOF(aperiodic_mode=aperiodic_mode,
               peak_width_limits=gamma_bw_limits,
               min_peak_height=min_peak_height,
               peak_threshold=1.0,
               max_n_peaks=max_n_peaks)
    try:
        fm.fit(freqs, psd, freq_range=fit_range)
    except Exception:
        return {"cf_gamma": None, "gamma_detectado": False,
                "erro_ajuste": None, "qualidade_ok": False,
                "n_picos": 0, "cf_gamma_alternativo": []}

    if not fm.has_model:
        return {"cf_gamma": None, "gamma_detectado": False,
                "erro_ajuste": None, "qualidade_ok": False,
                "n_picos": 0, "cf_gamma_alternativo": []}

    try:
        erro_ajuste = fm.get_params('error')
        todos_picos = fm.get_params('peak_params')
    except Exception:
        return {"cf_gamma": None, "gamma_detectado": False,
                "erro_ajuste": None, "qualidade_ok": False,
                "n_picos": 0, "cf_gamma_alternativo": []}

    if erro_ajuste is None:
        erro_ajuste = float('inf')

    cf_gamma = None
    gamma_detectado = False
    alternativas = []

    if todos_picos is not None and len(todos_picos) > 0:
        picos_arr = todos_picos if todos_picos.ndim > 1 else todos_picos.reshape(1, -1)
        mask = ((picos_arr[:, 0] >= gamma_cf_bounds[0]) &
                (picos_arr[:, 0] <= gamma_cf_bounds[1]))
        candidatos = picos_arr[mask]
        if len(candidatos) > 0:
            # Pico de maior amplitude (potência periódica)
            cf_gamma = float(candidatos[np.argmax(candidatos[:, 1]), 0])
            gamma_detectado = True
            alternativas = [float(c[0]) for c in candidatos]

    n_picos = int(len(todos_picos)) if todos_picos is not None else 0
    # Limiar de qualidade: erro < 0.20 (Gamma é mais ruidoso que Theta)
    qualidade_ok = gamma_detectado and erro_ajuste < 0.20

    return {
        "cf_gamma": cf_gamma,
        "gamma_detectado": gamma_detectado,
        "erro_ajuste": erro_ajuste,
        "qualidade_ok": qualidade_ok,
        "n_picos": n_picos,
        "cf_gamma_alternativo": alternativas,
    }


# ============================================================
# PLV Gamma → HFO
# ============================================================

def compute_plv_gamma_hfo(sinal, fs, f_gamma, f_hfo, n, bw_gamma=5.0, bw_hfo=5.0):
    """
    PLV entre n × phi_gamma e phi_hfo.
    Harmônico matemático: PLV ≈ 1 (fase travada ciclo a ciclo).
    Acoplamento genuíno: PLV mais baixo e variável.
    """
    nyq = fs * 0.5

    def narrow_band(sig, f, bw):
        lo = max(0.5, f - bw / 2)
        hi = min(nyq * 0.98, f + bw / 2)
        if lo >= hi:
            return sig
        b, a = butter(4, [lo / nyq, hi / nyq], btype='band')
        return filtfilt(b, a, sig)

    s_gamma = narrow_band(sinal, f_gamma, bw_gamma)
    s_hfo   = narrow_band(sinal, f_hfo,   bw_hfo)

    phi_gamma = np.angle(hilbert(s_gamma))
    phi_hfo   = np.angle(hilbert(s_hfo))

    diff = phi_hfo - (n * phi_gamma)
    return float(np.abs(np.mean(np.exp(1j * diff))))


# ============================================================
# Teste de razão harmônica
# ============================================================

def testa_razao_harmonica_hfo(cf_gamma, f_hfo, tol_rel=0.15, n_max=8):
    """
    Testa se f_hfo ≈ n × cf_gamma com tolerância relativa (15% de cf_gamma).
    Tolerância maior que em theta→gamma porque Gamma é mais variável.
    """
    for n in range(2, n_max + 1):
        tol_abs = tol_rel * cf_gamma
        if abs(f_hfo - n * cf_gamma) <= tol_abs:
            return True, n, abs(f_hfo - n * cf_gamma)
    return False, None, None


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
                           "veredito": "ERRO_CARGA", "motivo": str(e)})
            continue

        # FOOOF para cf_gamma
        res = extrai_cf_gamma_fooof(sinal_ctx, fs, f_linha=args.f_linha)

        # Potência ratio HFO/Gamma na janela candidata
        def pot_banda(sig, lo, hi, fs_):
            nyq = fs_ * 0.5
            lo = min(lo, nyq * 0.9)
            hi = min(hi, nyq * 0.98)
            if lo >= hi:
                return 1e-12
            b, a = butter(3, [lo / nyq, hi / nyq], btype='band')
            filt = filtfilt(b, a, sig)
            return float(np.mean(filt ** 2)) + 1e-12

        p_gamma = pot_banda(sinal_cand, 30, 80, fs)
        p_hfo   = pot_banda(sinal_cand, 150, min(250, fs * 0.45), fs)
        ratio = p_hfo / p_gamma

        if not res["qualidade_ok"]:
            veredito = "SEM_REFERENCIA_GAMMA"
            ordem_str = "n/a"
            plv_str = "n/a"
            linhas.append({
                "janela": f"{ini:.0f}-{fim:.0f}s", "cf_gamma": None,
                "erro_fooof": res["erro_ajuste"], "ordem": None, "plv": None,
                "ratio_hfo_gamma": round(ratio, 4), "veredito": veredito,
                "arquivo": arquivo, "canal": canal_str
            })
            print(f"  {ini:.0f}-{fim:.0f}s     | n/a       | n/a   | n/a    | "
                  f"{ratio:.3f}  | {veredito}")
            continue

        cf_gamma = res["cf_gamma"]
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

        suspeito, ordem, desvio = testa_razao_harmonica_hfo(
            cf_gamma, f_hfo_pico, args.tol_rel)

        plv = np.nan
        if suspeito and ordem is not None:
            try:
                plv = compute_plv_gamma_hfo(sinal_cand, fs, cf_gamma, f_hfo_pico, ordem)
            except Exception:
                plv = np.nan

        plv_alto = not np.isnan(plv) and plv > args.limiar_plv
        ratio_alto = ratio > args.limiar_ratio

        if suspeito and plv_alto and ratio_alto:
            veredito = f"SUSPEITO_HARMONICO_FORTE ({ordem}x, PLV={plv:.2f})"
        elif suspeito and plv_alto:
            veredito = f"REVISAR_FASE_TRAVADA ({ordem}x, PLV={plv:.2f})"
        elif suspeito:
            veredito = f"REVISAR_RAZAO_INTEIRA ({ordem}x)"
        else:
            veredito = "CLEAN"

        cf_str = f"{cf_gamma:.1f}"
        ord_str = str(ordem) if ordem is not None else "n/a"
        plv_str2 = f"{plv:.3f}" if not np.isnan(plv) else "n/a"
        print(f"  {ini:.0f}-{fim:.0f}s     | {cf_str:<9} | {ord_str:^5} | "
              f"{plv_str2:<6} | {ratio:.3f}  | {veredito}")

        linhas.append({
            "janela": f"{ini:.0f}-{fim:.0f}s",
            "arquivo": arquivo, "canal": canal_str,
            "cf_gamma_fooof": round(cf_gamma, 2),
            "erro_fooof": round(res["erro_ajuste"], 4) if res["erro_ajuste"] else None,
            "f_hfo_pico": round(f_hfo_pico, 1),
            "ordem_harmonico": ordem,
            "plv_gamma_hfo": round(plv, 4) if not np.isnan(plv) else None,
            "ratio_hfo_gamma": round(ratio, 4),
            "veredito": veredito,
        })

    out = pd.DataFrame(linhas)
    out.to_csv(args.saida, index=False)
    print(f"\nSalvo: {args.saida} ({len(out)} linhas)")
    print("Nota: SEM_REFERENCIA_GAMMA = FOOOF nao detectou pico de Gamma na janela de contexto.")
    print("CLEAN = nao e' harmonico de Gamma por nenhum dos tres criterios.")


if __name__ == "__main__":
    main()
