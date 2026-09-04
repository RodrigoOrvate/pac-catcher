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
from audita_skewness import theta_skewness_for_window  # REUSO, nao duplicacao


def extrai_cf_teta_fooof(sinal, fs, theta_range=(4, 12),
                          theta_cf_bounds=(5, 9.5), theta_bw_limits=(2, 5),
                          min_peak_height=0.05, nperseg_s=1.2,
                          aperiodic_mode='knee'):
    """
    Estima a frequencia central (cf) de teta via FOOOF (metodo modificado
    de Kuhn et al. 2026, LFP_FOOOF).

    Ajusta o componente aperiodico (1/f) e restringe a deteccao a UMA
    Gaussiana na banda teta (4-12 Hz), com cf limitado a 5-9.5 Hz.

    IMPORTANTE sobre aperiodic_mode:
        Default = 'knee' porque LFP real de CA1/DG tem 'knee frequency'
        real (~28 Hz em CA1, ~70 Hz em DG segundo Kuhn et al. 2026).
        'fixed' so deve ser usado em sinais sem componente 1/f ou em
        testes sinteticos com estrutura simples. A escolha de 'fixed'
        baseada em teste sintetico isolado NAO generaliza para LFP real.
    """
    nperseg = int(nperseg_s * fs)
    freqs, psd = welch(sinal, fs=fs, window='hann',
                        nperseg=nperseg, noverlap=nperseg // 2)

    fm = FOOOF(aperiodic_mode=aperiodic_mode, peak_width_limits=theta_bw_limits,
               min_peak_height=min_peak_height, peak_threshold=1.0,
               max_n_peaks=1)

    mask = (freqs >= theta_range[0]) & (freqs <= theta_range[1])
    try:
        fm.fit(freqs[mask], psd[mask], freq_range=theta_range)
    except Exception as e:
        # FOOOF pode falhar em sinais fracos/sem pico detectavel
        return {"cf_teta": None, "teta_detectado": False,
                "erro_ajuste": None, "qualidade_ok": False}

    if not fm.has_model:
        # Modelo nao foi ajustado (pico insuficiente)
        return {"cf_teta": None, "teta_detectado": False,
                "erro_ajuste": None, "qualidade_ok": False}

    try:
        erro_ajuste = fm.get_params('error')
        picos = fm.get_params('peak_params')
    except Exception:
        # Em casos raros, get_params pode falhar mesmo com has_model=True
        return {"cf_teta": None, "teta_detectado": False,
                "erro_ajuste": None, "qualidade_ok": False}

    if erro_ajuste is None:
        erro_ajuste = float('inf')

    cf_teta, teta_detectado = None, False
    if picos is not None and len(picos) > 0:
        cf_cand = picos[0] if picos.ndim == 1 else picos[0, 0]
        if theta_cf_bounds[0] <= cf_cand <= theta_cf_bounds[1]:
            cf_teta, teta_detectado = cf_cand, True

    # TODO: calibrar empiricamente
    qualidade_ok = teta_detectado and erro_ajuste < 0.15
    return {"cf_teta": cf_teta, "teta_detectado": teta_detectado,
            "erro_ajuste": erro_ajuste, "qualidade_ok": qualidade_ok}


def compute_plv_harmonico(sinal, fs, f_theta, f_gamma, n, t_ini, t_fim,
                           bw_theta=2.0, bw_gamma=2.0):
    """
    Calcula o PLV (Phase Locking Value) entre n*phi_theta e phi_gamma.

    Harmonico verdadeiro: PLV alto (~1), fase constante ciclo a ciclo.
    Acoplamento genuino: PLV mais baixo/variavel.

    Filtros narrow-band Butterworth ordem 4:
      - bw_theta=2.0Hz: extrai a fase do teta sem contaminacao de harmônicos
      - bw_gamma=2.0Hz: extrai a fase do gamma em torno de f_gamma

    Nota: Bandas muito estreitas (<1Hz) podem distorcer a fase em frequencias
    baixas (teta). 2Hz e um bom compromisso entre precisao e robustez.
    """
    # Filtragem estreita para extrair fase (Butterworth 4a ordem)
    def narrow_band(sig, f, bw):
        nyq = 0.5 * fs
        lo = max(0.1, f - bw / 2)
        hi = min(fs / 2 - 0.1, f + bw / 2)
        b, a = butter(4, [lo / nyq, hi / nyq], btype='band')
        return filtfilt(b, a, sig)

    # Extrair fase via Hilbert
    s_theta = narrow_band(sinal, f_theta, bw_theta)
    s_gamma = narrow_band(sinal, f_gamma, bw_gamma)

    phi_theta = np.angle(hilbert(s_theta))
    phi_gamma = np.angle(hilbert(s_gamma))

    # PLV n:1 -> |mean(exp(i * (phi_gamma - n * phi_theta)))|
    diff = phi_gamma - (n * phi_theta)
    plv = np.abs(np.mean(np.exp(1j * diff)))

    return plv


def testa_razao_harmonica(fase_hz, amp_hz, tol_rel=0.1, n_max=8):
    """
    Testa se amp_hz ~= n * fase_hz usando tolerancia relativa (10% de f_theta).
    """
    for n in range(2, n_max + 1):
        # Tolerancia relativa: 10% da freq de teta
        tol_abs = tol_rel * fase_hz
        desvio = abs(amp_hz - n * fase_hz)
        if desvio <= tol_abs:
            return True, n, desvio
    return False, None, None


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
    args = ap.parse_args()

    df = pd.read_csv(args.csv)
    linhas = []

    print(f"{'Canal':<10} | {'Janela':<14} | {'cf_teta':<8} | {'PLV':<6} | {'skew':<7} | Veredito")
    print("-" * 85)

    for _, r in df.iterrows():
        canal, arquivo = str(r["canal"]), str(r["arquivo"])
        path = os.path.join(args.pasta_ns2, arquivo)
        ini, fim = float(r["inicio_s"]), float(r["fim_s"])
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
            linhas.append({"rotulo": r.get("rotulo"), "canal": canal,
                            "janela": f"{ini:.0f}-{fim:.0f}s",
                            "cf_teta_fooof": None, "veredito": "ERROR", "motivo": str(e)})
            print(f"{canal:<10} | {ini:.0f}-{fim:.0f}s | ERROR: {e}")
            continue

        res_fooof = extrai_cf_teta_fooof(sinal_ctx, fs)

        try:
            skew, _ = theta_skewness_for_window(path, canal, ini, fim)
        except Exception:
            skew = None

        if not res_fooof["qualidade_ok"]:
            veredito = "SEM_REFERENCIA_TETA"
            plv_val = np.nan
            ordem = None
        else:
            # 1. Teste de Razao de Frequencia (Tolerancia Relativa)
            suspeito, ordem, desvio = testa_razao_harmonica(
                res_fooof["cf_teta"], amp_pico, args.tol_rel)

            # 2. Teste de Rigidez de Fase (PLV) - roda apenas se houver suspeita de razao
            plv_val = np.nan
            if suspeito:
                plv_val = compute_plv_harmonico(sinal_cand, fs,
                                               res_fooof["cf_teta"], amp_pico,
                                               ordem, ini, fim)

            skew_alto = skew is not None and abs(skew) > args.limiar_skew
            plv_alto = plv_val is not None and plv_val > args.limiar_plv

            if suspeito and skew_alto and plv_alto:
                veredito = f"SUSPEITO_HARMONICO_FORTE ({ordem}x, PLV={plv_val:.2f})"
            elif suspeito and plv_alto:
                veredito = f"REVISAR_FASE_TRAVADA ({ordem}x, PLV={plv_val:.2f})"
            elif suspeito:
                veredito = f"REVISAR_RAZAO_INTEIRA ({ordem}x)"
            elif skew_alto:
                veredito = "REVISAR_TETA_ASSIMETRICO"
            else:
                veredito = "CLEAN"

        cf_str = f"{res_fooof['cf_teta']:.2f}" if res_fooof['cf_teta'] is not None else "n/a"
        plv_str = f"{plv_val:.2f}" if not np.isnan(plv_val) else "n/a"
        skw_str = f"{skew:+.2f}" if skew is not None else "n/a"
        print(f"{canal:<10} | {ini:.0f}-{fim:.0f}s | {cf_str:<8} | {plv_str:<6} | {skw_str:<7} | {veredito}")

        linhas.append({
            "rotulo": r.get("rotulo"), "canal": canal,
            "janela": f"{ini:.0f}-{fim:.0f}s",
            "cf_teta_fooof": res_fooof["cf_teta"],
            "erro_ajuste_fooof": res_fooof["erro_ajuste"],
            "plv_harmonico": plv_val,
            "skewness": skew, "veredito": veredito,
        })

    pd.DataFrame(linhas).to_csv(args.saida, index=False)
    print(f"\nSalvo: {args.saida} ({len(linhas)} linhas)")
    print("Lembrete: este CSV e de AUDITORIA. Nao substitui vencedores.csv.")


if __name__ == "__main__":
    main()
