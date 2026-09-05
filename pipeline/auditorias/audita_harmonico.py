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
from linha_noise_kuhn import aplica_modo  # REUSO: limpeza de linha Kuhn (60Hz)


def extrai_cf_teta_fooof(sinal, fs, fit_range=(4, 100), theta_range=(4, 12),
                          theta_cf_bounds=(5, 9.5), theta_bw_limits=(2, 5),
                          min_peak_height=0.05, nperseg_s=1.2,
                          aperiodic_mode='knee', max_n_peaks=4,
                          preprocessar_linha=True, modo_preprocesso='hibrido',
                          f_linha=60.0):
    """
    Estima a frequencia central (cf) de teta via FOOOF (metodo modificado
    de Kuhn et al. 2026, LFP_FOOOF).

    ARQUITETURA DE DOIS PASSOS (segue o artigo, Eqs. 1-5 e Tabela 1):
        1. Fit amplo do modelo completo (1/f aperiodico + todos os picos
           periodicos) sobre uma faixa larga (default 4-100Hz). Isso da
           ao FOOOF espaco dinamico suficiente para ancorar a lei de
           potencia 1/f^n com confianca estatistica - 8Hz de largura
           (theta_range sozinho) e' pouquissimo para estimar expoente.
        2. Extracao do pico de teta por filtragem dos picos ja' ajustados
           via theta_cf_bounds (default 5-9.5Hz). NAO re-ajustamos um
           modelo novo dentro de uma janela estreita.

    Por que o fit amplo (e nao so' teta):
        O artigo reporta erros de "full model fit" (identico ao que
        fm.get_params('error') retorna, "mean absolute difference of full
        model fit" segundo a secao "Errors estimation") na faixa de
        0.014-0.048 mesmo em casos dificeis. Quando restringimos o fit
        a 4-12Hz, o erro sobe para ~0.18 no mesmo sinal sintetico
        realista: o FOOOF nao tem informacao suficiente para ancorar a
        curva aperiodica, entao o residuo explode.

    IMPORTANTE sobre aperiodic_mode:
        Default = 'knee' porque LFP real de CA1/DG tem 'knee frequency'
        real (~28 Hz em CA1, ~70 Hz em DG segundo Kuhn et al. 2026).
        'fixed' so deve ser usado em sinais sem componente 1/f ou em
        testes sinteticos com estrutura simples.

    Parametros:
        fit_range: tupla (f_min, f_max) para o fit amplo do FOOOF.
                   Default (4, 100) - cobre 1/f^slope bem abaixo do
                   knee (~28Hz) e a maior parte da banda gamma. Acima
                   de 100Hz, o modelo aperiodico de 1 knee comeca a
                   divergir e exige os modelos 2exp/3exp do artigo.
        theta_range, theta_cf_bounds: usados APENAS para filtrar picos
                   ja' ajustados (passo 2), NAO para restringir o fit.
        preprocessar_linha: se True (default), limpa o ruido de linha antes
                   do fit. modo_preprocesso: 'gaussiana' (subtrai a
                   gaussiana em log10 em 60/120/180) | 'cirurgica' (repoe a
                   banda ±2Hz pela 1/f so em 60Hz) | 'hibrido' (default;
                   repoe banda em 60/120/180 — vencedor da comparacao
                   compara_preprocesso_linha.py, menor erro sintetico sem
                   mutilar oscilacoes reais). f_linha: frequencia da rede
                   (BR=60).
    """
    nperseg = int(nperseg_s * fs)
    # nfft=4000 (zero-padding): segue Kuhn et al. 2026 (LFP_FOOOF).
    # Welch ainda janelado em 1.2s (resolucao estatistica real do espectro),
    # mas FFT em 4000 pontos interpola o espectro para grade fina
    # (~0.25 Hz/bin), dando ao FOOOF pontos suficientes para convergir
    # em Gaussiana de 2-5 Hz sem instabilidade numerica.
    nfft = 4000 if nperseg <= 4000 else nperseg
    freqs, psd = welch(sinal, fs=fs, window='hann',
                        nperseg=nperseg, noverlap=nperseg // 2,
                        nfft=nfft)

    # Limpeza de linha (receita Kuhn) ANTES do fit: a linha de 60Hz (BR)
    # cai dentro do gamma e, se o FOOOF a tratar como pico, infla o erro e
    # rouba vaga de pico. Modo default 'gaussiana' (subtrai a gaussiana em
    # log10 em 60/120/180); 'cirurgica'/'hibrido' disponiveis para comparar.
    if preprocessar_linha:
        psd = aplica_modo(freqs, psd, fs, modo_preprocesso,
                          f_linha=f_linha, verbose=False)

    # max_n_peaks=4: artigo detecta slow_gamma, fast_gamma, ripples etc.
    # no mesmo fit. Restringir a 1 so' faz sentido na hora de extrair por
    # banda (passo 2), nao no fit em si.
    fm = FOOOF(aperiodic_mode=aperiodic_mode, peak_width_limits=theta_bw_limits,
               min_peak_height=min_peak_height, peak_threshold=1.0,
               max_n_peaks=max_n_peaks)

    # PASSO 1: fit amplo sobre fit_range - dados de entrada NAO restringidos
    # a banda de teta. O FOOOF recebe o PSD inteiro (ou no max ate 100Hz)
    # e ajusta o modelo completo (aperiodico + periodicos) numa so' passada.
    try:
        fm.fit(freqs, psd, freq_range=fit_range)
    except Exception as e:
        # FOOOF pode falhar em sinais fracos/sem pico detectavel
        return {"cf_teta": None, "teta_detectado": False,
                "erro_ajuste": None, "qualidade_ok": False,
                "n_picos": 0}

    if not fm.has_model:
        # Modelo nao foi ajustado (pico insuficiente)
        return {"cf_teta": None, "teta_detectado": False,
                "erro_ajuste": None, "qualidade_ok": False,
                "n_picos": 0}

    try:
        erro_ajuste = fm.get_params('error')
        todos_picos = fm.get_params('peak_params')
    except Exception:
        # Em casos raros, get_params pode falhar mesmo com has_model=True
        return {"cf_teta": None, "teta_detectado": False,
                "erro_ajuste": None, "qualidade_ok": False,
                "n_picos": 0}

    if erro_ajuste is None:
        erro_ajuste = float('inf')

    # PASSO 2: extrair pico de teta por filtragem dos picos ja' ajustados.
    # Isso replica o passo "Detection range" da Tabela 1 do artigo:
    # o fit ja' foi feito na faixa ampla; agora restringimos o ROTULO
    # do pico, nao o modelo.
    cf_teta, teta_detectado = None, False
    if todos_picos is not None and len(todos_picos) > 0:
        # picos vem como (N, 3): cf, amp, bw. Garantir 2D
        picos_arr = todos_picos if todos_picos.ndim > 1 else todos_picos.reshape(1, -1)
        candidatos = picos_arr[(picos_arr[:, 0] >= theta_cf_bounds[0]) &
                                (picos_arr[:, 0] <= theta_cf_bounds[1])]
        if len(candidatos) > 0:
            # Se houver mais de um pico na banda, escolhemos o de maior
            # amplitude (potencia do pico periodico)
            cf_teta = float(candidatos[np.argmax(candidatos[:, 1]), 0])
            teta_detectado = True

    # Limiar 0.15: agora comparável aos valores do artigo (0.014-0.048
    # em casos bem-comportados; ate ~0.15 em casos "ruins" do modelo
    # mais simples 1exp). Com o fit amplo, este limiar vira conservador,
    # nao apertado. Calibracao empirica fina em LFP real segue pendente
    # mas NAO e' a unica coisa que falta.
    qualidade_ok = teta_detectado and erro_ajuste < 0.15
    n_picos = int(len(todos_picos)) if todos_picos is not None else 0
    return {"cf_teta": cf_teta, "teta_detectado": teta_detectado,
            "erro_ajuste": erro_ajuste, "qualidade_ok": qualidade_ok,
            "n_picos": n_picos}


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

        res_fooof = extrai_cf_teta_fooof(sinal_ctx, fs,
                                         modo_preprocesso=args.modo_preprocesso,
                                         f_linha=args.f_linha)

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
