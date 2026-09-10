"""
utils_harmonico.py
==========================================
Coleção de funções compartilhadas para auditoria de harmônicos e extração
de características via FOOOF (Kuhn et al. 2026).

Consolida lógica antes duplicada em audita_harmonico.py e audita_harmonico_hfo.py,
e introduz o teste de razão harmônica com detecção de ambiguidade de múltiplos
candidatos para bandas de alta frequência (HFO).
"""

import math
import numpy as np
from scipy.signal import welch, butter, filtfilt, hilbert

try:
    from fooof import FOOOF
except ImportError:
    pass  # O script principal lidará com isso

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from pipeline.auditorias.linha_noise_kuhn import aplica_modo


def calcula_n_max(f_fase, f_amp_max_banda, tolerancia_hz):
    """
    n_max = ceil((f_amp_max_banda + tolerancia_hz) / f_fase)
    Ex: teta=5Hz, banda HFO até 250Hz, tolerancia=1Hz -> n_max=51.
    """
    if f_fase <= 0:
        return 1
    return math.ceil((f_amp_max_banda + tolerancia_hz) / f_fase)


def testa_razao_harmonica(fase_hz, amp_hz, tol_rel=0.1, n_max=8):
    """
    Testa se amp_hz é aproximadamente n * fase_hz.
    
    Retorna (suspeito, ordem, desvio, ambiguo, candidatos).

    NOVO: em vez de retornar no primeiro n que satisfaz a tolerância,
    varre TODO o range [2, n_max] e coleta todos os candidatos. Se mais
    de um n cair dentro da tolerância (mais provável quanto maior o
    n_max, ex. para HFO com n_max~50), o veredito retorna ambiguo=True.
    """
    tol_abs = tol_rel * fase_hz
    candidatos = []
    
    for n in range(2, n_max + 1):
        desvio = abs(amp_hz - n * fase_hz)
        if desvio <= tol_abs:
            candidatos.append((n, desvio))
            
    if not candidatos:
        return False, None, None, False, []

    ambiguo = len(candidatos) > 1
    # melhor ajuste
    candidatos.sort(key=lambda x: x[1])
    ordem, desvio = candidatos[0]
    
    return True, ordem, desvio, ambiguo, candidatos


def calcula_ratio_banda(sinal, fs, banda_num, banda_den):
    """
    Retorna potência integrada em banda_num / potência integrada em banda_den.
    Usado como corroboração adicional: um HFO com razão harmônica suspeita E
    ratio HFO/Gamma muito alto é evidência forte de contaminação.
    """
    nyq = fs * 0.5
    
    def pot(lo, hi):
        hi = min(hi, nyq * 0.98)
        lo = min(lo, nyq * 0.9)
        if lo >= hi:
            return 1e-12
        b, a = butter(3, [lo / nyq, hi / nyq], btype='band')
        filt = filtfilt(b, a, sinal)
        return float(np.mean(filt ** 2)) + 1e-12
        
    p_num = pot(banda_num[0], banda_num[1])
    p_den = pot(banda_den[0], banda_den[1])
    return p_num / p_den


def _knee_valido(knee, exp, fit_range):
    """
    Diagnostico padrao de identificabilidade do parametro de knee no FOOOF
    (Donoghue et al. 2020): a frequencia de joelho derivada f = knee^(1/chi)
    deve cair DENTRO da faixa que foi de fato ajustada. Fora disso (ou knee<=0)
    e sinal de parametro nao-identificavel — o otimizador empurrou o knee pra
    um valor sem sustentacao nos dados (pode ir a zero, negativo, ou explodir
    pra valores astronomicos), normalmente absorvendo ruido/variancia residual
    em vez de capturar curvatura aperiodica real. Nesses casos o expoente
    (usado como proxy de E/I) tambem fica contaminado, porque knee e expoente
    trocam vies entre si quando o joelho nao esta bem restringido.
    """
    if knee is None or exp is None or knee <= 0 or exp <= 0:
        return False
    try:
        f_joelho = knee ** (1.0 / exp)
    except Exception:
        return False
    return fit_range[0] <= f_joelho <= fit_range[1]


def extrai_cf_teta_fooof(sinal, fs, fit_range=(2.0, 45.0), theta_range=(4, 12),
                          theta_cf_bounds=(5, 9.5), theta_bw_limits=(2, 5),
                          min_peak_height=0.05, nperseg_s=1.2,
                          aperiodic_mode='knee', max_n_peaks=4,
                          preprocessar_linha=True, modo_preprocesso='hibrido',
                          f_linha=60.0):
    nperseg = int(nperseg_s * fs)
    nfft = 4000 if nperseg <= 4000 else nperseg
    freqs, psd = welch(sinal, fs=fs, window='hann',
                        nperseg=nperseg, noverlap=nperseg // 2,
                        nfft=nfft)

    if preprocessar_linha:
        psd = aplica_modo(freqs, psd, fs, modo_preprocesso,
                          f_linha=f_linha, verbose=False)

    fm = FOOOF(aperiodic_mode=aperiodic_mode, peak_width_limits=theta_bw_limits,
               min_peak_height=min_peak_height, peak_threshold=1.0,
               max_n_peaks=max_n_peaks)

    dict_vazio = {
        "cf_teta": None, "teta_detectado": False,
        "erro_ajuste": None, "qualidade_ok": False, "n_picos": 0,
        "expoente_teta": None, "knee_teta": None, "offset_teta": None, "r2_teta": None,
        "knee_valido_teta": False
    }

    try:
        fm.fit(freqs, psd, freq_range=fit_range)
    except Exception:
        return dict_vazio

    if not fm.has_model:
        return dict_vazio

    try:
        erro_ajuste = fm.get_params('error')
        todos_picos = fm.get_params('peak_params')
        ap_params = fm.get_params('aperiodic_params')
        r2 = fm.get_params('r_squared')
    except Exception:
        return dict_vazio

    if erro_ajuste is None:
        erro_ajuste = float('inf')

    cf_teta, teta_detectado = None, False
    if todos_picos is not None and len(todos_picos) > 0:
        picos_arr = todos_picos if todos_picos.ndim > 1 else todos_picos.reshape(1, -1)
        candidatos = picos_arr[(picos_arr[:, 0] >= theta_cf_bounds[0]) &
                                (picos_arr[:, 0] <= theta_cf_bounds[1])]
        if len(candidatos) > 0:
            cf_teta = float(candidatos[np.argmax(candidatos[:, 1]), 0])
            teta_detectado = True

    qualidade_ok = teta_detectado and erro_ajuste < 0.15
    n_picos = int(len(todos_picos)) if todos_picos is not None else 0

    exp_val = float(ap_params[2]) if (ap_params is not None and len(ap_params) == 3) else (
        float(ap_params[1]) if (ap_params is not None and len(ap_params) == 2) else None
    )
    knee_val = float(ap_params[1]) if (ap_params is not None and len(ap_params) == 3) else None
    offset_val = float(ap_params[0]) if (ap_params is not None and len(ap_params) >= 1) else None
    r2_val = float(r2) if r2 is not None else None
    knee_valido = _knee_valido(knee_val, exp_val, fit_range)

    return {
        "cf_teta": cf_teta,
        "teta_detectado": teta_detectado,
        "erro_ajuste": erro_ajuste,
        "qualidade_ok": qualidade_ok,
        "n_picos": n_picos,
        "expoente_teta": exp_val,
        "knee_teta": knee_val,
        "offset_teta": offset_val,
        "r2_teta": r2_val,
        "knee_valido_teta": knee_valido
    }


def extrai_cf_gamma_fooof(sinal, fs, fit_range=None, gamma_cf_bounds=(25, 90),
                           gamma_bw_limits=(4, 30), min_peak_height=0.05,
                           nperseg_s=1.2, aperiodic_mode='knee', max_n_peaks=6,
                           preprocessar_linha=True, modo_preprocesso='hibrido',
                           f_linha=60.0):
    if fit_range is None:
        f_max = min(250.0, fs * 0.5 * 0.95)
        fit_range = (35.0, f_max)

    nperseg = int(nperseg_s * fs)
    nfft = max(nperseg, 4 * int(fs)) 

    freqs, psd = welch(sinal, fs=fs, window='hann',
                        nperseg=nperseg, noverlap=nperseg // 2, nfft=nfft)

    if preprocessar_linha:
        psd = aplica_modo(freqs, psd, fs, modo_preprocesso,
                          f_linha=f_linha, verbose=False)

    fm = FOOOF(aperiodic_mode=aperiodic_mode,
               peak_width_limits=gamma_bw_limits,
               min_peak_height=min_peak_height,
               peak_threshold=1.0,
               max_n_peaks=max_n_peaks)

    dict_vazio = {
        "cf_gamma": None, "gamma_detectado": False,
        "erro_ajuste": None, "qualidade_ok": False,
        "n_picos": 0, "cf_gamma_alternativo": [],
        "expoente_gamma": None, "knee_gamma": None, "offset_gamma": None, "r2_gamma": None,
        "knee_valido_gamma": False
    }

    try:
        fm.fit(freqs, psd, freq_range=fit_range)
    except Exception:
        return dict_vazio

    if not fm.has_model:
        return dict_vazio

    try:
        erro_ajuste = fm.get_params('error')
        todos_picos = fm.get_params('peak_params')
        ap_params = fm.get_params('aperiodic_params')
        r2 = fm.get_params('r_squared')
    except Exception:
        return dict_vazio

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
            cf_gamma = float(candidatos[np.argmax(candidatos[:, 1]), 0])
            gamma_detectado = True
            alternativas = [float(c[0]) for c in candidatos]

    n_picos = int(len(todos_picos)) if todos_picos is not None else 0
    qualidade_ok = gamma_detectado and erro_ajuste < 0.20

    exp_val = float(ap_params[2]) if (ap_params is not None and len(ap_params) == 3) else (
        float(ap_params[1]) if (ap_params is not None and len(ap_params) == 2) else None
    )
    knee_val = float(ap_params[1]) if (ap_params is not None and len(ap_params) == 3) else None
    offset_val = float(ap_params[0]) if (ap_params is not None and len(ap_params) >= 1) else None
    r2_val = float(r2) if r2 is not None else None

    return {
        "cf_gamma": cf_gamma,
        "gamma_detectado": gamma_detectado,
        "erro_ajuste": erro_ajuste,
        "qualidade_ok": qualidade_ok,
        "n_picos": n_picos,
        "cf_gamma_alternativo": alternativas,
        "expoente_gamma": exp_val,
        "knee_gamma": knee_val,
        "offset_gamma": offset_val,
        "r2_gamma": r2_val,
        "knee_valido_gamma": _knee_valido(knee_val, exp_val, fit_range)
    }


def compute_plv_harmonico(sinal, fs, f_fase, f_amp, n,
                           bw_fase=2.0, bw_amp=5.0):
    """
    Calcula o PLV entre n*phi_fase e phi_amp.
    Usado tanto para Theta->Gamma quanto para Gamma->HFO ou Theta->HFO.
    """
    nyq = 0.5 * fs
    def narrow_band(sig, f, bw):
        lo = max(0.5, f - bw / 2)
        hi = min(fs / 2 - 0.5, f + bw / 2)
        if lo >= hi:
            return sig
        b, a = butter(4, [lo / nyq, hi / nyq], btype='band')
        return filtfilt(b, a, sig)

    s_fase = narrow_band(sinal, f_fase, bw_fase)
    s_amp = narrow_band(sinal, f_amp, bw_amp)

    phi_fase = np.angle(hilbert(s_fase))
    phi_amp = np.angle(hilbert(s_amp))

    diff = phi_amp - (n * phi_fase)
    return float(np.abs(np.mean(np.exp(1j * diff))))
