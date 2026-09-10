"""
deteccao_ripple.py
==========================================
Módulo dedicado à detecção hierárquica de eventos transitórios (Ripples ⊆ HFO).

O pipeline espectral (Welch em janelas de 10s) não tem resolução temporal
para isolar eventos de alta frequência que duram apenas ~25-50 ms. Usar a
energia média da janela inteira na banda HFO como detector de Ripple resulta
num critério excessivamente permissivo ("paradoxo HFO constante").

Este módulo opera em resolução de amostra, calculando o envelope analítico
e aplicando limites de duração mínima e coocorrência (sharp-wave).
"""

import numpy as np
from scipy.signal import butter, filtfilt, hilbert
from pipeline.etapa1_triagem.preprocessa_referencia_diferencial import aplica_referencia_diferencial

def _filtra_banda(sinal, fs, lowcut, highcut, order=4):
    nyq = 0.5 * fs
    lo = max(lowcut / nyq, 0.001)
    hi = min(highcut / nyq, 0.999)
    if lo >= hi:
        return np.zeros_like(sinal)
    b, a = butter(order, [lo, hi], btype="bandpass")
    return filtfilt(b, a, sinal)

def detecta_eventos_ripple(sinal_lfp, fs, banda_ripple=(150, 250),
                            limiar_dp=3.0, duracao_min_ms=15.0, duracao_max_ms=200.0,
                            banda_sharp_wave=(1, 30), exigir_sharp_wave=False, sinal_referencia=None):
    """
    Detector hierárquico: um evento só conta como "ripple" (não apenas "HFO cru")
    se satisfizer TODOS os critérios abaixo simultaneamente (AND):

      1. Envelope (Hilbert) na banda_ripple ultrapassa limiar_dp desvios-padrão
         acima da mediana do envelope de toda a sessão/janela de referência.
      2. A ultrapassagem dura pelo menos duracao_min_ms de forma contígua.
      3. [Opcional, se exigir_sharp_wave=True] Há uma deflexão correspondente
         na banda_sharp_wave no mesmo intervalo de tempo.

    Retorna:
        - lista de eventos (inicio_s, fim_s, pico_amplitude, tem_sharp_wave)
        - série temporal booleana "hfo_cru" (somente critério 1)
        
    TODO: calibrar limiar_dp e duracao_min_ms empiricamente em sessão real,
    comparando visualmente eventos detectados contra o traço bruto.
    """
    # Aplica referencia diferencial se fornecida
    sinal = aplica_referencia_diferencial(sinal_lfp, sinal_referencia)
    sinal_np = np.asarray(sinal, dtype=np.float64)
    
    # 1. Filtra HFO/Ripple e obtém o envelope
    sinal_ripple = _filtra_banda(sinal_np, fs, *banda_ripple)
    envelope_ripple = np.abs(hilbert(sinal_ripple))
    
    # Referência robusta (mediana + MAD estimado ou desvio padrão clássico)
    # Como as caudas do HFO podem ser pesadas, mediana é mais segura que média.
    mediana = np.median(envelope_ripple)
    desvio = np.std(envelope_ripple)  # Ou usar 1.4826 * median(abs(x - median))
    if desvio == 0:
        desvio = 1e-12
        
    limiar_abs = mediana + (limiar_dp * desvio)
    
    # Critério 1: hfo_cru (simples cruzamento de limiar)
    hfo_cru_mask = envelope_ripple >= limiar_abs
    
    # Critério 2: duração mínima
    n_amostras_min = int(duracao_min_ms * 1e-3 * fs)
    
    eventos = []
    
    # Encontra os limites contíguos de True em hfo_cru_mask
    diffs = np.diff(hfo_cru_mask.astype(int))
    starts = np.where(diffs == 1)[0] + 1
    ends = np.where(diffs == -1)[0] + 1
    
    if hfo_cru_mask[0]:
        starts = np.insert(starts, 0, 0)
    if hfo_cru_mask[-1]:
        ends = np.append(ends, len(hfo_cru_mask))
        
    # Se formos exigir sharp wave
    if exigir_sharp_wave:
        sinal_sw = _filtra_banda(sinal_np, fs, *banda_sharp_wave)
        envelope_sw = np.abs(hilbert(sinal_sw))
        limiar_sw = np.median(envelope_sw) + (2.0 * np.std(envelope_sw)) # limiar frouxo para SW
    else:
        sinal_sw = None
        limiar_sw = None

    for ini, fim in zip(starts, ends):
        duracao = fim - ini
        if duracao >= n_amostras_min:
            pico_amp = np.max(envelope_ripple[ini:fim])
            
            tem_sharp_wave = False
            if exigir_sharp_wave and sinal_sw is not None:
                # Procura por aumento de energia na banda baixa no mesmo intervalo
                tem_sharp_wave = np.any(envelope_sw[ini:fim] > limiar_sw)
            else:
                tem_sharp_wave = True # Se não for exigido, consideramos satisfeito
                
            if (not exigir_sharp_wave) or tem_sharp_wave:
                eventos.append({
                    "inicio_s": ini / fs,
                    "fim_s": fim / fs,
                    "pico_amplitude": float(pico_amp),
                    "tem_sharp_wave": tem_sharp_wave
                })
                
    return eventos, hfo_cru_mask

def resume_eventos_por_janela(eventos, hfo_cru_mask, fs, janela_ini_s, janela_fim_s):
    """
    Dado o output de detecta_eventos_ripple e a máscara contínua, verifica
    se a janela [janela_ini_s, janela_fim_s] contém pelo menos um ripple
    totalmente ou parcialmente dentro dela, e se contém hfo_cru.
    
    Retorna (tem_ripple: bool, tem_hfo_cru: bool)
    """
    ini_amostra = int(janela_ini_s * fs)
    fim_amostra = int(janela_fim_s * fs)
    fim_amostra = min(fim_amostra, len(hfo_cru_mask))
    
    if ini_amostra >= len(hfo_cru_mask):
        return False, False
        
    tem_hfo_cru = bool(np.any(hfo_cru_mask[ini_amostra:fim_amostra]))
    
    tem_ripple = False
    for ev in eventos:
        # Se o evento cruza ou cai dentro da janela
        if (ev["inicio_s"] <= janela_fim_s) and (ev["fim_s"] >= janela_ini_s):
            tem_ripple = True
            break
            
    return tem_ripple, tem_hfo_cru
