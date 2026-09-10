"""
pac_core/filtering.py
==========================================
Filtros compartilhados do pipeline PAC: passa-faixa Butterworth e notch
multi-harmônico. Extraído das 4 cópias idênticas de `filtra_sinal` e das 2
cópias idênticas de `aplica_notch` que existiam em `comodulogram.py`,
`triagem_pac.py`, `refina_candidatos.py` e `auditorias/diagnostico_janela.py`.

Variantes de filtro com `order`/clamp diferentes (ex.: `deteccao_ripple.py`,
`auditorias/audita_held_out.py`, `auditorias/utils_harmonico.py`) NÃO foram
fundidas aqui — permanecem nos seus arquivos originais até uma decisão
científica explícita sobre unificá-las.
"""

import numpy as np
import scipy.signal as signal


def filtra_sinal(sinal_in, lowcut, highcut, fs, order=3):
    nyq = 0.5 * fs
    low = max(lowcut / nyq, 1e-6)
    high = min(highcut / nyq, 0.999)
    b, a = signal.butter(order, [low, high], btype="bandpass")
    return signal.filtfilt(b, a, sinal_in)


def aplica_notch(sinal_in, fs, freqs_notch=None, q_factor=30.0, linha_hz=None):
    """
    Rejeita a frequência da rede elétrica e seus harmônicos (iirnotch).
    Aceita uma lista de frequências (ex: [60, 120, 180, 240]).
    Necessário porque picos de linha dentro da banda de amplitude inflam
    o MI simulando acoplamento onde não há.

    `linha_hz` é um alias legado de `freqs_notch` (mantido por compatibilidade
    com call sites que ainda chamam pelo nome antigo).
    """
    if freqs_notch is None:
        freqs_notch = linha_hz
    if not freqs_notch:
        return sinal_in

    nyq = 0.5 * fs
    out = sinal_in

    # Suporta tanto um número (comportamento antigo, mas aqui usamos lista) quanto lista
    if not isinstance(freqs_notch, (list, tuple, np.ndarray)):
        freqs_notch = [freqs_notch]

    for f in freqs_notch:
        if f >= nyq * 0.98:
            continue
        b, a = signal.iirnotch(f / nyq, q_factor)
        out = signal.filtfilt(b, a, out)
    return out
