"""
linha_noise_kuhn.py - Limpeza de ruido de linha (60Hz e harmonicos) estilo
                      Kuhn et al. 2026 (LFP_FOOOF) para o FOOOF de producao.

Duas primitivas de remocao (escolha via `aplica_modo`):

  `remove_pico_gaussiana` - subtrai APENAS a(s) Gaussiana(s) do pico em log10,
    preservando o componente aperiodico local (nao cria buraco).
    Validada em test_fooof_linha_preprocess.py (erro cen. D: 0.265 -> 0.165).

  `remove_faixa_1f`       - port fiel do `rem_noise.m` do artigo: ajusta FOOOF
    'fixed' numa janela estreita e REPOE a banda estrita (±2 Hz) pela 1/f pura
    reconstruida dos parametros aperiodicos. Cirurgica (so toca a linha), mas
    so faz sentido na frequencia de linha (ou em seus harmonicos, via driver).

O driver `preprocessa_linha` aplica qualquer primitiva iterativamente em
`f_linha, 2*f_linha, ..., max_harmonicos*f_linha` (encadeando o PSD corrigido).
`aplica_modo` traduz os tres modos de producao concretos:
  - 'gaussiana' : remove_pico_gaussiana em 60/120/180 Hz
  - 'cirurgica' : remove_faixa_1f somente em 60 Hz (replica rem_noise.m)
  - 'hibrido'   : remove_faixa_1f em 60/120/180 Hz (cirurgica + rede inteira)

IMPORTANTE - ALCANCE: estas funcoes devem ser usadas SOMENTE para ruido de
linha (60+120+180 Hz no Brasil; 50+100+150 na Europa). NAO aplicar para
remover harmonicos de teta - seria circular no contexto do audita_harmonico.py,
que existe para DETECTAR se um pico em gamma e harmonico de teta. O artigo Kuhn
remove harmonicos de teta porque quer limpar o espectro; nosso script quer
ENCONTRAR relacoes harmonicas, nao apaga-las.

IMPORTANTE - ESCALA (remove_pico_gaussiana): FOOOF ajusta internamente em
log10. Subtrair em escala linear destrói o espectro (bug documentado
04/09/2026). Operamos em log10.
"""
import numpy as np
from scipy.ndimage import uniform_filter1d
from fooof import FOOOF


def remove_pico_gaussiana(freqs, psd, freq_centro, largura=7.0,
                          aperiodic_mode='fixed', verbose=False):
    """Subtrai APENAS a Gaussiana (pico) do PSD, em log10, em torno de freq_centro.

    Procedimento (Methods, 'Preprocessing of neural data' de Kuhn et al. 2026):
      1. Recorta janela estreita ao redor de freq_centro (±largura).
      2. Ajusta FOOOF 1exp nessa janela (dá aperiódico local + Gaussiana).
      3. Subtrai APENAS a(s) Gaussiana(s) ajustada(s) do PSD, em log10.
      4. Preserva o componente aperiódico local (não cria buraco).

    Retorna (psd_corrigido, ok) — ok=False se nada foi removido (sem pico
    detectavel / fit falhou).
    """
    lo, hi = freq_centro - largura, freq_centro + largura
    mask = (freqs >= lo) & (freqs <= hi)
    if mask.sum() < 5:
        return psd, False

    fm_local = FOOOF(aperiodic_mode=aperiodic_mode,
                     max_n_peaks=1, peak_width_limits=(1, 6),
                     min_peak_height=0.05, peak_threshold=0.5,
                     verbose=False)
    try:
        fm_local.fit(freqs[mask], psd[mask], freq_range=(lo, hi))
    except Exception:
        return psd, False
    if not fm_local.has_model or fm_local.n_peaks_ == 0:
        return psd, False

    # _peak_fit: componente de picos em log10 (shape (n_pts_janela,)).
    if fm_local._peak_fit is None:
        return psd, False

    psd_log = np.log10(psd[mask])
    psd_log_clean = psd_log - fm_local._peak_fit
    psd_corrigido = psd.copy()
    psd_corrigido[mask] = np.maximum(10 ** psd_log_clean, 1e-30)
    if verbose:
        print(f"    Kuhn: Gaussiana removida em {freq_centro:.0f}Hz")
    return psd_corrigido, True


def remove_faixa_1f(freqs, psd, freq_centro, largura_rep=2.0,
                    janela_fit=(-4.0, 5.0), aperiodic_mode='fixed',
                    verbose=False):
    """Port fiel do `rem_noise.m` do artigo Kuhn et al. 2026.

    Ajusta FOOOF 'fixed' na janela `freq_centro + janela_fit` (default
    [f-4, f+5], como o MATLAB usa [46,55] para 50Hz) com settings
    permissivas (min_peak_height=0, peak_threshold=0, largura de pico
    0.5-2 Hz — so a linha, nao gaussiana larga de gamma). Se achar pico,
    REPOE a banda `freq_centro ± largura_rep` (default ±2 Hz, como o
    MATLAB usa [48,52] para 50Hz) pela 1/f pura reconstruida dos
    parametros aperiodicos:  psd = 10**offset * f**(-exponente).

    Retorna (psd_corrigido, ok). So toca a banda estrita; nada fora dela
    e alterado (preserva o resto do espectro intacto).
    """
    lo_fit, hi_fit = freq_centro + janela_fit[0], freq_centro + janela_fit[1]
    mask_fit = (freqs >= lo_fit) & (freqs <= hi_fit)
    if mask_fit.sum() < 5:
        return psd, False

    # Settings identicas ao rem_noise.m (max_n_peaks=1, permissivo).
    fm_local = FOOOF(aperiodic_mode=aperiodic_mode,
                     max_n_peaks=1, peak_width_limits=(0.5, 2),
                     min_peak_height=0.0, peak_threshold=0.0,
                     verbose=False)

    # Estabilidade numerica: senoides puras (sintetico) produzem spike de
    # ~140 dB de faixa dinamica num bin. Numa janela de 9 Hz, o ajuste do
    # modelo aperiodico 'fixed' estoura em xs**exp (ydata vira NaN).
    # Suavizamos o PSD em log10 (5 bins) SO no fit local — nao altera o
    # PSD retornado, so da ao FOOOF numeros fitaveis.
    psd_log_fit = np.log10(psd[mask_fit])
    psd_fit_suave = 10 ** uniform_filter1d(psd_log_fit, size=5)
    try:
        fm_local.fit(freqs[mask_fit], psd_fit_suave, freq_range=(lo_fit, hi_fit))
    except Exception:
        return psd, False
    if not fm_local.has_model or fm_local.n_peaks_ == 0:
        return psd, False

    # Reconstitui a banda apenas pela componente aperiodica local.
    lo_rep, hi_rep = freq_centro - largura_rep, freq_centro + largura_rep
    mask_rep = (freqs >= lo_rep) & (freqs <= hi_rep)
    if mask_rep.sum() == 0:
        return psd, False

    offset, expoente = fm_local.aperiodic_params_[0], fm_local.aperiodic_params_[1]
    psd_corrigido = psd.copy()
    psd_corrigido[mask_rep] = (10.0 ** offset) * np.power(freqs[mask_rep], -expoente)
    if verbose:
        print(f"    Kuhn: banda [{lo_rep:.0f},{hi_rep:.0f}]Hz reposta pela 1/f")
    return psd_corrigido, True


def preprocessa_linha(freqs, psd, fs, func, f_linha=60.0, max_harmonicos=3,
                      verbose=False):
    """Aplica `func` iterativamente em `f_linha` e seus harmonicos.

    fs: taxa de amostragem (para nao ultrapassar Nyquist). f_linha=60Hz
    (rede BR) remove 60/120/180; f_linha=50 (EU) remove 50/100/150.
    O PSD corrigido de cada iteracao alimenta a proxima.
    """
    psd_out = psd.copy()
    nyq = fs / 2.0
    for h in range(1, max_harmonicos + 1):
        fc = f_linha * h
        if fc >= nyq:
            break
        psd_out, ok = func(freqs, psd_out, fc, verbose=verbose)
    return psd_out


MODOS = {
    'gaussiana': (remove_pico_gaussiana, 4),  # subtrai gaussiana em 60/120/180/240
    'cirurgica': (remove_faixa_1f, 1),        # repõe banda so em 60 (rem_noise.m)
    'hibrido': (remove_faixa_1f, 4),          # repõe banda em 60/120/180/240
}


def aplica_modo(freqs, psd, fs, modo, f_linha=60.0, verbose=False):
    """Aplica uma das tres configs de producao a um PSD ja calculado (Welch).

    modo: 'gaussiana' | 'cirurgica' | 'hibrido' (ver MODOS).
    Retorna o PSD limpo (mesmo shape de entrada).
    """
    if modo not in MODOS:
        raise ValueError(f"modo invalido: {modo!r}. Validos: {list(MODOS)}")
    func, n_harm = MODOS[modo]
    return preprocessa_linha(freqs, psd, fs, func,
                             f_linha=f_linha, max_harmonicos=n_harm,
                             verbose=verbose)