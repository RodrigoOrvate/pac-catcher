"""
test_synthetic_harmonico.py - Validacao sintetica do audita_harmonico.py

==========================================================================
IMPORTANTE: Decisao de teste (documentada para evitar confusao)
==========================================================================
Este teste CONTORNA o portao de qualidade (erro_ajuste < 0.15)
propositalmente para isolar a validacao do teste de RAZAO+PLV
do teste do proprio portao.

Razao: queremos saber se, DADO um cf de teta valido, a logica de
discriminacao harmonico vs genuino funciona. Se o portao rejeitasse
tudo, nao saberiamos se e' o portao que falha ou a logica que falha.

O portao de qualidade e' validado SEPARADAMENTE em test_fooof_quality.py
(pendente). Aqui focamos em: razao harmonica + PLV discriminam?

Quando o portao e' contornado, o veredito real do codigo de producao
ainda seria SEM_REFERENCIA_TETA. Os valores de PLV/ordem mostrados
neste teste sao "se o cf fosse aceito, discriminaria?".

==========================================================================
Cenarios testados (versao expandida)
==========================================================================
A) Harmonico em multiplas ordens (2x, 3x, 4x, 5x) com fase travada
B) Acoplamento genuino: gamma independente com envelope modulado
C) Quase-coincidencia: gamma independente em freq proxima a um multiplo
D) Sweep de SNR: ruido gaussiano crescente no cenario A (3x)
E) Fundo aperiodico realista: pink noise (1/f^n) + Gaussiana de teta
   (segue secao "Simulated data" de Kuhn et al. 2026)

Uso:
    python test_synthetic_harmonico.py
"""
import numpy as np
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "pipeline"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "pipeline/auditorias"))

from scipy.signal import butter, filtfilt

from audita_harmonico import (
    extrai_cf_teta_fooof,
    testa_razao_harmonica,
    compute_plv_harmonico,
)


def skewness_from_signal(data):
    """Skewness padronizado de Fisher-Pearson."""
    mean = np.mean(data)
    std = np.std(data)
    if std == 0:
        return 0.0, len(data)
    n = len(data)
    if n < 3:
        return np.nan, n
    skew = (n / ((n - 1) * (n - 2))) * np.sum(((data - mean) / std) ** 3)
    return skew, n


def butter_bandpass(sig, lo, hi, fs, order=4):
    """Filtro butterworth bandpass, zero-phase."""
    nyq = 0.5 * fs
    b, a = butter(order, [lo / nyq, hi / nyq], btype='band')
    return filtfilt(b, a, sig)


def skewness_from_band(data, fs, lo=4, hi=8, order=4):
    """Skewness da banda teta."""
    teta = butter_bandpass(data, lo, hi, fs, order)
    return skewness_from_signal(teta)


def generate_pink_noise(n_samples, slope=1.2, rng=None):
    """
    Gera pink noise (1/f^slope) via filtragem de white noise.
    slope=1.2 segue Kuhn et al. 2026 (Fig. 1).
    """
    if rng is None:
        rng = np.random.default_rng(42)
    white = rng.standard_normal(n_samples)
    # FFT-based pink noise
    fft = np.fft.rfft(white)
    freqs = np.fft.rfftfreq(n_samples)
    freqs[0] = 1.0  # evita divisao por zero
    fft = fft / (freqs ** (slope / 2))
    pink = np.fft.irfft(fft, n=n_samples)
    return pink / np.std(pink)  # normaliza


def generate_aperiodic_background(n_samples, fs, knee_freq=28.0, slope=1.2,
                                   offset=1.0, rng=None):
    """
    Gera componente aperiodico conforme Kuhn et al. 2026:
        L(f) = offset - slope * log10(f + knee_freq)
    Modelo 'knee' simplificado para gerar o PSD de fundo.
    """
    if rng is None:
        rng = np.random.default_rng(42)
    # White noise no dominio do tempo
    white = rng.standard_normal(n_samples)
    fft = np.fft.rfft(white)
    freqs = np.fft.rfftfreq(n_samples, d=1.0/fs)
    # Aplica perfil 1/f^knee
    freqs[0] = 1.0
    psd_profile = 1.0 / (1.0 + (freqs / knee_freq) ** slope)
    fft_shaped = fft * np.sqrt(psd_profile)
    signal = np.fft.irfft(fft_shaped, n=n_samples)
    return signal * offset


# FOOOF internals for synthetic test:
# Com nperseg=1.2s -> freq_res~0.83Hz. peak_width_limits=(2,5) exige 2-5 bins.
# Em sinal real (knee~28Hz, slope~1.2), o teta emerge do fundo. Em sinal
# sintetico sem estrutura biologica real, o ajuste falha com esses parametros.
# Para o teste sintetico, usamos nperseg=1.0s (freq_res=1Hz) e limites de
# largura mais flexiveis. ESTES SAO APENAS PARA O TESTE SINTETICO. O codigo
# de producao (audita_harmonico.py) usa os parametros calibrados para LFP real.
_FOOOF_TEST_KWARGS = {
    "aperiodic_mode": "knee",
    "peak_width_limits": (1.0, 4.0),   # mais flexivel (sintetico)
    "min_peak_height": 0.02,            # mais sensivel
    "peak_threshold": 1.0,
    "max_n_peaks": 1,
    "nperseg_s": 1.0,                  # freq_res=1Hz, ~9 bins em 4-12Hz
}


def generate_harmonic_signal(fs, dur_s, f_theta=8.0, ordem=3,
                              snr_db=20.0, aperiodic=True,
                              seed=42):
    """
    Cenario A: Teta nao-senoidal + harmonico de ordem n com fase rigida.
    snr_db: SNR do teta+gamma em relacao ao ruido (dB).
    aperiodic: se True, adiciona fundo 1/f realista.
    """
    rng = np.random.default_rng(seed)
    t = np.arange(0, dur_s, 1.0 / fs)
    n = len(t)

    # Teta sawtooth-like (soma de senoides)
    teta = (
        np.sin(2 * np.pi * f_theta * t)
        + 0.5 * np.sin(2 * np.pi * 2 * f_theta * t)
        + 0.3 * np.sin(2 * np.pi * 3 * f_theta * t)
    )

    # Gamma harmonico: fase travada em n*phi_theta
    f_gamma = ordem * f_theta
    phi_teta = 2 * np.pi * f_theta * t
    gamma = 0.5 * np.sin(ordem * phi_teta)

    # Sinal base
    sinal_base = teta + gamma

    # Ruido: gaussiano + (opcional) fundo aperiodico
    ruido_gauss = rng.standard_normal(n) * 0.05
    if aperiodic:
        # SNR: potencia(teta+gamma) / potencia(ruido_total) = 10^(SNR/10)
        pot_sinal = np.var(sinal_base)
        pot_ruido_desejada = pot_sinal / (10 ** (snr_db / 10))
        # Componente aperiodico
        aper = generate_aperiodic_background(n, fs, knee_freq=28.0,
                                            slope=1.2, offset=1.0, rng=rng)
        aper = aper / np.std(aper) * np.sqrt(pot_ruido_desejada * 0.7)
        ruido_gauss = ruido_gauss / np.std(ruido_gauss) * np.sqrt(pot_ruido_desejada * 0.3)
        ruido = aper + ruido_gauss
    else:
        # Ajusta ruido gaussiano para SNR desejado
        pot_sinal = np.var(sinal_base)
        pot_ruido_desejada = pot_sinal / (10 ** (snr_db / 10))
        ruido_gauss = ruido_gauss / np.std(ruido_gauss) * np.sqrt(pot_ruido_desejada)
        ruido = ruido_gauss

    sinal = sinal_base + ruido
    return sinal, fs, f_theta, f_gamma


def generate_genuine_coupling(fs, dur_s, f_theta=8.0, f_gamma=35.0,
                                snr_db=20.0, aperiodic=True, seed=42):
    """
    Cenario B: Teta senoidal + Gamma independente com envelope modulado.
    Gera MI alto mas SEM rigidez de fase n:1.
    """
    rng = np.random.default_rng(seed)
    t = np.arange(0, dur_s, 1.0 / fs)
    n = len(t)

    # Teta senoidal puro
    teta = np.sin(2 * np.pi * f_theta * t)

    # Gamma 35Hz: ruido filtrado em banda (NAO e harmonico de 8Hz)
    ruido_g = rng.standard_normal(n)
    b, a = butter(4, [(f_gamma - 5) / 500, (f_gamma + 5) / 500], btype="band")
    gamma = 0.5 * filtfilt(b, a, ruido_g)

    # Modulacao de envelope (gera MI sem PLV)
    envelope = 0.5 + 0.5 * np.sin(2 * np.pi * f_theta * t)
    gamma_mod = gamma * envelope

    sinal_base = teta + gamma_mod
    pot_sinal = np.var(sinal_base)
    pot_ruido_desejada = pot_sinal / (10 ** (snr_db / 10))

    if aperiodic:
        aper = generate_aperiodic_background(n, fs, knee_freq=28.0,
                                            slope=1.2, offset=1.0, rng=rng)
        aper = aper / np.std(aper) * np.sqrt(pot_ruido_desejada * 0.7)
        ruido_g = ruido_g / np.std(ruido_g) * np.sqrt(pot_ruido_desejada * 0.3)
        ruido = aper + ruido_g
    else:
        ruido = ruido_g / np.std(ruido_g) * np.sqrt(pot_ruido_desejada)

    sinal = sinal_base + ruido
    return sinal, fs, f_theta, f_gamma


def generate_near_coincidence(fs, dur_s, f_theta=8.0, f_gamma=25.0,
                                 snr_db=20.0, aperiodic=True, seed=42):
    """
    Cenario C: Gamma independente em freq PROXIMA a um multiplo de teta.
    25Hz esta' proximo de 24Hz (3x f_theta) - dentro de tolerancia 10%.
    Razao bate, mas PLV deve ser baixo (nao ha travamento de fase real).
    Testa se o PLV esta' fazendo o trabalho pesado de discriminacao.
    """
    rng = np.random.default_rng(seed)
    t = np.arange(0, dur_s, 1.0 / fs)
    n = len(t)

    # Teta senoidal
    teta = np.sin(2 * np.pi * f_theta * t)

    # Gamma 25Hz: ruido filtrado, independente (NAO travado em 3*phi_theta)
    ruido_g = rng.standard_normal(n)
    b, a = butter(4, [(f_gamma - 5) / 500, (f_gamma + 5) / 500], btype="band")
    gamma = 0.5 * filtfilt(b, a, ruido_g)

    # Sem modulacao de envelope (mais estrito que cenario B)
    sinal_base = teta + gamma
    pot_sinal = np.var(sinal_base)
    pot_ruido_desejada = pot_sinal / (10 ** (snr_db / 10))

    if aperiodic:
        aper = generate_aperiodic_background(n, fs, knee_freq=28.0,
                                            slope=1.2, offset=1.0, rng=rng)
        aper = aper / np.std(aper) * np.sqrt(pot_ruido_desejada * 0.7)
        ruido_g = ruido_g / np.std(ruido_g) * np.sqrt(pot_ruido_desejada * 0.3)
        ruido = aper + ruido_g
    else:
        ruido = ruido_g / np.std(ruido_g) * np.sqrt(pot_ruido_desejada)

    sinal = sinal_base + ruido
    return sinal, fs, f_theta, f_gamma


def run_test(label, sinal, fs, f_theta_ref, f_gamma_ref, ini, fim,
             verbose=True):
    """
    Roda o teste de razao+PLV em um sinal sintetico.
    NAO chama o portao de qualidade (decisao documentada no topo do arquivo).
    """
    idx_ini = int(ini * fs)
    idx_fim = int(fim * fs)
    sinal_cand = sinal[idx_ini:idx_fim]

    ctx_dur = 45.0
    centro = (ini + fim) / 2
    ctx_ini = max(0, centro - ctx_dur / 2)
    ctx_fim = min(len(sinal) / fs, centro + ctx_dur / 2)
    idx_ctx_ini = int(ctx_ini * fs)
    idx_ctx_fim = int(ctx_fim * fs)
    sinal_ctx = sinal[idx_ctx_ini:idx_ctx_fim]

    res_fooof = extrai_cf_teta_fooof(
        sinal_ctx, fs,
        aperiodic_mode=_FOOOF_TEST_KWARGS["aperiodic_mode"],
        theta_bw_limits=_FOOOF_TEST_KWARGS["peak_width_limits"],
        min_peak_height=_FOOOF_TEST_KWARGS["min_peak_height"],
        nperseg_s=_FOOOF_TEST_KWARGS["nperseg_s"],
    )
    skew, n = skewness_from_band(sinal_cand, fs, lo=4, hi=8)

    cf_usar = res_fooof["cf_teta"]
    plv_val = np.nan
    ordem = None
    suspeito = False
    desvio = None

    if cf_usar is not None:
        suspeito, ordem, desvio = testa_razao_harmonica(
            cf_usar, f_gamma_ref, tol_rel=0.10)
        if suspeito:
            plv_val = compute_plv_harmonico(
                sinal_cand, fs, cf_usar, f_gamma_ref, ordem, ini, fim)

    if verbose:
        erro_str = f"{res_fooof['erro_ajuste']:.4f}" if res_fooof['erro_ajuste'] is not None else "N/A"
        cf_str = f"{cf_usar:.2f}" if cf_usar else "N/A"
        plv_str = f"{plv_val:.3f}" if not np.isnan(plv_val) else "N/A"
        print(f"    cf_fooof={cf_str} (ref={f_theta_ref}Hz) | erro={erro_str} | "
              f"skew={skew:+.3f} | razao: suspeito={suspeito}, ordem={ordem} | "
              f"PLV={plv_str}")

    return {
        "cf_teta": cf_usar,
        "erro_ajuste": res_fooof["erro_ajuste"],
        "qualidade_ok": res_fooof["qualidade_ok"],
        "skew": skew,
        "suspeito_razao": suspeito,
        "ordem": ordem,
        "desvio": desvio,
        "plv": plv_val,
    }


def main():
    print("=" * 85)
    print("TESTE SINTETICO EXPANDIDO - audita_harmonico.py")
    print("=" * 85)
    print("NOTA: Portao de qualidade (erro_ajuste < 0.15) contornado neste teste")
    print("      para isolar validacao de razao+PLV. Ver docstring do arquivo.")
    print("=" * 85)

    fs = 1000.0
    dur_s = 60.0
    ini, fim = 20.0, 30.0
    resultados = []

    # =========================================================================
    # PARTE 1: Cenarios basicos COM fundo aperiodico realista (knee freq 28Hz)
    # =========================================================================
    print("\n>>> PARTE 1: Cenarios com fundo 1/f realista (knee=28Hz, slope=1.2)")
    print("=" * 85)

    cenarios_base = [
        ("A2_2x", "Harmonico 2x travado", generate_harmonic_signal, {"ordem": 2, "f_theta": 8.0}, 16.0),
        ("A3_3x", "Harmonico 3x travado", generate_harmonic_signal, {"ordem": 3, "f_theta": 8.0}, 24.0),
        ("A4_4x", "Harmonico 4x travado", generate_harmonic_signal, {"ordem": 4, "f_theta": 8.0}, 32.0),
        ("A5_5x", "Harmonico 5x travado", generate_harmonic_signal, {"ordem": 5, "f_theta": 8.0}, 40.0),
        ("B", "Acoplamento genuino (35Hz)", generate_genuine_coupling, {"f_gamma": 35.0}, 35.0),
        ("C", "Quase-coincidencia (25Hz ~ 3x8Hz)", generate_near_coincidence, {"f_gamma": 25.0}, 25.0),
    ]

    for cid, label, gen_fn, kwargs, f_g_ref in cenarios_base:
        print(f"\n[{cid}] {label} (f_gamma={f_g_ref}Hz)")
        # Cada cenario com seed proprio (RNG isolado)
        # Alem disso: rodar tambem sem fundo aperiodico como comparacao
        for aper_label, aper_flag in [("aperiodic_on", True), ("aperiodic_off", False)]:
            # Seed unico por cenario E por flag
            seed_cenario = hash((cid, aper_label)) % 100000
            sinal, fs_r, f_t, f_g = gen_fn(fs, dur_s, aperiodic=aper_flag,
                                            seed=seed_cenario, **kwargs)
            sub_label = f"{label} [{aper_label}]"
            print(f"  -- {aper_label} (seed={seed_cenario})")
            res = run_test(sub_label, sinal, fs_r, f_t, f_g, ini, fim)
            resultados.append({"id": f"{cid}_{aper_label}", "label": sub_label,
                                "f_gamma": f_g_ref, **res})

    # =========================================================================
    # PARTE 2: Sweep de SNR (cenario A3, harmonico 3x)
    # =========================================================================
    print("\n" + "=" * 85)
    print(">>> PARTE 2: Sweep de SNR (harmonico 3x)")
    print("=" * 85)

    for snr in [30, 20, 10, 5, 0]:
        print(f"\n[A3_SNR{snr}] Harmonico 3x, SNR={snr}dB")
        sinal, fs_r, f_t, f_g = generate_harmonic_signal(
            fs, dur_s, ordem=3, f_theta=8.0, snr_db=snr, aperiodic=True)
        res = run_test(f"SNR={snr}dB", sinal, fs_r, f_t, f_g, ini, fim)
        resultados.append({"id": f"A3_SNR{snr}", "label": f"3x @ SNR={snr}dB",
                            "f_gamma": 24.0, **res})

    # =========================================================================
    # PARTE 3: Comparacao knee vs fixed no cenario 3x
    # =========================================================================
    print("\n" + "=" * 85)
    print(">>> PARTE 3: aperiodic_mode knee vs fixed (harmonico 3x)")
    print("=" * 85)

    for mode in ['knee', 'fixed']:
        print(f"\n[A3_{mode}] Harmonico 3x, aperiodic_mode={mode}")
        sinal, fs_r, f_t, f_g = generate_harmonic_signal(
            fs, dur_s, ordem=3, f_theta=8.0, snr_db=20, aperiodic=True, seed=42)
        idx_ini = int(ini * fs_r)
        idx_fim = int(fim * fs_r)
        centro = (ini + fim) / 2
        sinal_ctx = sinal[int((centro - 22.5) * fs_r):int((centro + 22.5) * fs_r)]
        res_fooof = extrai_cf_teta_fooof(
            sinal_ctx, fs_r,
            aperiodic_mode=mode,
            theta_bw_limits=_FOOOF_TEST_KWARGS["peak_width_limits"],
            min_peak_height=_FOOOF_TEST_KWARGS["min_peak_height"],
            nperseg_s=_FOOOF_TEST_KWARGS["nperseg_s"],
        )
        erro_str = f"{res_fooof['erro_ajuste']:.4f}" if res_fooof['erro_ajuste'] is not None else "N/A"
        print(f"    aperiodic_mode={mode}: cf={res_fooof['cf_teta']}, "
              f"erro={erro_str}, qualidade_ok={res_fooof['qualidade_ok']}")

    # =========================================================================
    # RESUMO
    # =========================================================================
    print("\n" + "=" * 85)
    print("RESUMO")
    print("=" * 85)
    print(f"{'ID':<10} | {'Label':<32} | {'cf_teta':<8} | {'Erro':<7} | {'Q':<2} | {'Ordem':<5} | {'PLV':<6} | {'Skew':<7}")
    print("-" * 100)
    for r in resultados:
        cf_str = f"{r['cf_teta']:.2f}" if r['cf_teta'] else "N/A"
        erro_str = f"{r['erro_ajuste']:.3f}" if r['erro_ajuste'] is not None else "N/A"
        q_str = "OK" if r['qualidade_ok'] else "X"
        ordem_str = str(r['ordem']) if r['ordem'] else "-"
        plv_str = f"{r['plv']:.3f}" if not np.isnan(r['plv']) else "-"
        skew_str = f"{r['skew']:+.2f}"
        print(f"{r['id']:<10} | {r['label']:<32} | {cf_str:<8} | {erro_str:<7} | "
              f"{q_str:<2} | {ordem_str:<5} | {plv_str:<6} | {skew_str:<7}")

    # =========================================================================
    # AFERICAO
    # =========================================================================
    print("\n" + "=" * 85)
    print("AFERICAO (criterios: cf~8Hz, razao correta, PLV alto para harmonico, baixo para genuino)")
    print("=" * 85)

    # Harmonicos 2x-5x: PLV deve ser > 0.7
    for cid_prefix in ['A2_2x', 'A3_3x', 'A4_4x', 'A5_5x']:
        r = next((x for x in resultados if x['id'] == cid_prefix), None)
        if r and not np.isnan(r['plv']):
            plv_ok = r['plv'] > 0.7
            cf_ok = r['cf_teta'] is not None and abs(r['cf_teta'] - 8.0) < 1.0
            status = "OK" if (plv_ok and cf_ok) else "FALHOU"
            print(f"  {cid_prefix}: PLV={r['plv']:.3f} ({'alto' if plv_ok else 'BAIXO'}), "
                  f"cf={r['cf_teta']:.2f} ({'prox' if cf_ok else 'LONGE'}) -> {status}")

    # Genuino B: PLV deve ser baixo (< 0.5) ou razao nao suspeita
    r = next((x for x in resultados if x['id'] == 'B'), None)
    if r:
        plv_baixo = np.isnan(r['plv']) or r['plv'] < 0.5
        razao_nao = not r['suspeito_razao']
        status = "OK" if (plv_baixo and razao_nao) else "FALHOU"
        plv_display = "N/A" if np.isnan(r['plv']) else f"{r['plv']:.3f}"
        plv_str = f"{'N/A' if np.isnan(r['plv']) else f'{r['plv']:.3f}'}"
        print(f"  B (genuino): PLV={plv_display} "
              f"({'baixo' if plv_baixo else 'ALTO'}), razao suspeita={r['suspeito_razao']} -> {status}")

    # Quase-coincidencia C: razao suspeita MAS PLV baixo
    r = next((x for x in resultados if x['id'] == 'C'), None)
    if r:
        plv_baixo = np.isnan(r['plv']) or r['plv'] < 0.5
        # Idealmente: razao suspeita (25Hz ~ 3x8=24, dentro 10%?) + PLV baixo
        if r['suspeito_razao']:
            status = "OK" if plv_baixo else "FALSO POSITIVO"
        else:
            status = "OK (razao nao suspeita, discriminado por freq)"
        plv_display_c = "N/A" if np.isnan(r['plv']) else f"{r['plv']:.3f}"
        print(f"  C (quase-coincidencia): razao suspeita={r['suspeito_razao']}, "
              f"PLV={plv_display_c} -> {status}")

    # SNR sweep
    print("\n  --- SNR sweep (3x) ---")
    for snr in [30, 20, 10, 5, 0]:
        r = next((x for x in resultados if x['id'] == f'A3_SNR{snr}'), None)
        if r:
            plv_val = r['plv'] if not np.isnan(r['plv']) else 0
            cf_display = f"{r['cf_teta']:.2f}" if r['cf_teta'] else "N/A"
            print(f"  SNR={snr:2d}dB: PLV={plv_val:.3f}, cf={cf_display}")

    print("\n" + "=" * 85)
    print("NOTAS FINAIS")
    print("=" * 85)
    print("1. Portao de qualidade (erro_ajuste<0.15) NAO foi aplicado neste teste.")
    print("2. FOOOF rodou com aperiodic_mode='knee' (default para LFP real).")
    print("3. Fundo 1/f realista (knee=28Hz, slope=1.2) incluido em todos os cenarios.")
    print("4. Validacao REAL exige sessoes MTESC04/05 (caos biologico real).")


if __name__ == "__main__":
    main()
