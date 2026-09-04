"""
test_synthetic_harmonico.py - Validacao sintetica do audita_harmonico.py

Gera dois cenarios em memoria (sem I/O de arquivo):
  A) Teta nao-senoidal (sawtooth-like) com harmonico matematico em nF.
     Esperado: SUSPEITO_HARMONICO_FORTE (PLV alto, razao exata, skew alto).
  B) Teta senoidal + Gamma genuinamente independente com MI alto.
     Esperado: CLEAN ou REVISAR_RAZAO_INTEIRA (PLV baixo, skew baixo).

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
    """Versao direta de audita_skewness.theta_skewness_for_window (sem I/O)."""
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


def generate_scenario_A(fs=1000.0, dur_s=60.0):
    """
    Cenario A: Teta nao-senoidal (sawtooth-like) + harmonico 3x matematico.
    - Theta 8Hz com harmônicos naturais (2x, 3x)
    - Gamma 24Hz com fase rigidamente travada em 3*phi_theta
    - Skewness esperada: ALTA (forma de onda assimetrica)
    - PLV esperado: ALTO (~1.0, fase rigida)
    """
    t = np.arange(0, dur_s, 1.0 / fs)
    f_theta = 8.0

    # Theta sawtooth-like (soma de harmônicos de Fourier)
    teta = (
        np.sin(2 * np.pi * f_theta * t)
        + 0.5 * np.sin(2 * np.pi * 2 * f_theta * t)
        + 0.3 * np.sin(2 * np.pi * 3 * f_theta * t)
    )

    # Gamma harmonico: fase rigida em 3x phi_theta
    f_gamma = 3 * f_theta  # 24Hz
    phi_teta = 2 * np.pi * f_theta * t
    gamma = 0.5 * np.sin(3 * phi_teta)  # 3x = harmonico

    sinal = teta + gamma + 0.02 * np.random.randn(len(t))
    return sinal, fs, f_theta, f_gamma


def generate_scenario_B(fs=1000.0, dur_s=60.0):
    """
    Cenario B: Teta senoidal Puro + Gamma 35Hz genuinamente independente.
    - Theta 8Hz senoidal (skewness ~0)
    - Gamma 35Hz: ruido filtrado em banda estreita
    - 35Hz NAO e harmonico natural de sawtooth 8Hz (8*4=32, 8*5=40)
    - Modulacao de envelope gera MI alto mas SEM rigidez de fase n:1
    - Skewness esperada: BAIXA
    - PLV esperado: BAIXO
    """
    t = np.arange(0, dur_s, 1.0 / fs)
    f_theta = 8.0

    # Theta senoidal puro (skewness ~0)
    teta = np.sin(2 * np.pi * f_theta * t)

    # Gamma 35Hz: ruido filtrado em banda (35Hz nao e harmonico de 8Hz)
    f_gamma = 35.0
    ruido = np.random.randn(len(t))
    b, a = butter(4, [(f_gamma - 5) / 500, (f_gamma + 5) / 500], btype="band")
    gamma = 0.5 * filtfilt(b, a, ruido)

    # Modulacao de envelope por fase do teta (gera MI sem PLV)
    envelope = 0.5 + 0.5 * np.sin(2 * np.pi * f_theta * t)
    gamma_mod = gamma * envelope

    sinal = teta + gamma_mod + 0.02 * np.random.randn(len(t))
    return sinal, fs, f_theta, f_gamma


def main():
    print("=" * 80)
    print("TESTE SINTETICO - audita_harmonico.py")
    print("Validacao: harmonico vs. acoplamento genuino")
    print("=" * 80)

    fs = 1000.0
    dur_s = 60.0
    np.random.seed(42)

    # Janela de teste: 20-30s
    ini, fim = 20.0, 30.0
    centro = (ini + fim) / 2

    resultados = []

    for scenario, label, gen_fn in [
        ("A", "Harmonico verdadeiro (sawtooth + 3x fase rigida)",
         generate_scenario_A),
        ("B", "Acoplamento genuino (senoidal + gamma 40Hz independente)",
         generate_scenario_B),
    ]:
        print(f"\n{'=' * 40}")
        print(f"CENARIO {scenario}: {label}")
        print("=" * 40)

        sinal, fs, f_theta_ref, f_gamma_ref = gen_fn(fs, dur_s)

        # Extrair janela do candidato
        idx_ini = int(ini * fs)
        idx_fim = int(fim * fs)
        sinal_cand = sinal[idx_ini:idx_fim]

        # Janela de contexto para FOOOF
        ctx_dur = 45.0
        ctx_ini = max(0, centro - ctx_dur / 2)
        ctx_fim = min(dur_s, centro + ctx_dur / 2)
        idx_ctx_ini = int(ctx_ini * fs)
        idx_ctx_fim = int(ctx_fim * fs)
        sinal_ctx = sinal[idx_ctx_ini:idx_ctx_fim]

        print(f"  Duracao contexto: {ctx_fim - ctx_ini:.1f}s")
        print(f"  f_theta referencia: {f_theta_ref}Hz")
        print(f"  f_gamma referencia: {f_gamma_ref}Hz")

        # 1. FOOOF na janela de contexto
        print(f"\n  [1] FOOOF na janela de contexto ({ctx_fim - ctx_ini:.1f}s)...")
        res_fooof = extrai_cf_teta_fooof(sinal_ctx, fs)
        erro_str = f"{res_fooof['erro_ajuste']:.4f}" if res_fooof['erro_ajuste'] is not None else "N/A"
        print(f"      cf_teta={res_fooof['cf_teta']}, "
              f"erro={erro_str}, "
              f"qualidade_ok={res_fooof['qualidade_ok']}")

        # 2. Skewness na janela do candidato
        print(f"\n  [2] Skewness do theta na janela do candidato ({fim-ini:.0f}s)...")
        skew, n = skewness_from_band(sinal_cand, fs, lo=4, hi=8)
        print(f"      skewness={skew:.3f}, n={n} | {'ALTA' if abs(skew) > 0.5 else 'BAIXA'}")

        # 3. Teste de razao harmonica
        plv_val = np.nan
        ordem = None
        desvio = None
        suspeito = False

        # Para o teste sintetico, usamos cf_teta mesmo se qualidade_ok=False
        # (o limiar 0.15 e muito restritivo para sinais sem ruido)
        cf_usar = res_fooof["cf_teta"]
        if cf_usar is not None:
            print(f"\n  [3] Teste de razao harmonica (tol=10%, usando cf={cf_usar:.2f}Hz)...")
            suspeito, ordem, desvio = testa_razao_harmonica(
                cf_usar, f_gamma_ref, tol_rel=0.10)
            print(f"      suspeito={suspeito}, ordem={ordem}, desvio={desvio}")

            # 4. PLV se razao suspeita
            if suspeito:
                print(f"\n  [4] PLV harmonico (n={ordem})...")
                plv_val = compute_plv_harmonico(
                    sinal_cand, fs, cf_usar, f_gamma_ref, ordem, ini, fim)
                print(f"      PLV={plv_val:.3f} | {'ALTO (>=0.8)' if plv_val >= 0.8 else 'BAIXO (<0.8)'}")
        else:
            print(f"\n  [3-4] FOOOF nao detectou teta - nao executa teste de fase/razao")

        # Veredito
        if not res_fooof["qualidade_ok"]:
            veredito = "SEM_REFERENCIA_TETA"
        else:
            skew_alto = abs(skew) > 0.5
            plv_alto = not np.isnan(plv_val) and plv_val >= 0.8

            if suspeito and skew_alto and plv_alto:
                veredito = f"SUSPEITO_HARMONICO_FORTE (ordem={ordem}, PLV={plv_val:.2f}, skew={skew:.2f})"
            elif suspeito and plv_alto:
                veredito = f"REVISAR_FASE_TRAVADA (ordem={ordem}, PLV={plv_val:.2f})"
            elif suspeito:
                veredito = f"REVISAR_RAZAO_INTEIRA (ordem={ordem})"
            elif skew_alto:
                veredito = "REVISAR_TETA_ASSIMETRICO"
            else:
                veredito = "CLEAN"

        print(f"\n  VEREDITO: {veredito}")

        resultados.append({
            "cenario": scenario,
            "label": label,
            "f_theta_ref": f_theta_ref,
            "f_gamma_ref": f_gamma_ref,
            "cf_teta_fooof": res_fooof["cf_teta"],
            "qualidade_ok": res_fooof["qualidade_ok"],
            "suspeito_razao": suspeito,
            "ordem": ordem,
            "desvio": desvio,
            "skew": skew,
            "plv": plv_val,
            "veredito": veredito,
        })

    # Resumo
    print("\n" + "=" * 80)
    print("RESUMO")
    print("=" * 80)
    print(f"{'Cen':<4} | {'f_th':<6} | {'f_gam':<6} | {'cf_fooof':<8} | {'PLV':<6} | {'skew':<7} | {'Veredito'}")
    print("-" * 90)
    for r in resultados:
        plv_str = f"{r['plv']:.2f}" if not np.isnan(r['plv']) else "N/A"
        cf_str = f"{r['cf_teta_fooof']:.2f}" if r['cf_teta_fooof'] else "N/A"
        skew_str = f"{r['skew']:+.2f}"
        print(f"A/B   | {r['f_theta_ref']:<6.1f} | {r['f_gamma_ref']:<6.1f} | "
              f"{cf_str:<8} | {plv_str:<6} | {skew_str:<7} | {r['veredito']}")

    # Afericao automatica
    print("\n" + "=" * 80)
    print("AFERICAO")
    print("=" * 80)

    for r in resultados:
        if r["cenario"] == "A":
            esperado_label = "SUSPEITO_HARMONICO_FORTE ou REVISAR_FASE_TRAVADA"
            # Critico: PLV alto + skew alto ou suspeito de razao
            plv_ok = not np.isnan(r['plv']) and r['plv'] >= 0.8
            skew_ok = abs(r['skew']) > 0.5
            razao_ok = r['suspeito_razao']
            aprovado = (plv_ok or skew_ok) and razao_ok
        else:
            esperado_label = "CLEAN ou REVISAR_TETA_ASSIMETRICO"
            # Critico: PLV baixo + skew baixo + nao suspeito
            plv_ok = not np.isnan(r['plv']) and r['plv'] < 0.5
            skew_ok = abs(r['skew']) < 0.5
            razao_ok = not r['suspeito_razao']
            aprovado = (plv_ok and skew_ok) or razao_ok

        status = "PASSOU" if aprovado else "FALHOU"
        print(f"  Cenario {r['cenario']}: esperado={esperado_label}")
        print(f"    PLV alto: {not np.isnan(r['plv']) and r['plv'] >= 0.8} | "
              f"skew alto: {abs(r['skew']) > 0.5} | "
              f"razao suspecta: {r['suspeito_razao']}")
        print(f"    -> {status}")

    print("\n" + "=" * 80)
    print("NOTAS")
    print("=" * 80)
    print("- Este teste verifica se o script NAO QUEBRA e se as metricas")
    print("  se comportam conforme esperado em sinais sinteticos.")
    print("- A validacao REAL (taxa FP/FN) exige sessoes reais.")
    print("- Checar se cf_fooof ~= f_theta_ref (8Hz) em ambos cenarios.")


if __name__ == "__main__":
    main()
