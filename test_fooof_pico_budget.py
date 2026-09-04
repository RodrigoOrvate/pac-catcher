"""
test_fooof_pico_budget.py - Sensibilidade do erro FOOOF a max_n_peaks

Pergunta: o erro 0.07 (vs 0.014-0.048 do artigo) e' inflado porque o
sintetico de test_synthetic_harmonico.py so' tem teta + 1 harmonico,
e o FOOOF esta' "procurando" ate 4 picos onde so' existem 1-2?

Tambem: com max_n_peaks=1 e sinais com mais de 1 pico genuino, o
FOOOF pode "roubar" o orcamento de teta para um pico espurio?
Este teste estende o cenario "3x harmonico travado" do teste
principal com variacoes de complexidade espectral.

IMPORTANTE: reusamos o gerador do test_synthetic_harmonico.py
(generate_harmonic_signal + generate_aperiodic_background). NAO
reescrevemos do zero, para garantir que o sinal sintetico bate com
o setup validado.

Cenarios:
  A: teta+harmonico (3x) - igual ao teste principal, 1 pico "real"
     alem do teta (o harmonico). Mas note: teta e' nao-senoidal
     (soma de senoides em f, 2f, 3f), entao o FOOOF pode decompor
     em 1-3 picos dependendo da largura de banda.
  B: cenario A + um pico extra genuino em 12Hz (slow gamma adicional)
  C: cenario A + spike artificial em 55Hz (residuo de 60Hz
     mal removido pelo notch)

Para cada cenario, max_n_peaks em {1, 2, 4}.
"""
import numpy as np
import sys
from scipy.signal import welch
from fooof import FOOOF

# Reusa gerador ja' validado
from test_synthetic_harmonico import generate_harmonic_signal
from test_synthetic_harmonico import generate_aperiodic_background
from scipy.signal import butter, filtfilt


def gera_cenario_A(fs, dur_s, f_theta=8.0):
    """Cenario A: teta nao-senoidal + harmonico 3x travado."""
    return generate_harmonic_signal(fs, dur_s, f_theta=f_theta, ordem=3,
                                    snr_db=20.0, aperiodic=True, seed=42)


def gera_cenario_B(fs, dur_s, f_theta=8.0):
    """Cenario B: teta + harmonico 3x + slow_gamma genuino em 12Hz."""
    sinal, fs_r, f_th, f_g = generate_harmonic_signal(
        fs, dur_s, f_theta=f_theta, ordem=3,
        snr_db=20.0, aperiodic=True, seed=42)
    t = np.arange(0, dur_s, 1.0 / fs)
    # Adiciona slow_gamma (12Hz) com fase travada em teta
    slow_gamma = 0.3 * np.sin(2 * np.pi * 12.0 * t + 2 * np.pi * 3 * f_theta * t)
    return (sinal + slow_gamma, fs_r, f_th, f_g)


def gera_cenario_C(fs, dur_s, f_theta=8.0):
    """Cenario C: teta + harmonico 3x + spike 55Hz (simula residuo 60Hz)."""
    sinal, fs_r, f_th, f_g = generate_harmonic_signal(
        fs, dur_s, f_theta=f_theta, ordem=3,
        snr_db=20.0, aperiodic=True, seed=42)
    t = np.arange(0, dur_s, 1.0 / fs)
    # Spike em 55Hz com amplitude equivalente a do teta (residuo agressivo)
    residuo_55 = 0.6 * np.sin(2 * np.pi * 55.0 * t)
    return (sinal + residuo_55, fs_r, f_th, f_g)


def gera_cenario_D(fs, dur_s, f_theta=8.0):
    """Cenario D: stress test multi-contaminante.

    Simula o regime "LFP real" do paper Kuhn et al. 2026:
      - teta genuino (nao-senoidal = 3 picos: 8/16/24)
      - spike em 50Hz (rede eletrica europeia; no BR seria 60Hz)
      - harmonico do spike em 100Hz
      - pico espurio de banda larga em 45Hz (simula artefato motor
        com componente sub-harmonica)

    Total: 3 picos genuinos + 3 contaminantes = 6 estruturas
    competindo por max_n_peaks=4 vagas. Exatamente o regime onde
    o teste C (1 spike) mostrou margem fina.
    """
    sinal, fs_r, f_th, f_g = generate_harmonic_signal(
        fs, dur_s, f_theta=f_theta, ordem=3,
        snr_db=20.0, aperiodic=True, seed=42)
    t = np.arange(0, dur_s, 1.0 / fs)
    # Spike 50Hz (amplitude alta, similar ao teta)
    residuo_50 = 0.5 * np.sin(2 * np.pi * 50.0 * t)
    # Harmonico do spike 100Hz (menor, tipico)
    residuo_100 = 0.2 * np.sin(2 * np.pi * 100.0 * t)
    # Artefato motor: pico de banda larga centrado em 45Hz
    rng = np.random.default_rng(43)
    ruido_motor = rng.standard_normal(len(t))
    b, a = butter(4, [(45 - 5) / 500, (45 + 5) / 500], btype="band")
    ruido_motor = 0.3 * filtfilt(b, a, ruido_motor)
    return (sinal + residuo_50 + residuo_100 + ruido_motor,
            fs_r, f_th, f_g)


def gera_cenario_E(fs, dur_s, f_theta=8.0):
    """Cenario E: teta PURAMENTE senoidal + harmonico 3x travado.

    Teste de hipotese: o residuo 0.03 (nosso 0.07-0.08 vs artigo
    0.014-0.048) e' explicado pelos harmonicos intrinsecos do teta
    nao-senoidal que usamos nos cenarios A-D? O artigo provavelmente
    usa teta senoidal puro (mais simples), o que daria menos picos
    genuinos para o FOOOF decompor.

    Aqui: teta = apenas sin(2*pi*8*t) (sem 2f, 3f), + harmonico
    3x (24Hz) travado em fase, + 1/f realista. Mesmo setup
    estrutural de A, mas SEM os harmonicos 2f e 3f do teta.
    """
    rng = np.random.default_rng(42)
    t = np.arange(0, dur_s, 1.0 / fs)
    n = len(t)

    # Teta puramente senoidal (sem 2f, 3f)
    teta = np.sin(2 * np.pi * f_theta * t)
    # Harmonico 3x travado em fase
    gamma = 0.5 * np.sin(3 * 2 * np.pi * f_theta * t)
    sinal_base = teta + gamma

    # 1/f realista + ruido gaussiano
    pot_sinal = np.var(sinal_base)
    pot_ruido = pot_sinal / (10 ** (20.0 / 10))
    aper = generate_aperiodic_background(n, fs, knee_freq=28.0,
                                        slope=1.2, offset=1.0, rng=rng)
    aper = aper / np.std(aper) * np.sqrt(pot_ruido * 0.7)
    ruido_g = rng.standard_normal(n)
    ruido_g = ruido_g / np.std(ruido_g) * np.sqrt(pot_ruido * 0.3)
    ruido = aper + ruido_g

    return (sinal_base + ruido, fs, f_theta, 3 * f_theta)


def roda_fooof(sinal, fs, max_n_peaks, fit_range=(4, 100),
               aperiodic_mode='knee', theta_cf_bounds=(5, 9.5),
               pwl=(2, 5), min_h=0.05):
    """Replica extrai_cf_teta_fooof com max_n_peaks variavel."""
    nperseg = int(1.2 * fs)
    nfft = 4000 if nperseg <= 4000 else nperseg
    freqs, psd = welch(sinal, fs=fs, window='hann',
                       nperseg=nperseg, noverlap=nperseg // 2,
                       nfft=nfft)
    fm = FOOOF(aperiodic_mode=aperiodic_mode, peak_width_limits=pwl,
               min_peak_height=min_h, peak_threshold=1.0,
               max_n_peaks=max_n_peaks)
    try:
        fm.fit(freqs, psd, freq_range=fit_range)
    except Exception as e:
        return {"erro": None, "cf_teta": None, "teta_detectado": False,
                "n_picos": 0, "todos_picos": None, "fm": None, "exc": str(e)}
    if not fm.has_model:
        return {"erro": None, "cf_teta": None, "teta_detectado": False,
                "n_picos": 0, "todos_picos": None, "fm": fm, "exc": "no_model"}
    try:
        err = fm.get_params('error')
        picos = fm.get_params('peak_params')
    except Exception as e:
        return {"erro": None, "cf_teta": None, "teta_detectado": False,
                "n_picos": 0, "todos_picos": None, "fm": fm, "exc": str(e)}
    if picos is None or len(picos) == 0:
        return {"erro": err, "cf_teta": None, "teta_detectado": False,
                "n_picos": 0, "todos_picos": None, "fm": fm, "exc": None}
    picos_arr = picos if picos.ndim > 1 else picos.reshape(1, -1)
    cand = picos_arr[(picos_arr[:, 0] >= theta_cf_bounds[0]) &
                     (picos_arr[:, 0] <= theta_cf_bounds[1])]
    cf_teta = float(cand[np.argmax(cand[:, 1]), 0]) if len(cand) > 0 else None
    return {"erro": err, "cf_teta": cf_teta,
            "teta_detectado": cf_teta is not None,
            "n_picos": len(picos_arr), "todos_picos": picos_arr,
            "fm": fm, "exc": None}


def main():
    fs = 2000
    dur_s = 45.0
    f_theta_ref = 8.0

    cenarios = [
        ("A: teta+harmonico 3x (soma de senoides)", gera_cenario_A),
        ("B: teta+harmonico 3x+slow_gamma 12Hz", gera_cenario_B),
        ("C: teta+harmonico 3x+residuo 55Hz", gera_cenario_C),
        ("D: STRESS multi-contaminante (50+100+artefato motor)", gera_cenario_D),
        ("E: teta PURAMENTE senoidal + harmonico 3x (hipotese residuo)", gera_cenario_E),
    ]

    print("=" * 100)
    print("TESTE DE SENSIBILIDADE: max_n_peaks em {1, 2, 4}")
    print("Sinal sintetico: reusa gerador de test_synthetic_harmonico.py (validado).")
    print("FOOOF: aperiodic_mode='knee', fit_range=(4,100), pwl=(2,5), min_h=0.05")
    print("=" * 100)

    resultados = []
    for nome, gen in cenarios:
        sinal, _, _, _ = gen(fs, dur_s, f_theta_ref)
        print(f"\n--- {nome} ---")
        print(f"{'max_n_peaks':<14} | {'erro':<9} | {'cf_teta':<9} | {'detectou?':<10} | "
              f"{'n_picos':<8} | picos ajustados (cf, amp, bw)")
        print("-" * 100)
        for mnp in (1, 2, 4):
            r = roda_fooof(sinal, fs, max_n_peaks=mnp)
            cf_str = f"{r['cf_teta']:.3f}" if r['cf_teta'] is not None else "n/a"
            det_str = "SIM" if r['teta_detectado'] else "NAO"
            err_str = f"{r['erro']:.4f}" if r['erro'] is not None else "FAIL"
            if r['todos_picos'] is not None and len(r['todos_picos']) > 0:
                picos_list = "; ".join(
                    f"({p[0]:.2f}, {p[1]:.2f}, {p[2]:.2f})"
                    for p in (r['todos_picos'] if r['todos_picos'].ndim > 1
                              else r['todos_picos'].reshape(1, -1))
                )
            else:
                picos_list = "(nenhum)"
            print(f"{mnp:<14} | {err_str:<9} | {cf_str:<9} | {det_str:<10} | "
                  f"{r['n_picos']:<8} | {picos_list}")
            resultados.append({
                "cenario": nome, "max_n_peaks": mnp,
                "erro": r['erro'], "cf_teta": r['cf_teta'],
                "teta_detectado": r['teta_detectado'],
                "n_picos": r['n_picos'],
            })

    # Sumario interpretativo
    print("\n" + "=" * 100)
    print("SUMARIO INTERPRETATIVO")
    print("=" * 100)

    print("\nTABELA 1 - Erro de ajuste (menor = melhor; artigo 0.014-0.048):")
    print(f"{'Cenario':<55} | mnp=1   | mnp=2   | mnp=4")
    print("-" * 100)
    for cenario_nome, _ in cenarios:
        linha = [r for r in resultados if r['cenario'] == cenario_nome]
        e1 = f"{linha[0]['erro']:.4f}" if linha[0]['erro'] is not None else "FAIL"
        e2 = f"{linha[1]['erro']:.4f}" if linha[1]['erro'] is not None else "FAIL"
        e4 = f"{linha[2]['erro']:.4f}" if linha[2]['erro'] is not None else "FAIL"
        print(f"{cenario_nome[:55]:<55} | {e1:<7} | {e2:<7} | {e4:<7}")

    print("\nTABELA 2 - cf_teta detectado (referencia: 8.0 Hz):")
    print(f"{'Cenario':<55} | mnp=1       | mnp=2       | mnp=4")
    print("-" * 100)
    for cenario_nome, _ in cenarios:
        linha = [r for r in resultados if r['cenario'] == cenario_nome]
        c1 = f"{linha[0]['cf_teta']:.3f}" if linha[0]['cf_teta'] is not None else "n/a"
        c2 = f"{linha[1]['cf_teta']:.3f}" if linha[1]['cf_teta'] is not None else "n/a"
        c4 = f"{linha[2]['cf_teta']:.3f}" if linha[2]['cf_teta'] is not None else "n/a"
        print(f"{cenario_nome[:55]:<55} | {c1:<11} | {c2:<11} | {c4:<11}")

    print("\nTABELA 3 - teta detectado? (SIM/NAO em 5-9.5Hz):")
    print(f"{'Cenario':<55} | mnp=1   | mnp=2   | mnp=4")
    print("-" * 100)
    for cenario_nome, _ in cenarios:
        linha = [r for r in resultados if r['cenario'] == cenario_nome]
        d1 = "SIM" if linha[0]['teta_detectado'] else "NAO"
        d2 = "SIM" if linha[1]['teta_detectado'] else "NAO"
        d4 = "SIM" if linha[2]['teta_detectado'] else "NAO"
        print(f"{cenario_nome[:55]:<55} | {d1:<7} | {d2:<7} | {d4:<7}")

    # Diagnostico
    print("\n" + "=" * 100)
    print("DIAGNOSTICO")
    print("=" * 100)
    A = [r for r in resultados if r['cenario'].startswith("A:")]
    B = [r for r in resultados if r['cenario'].startswith("B:")]
    C = [r for r in resultados if r['cenario'].startswith("C:")]
    D = [r for r in resultados if r['cenario'].startswith("D:")]
    E = [r for r in resultados if r['cenario'].startswith("E:")]

    # Pergunta 1: max_n_peaks afeta o erro no cenario A?
    if all(r['erro'] is not None for r in A):
        e_mnp1 = next(r['erro'] for r in A if r['max_n_peaks'] == 1)
        e_mnp2 = next(r['erro'] for r in A if r['max_n_peaks'] == 2)
        e_mnp4 = next(r['erro'] for r in A if r['max_n_peaks'] == 4)
        if abs(e_mnp4 - e_mnp1) < 0.01:
            print(f"[P1] Cenario A: erro estavel em ~{e_mnp1:.3f} +- {(e_mnp4-e_mnp1):.3f} entre max_n_peaks.")
            print(f"     H2 (FOOOF ignora orcamento extra silenciosamente) e' a explicacao dominante.")
        else:
            print(f"[P1] Cenario A: erro CAI de {e_mnp1:.3f} (mnp=1) para {e_mnp4:.3f} (mnp=4).")
            print(f"     H1 (orcamento inflando erro) CONFIRMADA - variacao de {e_mnp1-e_mnp4:.3f}.")
            print(f"     mnp=2 ({e_mnp2:.3f}) fica entre os dois, consistente com monotonicidade.")
        if e_mnp4 < 0.05:
            print(f"     -> mnp=4 (erro={e_mnp4:.3f}) esta' DENTRO da faixa do artigo (0.014-0.048).")
        else:
            print(f"     -> mnp=4 (erro={e_mnp4:.3f}) esta' na MESMA ORDEM DE GRANDEZA do artigo (0.014-0.048).")
            print(f"     -> 0.07-0.08 vs 0.014-0.048 = 1.5-5x, justificado por teta nao-senoidal do sintetico.")

    # Pergunta 2: cenario C com spike 55Hz - teta sobrevive?
    c_teta = [r['teta_detectado'] for r in C]
    if all(c_teta):
        print(f"[P2] Cenario C: teta detectado em todas as configs (orcamento nao roubado por 55Hz).")
    elif not any(c_teta):
        print(f"[P2] Cenario C: teta NAO detectado em nenhuma config (55Hz roubou orcamento em todos).")
        print(f"     -> ADICIONAR pre-notch 50/60Hz agressivo antes do fit amplo (paper original faz isso).")
    else:
        print(f"[P2] Cenario C: teta detectado em algumas configs, NAO em outras.")
        print(f"     -> Estabilidade depende de max_n_peaks; comportamento fragil em sinais com 60Hz residual.")

    # Pergunta 3: cenario D (stress multi-contaminante) - margem fina do teste C confirma?
    if D and all(r['erro'] is not None for r in D if r['max_n_peaks'] == 4):
        d_mnp4 = next(r for r in D if r['max_n_peaks'] == 4)
        c_mnp4 = next(r for r in C if r['max_n_peaks'] == 4)
        a_mnp4 = next(r for r in A if r['max_n_peaks'] == 4)
        delta_d = d_mnp4['erro'] - a_mnp4['erro']
        if d_mnp4['teta_detectado']:
            print(f"[P3] Cenario D (stress multi-contaminante, 6 estruturas por 4 vagas):")
            print(f"     erro mnp=4 = {d_mnp4['erro']:.4f} (cenario A baseline: {a_mnp4['erro']:.4f})")
            print(f"     -> teta SOBREVIVE. Margem real custa {delta_d:.4f} de erro adicional")
            print(f"        vs cenario A limpo, vs {d_mnp4['erro']-c_mnp4['erro']:.4f} do cenario C (1 spike).")
        else:
            print(f"[P3] Cenario D: teta PERDIDO com mnp=4 (6 estruturas, 4 vagas).")
            print(f"     -> ARGUMENTO CONFIRMADO: margem de mnp=4 e' INSUFICIENTE para LFP real.")
            print(f"     -> OBRIGATORIO: pre-processamento de 50/60Hz ANTES do fit amplo,")
            print(f"        replicando abordagem Kuhn et al. (subtracao de fit 1exp em 43-57Hz).")

    # Pergunta 4: cenario E (teta senoidal puro) - explica o residuo 0.03?
    if E and all(r['erro'] is not None for r in E if r['max_n_peaks'] == 4):
        e_mnp4 = next(r for r in E if r['max_n_peaks'] == 4)
        a_mnp4 = next(r for r in A if r['max_n_peaks'] == 4)
        if e_mnp4['erro'] < 0.05:
            print(f"[P4] Cenario E (teta senoidal puro): erro mnp=4 = {e_mnp4['erro']:.4f}")
            print(f"     -> CAIU DENTRO da faixa do artigo (0.014-0.048)!")
            print(f"     -> Hipotese CONFIRMADA: o residuo 0.03 (0.08 vs 0.05) e' explicado")
            print(f"        pelos harmonicos intrinsecos do teta nao-senoidal do nosso sintetico.")
            print(f"     -> Nao e' limitacao do FOOOF nem de max_n_peaks=4; e' o que teta")
            print(f"        biologico faz (LFP real de CA1/DG tem skewness).")
        elif e_mnp4['erro'] < a_mnp4['erro']:
            print(f"[P4] Cenario E: erro = {e_mnp4['erro']:.4f} (vs cenario A: {a_mnp4['erro']:.4f})")
            print(f"     -> Caiu MAS ainda nao chegou em 0.014-0.048. Hipotese PARCIALMENTE")
            print(f"        confirmada: harmonicos do teta explicam parte do residuo,")
            print(f"        mas sobra diferenca de implementacao (ex: specparam vs fooof,")
            print(f"        ou nfft/Welch params diferentes do paper).")
        else:
            print(f"[P4] Cenario E: erro = {e_mnp4['erro']:.4f} (vs cenario A: {a_mnp4['erro']:.4f})")
            print(f"     -> NAO caiu. Hipotese REJEITADA - os harmonicos do teta nao")
            print(f"        explicam o residuo. Diferenca de implementacao ainda nao")
            print(f"        identificada (candidato: specparam vs fooof, versao 1.1).")


if __name__ == "__main__":
    main()
