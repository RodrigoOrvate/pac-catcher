"""
test_fooof_linha_preprocess.py - Validacao do pre-processamento de 50/60Hz
                                estilo Kuhn et al. 2026.

REUSO: a receita de limpeza de linha vive em production
`pipeline/auditorias/linha_noise_kuhn.py` (remove_pico_gaussiana /
preprocessa_linha / aplica_modo) -- este teste e' regressao desse codigo,
nao duplica a logica.

Tecnica (do paper Kuhn et al. 2026, LFP_FOOOF, Methods):
  - Subtrair APENAS a Gaussiana (pico) do PSD em log10, preservando o
    componente aperiodico local (subtrair o modelo completo cria buraco).
  - (Alternativa cirurgica em linha_noise_kuhn.remove_faixa_1f: repor a
    banda estrita pela 1/f, port do rem_noise.m do artigo.)

Aqui validamos nos cenarios sinteticos e contra a linha ela mesma:
  D original: spike 50Hz (amplitude alta) + 100Hz + artefato motor
  E estendido: D + teta genuino (nao-senoidal)

Para cada cenario, comparamos:
  - Sem preprocessamento (FOOOF amplo direto)
  - Com preprocessamento Kuhn (50Hz / 60Hz / hibrido - via aplica_modo)

Esperamos: com preprocessamento, o erro do cenario D cai de 0.265
para algo < 0.15 (dentro do limiar de producao).
"""
import numpy as np
from scipy.signal import welch
from fooof import FOOOF

from test_fooof_pico_budget import (
    gera_cenario_D, gera_cenario_E, roda_fooof
)
from pipeline.auditorias.linha_noise_kuhn import aplica_modo


def preprocessa_notch_sinal(sinal, fs, freq=60.0, bw=4.0, ordem=6):
    """Aplica notch IIR direto no SINAL (nao no PSD) antes do Welch.

    Esta e' a abordagem mais simples e a que o pipeline principal ja
    faz. Comparada com a abordagem Kuhn (subtracao de fit 1exp no PSD),
    notch no sinal e' menos elegante mas:
      - mais simples (1 linha de codigo)
      - nao introduz artefatos no componente 1/f
      - funciona bem para contaminantes estacionarios como 50/60Hz
    """
    from scipy.signal import iirnotch, filtfilt
    nyq = fs / 2
    Q = freq / bw  # fator de qualidade
    b, a = iirnotch(freq / nyq, Q, fs=fs)
    return filtfilt(b, a, sinal)


def preprocessa_notches_multiplos(sinal, fs, contaminantes):
    """Aplica multiplos notch IIR em cascata no sinal.

    contaminantes: lista de (freq, bw) tuples.
    Usado para testar o escenario D (50Hz + 100Hz + motor 45Hz).
    """
    from scipy.signal import iirnotch, filtfilt
    resultado = sinal.copy()
    for freq, bw in contaminantes:
        nyq = fs / 2
        Q = freq / bw
        b, a = iirnotch(freq / nyq, Q, fs=fs)
        resultado = filtfilt(b, a, resultado)
    return resultado


def roda_fooof_com_preprocess(sinal, fs, max_n_peaks=4, fit_range=(4, 100),
                               aperiodic_mode='knee', theta_cf_bounds=(5, 9.5),
                               pwl=(2, 5), min_h=0.05,
                               preprocessar_linha=True, f_linha=60.0,
                               modo_preprocesso='gaussiana'):
    """Replica extrai_cf_teta_fooof COM pre-processamento de 50/60Hz.

    Delega a limpeza de linha a `aplica_modo` (modulo de producao
    linha_noise_kuhn). `modo_preprocesso` seleciona gaussiana/cirurgica/
    hibrido; default gaussiana (subtrai a gaussiana em log10).
    """
    nperseg = int(1.2 * fs)
    nfft = 4000 if nperseg <= 4000 else nperseg
    freqs, psd = welch(sinal, fs=fs, window='hann',
                       nperseg=nperseg, noverlap=nperseg // 2,
                       nfft=nfft)
    if preprocessar_linha:
        print(f"  Aplicando preprocess Kuhn ({modo_preprocesso}, "
              f"f_linha={f_linha:.0f}Hz):")
        psd = aplica_modo(freqs, psd, fs, modo_preprocesso,
                          f_linha=f_linha, verbose=True)
    fm = FOOOF(aperiodic_mode=aperiodic_mode, peak_width_limits=pwl,
               min_peak_height=min_h, peak_threshold=1.0,
               max_n_peaks=max_n_peaks, verbose=False)
    try:
        fm.fit(freqs, psd, freq_range=fit_range)
    except Exception as e:
        return {"erro": None, "cf_teta": None, "teta_detectado": False,
                "n_picos": 0, "exc": str(e)}
    if not fm.has_model:
        return {"erro": None, "cf_teta": None, "teta_detectado": False,
                "n_picos": 0, "exc": "no_model"}
    try:
        err = fm.get_params('error')
        picos = fm.get_params('peak_params')
    except Exception as e:
        return {"erro": None, "cf_teta": None, "teta_detectado": False,
                "n_picos": 0, "exc": str(e)}
    if picos is None or len(picos) == 0:
        return {"erro": err, "cf_teta": None, "teta_detectado": False,
                "n_picos": 0}
    picos_arr = picos if picos.ndim > 1 else picos.reshape(1, -1)
    cand = picos_arr[(picos_arr[:, 0] >= theta_cf_bounds[0]) &
                     (picos_arr[:, 0] <= theta_cf_bounds[1])]
    cf_teta = float(cand[np.argmax(cand[:, 1]), 0]) if len(cand) > 0 else None
    return {"erro": err, "cf_teta": cf_teta,
            "teta_detectado": cf_teta is not None,
            "n_picos": len(picos_arr)}


def fmtErro(r):
    """Formata erro/ cf_teta de um resultado, tratando None."""
    e = f"{r['erro']:.4f}" if r['erro'] is not None else "n/a"
    c = f"{r['cf_teta']:.3f}" if r['cf_teta'] else "n/a"
    d = "SIM" if r['teta_detectado'] else "NAO"
    return e, c, d


def main():
    fs = 2000
    dur_s = 45.0
    f_theta = 8.0

    # Teste de sementes para D: margem de 0.015 pode ser instavel
    print("=" * 90)
    print("TESTE DE SENSIBILIDADE A SEMENTE (cenario D, Kuhn f=50Hz)")
    print("=" * 90)
    erros_sem = []
    erros_com = []
    n_sementes = 10
    for i in range(n_sementes):
        sinal, _, _, _ = gera_cenario_D(fs, dur_s, f_theta)
        # Troca a semente do rng dentro de gera_cenario_D nao da,
        # entao o cenario D e' deterministico (semente 43 fixa).
        # Para testar sementes, precisamos adicionar ruido branco
        # ao sinal D com sementes variaveis.
        rng = np.random.default_rng(100 + i)
        ruido = 0.01 * rng.standard_normal(len(sinal))
        sinal_r = sinal + ruido
        r_sem = roda_fooof(sinal_r, fs, max_n_peaks=4)
        r_com = roda_fooof_com_preprocess(sinal_r, fs, max_n_peaks=4,
                                          preprocessar_linha=True,
                                          f_linha=50.0)
        e_sem = r_sem['erro'] if r_sem['erro'] is not None else float('nan')
        e_com = r_com['erro'] if r_com['erro'] is not None else float('nan')
        erros_sem.append(e_sem)
        erros_com.append(e_com)
        print(f"  semente {i:2d}: sem_kuhn={e_sem:.4f}  com_kuhn={e_com:.4f}"
              f"  {'<0.15' if e_com < 0.15 else '>0.15'}")
    erros_sem = np.array(erros_sem)
    erros_com = np.array(erros_com)
    print(f"\n  media sem: {np.nanmean(erros_sem):.4f} ± {np.nanstd(erros_sem):.4f}")
    print(f"  media com: {np.nanmean(erros_com):.4f} ± {np.nanstd(erros_com):.4f}")
    print(f"  <0.15 (Kuhn): {(erros_com < 0.15).sum()}/{n_sementes} "
          f"({100*(erros_com < 0.15).sum()/n_sementes:.0f}%)")
    print(f"  >0.15 (Kuhn): {(erros_com >= 0.15).sum()}/{n_sementes} "
          f"({100*(erros_com >= 0.15).sum()/n_sementes:.0f}%)")

    cenarios = [
        ("A: teta+3x (soma de senoides)", gera_cenario_A),
        ("D: STRESS multi-contaminante (50+100+artefato motor)", gera_cenario_D),
    ]

    print("=" * 90)
    print("VALIDACAO: pre-processamento Kuhn (50/60Hz) antes do FOOOF amplo")
    print("=" * 90)

    for nome, gen in cenarios:
        sinal, _, _, _ = gen(fs, dur_s, f_theta)
        print(f"\n--- {nome} ---")

        # SEM preprocessamento
        print("SEM preprocessamento:")
        r_sem = roda_fooof(sinal, fs, max_n_peaks=4)
        e, c, d = fmtErro(r_sem)
        print(f"  erro={e}, cf_teta={c}, detectou={d}, n_picos={r_sem['n_picos']}")

        # COM preprocessamento Kuhn (Gauss-only, f=50Hz)
        print("COM preprocessamento Kuhn (Gauss-only, f=50Hz):")
        r_com = roda_fooof_com_preprocess(sinal, fs, max_n_peaks=4,
                                          preprocessar_linha=True,
                                          f_linha=50.0)
        e, c, d = fmtErro(r_com)
        print(f"  erro={e}, cf_teta={c}, detectou={d}, n_picos={r_com['n_picos']}")

        # Diagnostico
        if r_com['erro'] is not None and r_sem['erro'] is not None:
            delta = r_sem['erro'] - r_com['erro']
            if delta > 0.05:
                print(f"  -> Kuhn Gauss-only REDUZIU erro em {delta:.4f} "
                      f"({100*delta/r_sem['erro']:.0f}% de melhoria).")
            elif delta < -0.05:
                print(f"  -> Kuhn Gauss-only PIOROU erro em {-delta:.4f}.")
            else:
                print("  -> Kuhn Gauss-only SEM efeito significativo.")
            if r_com['erro'] < 0.15:
                print("  -> AGORA DENTRO do limiar 0.15.")
            else:
                print(f"  -> AINDA ACIMA do limiar 0.15 ({r_com['erro']:.4f}).")

        # Notch no sinal (alternativa simples, 50Hz = rede EU)
        print("COM notch no sinal em 50Hz (rede EU, banda 48-52Hz):")
        sinal_notch = preprocessa_notch_sinal(sinal, fs, freq=50.0, bw=4.0)
        r_notch = roda_fooof(sinal_notch, fs, max_n_peaks=4)
        e, c, d = fmtErro(r_notch)
        print(f"  erro={e}, cf_teta={c}, detectou={d}, n_picos={r_notch['n_picos']}")
        if r_notch['erro'] is not None and r_sem['erro'] is not None:
            dn = r_sem['erro'] - r_notch['erro']
            if dn > 0.05:
                print(f"  -> Notch REDUZIU erro em {dn:.4f}.")
            elif dn < -0.05:
                print(f"  -> Notch PIOROU erro em {-dn:.4f}.")
            print(f"  -> {'DENTRO' if r_notch['erro'] < 0.15 else 'ACIMA'}"
                  f" do limiar 0.15.")

        # So' no D: testes adicionais multi-contaminante
        if 'D' in nome:
            # Notch triplo (50+100+45Hz)
            print("COM notch TRIPLO (50Hz + 100Hz + 45Hz motor):")
            sinal_triplo = preprocessa_notches_multiplos(
                sinal, fs, [(50.0, 4.0), (100.0, 4.0), (45.0, 6.0)])
            r_triplo = roda_fooof(sinal_triplo, fs, max_n_peaks=4)
            e, c, d = fmtErro(r_triplo)
            print(f"  erro={e}, cf_teta={c}, detectou={d},"
                  f" n_picos={r_triplo['n_picos']}")
            if r_triplo['erro'] is not None and r_sem['erro'] is not None:
                dt = r_sem['erro'] - r_triplo['erro']
                if dt > 0.05:
                    print(f"  -> Notch triplo REDUZIU erro em {dt:.4f}.")
                elif dt < -0.05:
                    print(f"  -> Notch triplo PIOROU erro em {-dt:.4f}.")
                print(f"  -> {'DENTRO' if r_triplo['erro'] < 0.15 else 'ACIMA'}"
                      f" do limiar 0.15.")

            # Kuhn 60Hz (caso BR) — D nao tem contaminante em 60Hz
            print("COM preprocess Kuhn (Gauss-only, f=60Hz — caso BR):")
            r_br = roda_fooof_com_preprocess(sinal, fs, max_n_peaks=4,
                                             preprocessar_linha=True,
                                             f_linha=60.0)
            e, c, d = fmtErro(r_br)
            print(f"  erro={e}, cf_teta={c}, detectou={d},"
                  f" n_picos={r_br['n_picos']}")
            if r_br['erro'] is not None and r_sem['erro'] is not None:
                dbr = r_sem['erro'] - r_br['erro']
                if abs(dbr) > 0.05:
                    tag = "REDUZIU" if dbr > 0 else "PIOROU"
                    print(f"  -> Kuhn-60 {tag} erro em {abs(dbr):.4f}.")
                else:
                    print("  -> Kuhn-60 SEM efeito (esperado: D nao tem"
                          " contaminante em 60Hz).")

    print("\n" + "=" * 90)
    print("CONCLUSAO (atualizada com receita correta de Kuhn, 04/09/2026)")
    print("=" * 90)
    print("Receita correta (Kuhn et al. 2026): subtrair APENAS a Gaussiana")
    print("(pico) do PSD, preservando o componente aperiodico local.")
    print("Minha tentativa anterior subtraia o modelo COMPLETO (1/f+picos),")
    print("o que causava buraco no espectro (erro 0.265->0.700).")
    print("")
    print("Nota sobre 50Hz vs 60Hz:")
    print("  O artigo usa 50Hz (rede EU). Nosso pipeline usa 60Hz (rede BR).")
    print("  O cenario D tem contaminante em 50Hz (construido pra testar a")
    print("  receita do artigo). Testes com f=60Hz mostram comportamento")
    print("  do caso BR.")
    print("")
    print("Ver resultados acima para decidir se Kuhn Gauss-only resolve")
    print("o cenario D (erro precisa cair < 0.15 para portao ACEITAR).")


if __name__ == "__main__":
    from test_fooof_pico_budget import gera_cenario_A  # import local p/ main
    main()
