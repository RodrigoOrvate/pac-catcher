"""
triagem_pac.py
==========================================
Varredura automática de LFP para flagar janelas com suspeita de
Acoplamento Fase-Amplitude (PAC) em múltiplos pares de banda.

PARES SUPORTADOS (via --pares):
  theta_gamma  : Theta 4-8 Hz (fase) × Gamma  30-80  Hz (amplitude)  — exploração
  theta_hg     : Theta 4-8 Hz (fase) × HG     80-150 Hz (amplitude)  — misto
  theta_hfo    : Theta 4-8 Hz (fase) × HFO   150-250 Hz (amplitude)  — repouso/SWR

Para cada par calculado, o script reporta:
  mi_<par>    : KL-MI observado
  z_<par>     : z-score vs. 200 surrogates (deslocamento circular)
  p_<par>     : p-valor empírico

Colunas adicionais por janela:
  teta_ok      : 1 se theta (4-8 Hz) está presente na janela (SNR vs flancos)
                 CALCULADO POR JANELA — não reutilizado entre janelas.
  ratio_hfo_gamma : pot(HFO) / pot(Gamma) — proxy de harmônico de Gamma.
                    Se ≈ constante inteira e alto: HFO provavelmente é harmônico.

DIFERENÇA CRÍTICA vs. versão anterior:
  A versão anterior calculava apenas Theta-Gamma e não tinha teta_ok por janela.
  Esta versão calcula teta_ok DENTRO do loop de janelas (não cache externo)
  e suporta qualquer combinação dos três pares em uma única passagem.

Uso:
    # Um único par (compatibilidade com versão anterior)
    python triagem_pac.py --pasta /caminho/ns2 --saida resultados.csv

    # Todos os três pares
    python triagem_pac.py --pasta /caminho/ns2 --pares theta_gamma theta_hg theta_hfo \\
        --saida resultados_triplo.csv

    # Demo com dados sintéticos
    python triagem_pac.py --demo --pares theta_gamma theta_hg

Requer: neo, numpy, scipy, pandas
"""

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
import scipy.signal as signal
from scipy.signal import welch

# ==========================================
# CONFIGURAÇÃO DOS PARES DE BANDA
# ==========================================

BAND_PAIRS = {
    "theta_gamma": {
        "fase":  (4.0,   8.0),
        "amp":   (30.0,  80.0),
        "descricao": "Theta 4-8 Hz × Gamma 30-80 Hz (exploração locomotora)",
    },
    "theta_hg": {
        "fase":  (4.0,   8.0),
        "amp":   (80.0, 150.0),
        "descricao": "Theta 4-8 Hz × High-Gamma 80-150 Hz",
    },
    "theta_hfo": {
        "fase":  (4.0,   8.0),
        "amp":  (150.0, 250.0),
        "descricao": "Theta 4-8 Hz × HFO 150-250 Hz (repouso/SWR-associado)",
    },
}


# ==========================================
# FUNÇÕES DE PROCESSAMENTO
# ==========================================

def filtra_sinal(sinal, lowcut, highcut, fs, order=3):
    nyq = 0.5 * fs
    low  = max(lowcut  / nyq, 1e-6)
    high = min(highcut / nyq, 0.999)
    b, a = signal.butter(order, [low, high], btype="bandpass")
    return signal.filtfilt(b, a, sinal)


def calcula_mi(fase, envelope, n_bins=18):
    bins = np.linspace(-np.pi, np.pi, n_bins + 1)
    bin_idx = np.clip(np.digitize(fase, bins) - 1, 0, n_bins - 1)
    return _mi_de_bin_idx(bin_idx, envelope, n_bins)


def _mi_de_bin_idx(bin_idx, envelope, n_bins):
    """Núcleo vetorizado: soma/conta por bin via bincount em vez de loop Python."""
    soma_bins = np.bincount(bin_idx, weights=envelope, minlength=n_bins)
    cont_bins = np.bincount(bin_idx, minlength=n_bins)
    media_bins = np.divide(soma_bins, cont_bins, out=np.zeros(n_bins), where=cont_bins > 0)

    soma = np.sum(media_bins)
    if soma <= 0:
        return 0.0
    P = media_bins / soma
    H = -np.sum(P * np.log(P + 1e-10))
    return (np.log(n_bins) - H) / np.log(n_bins)


def mi_com_surrogates(fase, envelope, fs, n_surr=200, n_bins=18,
                       shift_min_s=1.0, rng=None):
    """
    Calcula MI observado e compara contra n_surr surrogates por deslocamento
    circular do envelope (quebra relação temporal fase-amplitude, preserva espectro).

    Retorna: mi_obs, z_score, p_empirico, mi_surr_media, mi_surr_dp
    """
    if rng is None:
        rng = np.random.default_rng()

    bins = np.linspace(-np.pi, np.pi, n_bins + 1)
    bin_idx = np.clip(np.digitize(fase, bins) - 1, 0, n_bins - 1)
    mi_obs = _mi_de_bin_idx(bin_idx, envelope, n_bins)

    n = len(envelope)
    shift_min = int(shift_min_s * fs)
    if n <= 2 * shift_min:
        shift_min = max(1, n // 10)

    mi_surr = np.empty(n_surr)
    deslocamentos = rng.integers(shift_min, n - shift_min, size=n_surr)
    for i, desloc in enumerate(deslocamentos):
        env_shift = np.roll(envelope, desloc)
        mi_surr[i] = _mi_de_bin_idx(bin_idx, env_shift, n_bins)

    media = np.mean(mi_surr)
    dp    = np.std(mi_surr)
    z     = (mi_obs - media) / dp if dp > 0 else 0.0
    p_emp = float(np.mean(mi_surr >= mi_obs))

    return mi_obs, z, p_emp, media, dp


def teta_ok_por_janela(trecho, fs, theta_band=(4.0, 8.0),
                        flancos=((2.0, 4.0), (8.0, 12.0)),
                        limiar_snr=1.5):
    """
    Verifica se theta está presente na janela.

    CALCULADO POR JANELA — não deve ser cacheado nem reutilizado entre janelas.

    Critério: SNR espectral = pot(theta) / mean(pot(flancos)).
    Flancos: 2-4 Hz e 8-12 Hz (imediatamente adjacentes ao theta 4-8 Hz).
    limiar_snr=1.5 → theta deve ter 50% mais potência que os flancos.

    Retorna: bool (True = theta presente)
    """
    nperseg = min(len(trecho), int(2 * fs))
    nfft = max(nperseg, int(4 * fs))
    f, psd = welch(trecho, fs=fs, nperseg=nperseg, nfft=nfft, window='hann')

    def pot_banda(lo, hi):
        mask = (f >= lo) & (f <= hi)
        return float(np.mean(psd[mask])) if mask.sum() > 0 else 1e-12

    p_theta  = pot_banda(*theta_band)
    p_flancos = np.mean([pot_banda(*fl) for fl in flancos])
    snr = p_theta / (p_flancos + 1e-12)
    return snr >= limiar_snr


def proxy_artefato_motor(sinal, fs, banda=(150, 450)):
    """
    Proxy heurístico de artefato muscular/movimento: potência relativa
    em banda larga de alta frequência (acima de Gamma), onde EMG
    tipicamente vaza no LFP. NÃO substitui EMG real.
    """
    nyq = 0.5 * fs
    alta = filtra_sinal(sinal, banda[0], min(banda[1], nyq * 0.98), fs, order=3)
    pot_alta  = np.mean(alta ** 2)
    pot_total = np.mean(sinal.astype(float) ** 2) + 1e-12
    return pot_alta / pot_total


def detecta_transiente(sinal, fs, limiar_diff=8.0, limiar_amp=8.0):
    """
    CAMADA 1 - FILTRO DE TRANSIENTE (Domínio do Tempo).
    Detecta artefatos de cabo/movimento abruptos.

    Retorna dict com:
      transiente_encontrado, frac_transiente, max_diff_z, max_amp_z
    """
    sinal = np.asarray(sinal, dtype=np.float64)

    diff_sinal = np.abs(np.diff(sinal))
    media_diff = np.mean(diff_sinal)
    std_diff   = np.std(diff_sinal)
    if std_diff > 0:
        diff_zscore = (diff_sinal - media_diff) / std_diff
        max_diff_z  = float(np.max(diff_zscore))
    else:
        diff_zscore = np.zeros(len(diff_sinal))
        max_diff_z  = 0.0

    media_amp = np.mean(np.abs(sinal))
    std_amp   = np.std(sinal)
    if std_amp > 0:
        amp_zscore = (np.abs(sinal) - media_amp) / std_amp
        max_amp_z  = float(np.max(amp_zscore))
    else:
        amp_zscore = np.zeros(len(sinal))
        max_amp_z  = 0.0

    mask_diff_trunc = np.zeros(len(amp_zscore), dtype=bool)
    mask_diff_trunc[:-1] = diff_zscore > limiar_diff
    mask_amp = amp_zscore > limiar_amp
    n_afetados = np.sum(mask_diff_trunc | mask_amp)
    frac_transiente = n_afetados / len(sinal)

    transiente_encontrado = (max_diff_z > limiar_diff) or (max_amp_z > limiar_amp)

    return {
        "transiente_encontrado": transiente_encontrado,
        "frac_transiente":       frac_transiente,
        "max_diff_z":            max_diff_z,
        "max_amp_z":             max_amp_z,
        "limiar_diff":           limiar_diff,
        "limiar_amp":            limiar_amp,
    }


def correlacao_gama_ruido(sinal, fs, theta_band=(4, 8), gamma_band=(30, 80),
                            ruido_band=(150, 250), limiar_r=0.6):
    """
    CAMADA 4 - PUNIÇÃO POR BANDA LARGA.
    Correlação de Pearson entre envelope do Gamma e banda de ruído.
    Se r > limiar_r: provavelmente transiente mecânico.
    Também retorna MVL (Mean Vector Length) para CAMADA 2.
    """
    lfp_theta = filtra_sinal(sinal, *theta_band, fs)
    lfp_gamma = filtra_sinal(sinal, *gamma_band, fs)
    nyq = 0.5 * fs
    lfp_ruido = filtra_sinal(sinal, ruido_band[0], min(ruido_band[1], nyq * 0.98), fs)

    env_gamma = np.abs(signal.hilbert(lfp_gamma))
    env_ruido = np.abs(signal.hilbert(lfp_ruido))

    if np.std(env_gamma) > 0 and np.std(env_ruido) > 0:
        correlacao = float(np.corrcoef(env_gamma, env_ruido)[0, 1])
    else:
        correlacao = 0.0

    fase = np.angle(signal.hilbert(lfp_theta))
    env_norm = (env_gamma - np.mean(env_gamma)) / (np.std(env_gamma) + 1e-12)
    complexo  = env_norm * np.exp(1j * fase)
    mvl = float(np.abs(np.mean(complexo)))

    return {
        "correlacao_ruido":    correlacao,
        "suspeito_banda_larga": correlacao > limiar_r,
        "mvl":                  mvl,
        "limiar_r":             limiar_r,
    }


def _ratio_hfo_gamma(trecho, fs):
    """Potência relativa HFO / Gamma. Valor alto + razão inteira → suspeita de harmônico."""
    nyq = 0.5 * fs
    def pot(lo, hi):
        hi = min(hi, nyq * 0.98)
        if lo >= hi:
            return 1e-12
        filt = filtra_sinal(trecho, lo, hi, fs)
        return float(np.mean(filt ** 2)) + 1e-12
    return pot(150, 250) / pot(30, 80)


# ==========================================
# VARREDURA MULTI-BANDA
# ==========================================

def varre_canal(sinal, fs, window_s=10.0, step_s=5.0, n_surr=200,
                 pares=None, rng=None, rotulo_progresso=None):
    """
    Varre um canal calculando z-score para N pares fase×amplitude.

    CORREÇÃO DO CACHE: teta_ok é calculado por janela via teta_ok_por_janela()
    — nunca reutilizado de uma janela para a próxima.

    pares: dict {nome: {"fase": (lo, hi), "amp": (lo, hi)}} ou None → theta_gamma.

    Colunas de saída: janela_ini_s, janela_fim_s, teta_ok,
      mi_<par>, z_<par>, p_<par> (para cada par),
      ratio_hfo_gamma (se theta_hfo presente),
      transiente_detectado, frac_transiente, max_diff_z, max_amp_z,
      mvl, correlacao_ruido, suspeito_banda_larga, proxy_artefato_motor.
    """
    if pares is None:
        pares = {"theta_gamma": BAND_PAIRS["theta_gamma"]}

    win    = int(window_s * fs)
    step   = int(step_s   * fs)
    inicios = list(range(0, len(sinal) - win + 1, step))
    total   = len(inicios)
    resultados = []

    calcula_hfo_ratio = "theta_hfo" in pares

    for n_janela, ini in enumerate(inicios):
        fim    = ini + win
        trecho = sinal[ini:fim].astype(float)

        # ── CORREÇÃO DO CACHE: teta_ok calculado AQUI, por janela ──
        t_ok = int(teta_ok_por_janela(trecho, fs))

        row = {
            "janela_ini_s": round(ini / fs, 3),
            "janela_fim_s": round(fim / fs, 3),
            "teta_ok":      t_ok,
        }

        # ── Calcula MI z-score para cada par ──
        for nome_par, cfg in pares.items():
            fase_band = cfg["fase"]
            amp_band  = cfg["amp"]
            nyq = 0.5 * fs
            amp_hi = min(amp_band[1], nyq * 0.98)
            if amp_band[0] >= amp_hi:
                # Banda não suportada por esta fs
                row[f"mi_{nome_par}"]  = float("nan")
                row[f"z_{nome_par}"]   = float("nan")
                row[f"p_{nome_par}"]   = float("nan")
                continue

            lfp_fase = filtra_sinal(trecho, *fase_band, fs)
            lfp_amp  = filtra_sinal(trecho, amp_band[0], amp_hi, fs)
            fase = np.angle(signal.hilbert(lfp_fase))
            env  = np.abs(signal.hilbert(lfp_amp))

            mi_obs, z, p_emp, mi_surr_m, mi_surr_dp = mi_com_surrogates(
                fase, env, fs, n_surr=n_surr, rng=rng)

            row[f"mi_{nome_par}"]  = round(mi_obs, 6)
            row[f"z_{nome_par}"]   = round(z, 4)
            row[f"p_{nome_par}"]   = round(p_emp, 4)

        # ── Proxy de harmônico HFO/Gamma ──
        if calcula_hfo_ratio:
            row["ratio_hfo_gamma"] = round(_ratio_hfo_gamma(trecho, fs), 4)

        # ── Proxy artefato motor (legado) ──
        row["proxy_artefato_motor"] = round(proxy_artefato_motor(trecho, fs), 6)

        # ── Camadas anti-falso-positivo ──
        trans = detecta_transiente(trecho, fs)
        row["transiente_detectado"] = trans["transiente_encontrado"]
        row["frac_transiente"]      = round(trans["frac_transiente"], 4)
        row["max_diff_z"]           = round(trans["max_diff_z"], 2)
        row["max_amp_z"]            = round(trans["max_amp_z"], 2)

        banda_info = correlacao_gama_ruido(trecho, fs)
        row["mvl"]                  = round(banda_info["mvl"], 4)
        row["correlacao_ruido"]     = round(banda_info["correlacao_ruido"], 3)
        row["suspeito_banda_larga"] = banda_info["suspeito_banda_larga"]

        resultados.append(row)

        if rotulo_progresso and (n_janela % 10 == 0 or n_janela == total - 1):
            print(f"  {rotulo_progresso}: janela {n_janela + 1}/{total}", end="\r")

    if rotulo_progresso:
        print()

    return pd.DataFrame(resultados)


# ==========================================
# LEITURA DE .ns2 VIA NEO (Blackrock)
# ==========================================

from ns2_utils import le_ns2  # leitura compartilhada


def varre_arquivo(caminho, canais=None, window_s=10.0, step_s=5.0,
                   n_surr=200, seed=None, pares=None):
    rng = np.random.default_rng(seed)
    dados, fs, nomes_canais = le_ns2(caminho)
    n_canais_total = dados.shape[1]

    if canais is None:
        canais_idx = range(n_canais_total)
    else:
        canais_idx = canais

    todos = []
    for ch in canais_idx:
        sinal = dados[:, ch]
        nome_canal = nomes_canais[ch] if ch < len(nomes_canais) else ch
        rotulo = f"{os.path.basename(caminho)} | canal {nome_canal}"
        df = varre_canal(sinal, fs, window_s, step_s, n_surr,
                          pares=pares, rng=rng, rotulo_progresso=rotulo)
        df.insert(0, "arquivo", os.path.basename(caminho))
        df.insert(1, "canal",   nome_canal)
        todos.append(df)

    return pd.concat(todos, ignore_index=True) if todos else pd.DataFrame()


# ==========================================
# MODO DEMO
# ==========================================

def roda_demo(pares=None):
    print("=== DEMO: validando critério de flagging em dados sintéticos ===\n")
    fs  = 1000.0
    dur = 60.0
    t   = np.arange(0, dur, 1 / fs)
    rng = np.random.default_rng(0)

    # Theta com jitter biológico
    freq_inst = 6.0 + np.cumsum(rng.normal(0, 0.002, len(t)))
    freq_inst = np.clip(freq_inst, 5.0, 7.0)
    fase_theta_acum = 2 * np.pi * np.cumsum(freq_inst) / fs
    theta = np.cos(fase_theta_acum)
    modulador = (-theta + 1) / 2.0

    # Parte 1: theta puro (sem gamma)
    # Parte 2: theta + gamma acoplado (30-80 Hz)
    gamma_acoplado = modulador * np.sin(2 * np.pi * 60.0 * t)
    hg_acoplado    = modulador * 0.3 * np.sin(2 * np.pi * 100.0 * t)

    meio = len(t) // 2
    sinal = np.zeros_like(t)
    sinal[:meio]  = theta[:meio] + rng.normal(0, 0.4, meio)
    sinal[meio:]  = (theta[meio:] + gamma_acoplado[meio:]
                     + hg_acoplado[meio:] + rng.normal(0, 0.4, len(t) - meio))

    if pares is None:
        pares = {"theta_gamma": BAND_PAIRS["theta_gamma"],
                 "theta_hg":    BAND_PAIRS["theta_hg"]}

    df = varre_canal(sinal, fs, window_s=10.0, step_s=5.0, n_surr=200,
                      pares=pares, rng=rng)
    df["esperado"] = np.where(df["janela_ini_s"] < dur / 2,
                               "ruído (sem acoplamento)", "acoplamento real")

    pd.set_option("display.width", 200)
    cols_z  = [f"z_{p}" for p in pares if f"z_{p}" in df.columns]
    cols_ok = ["janela_ini_s", "teta_ok"] + cols_z + ["transiente_detectado",
               "suspeito_banda_larga", "mvl", "esperado"]
    print(df[cols_ok].to_string(index=False))

    print("\nCheck: janelas 'ruído' → z baixo (~0-2); 'acoplamento real' → z alto (>>3).")
    if "ratio_hfo_gamma" in df.columns:
        print(f"ratio_hfo_gamma médio: {df['ratio_hfo_gamma'].mean():.3f}")
    print("teta_ok deve VARIAR (não ser constante em 58 janelas de sinal real).")


# ==========================================
# MAIN
# ==========================================

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pasta",  help="Pasta contendo arquivos .ns2")
    ap.add_argument("--canal",  type=int, default=None,
                     help="Índice de um canal específico (default: todos)")
    ap.add_argument("--janela", type=float, default=10.0, help="Tamanho da janela (s)")
    ap.add_argument("--passo",  type=float, default=5.0,  help="Passo entre janelas (s)")
    ap.add_argument("--n_surr", type=int,   default=200,  help="Número de surrogates")
    ap.add_argument("--pares",  nargs="+",
                     default=["theta_gamma"],
                     choices=list(BAND_PAIRS.keys()),
                     help="Pares a calcular (default: theta_gamma). "
                          "Ex: --pares theta_gamma theta_hg theta_hfo")
    ap.add_argument("--z_corte", type=float, default=3.0,
                     help="z-score mínimo para candidato (default 3.0)")
    ap.add_argument("--saida",  default="resultados_triagem.csv")
    ap.add_argument("--demo",   action="store_true",
                     help="Roda autoteste sintético")
    ap.add_argument("--limiar_teta_snr", type=float, default=1.5,
                     help="SNR mínimo para teta_ok=1 (default 1.5)")
    ap.add_argument("--limiar_transiente_diff", type=float, default=8.0)
    ap.add_argument("--limiar_transiente_amp",  type=float, default=8.0)
    ap.add_argument("--limiar_corr_ruido",      type=float, default=0.6)
    args = ap.parse_args()

    # Constrói dict de pares selecionados
    pares_sel = {p: BAND_PAIRS[p] for p in args.pares}

    if args.demo:
        roda_demo(pares=pares_sel)
        return

    if not args.pasta:
        print("Erro: forneça --pasta com os arquivos .ns2, ou use --demo.")
        sys.exit(1)

    arquivos = sorted(glob.glob(os.path.join(args.pasta, "*.ns2")))
    if not arquivos:
        print(f"Nenhum .ns2 encontrado em {args.pasta}")
        sys.exit(1)

    canais = [args.canal] if args.canal is not None else None

    print(f"Pares calculados: {list(pares_sel.keys())}")
    print(f"Janela {args.janela}s / passo {args.passo}s / {args.n_surr} surrogates")

    todos_resultados = []
    for i, arq in enumerate(arquivos):
        print(f"[{i+1}/{len(arquivos)}] Processando {os.path.basename(arq)} ...")
        try:
            df = varre_arquivo(arq, canais=canais, window_s=args.janela,
                                step_s=args.passo, n_surr=args.n_surr,
                                seed=i, pares=pares_sel)
            todos_resultados.append(df)
        except Exception as e:
            print(f"  -> ERRO ao processar {arq}: {e}")

    if not todos_resultados:
        print("Nenhum resultado gerado.")
        sys.exit(1)

    resultado_final = pd.concat(todos_resultados, ignore_index=True)

    # Ordena pelo z do primeiro par
    primeiro_par = list(pares_sel.keys())[0]
    z_col_ord = f"z_{primeiro_par}"
    if z_col_ord in resultado_final.columns:
        resultado_final = resultado_final.sort_values(z_col_ord, ascending=False)

    resultado_final.to_csv(args.saida, index=False)

    print(f"\nSalvo: {args.saida}")
    print(f"Total de janelas: {len(resultado_final)}")

    for par in pares_sel:
        z_col = f"z_{par}"
        if z_col in resultado_final.columns:
            n_cand = int((resultado_final[z_col] >= args.z_corte).sum())
            n_teta = int(((resultado_final[z_col] >= args.z_corte)
                 & resultado_final.get("teta_ok", pd.Series(True, index=resultado_final.index)).astype(bool)
                 ).sum()) if "teta_ok" in resultado_final.columns else n_cand
            print(f"  {par}: {n_cand} candidatos com z >= {args.z_corte} "
                  f"({n_teta} com teta_ok=1)")

    if "ratio_hfo_gamma" in resultado_final.columns:
        n_suspeito = int((resultado_final["ratio_hfo_gamma"] > 0.3).sum())
        print(f"\n  ratio_hfo_gamma > 0.3 (suspeita de harmônico Gamma→HFO): "
              f"{n_suspeito} janelas")
        print("  → Rode audita_harmonico_hfo.py sobre esses candidatos.")


if __name__ == "__main__":
    main()