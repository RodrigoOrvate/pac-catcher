"""
compara_preprocesso_linha.py - Comparador das 3 configs de limpeza de linha
                               (gaussiana / cirurgica / hibrido) vs. sem limpeza.

REUSO (principio do projeto): delega a extracao FOOOF a `extrai_cf_teta_fooof`
(de producao, audita_harmonico.py) e a limpeza de linha a `aplica_modo`
(linha_noise_kuhn.py). O gerador sintetico `gera_cenario_D` vem de
test_fooof_pico_budget.py. Nada de logica duplicada.

Objetivo (plano aprovado):
  - Sintetico EU  : cenario D (linha 50/100 Hz), f_linha=50
  - Sintetico BR  : cenario D com linha 60/120 Hz (caso Brasil), f_linha=60
  - Real .mat     : CA1_example.mat / DG_example.mat (tetrode 4 ch, 1000 Hz),
                    eletrodo de maior potencia (logica maxpower), f_linha
                    detectada no PSD (pico 45-65 Hz)
Para cada cenario: tabela erro_ajuste / cf_teta / n_picos por modo
{sem, gaussiana, cirurgica, hibrido} -> stdout + CSV. Figuras de PSD
comparativo por arquivo .mat (original + 3 limpos) para inspecao do
fast gamma 60-100 Hz e do knee (CA1 ~28 / DG ~70).

    python pipeline/auditorias/compara_preprocesso_linha.py
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.signal import welch
from scipy.io import loadmat
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # pipeline/auditorias/
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # pipeline/
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # SCRIPT/

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

from audita_harmonico import extrai_cf_teta_fooof
from linha_noise_kuhn import aplica_modo
from test_fooof_pico_budget import gera_cenario_D, generate_harmonic_signal

MODOS = ["sem", "gaussiana", "cirurgica", "hibrido"]
NPERSEG_S = 1.2


# ---------------------------------------------------------------------------
# Sinteticos
# ---------------------------------------------------------------------------
def gera_cenario_D_BR(fs, dur_s, f_theta=8.0):
    """Cenario D com linha BR (60/120 Hz) em vez de EU (50/100 Hz).

    Mesma estrutura do cenario D (teta nao-senoidal + linha + harmonico +
    artefato motor em 45 Hz), so que com contaminantes na rede brasileira.
    """
    from scipy.signal import butter, filtfilt
    sinal, fs_r, f_th, f_g = generate_harmonic_signal(
        fs, dur_s, f_theta=f_theta, ordem=3,
        snr_db=20.0, aperiodic=True, seed=42)
    t = np.arange(0, dur_s, 1.0 / fs)
    residuo_60 = 0.5 * np.sin(2 * np.pi * 60.0 * t)
    residuo_120 = 0.2 * np.sin(2 * np.pi * 120.0 * t)
    rng = np.random.default_rng(43)
    ruido_motor = rng.standard_normal(len(t))
    b, a = butter(4, [(45 - 5) / 500, (45 + 5) / 500], btype="band")
    ruido_motor = 0.3 * filtfilt(b, a, ruido_motor)
    return (sinal + residuo_60 + residuo_120 + ruido_motor,
            fs_r, f_th, f_g)


def roda_todos_modos(sinal, fs, f_linha):
    """Roda extrai_cf_teta_fooof nos 4 modos e devolve lista de dicts."""
    linhas = []
    for modo in MODOS:
        if modo == "sem":
            r = extrai_cf_teta_fooof(sinal, fs, preprocessar_linha=False)
        else:
            r = extrai_cf_teta_fooof(sinal, fs, preprocessar_linha=True,
                                     modo_preprocesso=modo, f_linha=f_linha)
        linhas.append({
            "modo": modo,
            "erro": r["erro_ajuste"],
            "cf_teta": r["cf_teta"],
            "n_picos": r["n_picos"],
            "qualidade": r["qualidade_ok"],
        })
    return linhas


def roda_sintetico(nome, gen, f_linha):
    fs, dur_s = 2000, 45.0
    sinal, fs_r, _, _ = gen(fs, dur_s)
    return nome, roda_todos_modos(sinal, fs_r, f_linha)


# ---------------------------------------------------------------------------
# Real .mat
# ---------------------------------------------------------------------------
def carrega_mat(path):
    """Carrega o .mat de Kuhn e devolve (lfp, fs, canal_maxpower).

    Usa o campo `lfp` (int16 cru, sem NaN), NAO o `mlfp`: o multi-unit
    filter de Kuhn zera janelas de spike com NaN, o que quebra o Welch/FOOOF.
    Logica maxpower: eletrodo com maior potencia media em [4, 200] Hz.
    """
    d = loadmat(path)
    lfp = d["lfp"].astype(float)
    fs = 1000
    potencias = []
    for ch in range(lfp.shape[1]):
        f, p = welch(lfp[:, ch], fs=fs, nperseg=int(NPERSEG_S * fs),
                     noverlap=int(NPERSEG_S * fs) // 2, nfft=4000)
        m = (f >= 4) & (f <= 200)
        potencias.append(np.log10(np.mean(p[m])))
    ch_max = int(np.argmax(potencias))
    return lfp, fs, ch_max


def detecta_f_linha(freqs, psd, faixa=(45, 65)):
    """Detecta a frequencia de linha real (pico na faixa 45-65 Hz)."""
    m = (freqs >= faixa[0]) & (freqs <= faixa[1])
    if m.sum() == 0:
        return None
    return float(freqs[m][np.argmax(psd[m])])


def roda_real(nome, path, outdir_figuras):
    mlfp, fs, ch = carrega_mat(path)
    sinal = mlfp[:, ch]

    # PSD original (para detectar linha e plotar)
    f_orig, p_orig = welch(sinal, fs=fs, nperseg=int(NPERSEG_S * fs),
                           noverlap=int(NPERSEG_S * fs) // 2, nfft=4000)
    f_linha = detecta_f_linha(f_orig, p_orig)
    if f_linha is None:
        f_linha = 60.0

    res = roda_todos_modos(sinal, fs, f_linha)

    # PSD de cada modo limpo (para a figura)
    psds = {"sem": p_orig}
    for modo in ["gaussiana", "cirurgica", "hibrido"]:
        psds[modo] = aplica_modo(f_orig, p_orig, fs, modo, f_linha=f_linha,
                                 verbose=False)

    figura_psd(nome, f_orig, psds, f_linha, outdir_figuras)
    return nome, ch, f_linha, res, (f_orig, psds)


def figura_psd(nome, freqs, psds, f_linha, outdir):
    os.makedirs(outdir, exist_ok=True)
    plt.figure(figsize=(12, 6))
    cores = {"sem": "k", "gaussiana": "tab:blue", "cirurgica": "tab:green",
             "hibrido": "tab:red"}
    rotulos = {"sem": "original", "gaussiana": "gaussiana",
               "cirurgica": "cirurgica", "hibrido": "hibrido"}
    for modo, p in psds.items():
        plt.semilogy(freqs, p, color=cores[modo], lw=1.2, label=rotulos[modo])
    plt.axvline(f_linha, color="gray", ls="--", lw=0.8,
                label=f"linha {f_linha:.0f} Hz")
    plt.axvspan(60, 100, color="orange", alpha=0.12, label="fast gamma 60-100")
    plt.xlabel("Frequencia (Hz)")
    plt.ylabel("PSD")
    plt.title(f"{nome} - PSD por modo de limpeza (ch. maior potencia, "
              f"f_linha={f_linha:.0f} Hz)")
    plt.legend()
    plt.xlim(0, 200)
    plt.tight_layout()
    path = os.path.join(outdir, f"comparacao_preprocesso_{nome}_psd.png")
    plt.savefig(path, dpi=110)
    plt.close()
    print(f"  Figura salva: {path}")


# ---------------------------------------------------------------------------
# Relatorio
# ---------------------------------------------------------------------------
def fmt_erro(e):
    return f"{e:.4f}" if e is not None else "FAIL"


def fmt_cf(c):
    return f"{c:.3f}" if c is not None else "n/a"


def imprime_tabela(nome, res, f_linha=None, extra=""):
    print(f"\n=== {nome} {extra} ===")
    print(f"{'modo':<10} | {'erro':<8} | {'cf_teta':<8} | {'n_picos':<7} | qual")
    print("-" * 55)
    for r in res:
        print(f"{r['modo']:<10} | {fmt_erro(r['erro']):<8} | "
              f"{fmt_cf(r['cf_teta']):<8} | {r['n_picos']:<7} | "
              f"{'OK' if r['qualidade'] else '--'}")
    return res


def decide_default(linhas_sint, linhas_real):
    """Sugere o modo default: menor erro medio sem mutilar o real.

    Prioridade: sintetico EU+BR (onde conhecemos a verdade) para derrubar o
    erro < 0.15; o real .mat serve para confirmar que o fast gamma e o knee
    nao sao mutilados (inspecao visual da figura) - nao ha ground truth de
    erro no real, apenas cf_teta plausivel.
    """
    melhor = None
    melhor_media = float("inf")
    for modo in MODOS[1:]:
        errs = []
        for _, res in linhas_sint:
            for r in res:
                if r["modo"] == modo and r["erro"] is not None:
                    errs.append(r["erro"])
        if not errs:
            continue
        media = float(np.mean(errs))
        if media < melhor_media:
            melhor_media = media
            melhor = modo
    return melhor, melhor_media


def main():
    outdir_figuras = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                  "figuras")
    csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "comparacao_preprocesso_linha.csv")

    print("=" * 80)
    print("COMPARADOR DE LIMPEZA DE LINHA (gaussiana / cirurgica / hibrido)")
    print("=" * 80)

    # 1. Sinteticos
    linhas_sint = []
    for nome, gen, fl in [("Sintetico_EU", gera_cenario_D, 50.0),
                          ("Sintetico_BR", gera_cenario_D_BR, 60.0)]:
        nome_c, res = roda_sintetico(nome, gen, fl)
        imprime_tabela(nome_c, res, f_linha=fl)
        linhas_sint.append((nome_c, res))

    # 2. Real .mat
    base = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))), "FOOOF")
    linhas_real = []
    for nome, arq in [("CA1", "CA1_example.mat"), ("DG", "DG_example.mat")]:
        path = os.path.join(base, arq)
        if not os.path.exists(path):
            print(f"\n[AVISO] {path} nao encontrado. Pulando {nome}.")
            continue
        nome_c, ch, fl, res, _ = roda_real(nome, path, outdir_figuras)
        imprime_tabela(nome_c, res, f_linha=fl,
                       extra=f"(ch {ch}, f_linha={fl:.1f} Hz)")
        linhas_real.append((nome_c, res))

    # 3. CSV (todos os cenarios)
    registros = []
    for nome, res in linhas_sint + linhas_real:
        for r in res:
            registros.append({"cenario": nome, **r})
    pd.DataFrame(registros).to_csv(csv_path, index=False)
    print(f"\nCSV salvo: {csv_path}")

    # 4. Resumo de decisao
    print("\n" + "=" * 80)
    print("RESUMO DE DECISAO (default de producao)")
    print("=" * 80)
    if linhas_sint:
        melhor, media = decide_default(linhas_sint, linhas_real)
        print(f"Menor erro medio no sintetico (EU+BR): '{melhor}' "
              f"(erro medio {media:.4f})")
        print("Confirmar visualmente nas figuras .mat que o fast gamma 60-100")
        print("e o knee (CA1 ~28 / DG ~70) sobrevivem no modo escolhido.")
    else:
        print("Sem dados sinteticos - nada a decidir.")


if __name__ == "__main__":
    main()
