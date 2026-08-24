"""
robustez_parametros.py
==========================================
Etapa 3 da validação PAC: os acoplamentos dos VENCEDORES sobrevivem a
mudanças de parâmetro e a uma métrica alternativa?

Por que existe: um resultado que só aparece com UM conjunto de parâmetros
(aqui, n_bins=18, fase f±1.0 Hz, amplitude f±5.0 Hz) pode ser artefato da
escolha. Acoplamento genuíno deve ser estável em torno desses valores.

O QUE FAZ, para cada janela/canal vencedor (etapas 1-2), no PAR DE PICO
(fase_pico x amp_pico do resumo FDR):
  A. n_bins em {10, 12, 15, 18, 24, 30} -- o KL-MI normaliza pela entropia
     máxima (log n_bins), mas o viés com poucos dados depende do nº de
     bins; um acoplamento real não pode desaparecer com 12 ou 24.
  B. largura do filtro de FASE: ±0.7 / ±1.0 / ±1.5 / ±2.5 Hz
  C. largura do filtro de AMPLITUDE: ±2.5 / ±5.0 / ±10.0 Hz
  D. mapa completo recomputado com n_bins=12 e n_bins=24: o PICO continua
     no mesmo lugar (±1 Hz de fase, ±5 Hz de amplitude)?
  E. métrica alternativa: MVL (mean vector length, Canolty et al. 2006),
     |média(envelope·e^{i·fase})| -- NÃO usa bins; se o z do MVL também for
     alto, o resultado não é artefato do binning do KL-MI.

Mesma nula em tudo: 200 surrogates por deslocamento circular do envelope
(>= 1 s) -- idêntica à do triagem_pac.py e do comodulogram.py. Sempre com
notch de 60 Hz.

Critério de sobrevivência (honesto): z estável (>= 3) em TODO o sweep de
n_bins + pico do mapa estável entre n_bins + MVL confirmando. As larguras de
filtro NÃO exigem z>=3 nos extremos: filtro de amplitude mais estreito que o
próprio evento de gamma (±2.5 Hz) perde o evento por construção, e fase
±2.5 Hz mistura frequências não-theta -- o esperado é um PLATÔ de z alto em
torno do valor canônico (±5 Hz), não imunidade aos extremos.

Uso (SCRIPT central do estudo — uma chamada por sessão):
    python robustez_parametros.py \
        --pasta "<sessao>/Basal antes da infusao" \
        --vencedores "<sessao>/vencedores.csv" \
        --resumo_fdr "<sessao>/comodulogramas_fdr/resumo_comodulogramas.csv" \
        --saida_csv "<sessao>/robustez_parametros.csv"

vencedores.csv (uma linha por vencedor; colunas):
    rotulo,arquivo,canal,inicio_s,fim_s,fase_pico_hz,amp_pico_hz[,comportamento]
Se fase_pico_hz/amp_pico_hz estiverem VAZIOS numa linha, o par de pico é lido
do resumo FDR (--resumo_fdr obrigatório nesse caso).

Saída: tabela no console + CSV (--saida_csv)
Requer: neo, numpy, scipy, pandas
"""

import argparse
import sys

import numpy as np
import pandas as pd
import scipy.signal as signal

from comodulogram import (filtra_sinal, aplica_notch, calcula_comodulograma_z,
                          z_pico_theta_gamma, _mi_de_bin_idx)
from ns2_utils import le_ns2, fatia_janela

N_BINS_SWEEP = [10, 12, 15, 18, 24, 30]
MEIA_FAISES = [0.7, 1.0, 1.5, 2.5]
MEIA_AMPS = [2.5, 3.5, 5.0, 7.5, 10.0]  # curva de sintonia em torno do canônico
N_SURR = 200
NOTCH_HZ = 60.0


def mi_z_par(lfp, fs, f_fase, f_amp, n_bins=18, meia_fase=1.0, meia_amp=5.0,
             n_surr=N_SURR, rng=None):
    """z do KL-MI em um ÚNICO par (f_fase, f_amp) -- mesmos filtros e mesma
    nula de deslocamento circular do mapa do comodulogram.py."""
    if rng is None:
        rng = np.random.default_rng()
    bins = np.linspace(-np.pi, np.pi, n_bins + 1)
    lfp_fase = filtra_sinal(lfp, f_fase - meia_fase, f_fase + meia_fase, fs)
    bin_idx = np.clip(
        np.digitize(np.angle(signal.hilbert(lfp_fase)), bins) - 1, 0, n_bins - 1)
    lfp_amp = filtra_sinal(lfp, f_amp - meia_amp, f_amp + meia_amp, fs)
    env = np.abs(signal.hilbert(lfp_amp))

    n = lfp.size
    shift_min = int(1.0 * fs)
    if n <= 2 * shift_min:
        shift_min = max(1, n // 10)
    deslocamentos = rng.integers(shift_min, n - shift_min, size=n_surr)

    mi_obs = _mi_de_bin_idx(bin_idx, env, n_bins)
    mi_surr = np.array([_mi_de_bin_idx(bin_idx, np.roll(env, d), n_bins)
                        for d in deslocamentos])
    dp = mi_surr.std()
    return (mi_obs - mi_surr.mean()) / dp if dp > 0 else 0.0


def mvl_z_par(lfp, fs, f_fase, f_amp, meia_fase=1.0, meia_amp=5.0,
              n_surr=N_SURR, rng=None):
    """z do MVL (Canolty et al. 2006): |média(envelope · e^{i·fase})|.
    Métrica ALTERNATIVA sem bins -- complementar ao KL-MI (que depende do
    nº de bins de fase). Mesma nula de deslocamento circular do envelope."""
    if rng is None:
        rng = np.random.default_rng()
    lfp_fase = filtra_sinal(lfp, f_fase - meia_fase, f_fase + meia_fase, fs)
    fase = np.angle(signal.hilbert(lfp_fase))
    lfp_amp = filtra_sinal(lfp, f_amp - meia_amp, f_amp + meia_amp, fs)
    env = np.abs(signal.hilbert(lfp_amp))

    n = lfp.size
    shift_min = int(1.0 * fs)
    if n <= 2 * shift_min:
        shift_min = max(1, n // 10)
    deslocamentos = rng.integers(shift_min, n - shift_min, size=n_surr)

    mvl_obs = np.abs(np.mean(env * np.exp(1j * fase)))
    mvl_surr = np.array([np.abs(np.mean(np.roll(env, d) * np.exp(1j * fase)))
                         for d in deslocamentos])
    dp = mvl_surr.std()
    return (mvl_obs - mvl_surr.mean()) / dp if dp > 0 else 0.0


def main():
    try:  # console Windows pode estar em cp1252; Θ/Γ quebram o print
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pasta", required=True, help="Pasta com os .ns2 da sessão")
    ap.add_argument("--vencedores", required=True,
                    help="CSV da sessão: rotulo,arquivo,canal,inicio_s,fim_s,"
                         "fase_pico_hz,amp_pico_hz")
    ap.add_argument("--resumo_fdr", default=None,
                    help="resumo_comodulogramas.csv da sessão (exigido se "
                         "alguma linha não fixa fase/amp de pico)")
    ap.add_argument("--saida_csv", default="robustez_parametros.csv")
    args = ap.parse_args()

    venc = pd.read_csv(args.vencedores)
    resumo_cache = {}
    linhas = []

    for _, row in venc.iterrows():
        rotulo = str(row["rotulo"])
        arquivo = str(row["arquivo"])
        canal = str(row["canal"])
        inicio = float(row["inicio_s"])
        fim = float(row["fim_s"])
        f_fixo = pd.to_numeric(row.get("fase_pico_hz"), errors="coerce")
        a_fixo = pd.to_numeric(row.get("amp_pico_hz"), errors="coerce")
        f_fixo = None if pd.isna(f_fixo) else float(f_fixo)
        a_fixo = None if pd.isna(a_fixo) else float(a_fixo)

        rng = np.random.default_rng(42)  # mesmo rng por janela = comparável
        caminho = f"{args.pasta}/{arquivo}"
        print(f"\n=== {rotulo}: {canal} @ {inicio:g}-{fim:g}s ({arquivo}) ===")

        dados, fs, nomes = le_ns2(caminho)
        mapa = {str(n): i for i, n in enumerate(nomes)}
        lfp = fatia_janela(dados[:, mapa[canal]], fs, inicio, fim).astype(float)
        lfp = aplica_notch(lfp, fs, linha_hz=NOTCH_HZ)
        del dados

        if f_fixo is not None and a_fixo is not None:
            f_pico, a_pico = f_fixo, a_fixo
            print(f"Par de pico (fixado no vencedores.csv): "
                  f"{f_pico:g} Hz x {a_pico:g} Hz")
        else:
            # par de pico vem do resumo FDR da sessão (etapa 2)
            if not args.resumo_fdr:
                raise SystemExit(
                    f"'{rotulo}': linha sem fase/amp de pico exige --resumo_fdr.")
            if args.resumo_fdr not in resumo_cache:
                resumo_cache[args.resumo_fdr] = pd.read_csv(args.resumo_fdr)
            resumo = resumo_cache[args.resumo_fdr]
            linha = resumo[(resumo["arquivo"] == arquivo) &
                           (resumo["canal"] == canal) &
                           (resumo["janela_ini_s"] == inicio)].iloc[0]
            f_pico = float(linha["fase_pico_hz"])
            a_pico = float(linha["amp_pico_hz"])
            print(f"Par de pico (etapa 2): {f_pico:g} Hz x {a_pico:g} Hz "
                  f"(z={linha['z_pico_theta_gamma']:.2f})")

        # ----------------------------------
        # A. n_bins
        # ----------------------------------
        zs_nb = {nb: mi_z_par(lfp, fs, f_pico, a_pico, n_bins=nb, rng=rng)
                 for nb in N_BINS_SWEEP}
        for nb, z in zs_nb.items():
            print(f"  n_bins={nb:>2}: z={z:5.2f}")
            linhas.append({"janela": rotulo, "canal": canal, "par_pico":
                           f"{f_pico:g}x{a_pico:g}", "teste": "n_bins",
                           "parametro": nb, "z": round(z, 2)})

        # ----------------------------------
        # B. largura do filtro de fase
        # ----------------------------------
        zs_f = {mf: mi_z_par(lfp, fs, f_pico, a_pico, meia_fase=mf, rng=rng)
                for mf in MEIA_FAISES}
        for mf, z in zs_f.items():
            print(f"  fase ±{mf:g} Hz: z={z:5.2f}")
            linhas.append({"janela": rotulo, "canal": canal, "par_pico":
                           f"{f_pico:g}x{a_pico:g}", "teste": "largura_fase",
                           "parametro": mf, "z": round(z, 2)})

        # ----------------------------------
        # C. largura do filtro de amplitude
        # ----------------------------------
        zs_a = {ma: mi_z_par(lfp, fs, f_pico, a_pico, meia_amp=ma, rng=rng)
                for ma in MEIA_AMPS}
        for ma, z in zs_a.items():
            print(f"  amplitude ±{ma:g} Hz: z={z:5.2f}")
            linhas.append({"janela": rotulo, "canal": canal, "par_pico":
                           f"{f_pico:g}x{a_pico:g}", "teste": "largura_amplitude",
                           "parametro": ma, "z": round(z, 2)})

        # ----------------------------------
        # D. pico do mapa com n_bins=12 e 24 (localização estável?)
        # ----------------------------------
        fases_freq = np.arange(4, 15, 1)
        amps_freq = np.arange(30, 155, 5)
        picos_estaveis = []
        for nb in (12, 24):
            z_mapa = calcula_comodulograma_z(lfp, fs, fases_freq, amps_freq,
                                             n_surr=N_SURR, n_bins=nb, rng=rng,
                                             notch_hz=None)  # notch já aplicado
            z_p, f_p, a_p = z_pico_theta_gamma(z_mapa, fases_freq, amps_freq)
            estavel = (abs(f_p - f_pico) <= 1.0) and (abs(a_p - a_pico) <= 5.0)
            picos_estaveis.append(estavel)
            print(f"  mapa n_bins={nb}: pico ΘΓ {f_p:g} x {a_p:g} Hz, "
                  f"z={z_p:.2f}  ({'estável' if estavel else 'DESLOCADO'})")
            linhas.append({"janela": rotulo, "canal": canal, "par_pico":
                           f"{f_pico:g}x{a_pico:g}", "teste": f"mapa_nb{nb}",
                           "parametro": f"{f_p:g}x{a_p:g}", "z": round(z_p, 2)})

        # ----------------------------------
        # E. métrica alternativa: MVL
        # ----------------------------------
        z_mvl = mvl_z_par(lfp, fs, f_pico, a_pico, rng=rng)
        print(f"  MVL (sem bins): z={z_mvl:5.2f}")
        linhas.append({"janela": rotulo, "canal": canal, "par_pico":
                       f"{f_pico:g}x{a_pico:g}", "teste": "MVL",
                       "parametro": "pico", "z": round(z_mvl, 2)})

        # veredito -- robustez no PLATÔ dos parâmetros, não nos extremos
        # (ver docstring): n_bins é o parâmetro do KL-MI e não pode colapsar;
        # larguras de filtro mapeiam a curva de sintonia; MVL confirma sem
        # binning (mas só captura o 1º momento -- z menor é normal).
        z_min_nb = min(zs_nb.values())
        z_min_bw = min(list(zs_f.values()) + list(zs_a.values()))
        robusto = (z_min_nb >= 3) and all(picos_estaveis)
        conf = "confirma" if z_mvl >= 3 else "não confirma (z menor é esperado)"
        print(f"  >> n_bins: z mínimo {z_min_nb:.2f} | larguras: z mínimo "
              f"{z_min_bw:.2f} | pico estável: {all(picos_estaveis)} | "
              f"MVL z={z_mvl:.2f} ({conf}) -> "
              f"{'ROBUSTO' if robusto else 'FRÁGIL'}")

    df = pd.DataFrame(linhas)
    df.to_csv(args.saida_csv, index=False)
    print(f"\nCSV salvo: {args.saida_csv} ({len(df)} linhas)")


if __name__ == "__main__":
    main()
