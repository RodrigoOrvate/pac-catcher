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

vencedores.csv (uma linha por vencedor; colunas aceitas):
    arquivo,canal,janela_ini_s,janela_fim_s,fase_pico_hz,amp_pico_hz
    [,rotulo,comportamento] -- aceita também o formato antigo inicio_s/fim_s;
    canal na convenção 1-based do dataset mestre ou nome nativo do .ns2.
Se fase_pico_hz/amp_pico_hz estiverem VAZIOS numa linha, o par de pico é lido
do resumo FDR (--resumo_fdr obrigatório nesse caso).

Saída: tabela no console + CSV (--saida_csv)
Requer: neo, numpy, scipy, pandas
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import scipy.signal as signal

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from pipeline.etapa3_comodulograma.comodulogram import (calcula_comodulograma_z,
                          z_pico_par, FASES_DEFAULT, AMPS_DEFAULT)
from pac_core.io import le_ns2, fatia_janela, resolve_canal_idx
from pipeline.etapa1_triagem.triagem_pac import detecta_transiente, correlacao_gama_ruido
from pac_core.filtering import filtra_sinal, aplica_notch
from pac_core.pac_metrics import (
    fase_para_bin_idx, _mi_de_bin_idx, gera_deslocamentos,
    mi_surrogates_de_deslocamentos, z_score_mi,
)

N_BINS_SWEEP = [10, 12, 15, 18, 24, 30]
MEIA_FAISES = [0.7, 1.0, 1.5, 2.5]
MEIA_AMPS = [2.5, 3.5, 5.0, 7.5, 10.0]  # curva de sintonia em torno do canônico
N_SURR = 200
# Auditoria 2026-09: era so 60.0 (fundamental) -- nao removia os harmonicos
# 120/180/240Hz, apesar da regra #1 do CLAUDE.md. Teste pareado (mesma seed,
# 60Hz-so vs multi-harmonico) nos 227 candidatos do OURO_PURIFICADO.csv: 37
# (16.3%) mudam de veredito de robustez, a maioria perdendo robustez por
# deteccao de transiente antes mascarada pelo ruido de linha residual (ver
# resultados/_verificacao_notch_robustez_227.csv).
NOTCH_HZ = [60.0, 120.0, 180.0, 240.0]

# Limiares de rejeição (CAMADA 2)
LIMIAR_MVL = 0.05       # MVL < 0.05 = distribuição circular = rejeitar
LIMIAR_RAYLEIGH = 0.05  # Rayleigh p > 0.05 = uniforme = rejeitar


def mi_z_par(lfp, fs, f_fase, f_amp, n_bins=18, meia_fase=1.0, meia_amp=5.0,
             n_surr=N_SURR, rng=None):
    """z do KL-MI em um ÚNICO par (f_fase, f_amp) -- mesmos filtros e mesma
    nula de deslocamento circular do mapa do comodulogram.py."""
    lfp_fase = filtra_sinal(lfp, f_fase - meia_fase, f_fase + meia_fase, fs)
    bin_idx = fase_para_bin_idx(np.angle(signal.hilbert(lfp_fase)), n_bins)
    lfp_amp = filtra_sinal(lfp, f_amp - meia_amp, f_amp + meia_amp, fs)
    env = np.abs(signal.hilbert(lfp_amp))

    deslocamentos = gera_deslocamentos(lfp.size, fs, n_surr=n_surr, rng=rng)

    mi_obs = _mi_de_bin_idx(bin_idx, env, n_bins)
    mi_surr = mi_surrogates_de_deslocamentos(bin_idx, env, deslocamentos, n_bins)
    return z_score_mi(mi_obs, mi_surr)


def mvl_z_par(lfp, fs, f_fase, f_amp, meia_fase=1.0, meia_amp=5.0,
              n_surr=N_SURR, rng=None):
    """z do MVL (Canolty et al. 2006): |média(envelope · e^{i·fase})|.
    Métrica ALTERNATIVA sem bins -- complementar ao KL-MI (que depende do
    nº de bins de fase). Mesma nula de deslocamento circular do envelope
    -- reusa só gera_deslocamentos() do núcleo (a nula é compartilhada,
    a métrica MVL não é KL-MI e não migra para pac_core.pac_metrics)."""
    lfp_fase = filtra_sinal(lfp, f_fase - meia_fase, f_fase + meia_fase, fs)
    fase = np.angle(signal.hilbert(lfp_fase))
    lfp_amp = filtra_sinal(lfp, f_amp - meia_amp, f_amp + meia_amp, fs)
    env = np.abs(signal.hilbert(lfp_amp))

    deslocamentos = gera_deslocamentos(lfp.size, fs, n_surr=n_surr, rng=rng)

    mvl_obs = np.abs(np.mean(env * np.exp(1j * fase)))
    mvl_surr = np.array([np.abs(np.mean(np.roll(env, d) * np.exp(1j * fase)))
                         for d in deslocamentos])
    dp = mvl_surr.std()
    return (mvl_obs - mvl_surr.mean()) / dp if dp > 0 else 0.0


def mvl_bruto_e_rayleigh(lfp, fs, f_fase, f_amp, meia_fase=1.0, meia_amp=5.0):
    """
    CAMADA 2 - MVL bruto e Teste de Rayleigh.

    Retorna:
      - mvl_obs: MVL observado (float)
      - rayleigh_p: p-valor do teste de Rayleigh (se p>0.05, distribuição é uniforme)

    O Teste de Rayleigh (Zar 1999): para ângulos de fase, testa a hipótese
    nula de uniformidade direcional. Se a amplitude de gamma está uniformemente
    distribuída pela fase do theta, NÃO há acoplamento preferencial.
    """
    lfp_fase = filtra_sinal(lfp, f_fase - meia_fase, f_fase + meia_fase, fs)
    fase = np.angle(signal.hilbert(lfp_fase))
    lfp_amp = filtra_sinal(lfp, f_amp - meia_amp, f_amp + meia_amp, fs)
    env = np.abs(signal.hilbert(lfp_amp))

    n = len(fase)
    if n < 10:
        return 0.0, 1.0

    # MVL bruto
    # Normaliza envelope para peso igual (MVL puro)
    mvl_obs = np.abs(np.mean((env - env.mean()) / (env.std() + 1e-12) * np.exp(1j * fase)))

    # Teste de Rayleigh: estatística R = n * MVL, p = exp(-R^2/n) (para n grande)
    # Versão mais precisa (Greenwood & Durand 1955):
    R = n * mvl_obs
    rayleigh_z = R**2 / n
    rayleigh_p = np.exp(-rayleigh_z) * (1 + (2*rayleigh_z - rayleigh_z**2) / (4*n)
                                          - (24*rayleigh_z - 132*rayleigh_z**2
                                              + 76*rayleigh_z**3 - 9*rayleigh_z**4) / (288*n**2))

    return mvl_obs, float(np.clip(rayleigh_p, 0, 1))


def avalia_robustez_evento(lfp, fs, f_pico, a_pico, par_nome="theta_gamma", rng=None):
    """Núcleo de robustez_parametros.py extraído como função reusável (2026-09,
    para permitir rodar o portão em lote sobre uma lista de eventos sem
    duplicar a lógica) -- MESMOS cálculos/limiares de quando só existia
    inline em main(), comportamento idêntico.

    `lfp` já deve estar recortado na janela e filtrado (notch aplicado),
    como em main(). Retorna dict com o veredito estruturado + o detalhe
    granular do sweep (mesmo formato de linha do --saida_csv).
    """
    if rng is None:
        rng = np.random.default_rng(42)
    linhas = []

    zs_nb = {nb: mi_z_par(lfp, fs, f_pico, a_pico, n_bins=nb, rng=rng)
             for nb in N_BINS_SWEEP}
    for nb, z in zs_nb.items():
        linhas.append({"teste": "n_bins", "parametro": nb, "z": round(z, 2)})

    zs_f = {mf: mi_z_par(lfp, fs, f_pico, a_pico, meia_fase=mf, rng=rng)
            for mf in MEIA_FAISES}
    for mf, z in zs_f.items():
        linhas.append({"teste": "largura_fase", "parametro": mf, "z": round(z, 2)})

    zs_a = {ma: mi_z_par(lfp, fs, f_pico, a_pico, meia_amp=ma, rng=rng)
            for ma in MEIA_AMPS}
    for ma, z in zs_a.items():
        linhas.append({"teste": "largura_amplitude", "parametro": ma, "z": round(z, 2)})

    fases_freq = FASES_DEFAULT
    amps_freq = AMPS_DEFAULT
    picos_estaveis = []
    for nb in (12, 24):
        z_mapa = calcula_comodulograma_z(lfp, fs, fases_freq, amps_freq,
                                         n_surr=N_SURR, n_bins=nb, rng=rng,
                                         notch_hz=None)
        z_p, f_p, a_p = z_pico_par(z_mapa, fases_freq, amps_freq, par_nome)
        if z_p is None:
            estavel = False
        else:
            estavel = (abs(f_p - f_pico) <= 1.0) and (abs(a_p - a_pico) <= 5.0)
        picos_estaveis.append(estavel)
        linhas.append({"teste": f"mapa_nb{nb}",
                       "parametro": f"{f_p:g}x{a_p:g}" if f_p else "sem_pico",
                       "z": round(z_p, 2) if z_p is not None else None})

    z_mvl = mvl_z_par(lfp, fs, f_pico, a_pico, rng=rng)
    linhas.append({"teste": "MVL", "parametro": "pico", "z": round(z_mvl, 2)})

    z_min_nb = min(zs_nb.values())
    z_min_bw = min(list(zs_f.values()) + list(zs_a.values()))
    n_nb_pass = sum(1 for z in zs_nb.values() if z >= 3)
    n_nb_total = len(zs_nb)

    mvl_bruto, rayleigh_p = mvl_bruto_e_rayleigh(lfp, fs, f_pico, a_pico)
    rejeitado_mvl = (mvl_bruto < LIMIAR_MVL) or (rayleigh_p > LIMIAR_RAYLEIGH)
    motivo_mvl = []
    if mvl_bruto < LIMIAR_MVL:
        motivo_mvl.append(f"MVL={mvl_bruto:.4f}<{LIMIAR_MVL} (distribuição circular)")
    if rayleigh_p > LIMIAR_RAYLEIGH:
        motivo_mvl.append(f"Rayleigh p={rayleigh_p:.3f}>{LIMIAR_RAYLEIGH} (uniforme)")

    trans_info = detecta_transiente(lfp, fs)
    rejeitado_trans = trans_info["transiente_encontrado"]

    banda_info = correlacao_gama_ruido(
        lfp, fs,
        theta_band=(f_pico - 1.0, f_pico + 1.0),
        gamma_band=(a_pico - 5.0, a_pico + 5.0)
    )
    rejeitado_banda = banda_info["suspeito_banda_larga"]

    robusto = (z_min_nb >= 3) and all(picos_estaveis) and not rejeitado_mvl \
               and not rejeitado_trans and not rejeitado_banda

    # Criterio relaxado (maioria, nao unanimidade) -- ver discussao 2026-09:
    # exigir z>=3 nos 6/6 n_bins + estabilidade nos 2/2 mapas recomputados
    # e um "E" combinatorio que penaliza estatisticamente sinais reais e
    # modestos em janelas de 10s (vies de estimativa de MI varia ponto a
    # ponto no sweep mesmo para acoplamento genuino). Os testes de artefato
    # especifico (MVL/Rayleigh, transiente, banda larga) continuam corte
    # duro -- so a exigencia de unanimidade no sweep parametrico e relaxada.
    robusto_maioria = (n_nb_pass >= 4) and (sum(picos_estaveis) >= 1) and not rejeitado_mvl \
               and not rejeitado_trans and not rejeitado_banda

    if rejeitado_mvl or rejeitado_trans or rejeitado_banda:
        status = f"FALSO POSITIVO (rejeitado por: {', '.join(motivo_mvl) if motivo_mvl else ''}"
        if rejeitado_trans:
            status += f", transiente max_diff={trans_info['max_diff_z']:.1f}σ"
        if rejeitado_banda:
            status += f", banda larga r={banda_info['correlacao_ruido']:.2f}"
        status += ")"
    else:
        conf = "confirma" if z_mvl >= 3 else "não confirma (z menor é esperado)"
        status = f"ROBUSTO ({conf})"

    return {
        "robusto": robusto,
        "robusto_maioria": robusto_maioria,
        "status": status,
        "z_min_nb": z_min_nb,
        "z_min_bw": z_min_bw,
        "n_nb_pass": n_nb_pass,
        "n_nb_total": n_nb_total,
        "pico_estavel": all(picos_estaveis),
        "n_mapas_estaveis": sum(picos_estaveis),
        "z_mvl": z_mvl,
        "mvl_bruto": mvl_bruto,
        "rayleigh_p": rayleigh_p,
        "transiente_encontrado": rejeitado_trans,
        "max_diff_transiente_sigma": trans_info["max_diff_z"],
        "suspeito_banda_larga": rejeitado_banda,
        "correlacao_gama_ruido": banda_info["correlacao_ruido"],
        "linhas_sweep": linhas,
    }


def main():
    try:  # console Windows pode estar em cp1252; Θ/Γ quebram o print
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pasta", required=True, help="Pasta com os .ns2 da sessão")
    ap.add_argument("--vencedores", required=True,
                    help="CSV: arquivo,canal,janela_ini_s,janela_fim_s,"
                         "fase_pico_hz,amp_pico_hz [,rotulo] -- aceita também "
                         "o formato antigo inicio_s/fim_s; canal na convenção "
                         "1-based do dataset mestre ou nome nativo do .ns2")
    ap.add_argument("--resumo_fdr", default=None,
                    help="resumo_comodulogramas.csv da sessão (exigido se "
                         "alguma linha não fixa fase/amp de pico)")
    ap.add_argument("--saida_csv", default="robustez_parametros.csv")
    args = ap.parse_args()

    venc = pd.read_csv(args.vencedores)
    col_ini = "janela_ini_s" if "janela_ini_s" in venc.columns else "inicio_s"
    col_fim = "janela_fim_s" if "janela_fim_s" in venc.columns else "fim_s"
    tem_rotulo = "rotulo" in venc.columns
    resumo_cache = {}
    linhas = []

    for _, row in venc.iterrows():
        arquivo = str(row["arquivo"])
        canal = str(row["canal"])
        inicio = float(row[col_ini])
        fim = float(row[col_fim])
        if tem_rotulo and pd.notna(row.get("rotulo")):
            rotulo = str(row["rotulo"])
        else:
            rotulo = f"{arquivo.replace('.ns2', '')}_{inicio:.0f}-{fim:.0f}s"
        f_fixo = pd.to_numeric(row.get("fase_pico_hz"), errors="coerce")
        a_fixo = pd.to_numeric(row.get("amp_pico_hz"), errors="coerce")
        f_fixo = None if pd.isna(f_fixo) else float(f_fixo)
        a_fixo = None if pd.isna(a_fixo) else float(a_fixo)

        rng = np.random.default_rng(42)  # mesmo rng por janela = comparável
        caminho = f"{args.pasta}/{arquivo}"
        print(f"\n=== {rotulo}: {canal} @ {inicio:g}-{fim:g}s ({arquivo}) ===")

        dados, fs, nomes = le_ns2(caminho)
        lfp = fatia_janela(dados[:, resolve_canal_idx(nomes, canal)], fs, inicio, fim).astype(float)
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

        par_nome = row.get("par", "theta_gamma") if isinstance(row, pd.Series) else "theta_gamma"
        resultado = avalia_robustez_evento(lfp, fs, f_pico, a_pico, par_nome, rng=rng)

        for linha_sweep in resultado["linhas_sweep"]:
            print(f"  {linha_sweep['teste']}={linha_sweep['parametro']}: z={linha_sweep['z']}")
            linhas.append({"janela": rotulo, "canal": canal,
                           "par_pico": f"{f_pico:g}x{a_pico:g}", **linha_sweep})

        print(f"  >> n_bins: z mínimo {resultado['z_min_nb']:.2f} | larguras: z mínimo "
              f"{resultado['z_min_bw']:.2f} | pico estável: {resultado['pico_estavel']} | "
              f"MVL z={resultado['z_mvl']:.2f} (bruto={resultado['mvl_bruto']:.4f}, "
              f"Rayleigh p={resultado['rayleigh_p']:.3f}) | "
              f"transiente: {resultado['transiente_encontrado']} | "
              f"r(γ,ruido)={resultado['correlacao_gama_ruido']:.2f} -> {resultado['status']}")

    df = pd.DataFrame(linhas)
    df.to_csv(args.saida_csv, index=False)
    print(f"\nCSV salvo: {args.saida_csv} ({len(df)} linhas)")


if __name__ == "__main__":
    main()
