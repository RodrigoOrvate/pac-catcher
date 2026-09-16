"""
pac_studio — API Python de 1 linha do ThetaGamma-Studio (backend do programa desktop, dashboard_desktop/).

Não reimplementa nada: só orquestra pac_core (io/filtering/pac_metrics) e a
config de bandas de pipeline.etapa1_triagem.triagem_pac (BAND_PAIRS), que hoje
só existiam via cópia manual da sequência notch -> filtro fase/amplitude ->
Hilbert -> MI+surrogates dentro de cada script (triagem_pac.py, comodulogram.py
etc.). Consumidor novo (dashboard, notebook, script) deve importar daqui, não
duplicar a sequência de novo.

Caminho-base do workspace vem de pac_core.workspace (nunca hardcode
C:\\ ou D:\\acoplamento_theta-gamma aqui).
"""
import numpy as np
from scipy import signal

from pac_core.io import carrega_dados, fatia_janela, resolve_canal_idx
from pac_core.filtering import filtra_sinal, aplica_notch
from pac_core.pac_metrics import (
    calcula_mi_com_surrogates, z_score_mi, p_empirico_mi,
    N_SURR_PADRAO, N_BINS_PADRAO, SHIFT_MIN_S_PADRAO,
)
from pac_core import workspace
from pipeline.etapa1_triagem.triagem_pac import BAND_PAIRS

NOTCH_PADRAO = (60.0, 120.0, 180.0, 240.0)

__all__ = ["analisa_janela", "BAND_PAIRS", "workspace", "transicao_estado", "CHAVE_TRANSICAO"]

# Mesma chave de agrupamento e limiar de gap de
# preditor/preditor_estado_comportamental.py::monta_dataset() -- importa a
# constante de lá em vez de retypar, pra nunca divergir se ela mudar.
from preditor.preditor_estado_comportamental import GAP_MAX_S as _GAP_MAX_S_PADRAO

CHAVE_TRANSICAO = ["sessao", "arquivo", "canal", "par"]


def transicao_estado(df_mestre, gap_max_s=None):
    """Reconstrói os pares consecutivos (estado comportamental N-1 -> janela
    N) do dataset mestre, MESMA lógica de setup usada em
    preditor_estado_comportamental.py::monta_dataset() (linhas antes do laço
    caro de extração de features por .ns2, que não é necessário pra esta
    análise). Valida H2 (docs/sintese_investigacao_preditor_pac.md): contra dataset_mestre_COM_COMPORTAMENTO.csv
    reproduz exatamente n=2715, chi2=41.26, p=2.57e-7 (ver
    tests/test_pac_studio.py).

    Retorna o df filtrado (só pares consecutivos, sem gap, sem 'Artefato /
    Cabo') com as colunas novas 'vencedor' (bool) e 'comport_anterior'.
    Quem chama decide o resto (crosstab, chi2_contingency, filtro por
    condicao/grupo etc.) -- essa função só entrega o dado limpo.
    """
    gap_max_s = _GAP_MAX_S_PADRAO if gap_max_s is None else gap_max_s

    df = df_mestre.sort_values(CHAVE_TRANSICAO + ["janela_ini_s"]).copy()
    df["vencedor"] = df["veredito_refino"] == "Candidato robusto"
    df["comport_anterior"] = df.groupby(CHAVE_TRANSICAO)["comportamento"].shift(1)
    janela_ini_anterior = df.groupby(CHAVE_TRANSICAO)["janela_ini_s"].shift(1)
    df["gap_s"] = df["janela_ini_s"] - janela_ini_anterior

    return df[(df["gap_s"] <= gap_max_s) & df["comport_anterior"].notna() &
              (df["comport_anterior"] != "Artefato / Cabo")].copy()


def analisa_janela(caminho_arquivo, canal, t_ini_s, t_fim_s, par="theta_gamma",
                    n_canais_bin=16, fs_bin=1000.0, n_surr=N_SURR_PADRAO,
                    n_bins=N_BINS_PADRAO, shift_min_s=SHIFT_MIN_S_PADRAO,
                    rng=42, notch=NOTCH_PADRAO):
    """Analisa PAC (MI observado + z-score + p-valor vs. surrogates) de uma
    janela de tempo de um canal, numa chamada única.

    Consolida a sequência hoje duplicada em triagem_pac.py/comodulogram.py:
    carrega -> fatia a janela -> notch -> filtra fase/amplitude -> Hilbert ->
    MI com surrogates. `par` escolhe a banda de amplitude via BAND_PAIRS
    (mesma config usada em produção -- "theta_gamma", "theta_hg" ou
    "theta_hfo"). `canal` aceita nome nativo do .ns2 ("chan5") ou índice
    1-based do dataset mestre (ver pac_core.io.resolve_canal_idx).

    `rng` default=42 segue a convenção do projeto (semente fixa p/
    reprodutibilidade, ver CLAUDE.md #5) -- passe outro valor ou None se
    precisar de surrogates diferentes.

    Retorna um dict; não decide "significativo" (isso é critério de cada
    consumidor -- ver refina_candidatos.py).
    """
    if par not in BAND_PAIRS:
        raise ValueError(f"par desconhecido: {par!r}. Use um de {list(BAND_PAIRS)}")
    banda_fase = BAND_PAIRS[par]["fase"]
    banda_amp = BAND_PAIRS[par]["amp"]

    dados, fs, canal_ids = carrega_dados(caminho_arquivo, n_canais_bin=n_canais_bin, fs_bin=fs_bin)
    idx = resolve_canal_idx(canal_ids, canal)
    trecho = fatia_janela(dados, fs, t_ini_s, t_fim_s)[:, idx]

    nyq = 0.5 * fs
    trecho = aplica_notch(trecho, fs, freqs_notch=list(notch))
    lfp_fase = filtra_sinal(trecho, banda_fase[0], min(banda_fase[1], nyq * 0.98), fs)
    lfp_amp = filtra_sinal(trecho, banda_amp[0], min(banda_amp[1], nyq * 0.98), fs)
    fase = np.angle(signal.hilbert(lfp_fase))
    envelope = np.abs(signal.hilbert(lfp_amp))

    mi_obs, mi_surr = calcula_mi_com_surrogates(
        fase, envelope, fs, n_surr=n_surr, n_bins=n_bins,
        shift_min_s=shift_min_s, rng=rng)

    return {
        "mi_observado": float(mi_obs),
        "z_score": float(z_score_mi(mi_obs, mi_surr)),
        "p_empirico": p_empirico_mi(mi_obs, mi_surr),
        "par": par,
        "canal": canal_ids[idx],
        "fs": fs,
        "banda_fase": banda_fase,
        "banda_amp": banda_amp,
        "n_amostras": len(trecho),
    }
