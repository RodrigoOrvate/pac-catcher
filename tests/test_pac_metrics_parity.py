"""
tests/test_pac_metrics_parity.py
==========================================
Teste de paridade numerica entre pac_core.pac_metrics e as implementacoes
legadas de KL-MI de Tort + surrogates que existiam em triagem_pac.py,
comodulogram.py, refina_candidatos.py e auditorias/diagnostico_janela.py
antes da migracao desses consumidores para o nucleo compartilhado.

Duas camadas:
  Camada A (estrutural): copias verbatim CONGELADAS das implementacoes
    legadas, mantidas como funcoes privadas neste arquivo (prefixo
    `_legado_`). So podem ser editadas junto com uma justificativa
    cientifica -- sao o contrato de paridade, nao um detalhe de
    implementacao. Comparadas ao nucleo novo com igualdade EXATA
    (np.testing.assert_array_equal / ==), nunca assert_allclose.
  Camada B (regressao literal): valores numericos capturados ANTES de
    qualquer consumidor ser migrado para pac_core.pac_metrics, com
    repr(float) de precisao total. Protege contra o caso em que nucleo
    novo e copia congelada da Camada A sao editados com o mesmo erro.

Convencao: script standalone (sem pytest), main() + sys.exit(1) em falha,
mesmo padrao dos demais arquivos em tests/.
"""
import os
import sys

_SCRIPT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _SCRIPT_ROOT)

import numpy as np
from scipy.signal import hilbert

from pac_core.filtering import filtra_sinal
from pac_core.pac_metrics import (
    _mi_de_bin_idx as mi_de_bin_idx_novo,
    fase_para_bin_idx,
    gera_deslocamentos,
    mi_surrogates_de_deslocamentos,
    z_score_mi,
    z_score_mi_mapa,
    p_empirico_mi,
    calcula_mi_com_surrogates,
)
from pipeline.etapa1_triagem.triagem_pac import gera_sinal_demo


# ==========================================
# CAMADA A: copias verbatim CONGELADAS das implementacoes legadas
# ==========================================

def _legado_mi_de_bin_idx(bin_idx, envelope, n_bins):
    """Copia verbatim (identica nas 5 implementacoes originais)."""
    soma_bins = np.bincount(bin_idx, weights=envelope, minlength=n_bins)
    cont_bins = np.bincount(bin_idx, minlength=n_bins)
    media_bins = np.divide(soma_bins, cont_bins, out=np.zeros(n_bins), where=cont_bins > 0)
    soma = np.sum(media_bins)
    if soma <= 0:
        return 0.0
    P = media_bins / soma
    H = -np.sum(P * np.log(P + 1e-10))
    return (np.log(n_bins) - H) / np.log(n_bins)


def _legado_mi_com_surrogates(fase, envelope, fs, n_surr=200, n_bins=18,
                               shift_min_s=1.0, rng=None):
    """Copia verbatim de triagem_pac.py::mi_com_surrogates (pre-migracao)."""
    if rng is None:
        rng = np.random.default_rng()

    bins = np.linspace(-np.pi, np.pi, n_bins + 1)
    bin_idx = np.clip(np.digitize(fase, bins) - 1, 0, n_bins - 1)
    mi_obs = _legado_mi_de_bin_idx(bin_idx, envelope, n_bins)

    n = len(envelope)
    shift_min = int(shift_min_s * fs)
    if n <= 2 * shift_min:
        shift_min = max(1, n // 10)

    mi_surr = np.empty(n_surr)
    deslocamentos = rng.integers(shift_min, n - shift_min, size=n_surr)
    for i, desloc in enumerate(deslocamentos):
        env_shift = np.roll(envelope, desloc)
        mi_surr[i] = _legado_mi_de_bin_idx(bin_idx, env_shift, n_bins)

    media = np.mean(mi_surr)
    dp    = np.std(mi_surr)
    z     = (mi_obs - media) / dp if dp > 0 else 0.0
    p_emp = float(np.mean(mi_surr >= mi_obs))

    return mi_obs, z, p_emp, media, dp


def _legado_nucleo_comodulograma(lfp_ativo, fs, fases_freq, amps_freq,
                                  n_surr=200, n_bins=18, rng=None):
    """Copia verbatim (reduzida: sem notch/retorna_mi, que nao afetam o
    nucleo de MI+surrogates) de comodulogram.py::calcula_comodulograma_z
    (pre-migracao)."""
    if rng is None:
        rng = np.random.default_rng()

    bins = np.linspace(-np.pi, np.pi, n_bins + 1)

    fases_por_freq = []
    for f_fase in fases_freq:
        lfp_fase = filtra_sinal(lfp_ativo, f_fase - 1.0, f_fase + 1.0, fs)
        fase = np.angle(hilbert(lfp_fase))
        fases_por_freq.append(np.clip(np.digitize(fase, bins) - 1, 0, n_bins - 1))

    envelopes_por_freq = []
    for f_amp in amps_freq:
        lfp_amp = filtra_sinal(lfp_ativo, f_amp - 5.0, f_amp + 5.0, fs)
        envelopes_por_freq.append(np.abs(hilbert(lfp_amp)))

    n = lfp_ativo.size
    shift_min = int(1.0 * fs)
    if n <= 2 * shift_min:
        shift_min = max(1, n // 10)
    deslocamentos = rng.integers(shift_min, n - shift_min, size=n_surr)

    mi_obs_mapa = np.zeros((len(amps_freq), len(fases_freq)))
    mi_surr_mapa = np.zeros((len(amps_freq), len(fases_freq), n_surr))
    for i, bin_idx in enumerate(fases_por_freq):
        for j, env in enumerate(envelopes_por_freq):
            mi_obs_mapa[j, i] = _legado_mi_de_bin_idx(bin_idx, env, n_bins)
            for k, desloc in enumerate(deslocamentos):
                mi_surr_mapa[j, i, k] = _legado_mi_de_bin_idx(bin_idx, np.roll(env, desloc), n_bins)

    media = mi_surr_mapa.mean(axis=2)
    dp = mi_surr_mapa.std(axis=2)
    z_mapa = np.divide(mi_obs_mapa - media, dp,
                        out=np.zeros_like(mi_obs_mapa), where=dp > 0)
    return z_mapa, mi_obs_mapa, mi_surr_mapa


# ==========================================
# FIXTURE: reusa gera_sinal_demo() de triagem_pac.py (extraida na Fase 1.0)
# ==========================================

def _fixture():
    sinal, fs, _ = gera_sinal_demo()
    trecho = sinal[30_000:40_000]  # segunda metade: teta+gamma+HG acoplados
    fase = np.angle(hilbert(filtra_sinal(trecho, 4, 8, fs)))
    env = np.abs(hilbert(filtra_sinal(trecho, 30, 80, fs)))
    return trecho, fs, fase, env


# ==========================================
# CASOS
# ==========================================

def caso_1_mi_de_bin_idx():
    _, fs, fase, env = _fixture()
    bin_idx = fase_para_bin_idx(fase, 18)

    novo = mi_de_bin_idx_novo(bin_idx, env, 18)
    legado = _legado_mi_de_bin_idx(bin_idx, env, 18)
    assert novo == legado, f"_mi_de_bin_idx: novo={novo!r} != legado={legado!r}"

    ESPERADO = 0.059218380378225785  # Camada B: capturado do legado vivo antes da migracao
    assert novo == ESPERADO, f"_mi_de_bin_idx: novo={novo!r} != ESPERADO={ESPERADO!r}"
    print(f"OK caso 1 (_mi_de_bin_idx isolado): {novo!r}")


def caso_2_mi_com_surrogates():
    _, fs, fase, env = _fixture()

    rng_novo = np.random.default_rng(12345)
    mi_obs, mi_surr = calcula_mi_com_surrogates(fase, env, fs, n_surr=200, n_bins=18, rng=rng_novo)
    z = z_score_mi(mi_obs, mi_surr)
    p_emp = p_empirico_mi(mi_obs, mi_surr)
    media = float(np.mean(mi_surr))
    dp = float(np.std(mi_surr))

    rng_legado = np.random.default_rng(12345)
    mi_obs_l, z_l, p_emp_l, media_l, dp_l = _legado_mi_com_surrogates(
        fase, env, fs, n_surr=200, n_bins=18, rng=rng_legado)

    assert mi_obs == mi_obs_l, f"mi_obs: {mi_obs!r} != {mi_obs_l!r}"
    assert z == z_l, f"z: {z!r} != {z_l!r}"
    assert p_emp == p_emp_l, f"p_emp: {p_emp!r} != {p_emp_l!r}"
    assert media == media_l, f"media: {media!r} != {media_l!r}"
    assert dp == dp_l, f"dp: {dp!r} != {dp_l!r}"

    # Camada B: capturado do legado vivo antes da migracao
    assert mi_obs == 0.059218380378225785
    assert z == 8.131957255232575
    assert p_emp == 0.0
    assert media == 0.010962269268132279
    assert dp == 0.005934132410625095
    print(f"OK caso 2 (mi_com_surrogates completo): mi_obs={mi_obs!r} z={z!r}")


def caso_3_comodulograma_grade():
    trecho, fs, _, _ = _fixture()
    fases_freq = np.array([4.0, 6.0, 8.0])
    amps_freq = np.array([30.0, 45.0, 60.0, 75.0])
    n_bins = 18
    n_surr = 20

    # Replica a logica de calcula_comodulograma_z usando as primitivas do
    # nucleo novo: UM conjunto de deslocamentos compartilhado por todas as
    # celulas, calculado fora do loop duplo i x j.
    bins = np.linspace(-np.pi, np.pi, n_bins + 1)
    fases_por_freq = []
    for f_fase in fases_freq:
        lfp_fase = filtra_sinal(trecho, f_fase - 1.0, f_fase + 1.0, fs)
        fase = np.angle(hilbert(lfp_fase))
        fases_por_freq.append(fase_para_bin_idx(fase, n_bins))

    envelopes_por_freq = []
    for f_amp in amps_freq:
        lfp_amp = filtra_sinal(trecho, f_amp - 5.0, f_amp + 5.0, fs)
        envelopes_por_freq.append(np.abs(hilbert(lfp_amp)))

    rng_novo = np.random.default_rng(999)
    deslocamentos = gera_deslocamentos(trecho.size, fs, n_surr=n_surr, rng=rng_novo)

    mi_obs_mapa = np.zeros((len(amps_freq), len(fases_freq)))
    mi_surr_mapa = np.zeros((len(amps_freq), len(fases_freq), n_surr))
    for i, bin_idx in enumerate(fases_por_freq):
        for j, env in enumerate(envelopes_por_freq):
            mi_obs_mapa[j, i] = mi_de_bin_idx_novo(bin_idx, env, n_bins)
            mi_surr_mapa[j, i, :] = mi_surrogates_de_deslocamentos(
                bin_idx, env, deslocamentos, n_bins)
    z_novo = z_score_mi_mapa(mi_obs_mapa, mi_surr_mapa)

    rng_legado = np.random.default_rng(999)
    z_legado, mi_obs_legado, mi_surr_legado = _legado_nucleo_comodulograma(
        trecho, fs, fases_freq, amps_freq, n_surr=n_surr, n_bins=n_bins, rng=rng_legado)

    np.testing.assert_array_equal(mi_obs_mapa, mi_obs_legado)
    np.testing.assert_array_equal(mi_surr_mapa, mi_surr_legado)
    np.testing.assert_array_equal(z_novo, z_legado)
    print(f"OK caso 3 (comodulograma grade, sorteio compartilhado): z max={z_novo.max():.6f}")


def caso_4_rng_compartilhado_avanca_e_reproduz():
    n, fs = 10_000, 1000.0

    rng_a = np.random.default_rng(55)
    d1 = gera_deslocamentos(n, fs, n_surr=10, rng=rng_a)
    d2 = gera_deslocamentos(n, fs, n_surr=10, rng=rng_a)  # mesmo objeto -> avanca
    assert not np.array_equal(d1, d2), (
        "duas chamadas com o MESMO rng deveriam produzir deslocamentos "
        "diferentes (o gerador avanca) -- isso e o que refina_candidatos.py "
        "e robustez_parametros.py dependem"
    )

    rng_b = np.random.default_rng(55)
    d1_reproduzido = gera_deslocamentos(n, fs, n_surr=10, rng=rng_b)
    np.testing.assert_array_equal(d1, d1_reproduzido)
    print("OK caso 4 (rng compartilhado avanca e e reproduzivel)")


def caso_5_rng_none_nao_levanta():
    n, fs = 10_000, 1000.0
    d = gera_deslocamentos(n, fs, n_surr=5, rng=None)
    assert len(d) == 5
    print("OK caso 5 (rng=None nao levanta, uso interativo preservado)")


def main():
    try:
        caso_1_mi_de_bin_idx()
        caso_2_mi_com_surrogates()
        caso_3_comodulograma_grade()
        caso_4_rng_compartilhado_avanca_e_reproduz()
        caso_5_rng_none_nao_levanta()
    except AssertionError as e:
        print(f"FALHOU: {e}")
        sys.exit(1)
    print("\nTODOS OS TESTES DE PARIDADE PASSARAM.")


if __name__ == "__main__":
    main()
