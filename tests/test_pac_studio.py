"""
tests/test_pac_studio.py
==========================================
Critério de aceite do backend do ThetaGamma-Studio (pac_studio): valida
que pac_studio.analisa_janela() é (1) uma API de 1 linha funcional sobre
dado real em disco (via formato .bin legado, sintético) e (2) tem paridade
BIT-A-BIT com a mesma sequência (notch -> filtra fase/amplitude -> Hilbert
-> MI+surrogates) chamada manualmente via pac_core -- nunca assert_allclose,
sempre comparação exata (mesmo padrão de tests/test_pac_metrics_parity.py).

Convencao: script standalone (sem pytest), main() + sys.exit(1) em falha,
mesmo padrao dos demais arquivos em tests/.
"""
import os
import sys
import tempfile

_SCRIPT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _SCRIPT_ROOT)

import numpy as np
import pandas as pd
from scipy.signal import hilbert
from scipy.stats import chi2_contingency

from pac_studio import analisa_janela, BAND_PAIRS, transicao_estado
from pac_core import workspace
from pac_core.io import carrega_dados, fatia_janela, resolve_canal_idx
from pac_core.filtering import filtra_sinal, aplica_notch
from pac_core.pac_metrics import calcula_mi_com_surrogates, z_score_mi, p_empirico_mi
from pipeline.etapa1_triagem.triagem_pac import gera_sinal_demo


def _escreve_bin_sintetico(sinal):
    """Grava `sinal` (1 canal) como .bin legado (int16), formato que
    pac_core.io.le_bin_legado já sabe ler -- sem inventar um formato novo
    só pro teste."""
    amp = np.max(np.abs(sinal)) or 1.0
    escala = 30000.0 / amp  # margem confortável dentro do range int16
    dados_int16 = (sinal * escala).astype(np.int16)
    f = tempfile.NamedTemporaryFile(suffix=".bin", delete=False)
    dados_int16.tofile(f.name)
    f.close()
    return f.name, escala


def caso_1_par_invalido_levanta():
    caminho, _ = _escreve_bin_sintetico(np.zeros(1000))
    try:
        try:
            analisa_janela(caminho, canal="Ch1", t_ini_s=0, t_fim_s=0.5,
                            par="banda_que_nao_existe", n_canais_bin=1, fs_bin=1000.0)
            raise AssertionError("Deveria ter levantado ValueError para par invalido")
        except ValueError:
            pass
    finally:
        os.remove(caminho)
    print("OK caso 1 (par invalido levanta ValueError)")


def caso_2_paridade_bit_exata():
    sinal, fs, _ = gera_sinal_demo(fs=1000.0, dur=60.0, seed=0)
    caminho, escala = _escreve_bin_sintetico(sinal)
    try:
        par = "theta_gamma"
        t_ini_s, t_fim_s = 30.0, 60.0  # metade com acoplamento real

        resultado = analisa_janela(
            caminho, canal="Ch1", t_ini_s=t_ini_s, t_fim_s=t_fim_s, par=par,
            n_canais_bin=1, fs_bin=fs, rng=42, notch=())

        # Replica a MESMA sequencia manualmente via pac_core, sem passar
        # por pac_studio, pra provar que o wrapper nao muda nenhum numero.
        dados, fs2, canal_ids = carrega_dados(caminho, n_canais_bin=1, fs_bin=fs)
        idx = resolve_canal_idx(canal_ids, "Ch1")
        trecho = fatia_janela(dados, fs2, t_ini_s, t_fim_s)[:, idx]
        trecho = aplica_notch(trecho, fs2, freqs_notch=[])  # notch=() no wrapper

        banda_fase = BAND_PAIRS[par]["fase"]
        banda_amp = BAND_PAIRS[par]["amp"]
        nyq = 0.5 * fs2
        lfp_fase = filtra_sinal(trecho, banda_fase[0], min(banda_fase[1], nyq * 0.98), fs2)
        lfp_amp = filtra_sinal(trecho, banda_amp[0], min(banda_amp[1], nyq * 0.98), fs2)
        fase = np.angle(hilbert(lfp_fase))
        envelope = np.abs(hilbert(lfp_amp))
        mi_obs, mi_surr = calcula_mi_com_surrogates(fase, envelope, fs2, rng=42)

        assert resultado["mi_observado"] == mi_obs, \
            f"mi_observado divergiu: {resultado['mi_observado']} != {mi_obs}"
        assert resultado["z_score"] == z_score_mi(mi_obs, mi_surr), "z_score divergiu"
        assert resultado["p_empirico"] == p_empirico_mi(mi_obs, mi_surr), "p_empirico divergiu"
        assert resultado["canal"] == "Ch1"
        assert resultado["banda_fase"] == banda_fase
        assert resultado["banda_amp"] == banda_amp
    finally:
        os.remove(caminho)
    print("OK caso 2 (paridade bit-exata com chamada manual via pac_core)")


def caso_3_contraste_acoplado_vs_nao_acoplado():
    # gera_sinal_demo: 1a metade e teta+ruido puro (sem acoplamento), 2a
    # metade tem gamma/HG travados na fase de teta (ver docstring da
    # funcao) -- z_score da 2a metade deve ser nitidamente maior.
    sinal, fs, _ = gera_sinal_demo(fs=1000.0, dur=60.0, seed=0)
    caminho, _ = _escreve_bin_sintetico(sinal)
    try:
        r_sem_acoplamento = analisa_janela(
            caminho, canal="Ch1", t_ini_s=0.0, t_fim_s=25.0, par="theta_gamma",
            n_canais_bin=1, fs_bin=fs, rng=42, notch=())
        r_com_acoplamento = analisa_janela(
            caminho, canal="Ch1", t_ini_s=35.0, t_fim_s=60.0, par="theta_gamma",
            n_canais_bin=1, fs_bin=fs, rng=42, notch=())

        assert r_com_acoplamento["z_score"] > r_sem_acoplamento["z_score"], (
            f"z_score com acoplamento ({r_com_acoplamento['z_score']:.2f}) deveria "
            f"ser maior que sem acoplamento ({r_sem_acoplamento['z_score']:.2f})"
        )
    finally:
        os.remove(caminho)
    print(f"OK caso 3 (contraste: sem acoplamento z={r_sem_acoplamento['z_score']:.2f} "
          f"< com acoplamento z={r_com_acoplamento['z_score']:.2f})")


def caso_4_transicao_estado_reproduz_h2():
    # Valida que transicao_estado() reproduz exatamente o resultado
    # documentado da Hipotese H2 (docs/sintese_investigacao_preditor_pac.md):
    # comportamento(N-1) -> acoplamento(N), chi2=41.3, p=2.6e-7, n=2715.
    caminho_mestre = os.path.join(workspace.BASE_SCRIPT, "resultados",
                                   "dataset_mestre_COM_COMPORTAMENTO.csv")
    if not os.path.exists(caminho_mestre):
        print("PULADO caso 4 (dataset_mestre_COM_COMPORTAMENTO.csv nao encontrado neste ambiente)")
        return

    df = pd.read_csv(caminho_mestre, encoding="utf-8-sig")
    consec = transicao_estado(df)
    assert len(consec) == 2715, f"esperado n=2715, veio {len(consec)}"

    tabela = pd.crosstab(consec["comport_anterior"], consec["vencedor"])
    chi2, p, dof, _ = chi2_contingency(tabela)
    assert abs(chi2 - 41.3) < 0.1, f"chi2 divergiu do documentado: {chi2:.2f}"
    assert abs(p - 2.6e-7) < 1e-8, f"p divergiu do documentado: {p:.3e}"
    print(f"OK caso 4 (transicao_estado reproduz H2: n={len(consec)}, "
          f"chi2={chi2:.2f}, p={p:.3e}, dof={dof})")


def main():
    try:
        caso_1_par_invalido_levanta()
        caso_2_paridade_bit_exata()
        caso_3_contraste_acoplado_vs_nao_acoplado()
        caso_4_transicao_estado_reproduz_h2()
    except AssertionError as e:
        print(f"FALHOU: {e}")
        sys.exit(1)
    print("\nTODOS OS TESTES DE pac_studio PASSARAM.")


if __name__ == "__main__":
    main()
