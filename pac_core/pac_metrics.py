"""
pac_core/pac_metrics.py
==========================================
Núcleo compartilhado de cálculo de acoplamento fase-amplitude (PAC): KL-MI
de Tort et al., binning de fase, e nula de surrogates por deslocamento
circular do envelope de amplitude.

Extraído das 5 implementações quase-idênticas que existiam em
`triagem_pac.py`, `comodulogram.py`, `refina_candidatos.py` e
`auditorias/diagnostico_janela.py`.

CONTRATO DE PARIDADE NUMÉRICA (não alterar sem atualizar
tests/test_pac_metrics_parity.py e justificar a mudança):
  1. bins = np.linspace(-pi, pi, n_bins+1); bin_idx = clip(digitize(fase,
     bins) - 1, 0, n_bins-1).
  2. MI: bincount + divide(..., where=cont>0); soma<=0 -> 0.0;
     H = -sum(P*log(P+1e-10)); (log(n_bins)-H)/log(n_bins).
  3. shift_min = int(shift_min_s*fs); fallback max(1, n//10) se
     n <= 2*shift_min.
  4. EXATAMENTE UMA chamada a rng.integers(shift_min, n-shift_min,
     size=n_surr) por invocação de gera_deslocamentos -- nunca uma
     chamada por surrogate. Isso é o que torna irrelevante a ordem em que
     os scripts legados calculam mi_obs vs. sorteiam (alguns fazem antes,
     outros depois) -- há um único consumo do gerador em qualquer ordem.
  5. dp = np.std(mi_surr) (ddof=0); z = (obs-mean)/dp if dp>0 else 0.0.
  6. p_emp = float(np.mean(mi_surr >= mi_obs)).

`rng` é sempre normalizado via `np.random.default_rng(rng)` na primeira
linha de cada função que precisa dele: aceita None (gerador novo,
não-determinístico -- necessário para o call site interativo de
comodulogram_interativo.py que não passa rng), int (semeado), ou um
Generator já existente (reaproveitado POR IDENTIDADE -- verificado que
np.random.default_rng(g) is g -- preservando o compartilhamento de estado
que refina_candidatos.py e robustez_parametros.py dependem entre chamadas
sucessivas do mesmo rng).

EXCLUSÕES DELIBERADAS (não migradas para este núcleo):
  - auditorias/audita_held_out.py: usa informação mútua via
    np.histogram2d (não é o KL-MI de Tort), e sorteia deslocamentos um a
    um dentro do loop (não vetorizado) -- consumiria o rng em ordem
    incompatível mesmo com a mesma seed. Testa uma pergunta metodológica
    distinta (double-dipping / reancoragem), autocontido.
  - robustez_parametros.py::mvl_z_par / mvl_bruto_e_rayleigh: métrica MVL
    (Canolty et al. 2006), não KL-MI. mvl_z_par reusa só
    gera_deslocamentos() daqui (a nula é compartilhada, a métrica não).
    mvl_bruto_e_rayleigh não usa rng/surrogates (Rayleigh analítico).
  - O ajuste de distribuição Gama em
    refina_candidatos.mi_com_surrogates_parametrico (p-valor analítico)
    permanece decisão científica local daquele script.
  - O default rng=np.random.default_rng(42) embutido no corpo de
    auditorias/audita_transientes.py::mi_z_par permanece lá -- é contrato
    deliberado daquele experimento (mesmos deslocamentos entre baseline e
    despike), não política do núcleo.
"""

import numpy as np

N_BINS_PADRAO = 18
N_SURR_PADRAO = 200
SHIFT_MIN_S_PADRAO = 1.0


def fase_para_bin_idx(fase, n_bins=N_BINS_PADRAO):
    bins = np.linspace(-np.pi, np.pi, n_bins + 1)
    return np.clip(np.digitize(fase, bins) - 1, 0, n_bins - 1)


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


def calcula_mi(fase, envelope, n_bins=N_BINS_PADRAO):
    bin_idx = fase_para_bin_idx(fase, n_bins)
    return _mi_de_bin_idx(bin_idx, envelope, n_bins)


def gera_deslocamentos(n_amostras, fs, n_surr=N_SURR_PADRAO,
                        shift_min_s=SHIFT_MIN_S_PADRAO, rng=None):
    """Sorteia os deslocamentos circulares da nula (>= shift_min_s de
    deslocamento). UMA única chamada a rng.integers() -- ver contrato no
    docstring do módulo."""
    rng = np.random.default_rng(rng)
    shift_min = int(shift_min_s * fs)
    if n_amostras <= 2 * shift_min:
        shift_min = max(1, n_amostras // 10)
    return rng.integers(shift_min, n_amostras - shift_min, size=n_surr)


def mi_surrogates_de_deslocamentos(bin_idx, envelope, deslocamentos, n_bins=N_BINS_PADRAO):
    """MI de cada surrogate, dado um array de deslocamentos já sorteado
    (ex.: por gera_deslocamentos). Permite reusar o MESMO conjunto de
    deslocamentos entre várias células (ex.: o mapa 2D do comodulograma)."""
    mi_surr = np.empty(len(deslocamentos))
    for i, desloc in enumerate(deslocamentos):
        env_shift = np.roll(envelope, desloc)
        mi_surr[i] = _mi_de_bin_idx(bin_idx, env_shift, n_bins)
    return mi_surr


def z_score_mi(mi_real, mi_surrogates):
    """z-score escalar de um único par fase-amplitude."""
    dp = np.std(mi_surrogates)
    if dp > 0:
        return (mi_real - np.mean(mi_surrogates)) / dp
    return 0.0


def z_score_mi_mapa(mi_obs_mapa, mi_surr_mapa, eixo=-1):
    """z-score vetorizado para um mapa (ex.: grade fase x amplitude do
    comodulograma), com mi_surr_mapa tendo o eixo de surrogates em `eixo`."""
    media = mi_surr_mapa.mean(axis=eixo)
    dp = mi_surr_mapa.std(axis=eixo)
    return np.divide(mi_obs_mapa - media, dp,
                      out=np.zeros_like(mi_obs_mapa), where=dp > 0)


def p_empirico_mi(mi_real, mi_surrogates):
    return float(np.mean(mi_surrogates >= mi_real))


def calcula_mi_com_surrogates(fase, envelope, fs, n_surr=N_SURR_PADRAO,
                               n_bins=N_BINS_PADRAO,
                               shift_min_s=SHIFT_MIN_S_PADRAO, rng=None):
    """Conveniência: calcula mi_obs e a distribuição completa de
    mi_surrogates para um par fase-amplitude. NÃO reduz a z -- quem quiser
    o z chama z_score_mi(mi_obs, mi_surr) depois; quem precisar da
    distribuição inteira (ajuste Gama, p-valor por célula do FDR) já a
    tem.

    Devolve: (mi_obs: float, mi_surrogates: np.ndarray de shape (n_surr,))
    """
    bin_idx = fase_para_bin_idx(fase, n_bins)
    mi_obs = _mi_de_bin_idx(bin_idx, envelope, n_bins)

    n = len(envelope)
    deslocamentos = gera_deslocamentos(n, fs, n_surr, shift_min_s, rng)
    mi_surr = mi_surrogates_de_deslocamentos(bin_idx, envelope, deslocamentos, n_bins)

    return mi_obs, mi_surr
