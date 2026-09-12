"""
comodulogram_interativo.py
==========================================
Estação de exploração visual de uma sessão antes da triagem automática.

Fluxo:
  1. PSD (linear + log) do canal, com 4 subplots.
  2. Comodulograma geral em janela de 1 min (pico ΘΓ).
  3. Time-series de MI (z do pico ΘΓ) e theta power (4-8 Hz).
  4. Click num pico do MI (ou --zoom_t_center headless) para zoom em janela
     de 10s: comodulograma da janela + ajuste FOOOF (knee) de teta e
     gama/HG sobre um contexto de 45s, mostrando PSD bruto, fundo
     aperiódico e picos periódicos isolados.

Uso:
  python pipeline/comodulogram_interativo.py \\
      --pasta_ns2 ".../Basal antes da infusao" \\
      --canal 5,16,22 \\
      --win_s 10 --passo_s 5 \\
      --sessao MTESC04_S1 \\
      --saida ".../RESULTADOS/diagnosticos_interativos/MTESC04_S1"

Requer: neo, numpy, scipy, matplotlib, fooof (ou specparam).
"""
import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize

# --- utilitários do pipeline --------------------------------------------------
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from pac_core.io import carrega_dados, fatia_janela
from pac_core.filtering import filtra_sinal, aplica_notch
from pac_core.pac_metrics import _mi_de_bin_idx
from pipeline.etapa1_triagem.triagem_pac import BAND_PAIRS
from pipeline.etapa3_comodulograma.comodulogram import (
    calcula_comodulograma_z,
    z_pico_par,
    bh_fdr_mapa,
    p_valores_por_celula,
    FASES_DEFAULT,
    AMPS_DEFAULT,
)

# fooof só é usado no zoom combinado (clique/--zoom_t_center); as demais
# features (PSD, comodulograma geral, timeline) funcionam sem ele.
try:
    from fooof import FOOOF
    from pipeline.auditorias.linha_noise_kuhn import aplica_modo as _aplica_modo_kuhn
    from pipeline.auditorias.utils_harmonico import _knee_valido
    _FOOOF_DISPONIVEL = True
except ImportError:
    FOOOF = None
    _FOOOF_DISPONIVEL = False

# --- constantes ---------------------------------------------------------------
BANDA_TETA = (4, 8)
FASES_FREQ = FASES_DEFAULT
AMPS_FREQ = AMPS_DEFAULT
N_BINS = 18
N_SURR = 200
FDR_Q = 0.05
# Limiar do z do pico ΘΓ para marcar como "evento" no time-series.
#
# ATENÇÃO: este limiar é INTENCIONALMENTE PERMISSIVO. Em canais com sinal
# fraco (ex.: referência, canal danificado), a distribuição de MI dos
# surrogates por deslocamento circular fica muito estreita (std → 0), e
# qualquer flutuação numérica pode dar z>3. Isso foi confirmado na
# validação com canal 0 da sessão 1 (silencioso, sem pico de teta/gama
# na PSD) — produziu 22 "eventos" espúrios com z entre 3.0 e 8.6.
#
# Por isso subimos o default para 5.0 (≈5σ sobre uma nula com std=1):
# - Canal 0 (controle silencioso): z máximo cai para ~1-2 → 0 eventos
# - Canais reais (canal 16 da sessão 1): z do pico ΘΓ ≈ 18 → mantém
# - Canal com sinal fraco mas genuíno: 5σ sobre a nula é uma barra
#   estatística razoável para triagem visual.
#
# O portão rigoroso do pipeline continua sendo o FDR no mapa 2D
# (n_sig_FDR da `calcula_comodulograma_z`) — não este limiar pontual.
LIMIAR_Z_EVENTO = 5.0   # z do pico ΘΓ para marcar como "evento"

# --- parâmetros FOOOF (réplica EXATA de utils_harmonico.py::extrai_cf_teta_fooof
# / extrai_cf_gamma_fooof -- a função REALMENTE usada por enriquece_dataset_mestre.py
# para gerar r2_teta_fooof_v2/knee_valido_teta/etc. no dataset mestre; NÃO
# gera_plot_fooof.py, que é só uma figura de apresentação sem limpeza de
# ruído de linha). Não recalibrar aqui: ver memória do projeto sobre
# max_n_peaks ser arquitetural (reduzir quebra detecção de harmônicos
# biológicos 2f/3f) e sobre Kuhn só remover ruído de linha, nunca teta.
FOOOF_CTX_S = 45.0  # janela de contexto p/ resolução espectral (mesmo valor de audita_harmonico.py)
FOOOF_MODO_PREPROCESSO = "hibrido"  # default de producao (linha_noise_kuhn)
FOOOF_TETA_FIT_RANGE = (2.0, 45.0)
FOOOF_TETA_KWARGS = dict(aperiodic_mode="knee", peak_width_limits=(2, 5),
                          min_peak_height=0.05, peak_threshold=1.0, max_n_peaks=4)
FOOOF_GAMMA_BW_LIMITS = (4, 30)
FOOOF_GAMMA_KWARGS = dict(aperiodic_mode="knee", peak_width_limits=(4, 30),
                           min_peak_height=0.05, peak_threshold=1.0, max_n_peaks=6)


# --- core: concatenação da sessão ---------------------------------------------

def concatena_sessao(pasta_basal, fs_esperado=1000.0):
    """Lê todos os .ns2 da pasta Basal e concatena em matriz contínua.

    Retorna:
      dados : ndarray (n_amostras_total, n_canais)
      fs    : float  (frequência de amostragem, em Hz)
      ids   : list of channel names (do primeiro arquivo)
      offsets : list of int (amostras de início de cada arquivo)
    """
    arquivos = sorted([
        os.path.join(pasta_basal, f)
        for f in os.listdir(pasta_basal)
        if f.lower().endswith((".ns2", ".bin"))
    ])
    if not arquivos:
        raise FileNotFoundError(f"Nenhum .ns2/.bin encontrado em {pasta_basal}")

    blocos = []
    fs = None
    ids = None
    offsets = [0]

    for caminho in arquivos:
        dados, fs_canal, ids_canal = carrega_dados(caminho)
        if fs is None:
            fs = fs_canal
            ids = ids_canal
        else:
            if fs_canal != fs:
                print(f"  AVISO: {os.path.basename(caminho)} tem fs={fs_canal},"
                      f" esperado {fs}. Usando {fs}.")
        blocos.append(dados.astype(np.float64))
        offsets.append(sum(b.shape[0] for b in blocos))

    dados = np.vstack(blocos)
    return dados, fs, ids, offsets


# --- core: time-series de MI + theta power ------------------------------------

def time_series_pac(dados, fs, canais,
                    win_s=10, passo_s=5,
                    fases_freq=None, amps_freq=None,
                    n_surr=N_SURR, n_bins=N_BINS,
                    notch_hz=60.0,
                    max_janelas=0,
                    par_ativo="theta_gamma"):
    """Computa MI (z do pico do par ativo) e theta power em janelas deslizantes.

    Para cada canal e cada janela deslizante:
      - recorta janela do sinal
      - aplica notch 60 Hz + harmônicos
      - computa comodulograma z (vs surrogates)
      - extrai z do pico no quadrante ΘΓ
      - computa theta power (Welch na banda 4-8 Hz)

    Args:
      dados    : ndarray (n_amostras, n_canais)
      fs       : float, Hz
      canais   : list of int
      win_s    : duração de cada janela em segundos
      passo_s  : passo da janela deslizante em segundos
      fases_freq, amps_freq : arrays de frequências (default: 4-14 e 30-80)
      n_surr, n_bins : parâmetros do comodulograma
      notch_hz : frequência da rede elétrica (60 para BR, 50 para EU)
      max_janelas : se >0, limita o número de janelas (validação rápida)

    Returns:
      t_centers    : ndarray (n_janelas,) — tempo em segundos do centro da janela
      mi_z_pico    : dict {canal: ndarray (n_janelas,)} — z do pico ΘΓ
      theta_power  : dict {canal: ndarray (n_janelas,)} — potência em 4-8 Hz
      win_inicios  : ndarray (n_janelas,) — índice de início de cada janela
    """
    if fases_freq is None:
        fases_freq = FASES_DEFAULT
    if amps_freq is None:
        amps_freq = AMPS_DEFAULT

    n_total = dados.shape[0]
    win_n = int(win_s * fs)
    passo_n = int(passo_s * fs)

    # centros das janelas
    inicios_all = np.arange(0, n_total - win_n + 1, passo_n)
    if max_janelas > 0:
        inicios = inicios_all[:max_janelas]
    else:
        inicios = inicios_all
    t_centers = (inicios + win_n // 2) / fs

    n_janelas = len(inicios)
    print(f"  {n_janelas} janelas ({win_s}s, passo {passo_s}s) sobre"
          f" {n_total/fs:.0f}s total")

    mi_z = {ch: np.full(n_janelas, np.nan) for ch in canais}
    th_pow = {ch: np.full(n_janelas, np.nan) for ch in canais}

    for i, ini in enumerate(inicios):
        fim = ini + win_n
        trecho = dados[ini:fim, :]

        for ch in canais:
            lfp = trecho[:, ch].astype(float)

            # notch 60 Hz
            lfp_n = aplica_notch(lfp, fs, linha_hz=notch_hz)

            # --- theta power (Welch na banda 4-8 Hz) -------------------
            from scipy.signal import welch
            nperseg_w = min(int(2.0 * fs), win_n)
            freqs_w, psd_w = welch(lfp_n, fs=fs, window='hann',
                                    nperseg=nperseg_w, noverlap=nperseg_w // 2)
            idx_teta = (freqs_w >= BANDA_TETA[0]) & (freqs_w <= BANDA_TETA[1])
            if idx_teta.any():
                # normalizar pela potência total (1-100 Hz) para comparabilidade
                idx_total = (freqs_w >= 1) & (freqs_w <= 100)
                th_pow[ch][i] = (psd_w[idx_teta].mean()
                                 / psd_w[idx_total & (freqs_w > 1)].mean()
                                 if idx_total.any() else np.nan)

            # --- MI z-scoredo -----------------------------------------
            # usa janela do comodulograma se for grande o suficiente
            if win_n < 2 * fs:
                # janela < 2s: comodulograma fica ruidoso; marca como nan
                mi_z[ch][i] = np.nan
                continue

            try:
                z_mapa = calcula_comodulograma_z(
                    lfp_n, fs, fases_freq, amps_freq,
                    n_surr=n_surr, n_bins=n_bins,
                    rng=np.random.default_rng(i), notch_hz=None,
                )
                z_pico, _, _ = z_pico_par(
                    z_mapa, fases_freq, amps_freq, par=par_ativo)
                mi_z[ch][i] = z_pico
            except Exception:
                mi_z[ch][i] = np.nan


        if (i + 1) % 50 == 0:
            print(f"    janela {i+1}/{n_janelas} ...")

    return t_centers, mi_z, th_pow, inicios


def _welch(x, fs, nperseg=1024):
    """Wrapper simples para scipy.signal.welch (evita re-importar)."""
    from scipy.signal import welch as _w
    return _w(x, fs=fs, window='hann', nperseg=nperseg, noverlap=nperseg // 2)


# --- core: FOOOF teta + gama/HG (zoom combinado) -------------------------------

def ajusta_fooof_teta_gamma(sinal_ctx, fs, modo_preprocesso=FOOOF_MODO_PREPROCESSO,
                            f_linha=60.0):
    """Ajusta FOOOF (aperiodic_mode='knee') separadamente em teta (2-45 Hz)
    e gama/HG (35 Hz-~0.95*Nyquist), replicando byte-a-byte
    utils_harmonico.py::extrai_cf_teta_fooof/extrai_cf_gamma_fooof -- os
    mesmos parâmetros já usados para gerar as colunas de qualidade FOOOF do
    dataset mestre (Welch independente por banda + limpeza de linha Kuhn
    ANTES do fit; sem isso, harmônicos de 60 Hz contaminam o ajuste).

    Retorna (freqs_teta, psd_teta, fm_teta, freqs_gamma, psd_gamma, fm_gamma).
    """
    from scipy.signal import welch

    # --- teta -----------------------------------------------------------
    nperseg_t = int(1.2 * fs)
    nfft_t = 4000 if nperseg_t <= 4000 else nperseg_t
    freqs_t, psd_t_bruto = welch(sinal_ctx, fs=fs, window="hann", nperseg=nperseg_t,
                                 noverlap=nperseg_t // 2, nfft=nfft_t)
    psd_t_limpo = _aplica_modo_kuhn(freqs_t, psd_t_bruto, fs, modo_preprocesso, f_linha=f_linha)

    fm_teta = FOOOF(**FOOOF_TETA_KWARGS)
    fm_teta.fit(freqs_t, psd_t_limpo, freq_range=list(FOOOF_TETA_FIT_RANGE))

    # --- gama/HG (faixa e nfft dependem de fs, igual extrai_cf_gamma_fooof) --
    f_max = min(250.0, fs * 0.5 * 0.95)
    fit_range_gamma = (35.0, f_max)
    nperseg_g = int(1.2 * fs)
    nfft_g = max(nperseg_g, 4 * int(fs))
    freqs_g, psd_g_bruto = welch(sinal_ctx, fs=fs, window="hann", nperseg=nperseg_g,
                                 noverlap=nperseg_g // 2, nfft=nfft_g)
    psd_g_limpo = _aplica_modo_kuhn(freqs_g, psd_g_bruto, fs, modo_preprocesso, f_linha=f_linha)

    fm_gamma = FOOOF(**FOOOF_GAMMA_KWARGS)
    fm_gamma.fit(freqs_g, psd_g_limpo, freq_range=list(fit_range_gamma))

    # psd "bruto" (pré-Kuhn) devolvido p/ o painel A (comparação honesta
    # do que a limpeza de linha removeu); fm_teta/fm_gamma já carregam o
    # espectro LIMPO usado no fit (fm.power_spectrum).
    return freqs_t, psd_t_bruto, fm_teta, freqs_g, psd_g_bruto, fm_gamma


def painel_fooof(ax, fm, cor_ap, cor_flat, cor_pico, faixa_pico, rotulo, fit_range):
    """Um painel: PSD real, fundo aperiódico (knee) e curva achatada com
    picos marcados na faixa de interesse. Subtração em log10 (fm já opera
    em log-potência) -- ver memória do projeto "fooof-log10-subtraction".
    """
    f = fm.freqs
    p = fm.power_spectrum
    ap = fm._ap_fit
    flat = p - ap
    offset_visual = 1.15  # só separa visualmente a curva achatada (igual gera_plot_fooof.py)

    ax.plot(f, p, color="#2c3e50", alpha=0.5, lw=1.3, label="PSD real")
    ax.plot(f, ap, color=cor_ap, linestyle="--", lw=2.0, label="Fundo aperiódico (knee)")
    ax.plot(f, flat + offset_visual, color=cor_flat, lw=1.8, label="Picos puros (achatado)")

    picos = fm.get_params("peak_params")
    if len(picos) > 0:
        ymin, ymax_dados = ax.get_ylim()
        picos_na_faixa = sorted(
            (float(cf) for cf, amp, bw in picos if faixa_pico[0] <= cf <= faixa_pico[1]))

        # empilha o rótulo em níveis verticais quando picos ficam a <15 Hz
        # um do outro (senão o texto se sobrepõe -- ex.: cluster de
        # harmônicos entre 150-210 Hz no painel gama/HG)
        DIST_MIN_HZ = 15.0
        cf_por_nivel = []  # último cf colocado em cada nível
        nivel_de_cada_pico = []
        for cf in picos_na_faixa:
            nivel = 0
            while nivel < len(cf_por_nivel) and abs(cf - cf_por_nivel[nivel]) < DIST_MIN_HZ:
                nivel += 1
            if nivel == len(cf_por_nivel):
                cf_por_nivel.append(cf)
            else:
                cf_por_nivel[nivel] = cf
            nivel_de_cada_pico.append(nivel)

        n_niveis = max(nivel_de_cada_pico) + 1 if nivel_de_cada_pico else 1
        passo = 0.08 * (ymax_dados - ymin)
        ax.set_ylim(ymin, ymax_dados + n_niveis * passo)

        for cf, nivel in zip(picos_na_faixa, nivel_de_cada_pico):
            ax.vlines(cf, ymin, ymax_dados, color=cor_pico, linestyle=":", lw=1.5)
            ax.text(cf, ymax_dados + nivel * passo, f"{cf:.2f} Hz", color=cor_pico,
                    ha="center", va="bottom", fontsize=8, fontweight="bold")

    exp_ = fm.get_params("aperiodic_params")[2]
    knee_ = fm.get_params("aperiodic_params")[1]
    r2_ = fm.get_params("r_squared")
    knee_ok = _knee_valido(knee_, exp_, fit_range)
    texto = (f"$R^2$={r2_:.3f}\nExpoente={exp_:.2f}\nKnee={knee_:.1f}\n"
             f"Knee válido={'Sim' if knee_ok else 'Não'}")
    ax.text(0.97, 0.95, texto, transform=ax.transAxes, fontsize=8, ha="right", va="top",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#f8f9fa", edgecolor="#bdc3c7", alpha=0.95))

    ax.set_title(f"FOOOF {rotulo} (knee)")
    ax.set_xlabel("Frequência (Hz)")
    ax.set_ylabel(r"Log$_{10}$ Potência")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, linestyle=":")


def plota_zoom_fooof_comodulograma(sinal_ctx, fs, z_mapa_plot, fases_freq, amps_freq,
                                   fp_pico, fa_pico, z_pico, par_ativo,
                                   canal, nomes_canais, t_center, nome_evento,
                                   saida_png, notch_hz=60.0):
    """4 painéis: PSD bruto | FOOOF teta (knee) | FOOOF gama/HG (knee) |
    comodulograma da janela de 10s. Sempre salva em disco (headless); se a
    figura for exibida via plt.show() (modo --interactive), funciona igual.
    """
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 11))
    titulo_canal = f"canal {canal+1}{f' ({nomes_canais[canal]})' if nomes_canais else ''}"
    fig.suptitle(f"Zoom FOOOF + Comodulograma — {titulo_canal} @ t={t_center:.1f}s", fontsize=13)

    erro_fooof = None
    if not _FOOOF_DISPONIVEL:
        erro_fooof = "pacote 'fooof' não instalado"
    else:
        try:
            freqs_t, psd_t_bruto, fm_teta, freqs_g, psd_g_bruto, fm_gamma = \
                ajusta_fooof_teta_gamma(sinal_ctx, fs, f_linha=notch_hz)
        except Exception as e:
            erro_fooof = str(e)

    if erro_fooof is not None:
        for ax, msg in [(ax1, "PSD bruto"), (ax2, "FOOOF teta"), (ax3, "FOOOF gama/HG")]:
            ax.text(0.5, 0.5, f"({msg} indisponível:\n{erro_fooof})",
                    ha="center", va="center", transform=ax.transAxes,
                    fontsize=9, color="gray")
    else:
        # --- painel A: PSD bruto sem FOOOF nem limpeza de linha -----------
        # (usa o Welch de teta, que cobre 0-fs/2 de qualquer forma; mostra
        # o sinal tal como chega, incluindo harmônicos de 60 Hz -- contraste
        # honesto com os painéis B/C, que já removem isso antes do fit)
        f_max_plot = min(150.0, fs * 0.5 * 0.95)
        mask_tot = (freqs_t >= 2.0) & (freqs_t <= f_max_plot)
        ax1.plot(freqs_t[mask_tot], np.log10(psd_t_bruto[mask_tot]), color="#2c3e50", lw=1.8,
                 label="PSD Original (Welch, sem limpeza)")
        ax1.axvspan(4, 8, color="#3498db", alpha=0.18, label="Teta (4-8 Hz)")
        ax1.axvspan(30, 80, color="#e74c3c", alpha=0.15, label="Gama (30-80 Hz)")
        ax1.axvspan(80, f_max_plot, color="#9b59b6", alpha=0.15, label="HG")
        ax1.set_title("PSD bruto (sem FOOOF, sem limpeza de linha)")
        ax1.set_xlabel("Frequência (Hz)")
        ax1.set_ylabel(r"Log$_{10}$ Potência")
        ax1.legend(fontsize=8)
        ax1.grid(alpha=0.3, linestyle=":")

        # --- painel B: FOOOF teta -------------------------------------------
        painel_fooof(ax2, fm_teta, cor_ap="#2980b9", cor_flat="#27ae60",
                     cor_pico="#e74c3c", faixa_pico=(4.0, 12.0), rotulo="Teta",
                     fit_range=FOOOF_TETA_FIT_RANGE)

        # --- painel C: FOOOF gama/HG -----------------------------------------
        fit_range_gamma = (35.0, min(250.0, fs * 0.5 * 0.95))
        painel_fooof(ax3, fm_gamma, cor_ap="#8e44ad", cor_flat="#e67e22",
                     cor_pico="#d35400", faixa_pico=fit_range_gamma, rotulo="Gama/HG",
                     fit_range=fit_range_gamma)

    # --- painel D: comodulograma da janela de 10s ---------------------------
    im = ax4.pcolormesh(fases_freq, amps_freq, z_mapa_plot, shading="auto", cmap="viridis")
    if fp_pico is not None and fa_pico is not None:
        ax4.scatter([fp_pico], [fa_pico], s=80, c="red", marker="X",
                   edgecolors="white", linewidths=1.5)
    ax4.set_xlabel("Fase (Hz)")
    ax4.set_ylabel("Amplitude (Hz)")
    if z_pico is not None and not (isinstance(z_pico, float) and np.isnan(z_pico)):
        titulo_pico = f"pico {par_ativo} z={z_pico:.2f} em {fp_pico:.1f}x{fa_pico:.1f}Hz"
    else:
        titulo_pico = f"(sem pico {par_ativo} nesta janela)"
    ax4.set_title(f"Comodulograma 10s\n{titulo_pico}")
    fig.colorbar(im, ax=ax4, label="z")

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(saida_png, dpi=150, bbox_inches="tight")
    print(f"  Zoom FOOOF+comodulograma salvo: {saida_png}")
    return fig


def gera_zoom_comodulograma_fooof(dados, fs, canal, t_center, nomes_canais, saida_base,
                                  par_ativo="theta_gamma", ctx_s=FOOOF_CTX_S,
                                  nome_evento=None, notch_hz=60.0):
    """Monta a figura combinada (FOOOF teta/gama + comodulograma) para a
    janela de 10s centrada em ``t_center`` (tempo GLOBAL da sessão
    concatenada). Reusado pelo click handler (--interactive) e por chamada
    direta headless (ex.: --zoom_t_center).
    """
    t0_global, t1_global = t_center - 5, t_center + 5
    i0 = max(0, int(t0_global * fs))
    i1 = min(dados.shape[0], int(t1_global * fs))
    if i1 - i0 < int(0.5 * fs):
        raise ValueError(f"janela de comodulograma inválida ({i1 - i0} amostras)")
    seg = dados[i0:i1, canal]

    z_mapa, mi_obs, mi_surr = calcula_comodulograma_z(
        seg, fs, FASES_FREQ, AMPS_FREQ, n_surr=N_SURR, n_bins=N_BINS, retorna_mi=True)
    p_mapa = p_valores_por_celula(mi_obs, mi_surr)
    sig = bh_fdr_mapa(p_mapa, alpha=FDR_Q)
    z_mapa_plot = np.where(sig, z_mapa, np.nan)

    z_pico, fp_pico, fa_pico = z_pico_par(z_mapa, FASES_FREQ, AMPS_FREQ, par=par_ativo)
    if z_pico is None or (isinstance(z_pico, float) and np.isnan(z_pico)):
        fp_pico, fa_pico = None, None

    t_total_s = dados.shape[0] / fs
    ctx_ini = max(0.0, t_center - ctx_s / 2)
    ctx_fim = min(t_total_s, t_center + ctx_s / 2)
    ic0, ic1 = int(ctx_ini * fs), int(ctx_fim * fs)
    sinal_ctx = dados[ic0:ic1, canal].astype(float)

    if nome_evento is None:
        nome_evento = (f"canal {canal+1}"
                       f"{' (' + nomes_canais[canal] + ')' if nomes_canais else ''}"
                       f" @ t_global={t_center:.1f}s")

    nome_zoom = f"zoom_fooof_ch{canal+1}_{t_center:.0f}s.png"
    caminho_zoom = os.path.join(saida_base, nome_zoom)
    fig = plota_zoom_fooof_comodulograma(
        sinal_ctx, fs, z_mapa_plot, FASES_FREQ, AMPS_FREQ,
        fp_pico, fa_pico, z_pico, par_ativo, canal, nomes_canais,
        t_center, nome_evento, caminho_zoom, notch_hz=notch_hz)
    return fig, z_pico


# --- core: PSD dupla escala ---------------------------------------------------

def plota_psd_dupla(lfp, fs, canal_id, saida_png, notch_hz=60.0):
    """PSD linear (topo) e log (embaixo), 4 subplots cada.

    Subplots: 0-200 Hz | 0-50 Hz | 50-150 Hz | 4-15 Hz (theta zoom)
    """
    lfp_n = aplica_notch(lfp, fs, linha_hz=notch_hz) if notch_hz else lfp

    nperseg = int(2.0 * fs)
    nfft = 4000
    freqs, psd = _welch(lfp_n, fs, nperseg=nperseg)
    psd = psd / fs  # normalize to power density

    fig, axes = plt.subplots(2, 2, figsize=(14, 8))
    fig.suptitle(f"PSD — canal {canal_id}", fontsize=14)

    ranges = [(0.2, 200), (0.2, 50), (50, 150), (3, 20)]
    titles = ["0–200 Hz (vista geral)",
              "0–50 Hz (baixas)",
              "50–150 Hz (gamma)",
              "3–20 Hz (theta zoom)"]
    row = 0
    for idx, ((lo, hi), tit) in enumerate(zip(ranges, titles)):
        ax = axes.flat[idx]
        mask = (freqs >= lo) & (freqs <= hi)
        ax.plot(freqs[mask], psd[mask], color="steelblue", lw=0.8)
        ax.set_xlim(lo, hi)
        ax.set_xlabel("Frequência (Hz)")
        ax.set_title(tit)
        if lo < 20:
            # marca bandas
            for (blo, bhi), cor, lab in [
                ((4, 8), "gold", "θ"),
                ((30, 80), "salmon", "γ"),
            ]:
                if lo <= blo and hi >= bhi:
                    ax.axvspan(blo, bhi, color=cor, alpha=0.15, label=lab)
            ax.legend(fontsize=7)
        ax.grid(alpha=0.3)

    fig.supxlabel("Frequência (Hz)")
    fig.supylabel("Densidade Espectral de Potência (μV²/Hz)")
    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    fig.savefig(saida_png, dpi=150, bbox_inches="tight")
    print(f"  PSD salva: {saida_png}")
    plt.close(fig)


# --- core: 3 painéis da sessão -----------------------------------------------

def plota_sessao_3painel(t_centers, mi_z, theta_power,
                          canal, saida_png,
                          dados, fs, t_inicio_sessao=0.0,
                          comod_geral=None,
                          nomes_canais=None,
                          par_ativo="theta_gamma",
                          notch_hz=60.0):
    """3 painéis (cima→baixo): comodulograma geral | MI time-series | theta power.

    Os painéis MI e theta power compartilham o mesmo eixo x (tempo).
    Clique no painel MI abre zoom de 10s.

    Parâmetros
    ----------
    dados : np.ndarray, shape (n_amostras, n_canais)
        Sinal completo concatenado da sessão (necessário para o click handler
        recortar a janela 10s).
    fs : float
        Frequência de amostragem.
    t_inicio_sessao : float
        Offset (s) entre o ``t=0`` do sinal concatenado e o início "real" do
        registro. Mantido só para log; o comodulograma geral usa tempo global
        do ``dados`` (t=0 no início do .ns2 #1).
    comod_geral : tuple (z_mapa_plot, z_mapa_full, fases_freq, amps_freq,
                t_ini_global, t_fim_global) ou None
                (z_mapa_plot = mascarado por FDR, z_mapa_full = original)
        Resultado de ``calcula_comodulograma_z`` aplicado a uma janela de
        1 min do registro. Se ``None``, o painel topo fica com placeholder.
    """
    ch = canal
    mi = mi_z.get(ch, np.full(len(t_centers), np.nan))
    th = theta_power.get(ch, np.full(len(t_centers), np.nan))

    # máscara de nan
    mask_ok = ~(np.isnan(mi) | np.isnan(th))

    fig, axes = plt.subplots(3, 1, figsize=(14, 10),
                              gridspec_kw={"height_ratios": [2, 1, 1]})
    # NÃO usamos sharex aqui: o painel topo é 2D (fase×amplitude), os de
    # baixo são 1D (tempo). Compartilhar x bagunçaria o comodulograma.
    # Alinhamento visual entre MI e theta_power é feito pelos eixos idênticos
    # construídos a partir do mesmo t_centers.
    fig.suptitle(f"Exploração — canal {ch+1}"
                 f"{f' ({nomes_canais[ch]})' if nomes_canais else ''}",
                 fontsize=13)

    # --- painel 1 (topo): comodulograma geral em janela de 1 min ------
    ax_cz = axes[0]
    if comod_geral is not None:
        z_mapa_plot, z_mapa_full, fases_freq_c, amps_freq_c, \
            t_1min_ini, t_1min_fim = comod_geral
        im = ax_cz.pcolormesh(fases_freq_c, amps_freq_c, z_mapa_plot,
                               shading="auto", cmap="viridis")
        # marcar pico do par ativo (sempre no z_mapa completo, não mascarado)
        try:
            z_pico, fp_pico, fa_pico = z_pico_par(
                z_mapa_full, fases_freq_c, amps_freq_c, par=par_ativo)
            if z_pico is not None and not np.isnan(z_pico):
                ax_cz.scatter([fp_pico], [fa_pico], s=80, c="red",
                              marker="X", edgecolors="white", linewidths=1.5,
                              label=f"pico z={z_pico:.2f}")
                ax_cz.legend(fontsize=8, loc="upper right")
        except Exception:
            pass
        ax_cz.set_title(
            f"Comodulograma geral (1 min, {t_1min_ini:.0f}–{t_1min_fim:.0f}s, FDR q={FDR_Q})"
            f"  [clique no painel MI para zoom 10s]")
        ax_cz.set_ylabel("Frequência de amplitude (Hz)")
        fig.colorbar(im, ax=ax_cz, label="z")
    else:
        t_total = t_centers[-1] if len(t_centers) > 0 else 0
        t_1min_ini = max(0, t_total / 2 - 30)
        t_1min_fim = min(t_total, t_total / 2 + 30)
        ax_cz.set_title(f"Comodulograma geral (1 min, {t_1min_ini:.0f}–{t_1min_fim:.0f}s)"
                        f"  [INDISPONÍVEL]")
        ax_cz.set_ylabel("Frequência de amplitude (Hz)")
        ax_cz.text(0.5, 0.5, "(comodulograma geral não calculado)",
                   ha="center", va="center", transform=ax_cz.transAxes,
                   fontsize=10, color="gray")

    # --- painel 2: MI z-scoredo ----------------------------------------
    ax_mi = axes[1]
    ax_mi.fill_between(t_centers, mi, alpha=0.4, color="coral")
    ax_mi.plot(t_centers, mi, color="firebrick", lw=0.7, label=f"MI z ({par_ativo})")
    ax_mi.axhline(LIMIAR_Z_EVENTO, color="red", ls="--",
                  lw=1, label=f"limiar={LIMIAR_Z_EVENTO}")
    ax_mi.set_ylabel("MI z-scoredo")
    ax_mi.grid(alpha=0.3)
    ax_mi.legend(fontsize=8)

    # --- painel 3 (embaixo): theta power ------------------------------
    ax_th = axes[2]
    ax_th.fill_between(t_centers, th, alpha=0.4, color="gold")
    ax_th.plot(t_centers, th, color="darkorange", lw=0.7, label="θ 4–8 Hz")
    ax_th.set_ylabel("Theta power\n(normalizado)")
    ax_th.set_xlabel("Tempo (s)")
    ax_th.grid(alpha=0.3)
    ax_th.legend(fontsize=8)

    plt.tight_layout(rect=[0, 0.03, 1, 0.96])

    # --- click handler (só funciona se fig for mostrada interativamente) --
    cid = fig.canvas.mpl_connect(
        "button_press_event",
        lambda event: _on_click(event, fig, t_centers, mi,
                                theta_power, ch, fs, dados, t_inicio_sessao,
                                nomes_canais=nomes_canais,
                                saida_base=os.path.dirname(saida_png),
                                par_ativo=par_ativo, notch_hz=notch_hz),
    )

    fig.savefig(saida_png, dpi=150, bbox_inches="tight")
    print(f"  3-painel salvo: {saida_png}")

    return fig, cid


def _on_click(event, fig, t_centers, mi_z, theta_pow,
              canal, fs, dados, t_inicio_sessao, nomes_canais, saida_base,
              par_ativo="theta_gamma", notch_hz=60.0):
    """Callback de click: recorta 10s em torno do clique e mostra o
    comodulograma detalhado da janela.

    Requer que ``fig`` esteja sendo mostrada com ``plt.show()`` — clicar em
    PNGs salvos não dispara callback.
    """
    if event.inaxes is None:
        return
    ax = event.inaxes
    # só dispara se clicar no painel MI (o do meio)
    if ax.get_ylabel() != "MI z-scoredo":
        return
    if event.xdata is None:
        return

    t_click = event.xdata
    idx = np.argmin(np.abs(t_centers - t_click))
    t_center = t_centers[idx]
    z_val = mi_z[idx] if not np.isnan(mi_z[idx]) else None

    print(f"\n>>> Click em t={t_click:.1f}s →"
          f" zoom na janela {t_center-5:.0f}–{t_center+5:.0f}s")

    nome_evento = (f"canal {canal+1}"
                   f"{' (' + nomes_canais[canal] + ')' if nomes_canais else ''}"
                   f" @ t_global={t_center:.1f}s"
                   f"{f', z={z_val:.1f}' if z_val is not None else ''}")

    try:
        fig_zoom, z_pico = gera_zoom_comodulograma_fooof(
            dados, fs, canal, t_center, nomes_canais, saida_base,
            par_ativo=par_ativo, nome_evento=nome_evento, notch_hz=notch_hz)
    except Exception as e:
        print(f"  [zoom] erro: {e}")
        return

    print(f"  (pico z={z_pico:.2f})" if z_pico is not None else "  (sem pico)")
    plt.show(block=False)  # mostra sem bloquear o callback


# --- main ---------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Exploração interativa de sessão — "
                    "PSD + comodulograma + time-series MI/theta",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--pasta_ns2", required=True,
                   help="Pasta 'Basal antes da infusao' (contém 3 .ns2)")
    p.add_argument("--canal", required=True,
                   help="Lista de canais NA CONVENÇÃO DO DATASET MESTRE "
                        "(1-based: mesmo número da coluna 'canal' do CSV, "
                        "ex.: candidatos_vencedores_consolidados.csv -- "
                        "processa_sessao.py nomeia as pastas de saída como "
                        "chan{indice_array+1}). Ex.: 5,16,22. Convertido "
                        "internamente para índice 0-based do array.")
    p.add_argument("--sessao", default="sessao",
                   help="Nome da sessão (usado no prefixo dos arquivos)")
    p.add_argument("--saida", required=True,
                   help="Pasta de saída (ex.: .../RESULTADOS/diagnosticos_interativos/MTESC04_S1)")
    p.add_argument("--win_s", type=float, default=10,
                   help="Janela para MI (default: 10 s)")
    p.add_argument("--passo_s", type=float, default=5,
                   help="Passo da janela deslizante (default: 5 s)")
    p.add_argument("--notch_hz", type=float, default=60.0,
                   help="Frequência da rede elétrica (default: 60 Hz para BR)")
    p.add_argument("--n_surr", type=int, default=N_SURR,
                   help=f"Número de surrogates (default: {N_SURR})")
    p.add_argument("--n_bins", type=int, default=N_BINS,
                   help=f"Número de bins de fase (default: {N_BINS})")
    p.add_argument("--fdr_q", type=float, default=FDR_Q,
                   help=f"Q para FDR (default: {FDR_Q})")
    p.add_argument("--limiar_z", type=float, default=LIMIAR_Z_EVENTO,
                   help=f"Limiar z do pico ΘΓ para marcar evento no time-series "
                        f"(default: {LIMIAR_Z_EVENTO}; permissivo p/ triagem visual — "
                        f"portão rigoroso é o FDR do mapa 2D)")
    p.add_argument("--dpi", type=int, default=150,
                   help="DPI das figuras salvas (default: 150)")
    p.add_argument("--fs", type=float, default=1000.0,
                   help="Frequência de amostragem (default: 1000 Hz)")
    p.add_argument("--max_janelas", type=int, default=0,
                   help="Se >0, trunca o número de janelas do time-series "
                        "(útil para validação rápida). Default: 0 (todas)")
    p.add_argument("--interactive", action="store_true",
                   help="Abre matplotlib GUI no fim: a janela 3-painel fica "
                        "viva e o clique no painel MI abre zoom 10s + salva "
                        "PNG. Use este flag se quiser 'mexer' na figura em "
                        "tempo real. Sem o flag, o script apenas gera PNGs "
                        "estáticos (modo headless).")
    p.add_argument("--par", default="theta_gamma",
                   choices=list(BAND_PAIRS.keys()),
                   help="Qual par usar para rastrear o z_pico no time-series (default: theta_gamma)")
    p.add_argument("--zoom_t_center", type=float, default=None,
                   help="Se definido, gera direto (sem clique/GUI) o zoom "
                        "FOOOF+comodulograma centrado neste tempo (s, global "
                        "da sessão concatenada) para os canais de --canal, e "
                        "encerra sem rodar PSD/timeline completos. Útil para "
                        "testes headless de uma janela específica.")
    return p.parse_args()


def main():
    args = parse_args()

    # parse canais: entrada na convenção do dataset mestre (1-based, ver
    # processa_sessao.py -- chan{c+1}); convertido aqui p/ índice 0-based
    # real do array `dados`. Rótulos/nomes de arquivo abaixo sempre
    # exibem de volta o número 1-based (ch+1) p/ bater com a coluna
    # "canal" do CSV.
    canais_1based = [int(c.strip()) for c in args.canal.split(",")]
    if any(c < 1 for c in canais_1based):
        raise ValueError(f"--canal deve usar a convenção 1-based do dataset "
                         f"(>=1); recebido: {canais_1based}")
    canais = [c - 1 for c in canais_1based]

    # cria pasta de saída
    os.makedirs(args.saida, exist_ok=True)
    prefixo = args.sessao

    print(f"Carregando sessão: {args.pasta_ns2}")
    dados, fs, ids_canais, offsets = concatena_sessao(args.pasta_ns2,
                                                      fs_esperado=args.fs)
    print(f"  {dados.shape[1]} canais, {dados.shape[0]/fs:.0f}s total"
          f" ({dados.shape[0]/fs/60:.1f} min), fs={fs} Hz")

    # sanity check: canais existem
    for ch in canais:
        if ch >= dados.shape[1]:
            print(f"  AVISO: canal {ch+1} não existe (só {dados.shape[1]} canais)."
                  " Pulando.")
            canais.remove(ch)

    if args.zoom_t_center is not None:
        print(f"\n=== Teste headless: zoom FOOOF+comodulograma @ t={args.zoom_t_center}s ===")
        for ch in canais:
            try:
                gera_zoom_comodulograma_fooof(
                    dados, fs, ch, args.zoom_t_center, ids_canais, args.saida,
                    par_ativo=args.par, notch_hz=args.notch_hz)
            except Exception as e:
                print(f"  ERRO no canal {ch+1}: {e}")
        return

    for ch in canais:
        print(f"\n=== Canal {ch+1} ({ids_canais[ch] if ch < len(ids_canais) else '?'}) ===")

        lfp_ch = dados[:, ch].astype(float)

        # 1. PSD
        saida_psd = os.path.join(args.saida, f"psd_{prefixo}_ch{ch+1}.png")
        plota_psd_dupla(lfp_ch, fs, ids_canais[ch]
                        if ch < len(ids_canais) else str(ch+1),
                        saida_psd, notch_hz=args.notch_hz)

    # 2. Time-series de MI + theta power (todos os canais juntos)
    print("\nComputando time-series de MI e theta power ...")

    t_centers, mi_z, th_pow, win_inicios = time_series_pac(
        dados, fs, canais,
        win_s=args.win_s, passo_s=args.passo_s,
        fases_freq=FASES_FREQ, amps_freq=AMPS_FREQ,
        n_surr=args.n_surr, n_bins=args.n_bins,
        notch_hz=args.notch_hz,
        max_janelas=args.max_janelas,
        par_ativo=args.par,
    )

    # 3. Comodulograma geral em 1 min (janela central) — por canal
    t_total = t_centers[-1] if len(t_centers) > 0 else 0
    t_mid = t_total / 2
    # pegar 60s ao redor do ponto médio do registro
    t_ini_1min = max(0, t_mid - 30)
    t_fim_1min = min(t_total, t_mid + 30)
    ini_1min = int(t_ini_1min * fs)
    fim_1min = int(t_fim_1min * fs)

    # comod_geral POR CANAL: o painel topo do 3-painel deve corresponder
    # ao canal cujos MI/theta estão plotados abaixo.
    comod_geral_por_canal = {}
    for ch in canais:
        lfp_1min = dados[ini_1min:fim_1min, ch].astype(float)
        lfp_1min = aplica_notch(lfp_1min, fs, linha_hz=args.notch_hz)
        print(f"  Comodulograma geral (canal {ch+1}) "
              f"em {t_ini_1min:.0f}–{t_fim_1min:.0f}s "
              f"({(t_fim_1min-t_ini_1min):.0f}s) ...")
        try:
            z_mapa_1min, mi_obs, mi_surr = calcula_comodulograma_z(
                lfp_1min, fs, FASES_FREQ, AMPS_FREQ,
                n_surr=args.n_surr, n_bins=args.n_bins, retorna_mi=True)
            p_mapa_1min = p_valores_por_celula(mi_obs, mi_surr)
            sig_1min = bh_fdr_mapa(p_mapa_1min, alpha=args.fdr_q)
            z_mapa_1min_plot = np.where(sig_1min, z_mapa_1min, np.nan)
            comod_geral_por_canal[ch] = (z_mapa_1min_plot, z_mapa_1min,
                                          FASES_FREQ, AMPS_FREQ,
                                          t_ini_1min, t_fim_1min)
            print(f"    OK (z_pico={np.nanmax(z_mapa_1min):.2f},"
                  f" n_sig_FDR={sig_1min.sum()})")
        except Exception as e:
            print(f"    ERRO no comodulograma geral (canal {ch+1}): {e}")
            comod_geral_por_canal[ch] = None

    # salva CSV de time-series
    import csv as _csv
    csv_ts = os.path.join(args.saida, f"timeseries_{prefixo}.csv")
    with open(csv_ts, "w", newline="") as f:
        w = _csv.writer(f)
        w.writerow(["t_center_s"] + [f"mi_z_ch{ch+1}" for ch in canais]
                   + [f"theta_pow_ch{ch+1}" for ch in canais])
        for i, t in enumerate(t_centers):
            row = [f"{t:.2f}"]
            row += [f"{mi_z[ch][i]:.4f}" if not np.isnan(mi_z[ch][i]) else "nan"
                    for ch in canais]
            row += [f"{th_pow[ch][i]:.6f}" if not np.isnan(th_pow[ch][i]) else "nan"
                    for ch in canais]
            w.writerow(row)
    print(f"  CSV time-series salvo: {csv_ts}")

    # 4. Figura 3-painéis por canal
    for ch in canais:
        saida_3p = os.path.join(args.saida,
                                f"sessao_3paineis_{prefixo}_ch{ch+1}.png")
        plota_sessao_3painel(t_centers, mi_z, th_pow, ch,
                             saida_3p,
                             dados=dados, fs=fs,
                             t_inicio_sessao=0.0,
                             comod_geral=comod_geral_por_canal.get(ch),
                             nomes_canais=ids_canais,
                             par_ativo=args.par,
                             notch_hz=args.notch_hz)

    # 5. CSV de eventos (z > limiar)
    csv_ev = os.path.join(args.saida, f"eventos_{prefixo}.csv")
    with open(csv_ev, "w", newline="") as f:
        w = _csv.writer(f)
        w.writerow(["t_center_s", "canal", "mi_z", "theta_power"])
        for ch in canais:
            for i, t in enumerate(t_centers):
                z = mi_z[ch][i]
                th = th_pow[ch][i]
                if not np.isnan(z) and z >= args.limiar_z:
                    w.writerow([f"{t:.2f}", ch+1,
                               f"{z:.4f}",
                               f"{th:.6f}" if not np.isnan(th) else "nan"])
    print(f"  CSV eventos salvo: {csv_ev}")

    n_eventos = 0
    for ch in canais:
        n_eventos += (mi_z[ch] >= args.limiar_z).sum()
    print(f"\nTotal de eventos (z > {args.limiar_z}): {n_eventos}")

    print(f"\nFiguras salvas em: {args.saida}")

    if args.interactive:
        # Modo interativo: reabre o último 3-painel com plt.show()
        # para que o clique no painel MI funcione de verdade.
        plt.ion()
        print("\n>>> Modo interativo ativo <<<")
        print(f"  Figura reaberta: {saida_3p}")
        print("  Clique no painel MI (meio) para zoom 10s → PNG salvo automaticamente.")
        print("  Cada zoom abre nova janela. Feche todas para sair.")
        print("  Pressione ENTER aqui no terminal para finalizar.")
        plt.show()
        input("Pressione ENTER para encerrar...")
    else:
        print("Abra os PNGs para exploração visual.")
        print("Click nos painel MI das figuras 3-painel abre zoom 10s.")


if __name__ == "__main__":
    main()
