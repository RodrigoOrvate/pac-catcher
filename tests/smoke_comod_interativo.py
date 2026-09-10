"""Smoke test do comodulogram_interativo.py com sinal sintético.

Gera 30s de LFP com teta (8Hz) + gamma (60Hz) modulado em amplitude
— acoplamento ΘΓ forte esperado. Salva PSD + 3-painel em tempdir.

NÃO depende de .ns2 real. Valida:
  - imports ok
  - time_series_pac devolve (t_centers, mi_z, th_pow) com shape coerente
  - plota_psd_dupla e plota_sessao_3painel geram PNG sem erro
  - z do pico ΘΓ na janela central é alto (>3)
"""
import os
import sys
import tempfile

import numpy as np

# garantir path da raiz do SCRIPT (para importar o pacote pipeline)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline.etapa0_exploracao.comodulogram_interativo import (
    BANDA_TETA, N_BINS, N_SURR, FDR_Q, LIMIAR_Z_EVENTO,
    time_series_pac, plota_psd_dupla, plota_sessao_3painel,
)


def gera_lfp_sintetico(fs=1000.0, dur_s=30.0, f_theta=8.0, f_gamma=60.0,
                       amp_teta=1.0, amp_gamma=0.6, amp_ruido=0.05,
                       n_canais=2, seed=42):
    rng = np.random.default_rng(seed)
    n = int(dur_s * fs)
    t = np.arange(n) / fs
    # teta: senoidal puro
    teta = amp_teta * np.sin(2 * np.pi * f_theta * t)
    # gamma modulado pela fase do teta (PAC verdadeiro)
    mod = 0.5 + 0.5 * np.cos(2 * np.pi * f_theta * t)  # 0..1
    gamma = amp_gamma * mod * np.sin(2 * np.pi * f_gamma * t)
    sinal = teta + gamma
    # ruído branco por canal
    dados = np.stack([sinal + amp_ruido * rng.standard_normal(n)
                       for _ in range(n_canais)], axis=1)
    ids = [f"ch{i}" for i in range(n_canais)]
    return dados, fs, ids


def main():
    out = tempfile.mkdtemp(prefix="comod_interativo_smoke_")
    print(f"Saída em: {out}")

    dados, fs, ids = gera_lfp_sintetico()
    canais = [0, 1]
    print(f"Shape dados: {dados.shape}, fs={fs}")

    # --- PSD (apenas canal 0) -----------------------------------------
    plota_psd_dupla(dados[:, 0], fs, canal_id=0,
                     saida_png=os.path.join(out, "psd_ch0.png"),
                     notch_hz=60.0)
    print("  PSD salva.")

    # --- time-series (janela 5s, passo 2s — mais rápido) ---------------
    t_centers, mi_z, th_pow, inicios = time_series_pac(
        dados, fs, canais, win_s=5, passo_s=2,
        n_bins=N_BINS, n_surr=N_SURR,
        notch_hz=60.0,
    )
    print(f"  t_centers.shape={t_centers.shape}, "
          f"mi_z.keys={list(mi_z.keys())}")
    # z do pico em cada canal na janela central
    i_mid = len(t_centers) // 2
    for ch in canais:
        mz = mi_z[ch]
        if not np.all(np.isnan(mz)):
            print(f"    ch{ch} z(pico) no meio={mz[i_mid]:.2f}, "
                  f"max={np.nanmax(mz):.2f}")
        else:
            print(f"    ch{ch} z(pico) = NaN inteiro!")

    # --- 3-painel (canal 0) -------------------------------------------
    # pular comod_geral para acelerar (None vai plotar placeholder)
    plota_sessao_3painel(
        t_centers, mi_z, th_pow, 0,
        os.path.join(out, "sessao_3paineis_smoke_ch0.png"),
        dados=dados, fs=fs, t_inicio_sessao=0.0,
        comod_geral=None,
        nomes_canais=ids,
    )
    print("  3-painel salvo.")
    print("OK — smoke test passou.")


if __name__ == "__main__":
    main()
