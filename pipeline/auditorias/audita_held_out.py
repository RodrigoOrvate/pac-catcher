"""
audita_held_out.py — Validação anti-double-dipping para janelas reancoradas

Testa se o MI de um candidato se concentra de fato na ilha reportada, ou se
era inflado por double-dipping (ilha escolhida por energia de theta ap´os ver
onde o MI era alto). Regiões comparadas:
  (a) ILHA: o trecho fino reportado.
  (b) RESTO-DA-JANELA: antes+depois da ilha dentro da janela de 10s.
  (c) JANELA ADJACENTE (opcional, teste mais forte): trecho de mesma duração
      imediatamente anterior, fora do cálculo de energia usado na âncora.
Cada região é convertida em z contra 200 surrogates (deslocamento circular,
semente 42, mínimo 1s), mesmo critério do pipeline. Só é chamado de "validado"
um held-out com z>=3. OBRIGATÓRIO apenas quando a janela foi reancorada (a
ilha está no vencedores.csv via colunas inicio_ilha_s/fim_ilha_s).

Refactor de test_split_half.py (lógica preservada) → por CLI, lendo o
vencedores.csv da sessão (regra de ouro).

    python audita_held_out.py --csv "<sessao>/RESULTADOS/vencedores.csv" \
        --pasta_ns2 "<sessao>/<BASAL>" --saida "<sessao>/RESULTADOS/held_out.csv"
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import signal

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

from pac_core.io import carrega_dados

N_SURR = 200
SEED = 42
Z_THRESHOLD = 3.0


def bandpass_filter(data, fs, low, high, order=4):
    nyq = 0.5 * fs
    low_n = max(low / nyq, 1e-3)
    high_n = min(high / nyq, 0.999)
    b, a = signal.butter(order, [low_n, high_n], btype='band')
    return signal.filtfilt(b, a, data)


def calculate_mi(phase, amp, n_bins=18):
    hist_2d = np.histogram2d(phase, amp, bins=n_bins)[0]
    p_joint = hist_2d / np.sum(hist_2d)
    p_phase = np.sum(p_joint, axis=1)
    p_amp = np.sum(p_joint, axis=0)
    mask = p_joint > 0
    return np.sum(p_joint[mask] * np.log(
        p_joint[mask] / (p_phase[:, None] * p_amp[None, :])[mask]
    ))


def get_phase_amp(segment, fs, f_phase, f_amp):
    phase_sig = bandpass_filter(segment, fs, f_phase - 1, f_phase + 1)
    phase = np.angle(signal.hilbert(phase_sig))
    amp_sig = bandpass_filter(segment, fs, f_amp - 5, f_amp + 5)
    amp = np.abs(signal.hilbert(amp_sig))
    return phase, amp


def mi_with_surrogate_z(segment, fs, f_phase, f_amp, n_surr=N_SURR, seed=SEED):
    if len(segment) < fs * 1:  # menos de 1s: MI nao confiavel
        return np.nan, np.nan, len(segment)

    phase, amp = get_phase_amp(segment, fs, f_phase, f_amp)
    mi_obs = calculate_mi(phase, amp)

    rng = np.random.default_rng(seed)
    n = len(amp)
    min_shift = int(fs * 1)
    surr_mis = np.empty(n_surr)
    for i in range(n_surr):
        shift = rng.integers(min_shift, n - min_shift) if n > 2 * min_shift else rng.integers(1, n)
        amp_shifted = np.roll(amp, shift)
        surr_mis[i] = calculate_mi(phase, amp_shifted)

    mu, sigma = surr_mis.mean(), surr_mis.std()
    z = (mi_obs - mu) / sigma if sigma > 0 else np.nan
    return mi_obs, z, len(segment)


def load_channel(file_path, channel_name):
    dados, fs, canal_ids = carrega_dados(file_path)
    chan_idx = canal_ids.index(channel_name)
    return dados[:, chan_idx], fs


def held_out_validation(file_path, channel_name, full_window, island_window,
                        f_phase, f_amp):
    data, fs = load_channel(file_path, channel_name)
    fs_start, fs_end = full_window
    is_start, is_end = island_window
    full = data[int(fs_start * fs):int(fs_end * fs)]

    before = data[int(fs_start * fs):int(is_start * fs)]
    after = data[int(is_end * fs):int(fs_end * fs)]
    rest = np.concatenate([before, after])
    island = data[int(is_start * fs):int(is_end * fs)]

    mi_island, z_island, n_island = mi_with_surrogate_z(island, fs, f_phase, f_amp)
    mi_rest, z_rest, n_rest = mi_with_surrogate_z(rest, fs, f_phase, f_amp)
    return {
        'island': (mi_island, z_island, n_island),
        'rest_of_window': (mi_rest, z_rest, n_rest),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True, help="vencedores.csv da sessão")
    ap.add_argument("--pasta_ns2", required=True, help="pasta com os .ns2")
    ap.add_argument("--saida", default="held_out.csv")
    ap.add_argument("--z", type=float, default=Z_THRESHOLD,
                    help="limiar z de significância (default 3.0)")
    args = ap.parse_args()

    df = pd.read_csv(args.csv)
    if not {"inicio_ilha_s", "fim_ilha_s"}.issubset(df.columns):
        print("Aviso: vencedores.csv sem colunas inicio_ilha_s/fim_ilha_s — "
              "nenhuma janela foi reancorada. Held-out não é obrigatório. "
              "Grava apenas linhas N/A.")
        df["inicio_ilha_s"] = np.nan
        df["fim_ilha_s"] = np.nan

    linhas = []
    print(f"{'Canal':<10} | {'Regiao':<16} | {'N':<7} | {'MI':<10} | {'z':<8} | Veredito")
    print("-" * 75)

    for _, r in df.iterrows():
        canal = str(r["canal"])
        arquivo = str(r["arquivo"])
        path = os.path.join(args.pasta_ns2, arquivo)
        f_phase = float(r["fase_pico_hz"])
        f_amp = float(r["amp_pico_hz"])
        full = (float(r["inicio_s"]), float(r["fim_s"]))

        if pd.isna(r["inicio_ilha_s"]):
            linhas.append({"rotulo": r.get("rotulo"), "canal": canal,
                           "regiao": "nao_reancorado", "n": np.nan,
                           "mi": np.nan, "z": np.nan,
                           "veredito": "N/A (nao reancorado)"})
            print(f"{canal:<10} | N/A (nao reancorado)")
            continue

        island = (float(r["inicio_ilha_s"]), float(r["fim_ilha_s"]))
        res = held_out_validation(path, canal, full, island, f_phase, f_amp)

        for regiao, (mi, z, n) in res.items():
            if np.isnan(z):
                veredito = "N/A (segmento curto demais)"
            else:
                veredito = "SIGNIFICANT (z>=3)" if z >= args.z else "NOT significant"
            print(f"{canal:<10} | {regiao:<16} | {n:<7} | {mi:<10.6f} | {z:<8.2f} | {veredito}")
            linhas.append({"rotulo": r.get("rotulo"), "canal": canal,
                           "regiao": regiao, "n": n, "mi": round(mi, 6),
                           "z": round(z, 2) if not np.isnan(z) else np.nan,
                           "veredito": veredito})

        z_isl, z_rest = res['island'][1], res['rest_of_window'][1]
        if not np.isnan(z_isl) and not np.isnan(z_rest):
            if z_isl >= args.z and z_rest < args.z:
                print(f"{'':<10} |   -> Padrao FOCAL genuino (sig. na ilha, nao no resto).")
            elif z_isl >= args.z and z_rest >= args.z:
                print(f"{'':<10} |   -> MI tambem sig. FORA da ilha -> efeito menos localizado.")
            else:
                print(f"{'':<10} |   -> Ilha NAO sobrevive ao held-out -> NAO validado.")

    out = pd.DataFrame(linhas)
    out.to_csv(args.saida, index=False)
    print(f"\nSalvo: {args.saida} ({len(out)} linhas)")


if __name__ == "__main__":
    main()
