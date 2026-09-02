"""
testa_metodos_matlab.py
==========================================
Sanity check dos métodos de PAC (modulation index / PLV) usando dados de exemplo
fornecidos por terceiros (.mat), que ficam FORA do SCRIPT.

Testes:
  1. Sinal SINTÉTICO com acoplamento θ-γ conhecido  -> método deve detectar MI alto
     no par correto e ~0 no par independente (controle).
  2. Upload do `Example_plv_modindex.mat` (lfp1, lfp2, lfp_sw) -> rodar o KL-MI do
     projeto e PLV, e reportar se há acoplamento detectável entre os pares.

Isso valida QUE a implementação do pipeline funciona (e como ela se comporta em
dados de terceiros) — NÃO é validação das descobertas de MTESC04/05.

Reutiliza as funções reais do pipeline (filtra_sinal, _mi_de_bin_idx de
comodulogram.py), para testar exatamente o que o projeto usa.

Uso:
    python testa_metodos_matlab.py                 # sintético + exemplo
    python testa_metodos_matlab.py --mat C:/acoplamento_theta-gamma/Example_plv_modindex.mat
"""

import os
import sys
import argparse
import numpy as np
import scipy.io as sio
from scipy import signal

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PIPELINE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
sys.path.insert(0, PIPELINE_DIR)

# Importa as funções REAIS do pipeline
from comodulogram import filtra_sinal, _mi_de_bin_idx


def mi_pareado(sinal, fs, f_theta, f_gamma, n_bins=18, largura_gamma=4.0):
    """KL-MI entre fase da theta e envelope da gamma numa janela (1 canal)."""
    lfp_fase = filtra_sinal(sinal, f_theta - 1.0, f_theta + 1.0, fs)
    lfp_amp = filtra_sinal(sinal, f_gamma - largura_gamma, f_gamma + largura_gamma, fs)
    fase = np.angle(signal.hilbert(lfp_fase))
    env = np.abs(signal.hilbert(lfp_amp))
    bins = np.linspace(-np.pi, np.pi, n_bins + 1)
    bin_idx = np.clip(np.digitize(fase, bins) - 1, 0, n_bins - 1)
    return _mi_de_bin_idx(bin_idx, env, n_bins)


def plv_pareado(sinal, fs, f_theta, f_gamma, n_trial=5):
    """PLV (phase-locking value) entre a fase da theta e a fase do sinal em gamma."""
    lfp_fase = filtra_sinal(sinal, f_theta - 1.0, f_theta + 1.0, fs)
    lfp_gama = filtra_sinal(sinal, f_gamma - 4.0, f_gamma + 4.0, fs)
    f1 = np.angle(signal.hilbert(lfp_fase))
    f2 = np.angle(signal.hilbert(lfp_gama))
    # PLV padrão (fase1 - fase2), com janelamento por trial
    n = len(f1)
    ventana = n // n_trial
    plvs = []
    for i in range(n_trial):
        a = f1[i * ventana:(i + 1) * ventana]
        b = f2[i * ventana:(i + 1) * ventana]
        plvs.append(np.abs(np.mean(np.exp(1j * (a - b)))))
    return float(np.mean(plvs))


def testa_sintetico():
    print("=== Teste 1: sinal sintético com acoplamento conhecido ===")
    fs = 1024.0
    t = np.arange(0, 20, 1 / fs)
    f_theta, f_gamma = 6.0, 45.0

    # Sinal A: theta 6 Hz (onda real, dá fase bem definida) modulando a
    # amplitude de gamma 45 Hz (acoplamento real)
    fase = 2 * np.pi * f_theta * t
    env_mod = 0.5 + 0.5 * np.cos(fase)             # envoltória modulada pela fase theta
    sinal_theta = 1.0 * np.sin(fase)               # oscilação theta real no sinal
    sinal_mod = np.sin(2 * np.pi * f_gamma * t) * (0.5 + 0.8 * env_mod)
    sinal_a = sinal_theta + sinal_mod + 0.1 * np.random.default_rng(0).standard_normal(t.size)

    # Sinal B: theta 6 Hz (igual) + gamma 45 Hz com envelope CONSTANTE
    # (sem modulação pela theta) -> deve dar MI ~ 0
    sinal_b = 1.0 * np.sin(fase) + np.sin(2 * np.pi * f_gamma * t) \
        + 0.1 * np.random.default_rng(1).standard_normal(t.size)

    # Banda gamma LARGA (f_gamma ± 7 Hz): a modulação a 6 Hz num portador de
    # 45 Hz cria bandas laterais em 39 e 51 Hz; banda estreita (±4) as corta.
    mi_mod = mi_pareado(sinal_a, fs, f_theta, f_gamma, largura_gamma=7.0)
    mi_ctrl = mi_pareado(sinal_b, fs, f_theta, f_gamma, largura_gamma=7.0)
    print(f"  MI (theta 6 mod. gamma 45, acoplado)  = {mi_mod:.4f}")
    print(f"  MI (gamma 45 sem modulacao)           = {mi_ctrl:.4f}")

    # KL-MI normalizado (Tort) é de escala pequena: acoplamento forte ≈ 0.01–0.02,
    # controle ≈ 0. Critério RELATIVO, não absoluto.
    ok = (mi_mod > 0.005) and (mi_mod > 3 * mi_ctrl)
    print(f"  -> {'PASSOU: MI distingue acoplado de controle' if ok else 'FALHOU: MI nao distingue'}"
          " (critério: MI_coupled > 0.005 e > 3x controle)")
    return ok


def testa_exemplo(caminho_mat):
    print("\n=== Teste 2: dados de exemplo treinamento (.mat) ===")
    if not os.path.exists(caminho_mat):
        print(f"  Arquivo nao encontrado: {caminho_mat}")
        return False

    m = sio.loadmat(caminho_mat)
    srate = float(m['srate'].ravel()[0])
    lfp1 = m['lfp1'].ravel()
    lfp2 = m['lfp2'].ravel()
    lfp_sw = m['lfp_sw'].ravel()
    print(f"  srate={srate} | lfp1 n={len(lfp1)} ({len(lfp1)/srate:.1f}s)"
          f" | lfp2 n={len(lfp2)} | lfp_sw n={len(lfp_sw)}")

    # Espectro do exemplo: theta forte em 6–7 Hz em todos; lfp1 tem pico gamma
    # em ~51 Hz (os outros não). Testar PAC justamente nesse par.
    print("  Calculando MI (theta 6 Hz) x bandas gamma em cada sinal...")
    for nome, sinal, f_gamma in [("lfp1", lfp1, 51), ("lfp2", lfp2, 51), ("lfp_sw", lfp_sw, 51)]:
        mi = mi_pareado(sinal, srate, 6, f_gamma, largura_gamma=6.0)
        print(f"    {nome}: MI(theta6 x gamma{f_gamma}) = {mi:.4f}")

    # PLV entre lfp1 e lfp_sw (para ver se há bloqueio de fase na banda comum)
    plv = plv_pareado(lfp1, srate, 6, 6)  # PLV entre a fase theta dos dois lfps
    print(f"    PLV(fase theta) lfp1 vs lfp_sw = {plv:.4f}")
    return True


def main():
    ap = argparse.ArgumentParser(description="Sanity check de métodos PAC com dados .mat de exemplo")
    ap.add_argument("--mat", default=r"C:/acoplamento_theta-gamma/DADOS_EXEMPLO/Example_plv_modindex.mat")
    args = ap.parse_args()

    ok1 = testa_sintetico()
    ok2 = testa_exemplo(args.mat)

    print("\n=== Resumo ===")
    print(f"  Sintético distingue acoplado vs controle: {'SIM' if ok1 else 'NAO'}")
    print(f"  Leitura/processamento do exemplo .mat:    {'OK' if ok2 else 'FALHOU'}")
    if ok1:
        print("  -> Método de MI do pipeline está funcionando. Pode aplicar "
              "a mesma lógica nos seus dados (e no preditor).")
    else:
        print("  -> Revisar implementação de MI.")


if __name__ == "__main__":
    main()
