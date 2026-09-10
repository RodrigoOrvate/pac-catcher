"""
prever_pac_tempo_real.py
==========================================
Predictor "definitivo" de acoplamento theta-gamma em tempo real.

Lê continuamente o sinal de um arquivo .ns2 (ou de um fluxo), recorta uma janela
deslizante, extrai as mesmas features do treinamento (treinar_preditor.py) e
aplica o modelo treinado (modelo_pac.pkl). Quando a probabilidade de PAC iminente
ultrapassa o limiar, dispara um pulso TTL para optogenética.

O disparo TTL é feito por `dispara_ttl()` — uma função plugável. Adapte o corpo
dela ao seu hardware (DAQ National Instruments, Arduino, porta serial, picoDAC,
etc.). Aqui está como STAND-IN que só imprime/loga.

Uso:
    # Replay offline sobre um .ns2 (demonstra janelas que disparariam):
    python prever_pac_tempo_real.py --ns2 <arquivo.ns2> --canal chan16 \
        --modelo pipeline/modelo_pac.pkl --replay

    # Modo streaming abstrato (gera dados sintéticos p/ testar o loop):
    python prever_pac_tempo_real.py --simular
"""

import os
import sys
import argparse
import numpy as np
import joblib
from scipy.signal import welch

# Resolução de caminhos: este script fica em SCRIPT/preditor/, e pac_core
# compartilhado fica em SCRIPT/
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(SCRIPT_DIR, '..'))
from pac_core.io import carrega_dados

FS = 1000          # taxa de amostragem (Hz)
JANELA_S = 10      # duração da janela (mesma do treinamento)
PASSO_S = 1        # avanço da janela deslizante (s) a cada avaliação
LIMIAR = 0.7       # probabilidade mínima p/ disparar TTL

# Colunas de features na MESMA ordem do treinamento
FEATURES = ['Energia', 'Theta_Energy', 'Gamma_Energy', 'Ratio_TG', 'Theta_Peak']


def extrair_features(window_data):
    """Features de uma janela (n_amostras x n_canais) — espelha treinar_preditor."""
    sig = np.mean(window_data, axis=1)
    energia = np.var(sig)
    f, psd = welch(sig, FS, nperseg=min(1000, len(sig)))

    theta_mask = (f >= 4) & (f <= 8)
    theta_energia = np.mean(psd[theta_mask]) if np.any(theta_mask) else 0.0
    theta_pico = f[theta_mask][np.argmax(psd[theta_mask])] if np.any(theta_mask) else 0.0

    gamma_mask = (f >= 30) & (f <= 80)
    gamma_energia = np.mean(psd[gamma_mask]) if np.any(gamma_mask) else 0.0

    razao_tg = theta_energia / (gamma_energia + 1e-6)
    return np.array([energia, theta_energia, gamma_energia, razao_tg, theta_pico])


def dispara_ttl(t_seg):
    """PLUGUE SEU HARDWARE AQUI.

    Este é o momento de acionar o pulso TTL (ex.: 5V, 10 ms) que liga o laser de
    optogenética. Mantenha o pulso curto e debounced para não ficar refletindo.

    Exemplos de backend (descomente o seu):
      - NI-DAQ (nidaqmx):
          import nidaqmx
          with nidaqmx.Task() as task:
              task.do_channels.add_do_chan("Dev1/port0/line0")
              task.write(True); time.sleep(0.01); task.write(False)
      - Arduino via serial:
          import serial, time
          ser = serial.Serial('COM3', 9600)
          ser.write(b'H'); time.sleep(0.01); ser.write(b'L')
    """
    # STAND-IN: só registra. Troque pelo disparo real quando houver hardware.
    print(f"    >>> [TTL] Disparo de pulso em t={t_seg:.1f}s (optogenética) <<<")


def carrega_modelo(caminho_modelo):
    if not caminho_modelo or not os.path.exists(caminho_modelo):
        # fallback para o modelo gerado pelo treinamento
        cand = os.path.join(SCRIPT_DIR, 'modelo_pac.pkl')
        if os.path.exists(cand):
            caminho_modelo = cand
        else:
            raise FileNotFoundError(
                "Modelo não encontrado. Rode primeiro: python treinar_preditor.py")
    return joblib.load(caminho_modelo)


def modulo_replay(ns2_path, canal_nome, modelo, limiar=LIMIAR):
    """Roda a janela deslizante sobre um .ns2 inteiro e marca onde dispararia."""
    dados, fs, canal_ids = carrega_dados(ns2_path)
    if canal_nome in canal_ids:
        ch = canal_ids.index(canal_nome)
    else:
        ch = 0  # fallback: 1º canal

    n_amostras = JANELA_S * int(fs)
    passo = PASSO_S * int(fs)
    print(f"Replay: {os.path.basename(ns2_path)} | canal {canal_nome} | "
          f"{dados.shape[0]/fs:.0f}s | janela={JANELA_S}s passo={PASSO_S}s")

    for inicio in range(0, dados.shape[0] - n_amostras + 1, passo):
        janela = dados[inicio:inicio + n_amostras]
        props = modelo.predict_proba(extrair_features(janela).reshape(1, -1))[0]
        p_pac = props[1] if len(props) == 2 else props[0]
        t_seg = inicio / fs
        marca = "  <-- DISPARARIA" if p_pac >= limiar else ""
        print(f"  t={t_seg:6.1f}s  P(PAC)={p_pac:.3f}{marca}")
        if p_pac >= limiar:
            dispara_ttl(t_seg)


def modulo_simulado(modelo, limiar=LIMIAR):
    """Modo demo: gera streaming sintético em tempo real (theta forte) p/ testar
    o loop de avaliação + disparo sem hardware."""
    import time
    rng = np.random.default_rng(0)
    print("Modo simulado: gerando sinal theta forte para forçar um disparo "
          "demonstrativo (Ctrl+C para encerrar).")
    t = 0.0
    while True:
        # Sinal sintético: theta 6 Hz + ruído
        n = JANELA_S * FS
        tt = np.linspace(t, t + JANELA_S, n, endpoint=False)
        janela = np.sin(2 * np.pi * 6 * tt).reshape(-1, 1) + 0.3 * rng.standard_normal((n, 1))
        props = modelo.predict_proba(extrair_features(janela).reshape(1, -1))[0]
        p_pac = props[1] if len(props) == 2 else props[0]
        if p_pac >= limiar:
            print(f"  t={t:6.1f}s P(PAC)={p_pac:.3f} >>> [TTL] DISPARO <<<")
            dispara_ttl(t)
        else:
            print(f"  t={t:6.1f}s P(PAC)={p_pac:.3f}")
        t += PASSO_S
        time.sleep(0.05)


def main():
    ap = argparse.ArgumentParser(description="Predictor em tempo real de PAC theta-gamma (theta x gamma)")
    ap.add_argument("--ns2", help="caminho do arquivo .ns2 a varrer (replay)")
    ap.add_argument("--canal", default="chan1", help="canal de interesse (ex.: chan16)")
    ap.add_argument("--modelo", help="caminho do modelo .pkl (default: pipeline/modelo_pac.pkl)")
    ap.add_argument("--limiar", type=float, default=LIMIAR, help="prob. mínima p/ disparo")
    ap.add_argument("--replay", action="store_true", help="modo replay sobre --ns2")
    ap.add_argument("--simular", action="store_true", help="modo demo com sinal sintético")
    args = ap.parse_args()

    modelo = carrega_modelo(args.modelo)

    if args.simular:
        modulo_simulado(modelo, args.limiar)
    elif args.ns2:
        modulo_replay(args.ns2, args.canal, modelo, args.limiar)
    else:
        ap.print_help()


if __name__ == "__main__":
    main()
