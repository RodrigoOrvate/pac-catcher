"""
prever_pac_tempo_real.py
==================================================
Módulo de tempo real / replay para PAC teta-gama em CA1.

Baseado no modelo validado de estado comportamental + potência teta
(modelo_estado_comportamental.pkl, AUC-ROC 0.618), com limiar calibrado
por F1 (0.431 -> recall ~90% do vencedor com precisão ~45%).

Papel no Circuito Fechado (Closed-Loop):
Atua como braço de GATING (habilitação) de estado para optogenética:
monitora a dinâmica de teta rápido (tipo 1, 7-10 Hz) e estado comportamental
recente. Quando P(PAC) >= limiar, dispara o sinal TTL para armar/acionar
a estimulação em malha fechada.

Uso:
    # Replay offline sobre um .ns2 (demonstra janelas que disparariam):
    python prever_pac_tempo_real.py --ns2 <arquivo.ns2> --canal chan16 --replay

    # Modo streaming demonstrativo com sinal sintético:
    python prever_pac_tempo_real.py --simular
"""

import os
import sys
import argparse
import time
import numpy as np
import pandas as pd
import joblib
from scipy.signal import welch

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(SCRIPT_DIR, '..'))
from pac_core.io import carrega_dados
from pac_core.filtering import aplica_notch

FS = 1000               # taxa de amostragem padrão (Hz)
JANELA_S = 10           # duração da janela anterior N-1 (s)
PASSO_S = 1             # avanço da janela deslizante (s)
N_SUBJANELAS = 5        # fatias temporais para sub-janela próxima e tendência
BANDA_LENTA = (4.0, 7.0)    # teta tipo 2 (sniffing / atenção)
BANDA_RAPIDA = (7.0, 10.0)  # teta tipo 1 (locomoção)
NOTCH_HZ = [60.0, 120.0, 180.0, 240.0]

CATEGORIAS_COMPORTAMENTO = [
    "Exploração / Locomoção",
    "Grooming / Limpeza",
    "Imóvel / Descanso",
    "Movimento de Cabeça",
    "Rearing / Em pé",
    "Sniffing / Farejando",
    "Sono"
]


def banda_power(sig, fs, banda):
    if len(sig) < 20:
        return 1e-6
    f, psd = welch(sig, fs, nperseg=min(500, len(sig)))
    m = (f >= banda[0]) & (f <= banda[1])
    return float(np.mean(psd[m])) if np.any(m) else 1e-6


def clean_str(s):
    import unicodedata
    nfkd = unicodedata.normalize('NFKD', str(s))
    return ''.join(c for c in nfkd if not unicodedata.combining(c) and c.isalnum()).lower()


def extrair_features_estado(sig_full, fs, comportamento_str, feature_names):
    """Extrai exatamente as features esperadas por modelo_estado_comportamental.pkl."""
    sig_limpo = aplica_notch(sig_full, fs, NOTCH_HZ)
    lenta_whole = banda_power(sig_limpo, fs, BANDA_LENTA)
    rapida_whole = banda_power(sig_limpo, fs, BANDA_RAPIDA)

    # 5 sub-janelas para "close" (última) e "tendência" (inclinação)
    bordas = np.linspace(0, len(sig_limpo), N_SUBJANELAS + 1).astype(int)
    lentas_sub, rapidas_sub = [], []
    for j in range(N_SUBJANELAS):
        sub = sig_limpo[bordas[j]:bordas[j + 1]]
        lentas_sub.append(banda_power(sub, fs, BANDA_LENTA))
        rapidas_sub.append(banda_power(sub, fs, BANDA_RAPIDA))

    x_sub = np.arange(N_SUBJANELAS)
    incl_lenta = float(np.polyfit(x_sub, lentas_sub, 1)[0])
    incl_rapida = float(np.polyfit(x_sub, rapidas_sub, 1)[0])

    close_lenta = lentas_sub[-1]
    close_rapida = rapidas_sub[-1]

    valores_limpos = {
        clean_str("log_teta_lenta_whole"): np.log(lenta_whole + 1e-6),
        clean_str("log_teta_rapida_whole"): np.log(rapida_whole + 1e-6),
        clean_str("teta_lenta_tendencia"): incl_lenta,
        clean_str("teta_rapida_tendencia"): incl_rapida,
        clean_str("log_teta_lenta_close"): np.log(close_lenta + 1e-6),
        clean_str("log_teta_rapida_close"): np.log(close_rapida + 1e-6),
    }

    # One-hot de comportamento
    comport_alvo = clean_str(comportamento_str)
    for c in CATEGORIAS_COMPORTAMENTO:
        chave_col = clean_str(f"comport_{c}")
        valores_limpos[chave_col] = 1.0 if clean_str(c) in comport_alvo or comport_alvo in clean_str(c) else 0.0

    # Monta vetor na mesma ordem exigida pelo modelo
    vetor = []
    for f in feature_names:
        chave_f = clean_str(f)
        vetor.append(valores_limpos.get(chave_f, 0.0))

    return np.array(vetor).reshape(1, -1)


def dispara_ttl(t_seg):
    """PLUGUE SEU HARDWARE DE OPTOGENÉTICA AQUI.

    Aciona o pulso TTL (ex.: 5V, 10 ms) para habilitar o laser de optogenética.
    Exemplos de backend:
      - NI-DAQ (nidaqmx):
          import nidaqmx
          with nidaqmx.Task() as task:
              task.do_channels.add_do_chan("Dev1/port0/line0")
              task.write(True); time.sleep(0.01); task.write(False)
      - Arduino via serial:
          import serial
          ser = serial.Serial('COM3', 9600)
          ser.write(b'H'); time.sleep(0.01); ser.write(b'L')
    """
    print(f"    >>> [TTL] Disparo de pulso em t={t_seg:.1f}s (Optogenética / Closed-Loop) <<<")


def carrega_modelo(caminho_modelo=None):
    if not caminho_modelo or not os.path.exists(caminho_modelo):
        padrao = os.path.join(SCRIPT_DIR, 'modelo_estado_comportamental.pkl')
        if os.path.exists(padrao):
            caminho_modelo = padrao
        else:
            legado = os.path.join(SCRIPT_DIR, 'modelo_pac.pkl')
            if os.path.exists(legado):
                caminho_modelo = legado
            else:
                raise FileNotFoundError("Nenhum modelo encontrado em SCRIPT/preditor/.")

    obj = joblib.load(caminho_modelo)
    if isinstance(obj, dict) and "modelo" in obj:
        return obj["modelo"], obj.get("features", []), obj.get("limiar_f1", 0.431)
    else:
        # Modelo legado (sklearn puro)
        return obj, [], 0.50


def resolve_canal(canal_arg, canal_ids):
    """Aceita nome nativo ('chan16') OU indice 1-based numerico ('16'/16),
    igual a convencao usada no resto do pipeline (dataset mestre). Falha
    alto em vez de cair silenciosamente no canal 0 quando nao resolve --
    mesmo tipo de bug ja corrigido em comodulogram_interativo.py/
    figura_apresentacao.py (resolucao por nome nativo cega ao formato
    1-based do dataset mestre)."""
    s = str(canal_arg)
    if s in canal_ids:
        return canal_ids.index(s)
    try:
        idx = int(s.lower().replace("chan", "")) - 1
    except ValueError:
        raise ValueError(f"Canal {canal_arg!r} nao reconhecido. Canais nativos: {canal_ids}")
    if not (0 <= idx < len(canal_ids)):
        raise ValueError(f"Canal {canal_arg!r} fora do range (1-{len(canal_ids)}).")
    return idx


def modulo_replay(ns2_path, canal_nome, modelo, features, limiar, comportamento="Exploração / Locomoção"):
    """Roda a janela deslizante sobre um .ns2 inteiro e avalia decisões."""
    dados, fs, canal_ids = carrega_dados(ns2_path)
    ch = resolve_canal(canal_nome, canal_ids)
    n_amostras = JANELA_S * int(fs)
    passo = PASSO_S * int(fs)
    print(f"Replay: {os.path.basename(ns2_path)} | canal {canal_nome} | "
          f"{dados.shape[0]/fs:.0f}s | janela={JANELA_S}s passo={PASSO_S}s | limiar={limiar:.3f}")

    for inicio in range(0, dados.shape[0] - n_amostras + 1, passo):
        janela = dados[inicio:inicio + n_amostras, ch].astype(float)
        X = extrair_features_estado(janela, fs, comportamento, features)
        p_pac = modelo.predict_proba(X)[0, 1]
        t_seg = inicio / fs
        marca = "  <-- DISPARARIA TTL" if p_pac >= limiar else ""
        print(f"  t={t_seg:6.1f}s  P(PAC)={p_pac:.3f}{marca}")
        if p_pac >= limiar:
            dispara_ttl(t_seg)


def modulo_simulado(modelo, features, limiar, max_passos=15):
    """Modo demo: simula streaming em tempo real com transição de repouso para locomoção."""
    rng = np.random.default_rng(42)
    print("=" * 70)
    print(f"Modo simulado: demonstrando streaming com transição comportamental")
    print(f"Limiar de disparo calibrado: {limiar:.3f}")
    print("=" * 70)

    t = 0.0
    for passo in range(max_passos):
        n = JANELA_S * FS
        tt = np.linspace(t, t + JANELA_S, n, endpoint=False)

        # Simula transição: primeiros passos em Sono (teta baixo, inativo), depois locomoção (teta 8Hz ativo)
        if passo < 5:
            comportamento = "Sono"
            sinal = 0.05 * np.sin(2 * np.pi * 1.5 * tt) + 0.1 * rng.standard_normal(n)
        else:
            comportamento = "Exploração / Locomoção"
            # Teta rápido acelerando na locomoção voluntária
            freq = 7.5 + min(1.5, (passo - 5) * 0.2)
            sinal = 1.8 * np.sin(2 * np.pi * freq * tt) + 0.3 * rng.standard_normal(n)

        X = extrair_features_estado(sinal, FS, comportamento, features)
        p_pac = modelo.predict_proba(X)[0, 1]

        marca = " >>> [TTL] DISPARO <<<" if p_pac >= limiar else ""
        print(f"  t={t:5.1f}s | Estado: {comportamento:<22} | P(PAC)={p_pac:.3f}{marca}")
        if p_pac >= limiar:
            dispara_ttl(t)

        t += PASSO_S
        time.sleep(0.08)

    print("\nDemonstração de streaming concluída com sucesso.")


def main():
    ap = argparse.ArgumentParser(description="Predictor de PAC theta-gamma em tempo real / replay")
    ap.add_argument("--ns2", help="caminho do arquivo .ns2 a varrer (replay)")
    ap.add_argument("--canal", default="chan1", help="canal de interesse (ex.: chan16)")
    ap.add_argument("--modelo", help="caminho do modelo .pkl")
    ap.add_argument("--limiar", type=float, default=None, help="prob. mínima p/ disparo")
    ap.add_argument("--comportamento", default="Exploração / Locomoção", help="estado comportamental da janela")
    ap.add_argument("--replay", action="store_true", help="modo replay sobre --ns2")
    ap.add_argument("--simular", action="store_true", help="modo demo com sinal sintético")
    args = ap.parse_args()

    modelo, features, limiar_calibrado = carrega_modelo(args.modelo)
    limiar = args.limiar if args.limiar is not None else limiar_calibrado

    if args.simular:
        modulo_simulado(modelo, features, limiar)
    elif args.ns2:
        modulo_replay(args.ns2, args.canal, modelo, features, limiar, args.comportamento)
    else:
        ap.print_help()


if __name__ == "__main__":
    main()
