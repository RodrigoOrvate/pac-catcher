"""
inspeciona_evento.py -- plota o 1o candidato a ripple de um canal (bruto,
150-250 Hz e sharp-wave 1-30 Hz, +-200 ms em torno do pico).

    python pipeline/etapa1_triagem/inspeciona_evento.py
    python pipeline/etapa1_triagem/inspeciona_evento.py --arquivo <x.ns2> --canal chan5 --limiar_dp 2.5
"""
import argparse
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from pipeline.etapa1_triagem.deteccao_ripple import detecta_eventos_ripple
from pipeline.etapa1_triagem.triagem_coocorrencia import carrega_sinal
from pac_core import workspace

ARQUIVO_PADRAO = os.path.join(workspace.BASE_LAC_NOCI, "MTESC04_NOCI", "MTESC04 -- 1 - infusao - 08-07-2024",
                              "Basal antes da infusao", "20240708-123605-002.ns2")

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    y = filtfilt(b, a, data)
    return y

def butter_lowpass_filter(data, cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y

def inspeciona_candidato(arquivo=ARQUIVO_PADRAO, canal="chan1", limiar_dp=3.0, duracao_min_ms=15.0,
                         saida_dir=None):
    print("Carregando sinal...")
    sinal, fs, _ = carrega_sinal(arquivo, canal)
    print(f"Sinal carregado. fs={fs}")

    print(f"Detectando ripple ({duracao_min_ms:g}ms, DP={limiar_dp:g})...")
    df_ripples, limiar = detecta_eventos_ripple(sinal, fs, limiar_dp=limiar_dp, duracao_min_ms=duracao_min_ms)
    
    if not df_ripples:
        print("Nenhum candidato encontrado!")
        return
        
    print(f"Encontrados {len(df_ripples)} candidatos.")
    
    candidato = df_ripples[0]
    ini_s = candidato['inicio_s']
    fim_s = candidato['fim_s']
    pico_s = candidato['pico_s']
    
    print(f"Candidato: ini={ini_s:.3f}s, fim={fim_s:.3f}s, pico={pico_s:.3f}s, dur={(fim_s-ini_s)*1000:.1f}ms")
    
    # Extrair +- 200ms
    margem = 0.2
    idx_ini = max(0, int((pico_s - margem) * fs))
    idx_fim = min(len(sinal), int((pico_s + margem) * fs))
    
    t_plot = np.arange(idx_ini, idx_fim) / fs
    sinal_recorte = sinal[idx_ini:idx_fim]
    
    # Filtros
    sinal_ripple = butter_bandpass_filter(sinal_recorte, 150, 250, fs, order=4)
    sinal_sw = butter_bandpass_filter(sinal_recorte, 1, 30, fs, order=2) # Sharp-wave (low freq)
    
    plt.figure(figsize=(10, 8))
    
    plt.subplot(3, 1, 1)
    plt.plot(t_plot, sinal_recorte, 'k')
    plt.axvline(pico_s, color='r', linestyle='--', alpha=0.5)
    plt.axvspan(ini_s, fim_s, color='yellow', alpha=0.3)
    plt.title("Sinal Bruto")
    plt.ylabel("uV")
    
    plt.subplot(3, 1, 2)
    plt.plot(t_plot, sinal_ripple, 'b')
    plt.axvline(pico_s, color='r', linestyle='--', alpha=0.5)
    plt.axvspan(ini_s, fim_s, color='yellow', alpha=0.3)
    plt.title("Banda Ripple (150-250Hz)")
    plt.ylabel("uV")
    
    plt.subplot(3, 1, 3)
    plt.plot(t_plot, sinal_sw, 'g')
    plt.axvline(pico_s, color='r', linestyle='--', alpha=0.5)
    plt.axvspan(ini_s, fim_s, color='yellow', alpha=0.3)
    plt.title("Banda Sharp-Wave (1-30Hz)")
    plt.ylabel("uV")
    plt.xlabel("Tempo (s)")
    
    plt.tight_layout()
    
    saida_dir = saida_dir or workspace.figuras("inspecao_ripple")
    os.makedirs(saida_dir, exist_ok=True)
    out_path = os.path.join(saida_dir, f"candidato_ripple_{os.path.splitext(os.path.basename(arquivo))[0]}_{canal}.png")
    plt.savefig(out_path)
    plt.close()
    print(f"Plot salvo em {out_path}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arquivo", default=ARQUIVO_PADRAO, help=".ns2 a inspecionar")
    ap.add_argument("--canal", default="chan1", help="nome do canal (default chan1)")
    ap.add_argument("--limiar_dp", type=float, default=3.0)
    ap.add_argument("--duracao_min_ms", type=float, default=15.0)
    ap.add_argument("--saida_dir", default=None, help="default: SCRIPT/figuras/inspecao_ripple")
    a = ap.parse_args()
    inspeciona_candidato(a.arquivo, a.canal, a.limiar_dp, a.duracao_min_ms, a.saida_dir)
