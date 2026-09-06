import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from deteccao_ripple import detecta_eventos_ripple
from triagem_coocorrencia import carrega_sinal

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

def inspeciona_candidato():
    arquivo = r"C:\acoplamento_theta-gamma\MTESC04_NOCI\MTESC04 -- 1 - infusao - 08-07-2024\Basal antes da infusao\20240708-123605-002.ns2"
    canal = "chan1"
    
    print("Carregando sinal...")
    sinal, fs, _ = carrega_sinal(arquivo, canal)
    print(f"Sinal carregado. fs={fs}")
    
    print("Detectando ripple (15ms, DP=3.0)...")
    df_ripples, limiar = detecta_eventos_ripple(sinal, fs, limiar_dp=3.0, duracao_min_ms=15.0)
    
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
    
    out_path = r"C:\Users\rodri\.gemini\antigravity-ide\brain\4887ffb1-11dd-4c40-9405-5bdadf2068d3\candidato_ripple.png"
    plt.savefig(out_path)
    print(f"Plot salvo em {out_path}")

if __name__ == "__main__":
    inspeciona_candidato()
