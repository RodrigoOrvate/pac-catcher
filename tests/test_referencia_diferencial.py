import sys
import os
import numpy as np
from scipy.signal import welch
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "pipeline"))
from pipeline.ns2_utils import carrega_dados
from pipeline.preprocessa_referencia_diferencial import seleciona_pool_referencia, constroi_referencia
from pipeline.deteccao_ripple import detecta_eventos_ripple

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    y = filtfilt(b, a, data)
    return y

def inspeciona_candidato(sinal_bruto, fs, candidato, plot_name):
    ini_s = candidato['inicio_s']
    fim_s = candidato['fim_s']
    pico_s = (ini_s + fim_s) / 2.0
    
    margem = 0.2
    idx_ini = max(0, int((pico_s - margem) * fs))
    idx_fim = min(len(sinal_bruto), int((pico_s + margem) * fs))
    
    t_plot = np.arange(idx_ini, idx_fim) / fs
    sinal_recorte = sinal_bruto[idx_ini:idx_fim]
    
    sinal_ripple = butter_bandpass_filter(sinal_recorte, 150, 250, fs, order=4)
    sinal_sw = butter_bandpass_filter(sinal_recorte, 1, 30, fs, order=2)
    
    plt.figure(figsize=(10, 8))
    
    plt.subplot(3, 1, 1)
    plt.plot(t_plot, sinal_recorte, 'k')
    plt.axvline(pico_s, color='r', linestyle='--', alpha=0.5)
    plt.axvspan(ini_s, fim_s, color='yellow', alpha=0.3)
    plt.title(f"Sinal Bruto Diferencial (Tamanho {fim_s - ini_s:.3f}s)")
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
    
    out_path = os.path.join(r"C:\Users\rodri\.gemini\antigravity-ide\brain\4887ffb1-11dd-4c40-9405-5bdadf2068d3", plot_name)
    plt.savefig(out_path)
    plt.close()

def detecta_flags_simples(dados, fs):
    n_canais = dados.shape[1]
    flags = {}
    p_gamas = []
    p_hgs = []
    
    for c in range(n_canais):
        f, psd = welch(dados[:, c], fs=fs, nperseg=int(2*fs), nfft=int(4*fs))
        p_gamas.append(np.mean(psd[(f >= 30) & (f <= 80)]))
        p_hgs.append(np.mean(psd[(f >= 80) & (f <= 150)]))
        
    lim_g = np.median(p_gamas) * 1.5
    lim_h = np.median(p_hgs) * 1.5
    
    for c in range(n_canais):
        f_dict = {}
        if p_gamas[c] > lim_g:
            f_dict['teta_gama'] = True
        if p_hgs[c] > lim_h:
            f_dict['teta_hg'] = True
        flags[f'chan{c+1}'] = f_dict
        
    return flags

def roda_testes_diferencial_loo(arquivo, limiar_dp=3.0, duracao_ms=15.0):
    print(f"\n======================================")
    print(f"Carregando {os.path.basename(arquivo)}...")
    dados, fs, canais = carrega_dados(arquivo)
    
    flags = detecta_flags_simples(dados, fs)
    canal_alvo_idx = 0 # chan1
    sinal_alvo = dados[:, canal_alvo_idx]
    
    # Baseline sem diferencial
    df_base, _ = detecta_eventos_ripple(sinal_alvo, fs, limiar_dp=limiar_dp, duracao_min_ms=duracao_ms)
    cont_base = len(df_base) if df_base else 0
    print(f"[Baseline Sem Referencial]: {cont_base} ripples")
    
    for tamanho in [4, 8, 16]:
        print(f"\n--- POOL = {tamanho} ---")
        pool = seleciona_pool_referencia(dados, fs, flags, canal_alvo_idx=canal_alvo_idx, n_canais_pool=tamanho)
        sinal_ref = constroi_referencia(dados, pool)
        df_dif, _ = detecta_eventos_ripple(sinal_alvo, fs, limiar_dp=limiar_dp, duracao_min_ms=duracao_ms, sinal_referencia=sinal_ref)
        cont = len(df_dif) if df_dif else 0
        print(f"[Diferencial Completo Pool={tamanho}]: {cont} ripples (Pool: {pool})")
        
        if cont > 0:
            inspeciona_candidato(sinal_alvo - sinal_ref, fs, df_dif[0], f"cand_pool{tamanho}_{os.path.basename(arquivo)}.png")
            
        print("  Rodando Leave-One-Out (removendo 1 por vez):")
        for p in pool:
            pool_loo = [c for c in pool if c != p]
            sinal_ref_loo = constroi_referencia(dados, pool_loo)
            df_loo, _ = detecta_eventos_ripple(sinal_alvo, fs, limiar_dp=limiar_dp, duracao_min_ms=duracao_ms, sinal_referencia=sinal_ref_loo)
            cont_loo = len(df_loo) if df_loo else 0
            print(f"  [-{p}]: {cont_loo} ripples")
            if cont_loo > 0:
                inspeciona_candidato(sinal_alvo - sinal_ref_loo, fs, df_loo[0], f"cand_pool{tamanho}_sem_{p}_{os.path.basename(arquivo)}.png")

if __name__ == "__main__":
    base_dir = r"C:\acoplamento_theta-gamma\MTESC04_NOCI\MTESC04 -- 1 - infusao - 08-07-2024\Basal antes da infusao"
    
    roda_testes_diferencial_loo(os.path.join(base_dir, "20240708-123605-002.ns2"), limiar_dp=3.0, duracao_ms=15.0)
    roda_testes_diferencial_loo(os.path.join(base_dir, "20240708-123605-003.ns2"), limiar_dp=3.0, duracao_ms=15.0)
