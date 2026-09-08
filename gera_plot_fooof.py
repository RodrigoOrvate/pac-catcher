import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
from neo.rawio import BlackrockRawIO
from fooof import FOOOF

# Caminho do arquivo e parametros
ns2_path = r"C:\acoplamento_theta-gamma\LAC_NOCI\MTESC04_NOCI\MTESC04 -- 1 - infusao - 08-07-2024\Basal antes da infusao\20240708-123605-003.ns2"
chan_idx = 10 - 1  # canal 10
fs = 1000.0
janela_ini = 40.0
janela_fim = 50.0

# 1. Carregar o sinal
reader = BlackrockRawIO(filename=ns2_path)
reader.parse_header()
start_sample = int(janela_ini * fs)
end_sample = int(janela_fim * fs)
raw_sigs = reader.get_analogsignal_chunk(block_index=0, seg_index=0,
                                        i_start=start_sample, i_stop=end_sample,
                                        channel_indexes=[chan_idx])
sig = raw_sigs[:, 0].astype(np.float64)

# 2. Calcular o PSD (Welch)
freqs, psd = signal.welch(sig, fs=fs, nperseg=int(fs*2)) # Resolucao 0.5Hz

# 3. Ajustar o FOOOF
fm = FOOOF(peak_width_limits=[1, 8], max_n_peaks=6, min_peak_height=0.3, aperiodic_mode='fixed')
freq_range = [2, 100] # Para mostrar o Teta e o Gama
fm.fit(freqs, psd, freq_range)

# Extrair os componentes do FOOOF para plotagem manual
f = fm.freqs
p = fm.power_spectrum
ap = fm._ap_fit
flat = p - ap

# 4. Criar a figura de comparacao
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# --- Painel A: Espectro Bruto (Sem FOOOF) ---
ax1.plot(f, p, color='black', linewidth=1.5, label='PSD Original')
ax1.set_title("A. Sem FOOOF (Espectro de Potência Bruto)", fontsize=14)
ax1.set_xlabel("Frequência (Hz)", fontsize=12)
ax1.set_ylabel("Log Potência", fontsize=12)
# Marcar onde o Theta deveria estar
ax1.axvspan(4, 10, color='blue', alpha=0.1, label='Banda Theta (4-10 Hz)')
ax1.axvspan(30, 80, color='red', alpha=0.1, label='Banda Gama (30-80 Hz)')
ax1.legend(loc='upper right')
ax1.grid(True, alpha=0.3)

# --- Painel B: Com FOOOF ---
ax2.plot(f, p, color='black', alpha=0.4, linewidth=1.5, label='PSD Original')
ax2.plot(f, ap, color='blue', linestyle='--', linewidth=2, label='Fundo Aperiódico (1/f)')
ax2.plot(f, flat + np.min(p), color='green', linewidth=2, label='Sinal "Achatado" (Picos Reais)')
ax2.set_title("B. Com FOOOF (Separação Periódico/Aperiódico)", fontsize=14)
ax2.set_xlabel("Frequência (Hz)", fontsize=12)
ax2.set_ylabel("Log Potência", fontsize=12)
ax2.legend(loc='upper right')
ax2.grid(True, alpha=0.3)

# Marcar os picos reais encontrados pelo FOOOF
for peak in fm.get_params('peak_params'):
    cf, amp, bw = peak
    ax2.axvline(cf, color='red', linestyle=':', alpha=0.6)
    ax2.text(cf, np.max(p), f'{cf:.1f} Hz', color='red', ha='center', va='bottom', fontsize=10)

plt.tight_layout()
fig.savefig(r"C:\acoplamento_theta-gamma\comparacao_fooof.png", dpi=150)
print("Imagem salva em C:\\acoplamento_theta-gamma\\comparacao_fooof.png")
