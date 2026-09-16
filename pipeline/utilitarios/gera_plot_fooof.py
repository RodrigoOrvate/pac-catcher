"""
gera_plot_fooof.py
==========================================
Gera a figura de 3 painéis comparando o espectro bruto (sem FOOOF) contra
o ajuste FOOOF v2 (knee + fit particionado teta/gama) sobre um trecho de
LFP real. Uso pontual para apresentação/lab meeting -- não é parte do
pipeline automático (`processa_sessao.py` não chama este script).

Os defaults reproduzem o caso da figura original (MTESC04 sessão 1,
003.ns2, chan10, 40-50 s); qualquer outro caso via CLI:

    python pipeline/utilitarios/gera_plot_fooof.py
    python pipeline/utilitarios/gera_plot_fooof.py --arquivo <x.ns2> --canal 13 --ini 100 --fim 110
"""
import argparse
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal

# Configurar caminhos
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from pac_core.io import carrega_dados, fatia_janela, resolve_canal_idx
from pac_core import workspace
from pipeline.auditorias.utils_harmonico import _knee_valido
from pipeline.auditorias.linha_noise_kuhn import MODOS, aplica_modo
from fooof import FOOOF

F_LINHA = 60.0
LARGURA_LINHA_HZ = 2.0  # faixa que o modo 'hibrido'/'cirurgica' repõe pela 1/f (remove_faixa_1f)

# 1. Parâmetros e Arquivo
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--arquivo", default=os.path.join(workspace.BASE_LAC_NOCI, "MTESC04_NOCI",
                                                  "MTESC04 -- 1 - infusao - 08-07-2024",
                                                  "Basal antes da infusao", "20240708-123605-003.ns2"))
ap.add_argument("--canal", default="chan10", help="nome nativo (chan10) ou número 1-based do dataset mestre")
ap.add_argument("--ini", type=float, default=40.0, help="início da janela (s)")
ap.add_argument("--fim", type=float, default=50.0, help="fim da janela (s)")
ap.add_argument("--ctx_s", type=float, default=45.0, help="contexto em torno do centro da janela (s)")
ap.add_argument("--saida", default=None, help="PNG de saída (default: SCRIPT/figuras/fooof/comparacao_fooof_v2.png)")
ap.add_argument("--modo_linha", default="hibrido", choices=list(MODOS) + ["nenhum"],
                help="limpeza do ruído de rede no PSD antes do FOOOF (linha_noise_kuhn, mesmo default "
                     "de audita_harmonico.py). 'hibrido' só repõe ±2 Hz em 60/120/180/240 Hz pela 1/f "
                     "local, e só se houver pico ESTREITO (0,5-2 Hz) -- gama largo perto de 120 Hz fica. "
                     "'nenhum' = PSD cru (figura antiga).")
args = ap.parse_args()

ns2_path = args.arquivo
centro = (args.ini + args.fim) / 2.0
ctx_s = args.ctx_s

print("Carregando LFP...")
dados, fs, canal_ids = carrega_dados(ns2_path)
chan_idx = resolve_canal_idx(canal_ids, args.canal)
sinal_ctx = fatia_janela(dados, fs, centro - ctx_s / 2.0, centro + ctx_s / 2.0)[:, chan_idx]

# 2. Executar FOOOF V2 (Piecewise + Knee)
print("Ajustando FOOOF V2 (Baixa e Alta frequência)...")
nperseg = int(1.2 * fs)
nfft = 4000 if nperseg <= 4000 else nperseg
freqs, psd = signal.welch(sinal_ctx, fs=fs, window='hann', nperseg=nperseg, noverlap=nperseg//2, nfft=nfft)
psd_cru = psd
harmonicos = []
if args.modo_linha != "nenhum":
    psd = aplica_modo(freqs, psd_cru, fs, args.modo_linha, f_linha=F_LINHA)
    harmonicos = [F_LINHA * h for h in range(1, MODOS[args.modo_linha][1] + 1) if F_LINHA * h < fs / 2]
    print(f"Ruído de rede limpo no PSD (modo '{args.modo_linha}': {', '.join(f'{h:.0f}' for h in harmonicos)} Hz)")


def perto_da_rede(cf, tol=3.0):
    """Pico a menos de `tol` Hz de um harmônico da rede: não dá pra separar
    gama real de resíduo de linha só pelo espectro -- rotular com cautela."""
    return any(abs(cf - F_LINHA * h) <= tol for h in range(1, int(fs / 2 // F_LINHA) + 1))


def marca_linha(ax, f_max, rotular=False):
    for i, h in enumerate(h for h in harmonicos if h <= f_max):
        ax.axvspan(h - LARGURA_LINHA_HZ, h + LARGURA_LINHA_HZ, facecolor='none', edgecolor='#7f8c8d',
                   hatch='///', linewidth=0, alpha=0.6,
                   label=f'Rede limpa (±{LARGURA_LINHA_HZ:g} Hz, {args.modo_linha})' if rotular and i == 0 else None)

# FOOOF Teta (2-45 Hz com knee)
fm_teta = FOOOF(aperiodic_mode='knee', peak_width_limits=(2, 5), min_peak_height=0.05, max_n_peaks=4)
fm_teta.fit(freqs, psd, freq_range=[2.0, 45.0])

# FOOOF Gama (35-150 Hz com knee)
fm_gamma = FOOOF(aperiodic_mode='knee', peak_width_limits=(4, 30), min_peak_height=0.05, max_n_peaks=5)
fm_gamma.fit(freqs, psd, freq_range=[35.0, 150.0])

# 3. Criar Figura Elegante (3 Painéis)
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5.5))

# Estilos globais
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'

# --- PAINEL A: Espectro Bruto Sem FOOOF (2 - 150 Hz) ---
mask_total = (freqs >= 2.0) & (freqs <= 150.0)
f_tot = freqs[mask_total]
p_tot = np.log10(psd[mask_total])

if harmonicos:
    ax1.plot(f_tot, np.log10(psd_cru[mask_total]), color='#95a5a6', linewidth=1.0, label='PSD cru (com rede)')
ax1.plot(f_tot, p_tot, color='#2c3e50', linewidth=2.0,
         label='PSD (Welch, rede limpa)' if harmonicos else 'PSD Original (Welch)')
marca_linha(ax1, 150.0, rotular=True)
ax1.axvspan(4, 10, color='#3498db', alpha=0.18, label='Banda Teta (4–10 Hz)')
ax1.axvspan(30, 80, color='#e74c3c', alpha=0.15, label='Banda Gama (30–80 Hz)')
ax1.axvspan(80, 150, color='#9b59b6', alpha=0.15, label='Banda HG (80–150 Hz)')

ax1.set_title("A. Sem FOOOF (Espectro de Potência Bruto)\nIlusão 1/f: potências altas mascaram oscilações rápidas", fontsize=12, fontweight='bold')
ax1.set_xlabel("Frequência (Hz)", fontsize=11)
ax1.set_ylabel(r"Log$_{10}$ Potência ($\mu$V$^2$/Hz)", fontsize=11)
ax1.legend(loc='upper right', frameon=True, fontsize=9)
ax1.grid(True, linestyle=':', alpha=0.5)

# Limites unificados do eixo Y para os 3 painéis (escala de potência homogênea),
# tirados dos dados (antes eram fixos em 1.0-4.85, calibrados só pro caso default)
_curvas = [p_tot, np.log10(psd_cru[mask_total]), fm_teta.power_spectrum, fm_teta.power_spectrum - fm_teta._ap_fit + 1.15,
           fm_gamma.power_spectrum, fm_gamma.power_spectrum - fm_gamma._ap_fit + 1.15]
_lo = min(float(np.min(c)) for c in _curvas)
_hi = max(float(np.max(c)) for c in _curvas)
Y_LIM = (_lo - 0.1, _hi + 0.6)

ax1.set_ylim(Y_LIM)
ax2.set_ylim(Y_LIM)
ax3.set_ylim(Y_LIM)

# --- PAINEL B: FOOOF V2 Banda Baixa (2 - 45 Hz, Knee Mode) ---
f_t = fm_teta.freqs
p_t = fm_teta.power_spectrum
ap_t = fm_teta._ap_fit
flat_t = p_t - ap_t

ax2.plot(f_t, p_t, color='#2c3e50', alpha=0.5, linewidth=1.5, label='PSD Real')
ax2.plot(f_t, ap_t, color='#2980b9', linestyle='--', linewidth=2.2, label='Fundo Aperiódico (Knee)')
ax2.plot(f_t, flat_t + 1.15, color='#27ae60', linewidth=2.0, label='Picos Puros (Achatado)')

# Marcar picos detectados
picos_t = fm_teta.get_params('peak_params')
if len(picos_t) > 0:
    for p in picos_t:
        cf, amp, bw = p
        if 4.0 <= cf <= 12.0:
            ax2.vlines(cf, Y_LIM[0], Y_LIM[1] - 0.25, color='#e74c3c', linestyle=':', linewidth=1.8)
            ax2.text(cf, Y_LIM[1] - 0.2, f"Pico Teta: {cf:.2f} Hz", color='#c0392b', ha='center', fontsize=9, fontweight='bold')

exp_t = fm_teta.get_params('aperiodic_params')[2]
knee_t = fm_teta.get_params('aperiodic_params')[1]
r2_t = fm_teta.get_params('r_squared')

ax2.set_title("B. FOOOF V2: Teta (2–45 Hz com Knee)\nModela curvatura temporal hipocampal sem distorção", fontsize=12, fontweight='bold')
ax2.set_xlabel("Frequência (Hz)", fontsize=11)
ax2.set_ylabel(r"Log$_{10}$ Potência", fontsize=11)
ax2.legend(loc='upper right', frameon=True, fontsize=9)
ax2.grid(True, linestyle=':', alpha=0.5)

# Box informativo posicionado diretamente abaixo da legenda no canto superior direito
knee_ok_t = "Sim" if _knee_valido(knee_t, exp_t, (2.0, 45.0)) else "Não"
texto_t = f"FOOOF V2 (Baixa Freq)\n$R^2$ = {r2_t:.4f}\nExpoente E/I ($\\chi$) = {exp_t:.2f}\nKnee ($k$) = {knee_t:.1f}\nKnee Válido = {knee_ok_t}"
ax2.text(0.97, 0.74, texto_t, transform=ax2.transAxes, fontsize=8.5, ha='right', va='top',
         bbox=dict(boxstyle='round,pad=0.45', facecolor='#f8f9fa', edgecolor='#bdc3c7', alpha=0.95))

# --- PAINEL C: FOOOF V2 Banda Alta (35 - 150 Hz) ---
f_g = fm_gamma.freqs
p_g = fm_gamma.power_spectrum
ap_g = fm_gamma._ap_fit
flat_g = p_g - ap_g

ax3.plot(f_g, p_g, color='#2c3e50', alpha=0.5, linewidth=1.5, label='PSD Real')
ax3.plot(f_g, ap_g, color='#8e44ad', linestyle='--', linewidth=2.2, label='Fundo Aperiódico Alta-Freq')
ax3.plot(f_g, flat_g + 1.15, color='#e67e22', linewidth=2.0, label='Picos Gama Isolados')
marca_linha(ax3, 150.0, rotular=True)

picos_g = fm_gamma.get_params('peak_params')
if len(picos_g) > 0:
    for n_pico, p in enumerate(sorted(picos_g, key=lambda p: p[0])):
        cf, amp, bw = p
        if cf >= 35.0:
            # Pico colado num harmônico da rede (mesmo após a limpeza): pode ser gama
            # real na mesma frequência ou resíduo de linha -- o espectro sozinho não
            # separa os dois, então fica marcado em cinza pra conferência, não como gama.
            cor, rot = (('#7f8c8d', f"{cf:.1f} Hz\n(rede? conferir)") if perto_da_rede(cf)
                        else ('#d35400', f"{cf:.1f} Hz"))
            topo = float(np.max(flat_g + 1.15)) + 0.1 + 0.18 * (n_pico % 2)  # alterna altura: picos vizinhos não se sobrepõem
            ax3.vlines(cf, Y_LIM[0], topo, color=cor, linestyle=':', linewidth=1.5)
            ax3.text(cf, topo + 0.1, rot, color=cor, ha='center', fontsize=8, fontweight='bold')

exp_g = fm_gamma.get_params('aperiodic_params')[2]
r2_g = fm_gamma.get_params('r_squared')

ax3.set_title("C. FOOOF V2: Gama / HG (35–150 Hz)\nDesacoplamento espectral: oscilações rápidas puras", fontsize=12, fontweight='bold')
ax3.set_xlabel("Frequência (Hz)", fontsize=11)
ax3.set_ylabel(r"Log$_{10}$ Potência", fontsize=11)
ax3.legend(loc='upper right', frameon=True, fontsize=9)
ax3.grid(True, linestyle=':', alpha=0.5)

# Box informativo posicionado diretamente abaixo da legenda no canto superior direito
texto_g = f"FOOOF V2 (Alta Freq)\n$R^2$ = {r2_g:.4f}\nExpoente Alta-Freq = {exp_g:.2f}\nDesacoplado de Teta"
ax3.text(0.97, 0.74, texto_g, transform=ax3.transAxes, fontsize=8.5, ha='right', va='top',
         bbox=dict(boxstyle='round,pad=0.45', facecolor='#f8f9fa', edgecolor='#bdc3c7', alpha=0.95))

plt.tight_layout()
saida_png = args.saida or os.path.join(workspace.figuras("fooof"), "comparacao_fooof_v2.png")
fig.savefig(saida_png, dpi=200)
print(f"Figura gerada e salva com sucesso em: {saida_png}")
