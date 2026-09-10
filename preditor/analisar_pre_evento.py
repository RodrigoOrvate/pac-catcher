"""
analisar_pre_evento.py
Versão 3: Extração robusta de pré-eventos com logging detalhado.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import os
import glob
import sys

# Adicionar pasta pipeline ao path para importar ns2_utils (SCRIPT/preditor/ -> SCRIPT/pipeline/)
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from pac_core.io import carrega_dados

PRE_WINDOW = 10

def load_signal_chunk(file_path, start_s, duration_s, channel_name):
    try:
        dados, fs, canal_ids = carrega_dados(file_path)

        # Encontrar o índice do canal correto
        if channel_name in canal_ids:
            canal_idx = canal_ids.index(channel_name)
        else:
            # Tentar correspondência numérica se o nome não for encontrado exatamente
            try:
                num = int(str(channel_name).replace('chan', ''))
                # Blackrock RawIO frequently names channels as 'Chan1', 'Chan2', etc.
                # or '1', '2', etc. Let's try to find a match.
                for i, cid in enumerate(canal_ids):
                    if str(num) in cid:
                        canal_idx = i
                        break
                else:
                    print(f"    [ERRO] Canal {channel_name} não encontrado no arquivo. Canais disponíveis: {canal_ids}")
                    return None, None, None
            except ValueError:
                print(f"    [ERRO] Canal {channel_name} não encontrado no arquivo. Canais disponíveis: {canal_ids}")
                return None, None, None

        idx_start = int(start_s * fs)
        idx_end = idx_start + int(duration_s * fs)

        print(f"    [DEBUG] Arquivo: {os.path.basename(file_path)}, Duração: {dados.shape[0]/fs:.2f}s, Início: {start_s:.2f}s, Fim: {start_s+duration_s:.2f}s")

        if idx_end <= dados.shape[0] and idx_start >= 0:
            # Retornar janela (amostras x canais), a taxa de amostragem e o índice do canal
            return dados[idx_start:idx_end], fs, canal_idx
        else:
            print(f"    [ERRO] Índices fora dos limites: {idx_start}-{idx_end} (máx {dados.shape[0]})")
            return None, None, None
    except Exception as e:
        print(f"    [ERRO NEO] {e}")
        return None, None, None

def process_session_vencedores():
    base = r"C:\acoplamento_theta-gamma"
    output_dir = os.path.join(base, "ANALISE_PRE_EVENTO")
    os.makedirs(output_dir, exist_ok=True)

    vencedores_files = glob.glob(os.path.join(base, "**", "RESULTADOS", "vencedores.csv"), recursive=True)

    print(f"Encontrados {len(vencedores_files)} arquivos vencedores.csv")

    for v_csv in vencedores_files:
        sessao_dir = os.path.dirname(os.path.dirname(v_csv))
        sessao_name = os.path.basename(sessao_dir)
        estudo_name = os.path.basename(os.path.dirname(sessao_dir))

        print(f"\n--- Sessão: {estudo_name} / {sessao_name} ---")
        try:
            df = pd.read_csv(v_csv)
        except Exception as e:
            print(f"  Erro ao ler CSV: {e}")
            continue

        if 'inicio_s' not in df.columns:
            print("  Coluna 'inicio_s' não encontrada.")
            continue

        for idx, row in df.iterrows():
            canal_name = row['canal']
            inicio = float(row['inicio_s'])
            filename = row['arquivo']

            ns2_path = os.path.join(sessao_dir, "Basal antes da infusao", filename)
            if not os.path.exists(ns2_path):
                ns2_path = os.path.join(sessao_dir, filename)

            print(f"  Processando {canal_name} @ {inicio}s...")

            if not os.path.exists(ns2_path):
                print(f"    Arquivo {filename} não encontrado em {sessao_dir}")
                continue

            sig_pre, fs, canal_idx = load_signal_chunk(ns2_path, inicio - PRE_WINDOW, PRE_WINDOW, canal_name)

            if sig_pre is not None:
                print(f"    Sinal extraído com sucesso. Formato: {sig_pre.shape}")

                # 1. Salvar CSV
                raw_name = f"raw_{estudo_name}_{sessao_name}_{canal_name}.csv"
                raw_path = os.path.join(output_dir, raw_name)
                pd.DataFrame(sig_pre).to_csv(raw_path, index=False, header=[f"ch{i}" for i in range(sig_pre.shape[1])])
                print(f"    CSV salvo: {raw_name}")

                # 2. Salvar Espectrograma
                f, t, Sxx = signal.spectrogram(sig_pre[:, canal_idx], fs, nperseg=256, noverlap=200)
                plt.figure(figsize=(10, 6))
                plt.pcolormesh(t - PRE_WINDOW, f, 10 * np.log10(Sxx), shading='gouraud')
                plt.title(f"Pré-Evento: {estudo_name} {sessao_name} - {canal_name}")
                plt.ylim(0, 100)
                plt.colorbar()
                fig_name = f"pre_{estudo_name}_{sessao_name}_{canal_name}.png"
                plt.savefig(os.path.join(output_dir, fig_name))
                plt.close()
                print(f"    Figura salva: {fig_name}")

if __name__ == "__main__":
    process_session_vencedores()
