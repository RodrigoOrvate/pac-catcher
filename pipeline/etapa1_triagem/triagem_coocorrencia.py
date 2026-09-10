"""
triagem_coocorrencia.py
==========================================
Triagem de coocorrência (Passo 0.5) para evitar rodar PAC em janelas inúteis.
Usa resolução temporal dupla: Welch (10s) para fenômenos contínuos (teta, gama, HG)
e deteccao_ripple (resolução amostral) para eventos rápidos (HFO/Ripple).

Uso:
  python triagem_coocorrencia.py --origem DADOS_EXEMPLO/LFP_HG_HFO.mat --canal lfpBruto
  python triagem_coocorrencia.py --origem dados.ns2 --canal 17
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import scipy.io
from scipy.signal import welch

# Adiciona a raiz do SCRIPT ao path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from pipeline.etapa1_triagem.deteccao_ripple import detecta_eventos_ripple, resume_eventos_por_janela
from pipeline.etapa1_triagem.triagem_pac import teta_ok_por_janela

def carrega_sinal(origem, canal):
    if origem.endswith('.mat'):
        m = scipy.io.loadmat(origem)
        for chave in [canal, 'lfp', 'LFP', 'lfpBruto', 'lfpHG', 'lfpHFO']:
            if chave in m:
                arr = np.squeeze(m[chave]).astype(float)
                fs = float(m.get('fs', m.get('Fs', [[1000.0]]))[0][0])
                return arr, fs, chave
        raise KeyError(f"Nenhuma variável de sinal encontrada em {origem}")
    elif origem.endswith('.ns2'):
        try:
            from pac_core.io import carrega_dados
            dados, fs, canal_ids = carrega_dados(origem)
            idx = canal_ids.index(str(canal)) if str(canal) in canal_ids else int(canal)
            nome_canal = canal_ids[idx]
            return dados[:, idx].astype(float), fs, nome_canal
        except Exception as e:
            raise RuntimeError(f"Erro ao carregar .ns2: {e}")
    else:
        raise ValueError("Formato de origem não suportado. Use .mat ou .ns2")


def triagem_coocorrencia(sinal, fs, limiar_dp=4.0, duracao_ms=25.0, janela_s=10.0, passo_s=5.0):
    duracao_total_s = len(sinal) / fs
    win = int(janela_s * fs)
    step = int(passo_s * fs)
    
    inicios = list(range(0, len(sinal) - win + 1, step))
    n_janelas = len(inicios)
    
    # 1. Resolução temporal fina (Amostras) para HFO/Ripple
    print("Detectando ripples em todo o sinal (resolução de amostra)...")
    eventos_ripple, hfo_cru_mask = detecta_eventos_ripple(
        sinal, fs,
        banda_ripple=(150, 250),
        limiar_dp=limiar_dp,
        duracao_min_ms=duracao_ms,
        banda_sharp_wave=(1, 30),
        exigir_sharp_wave=False  # Inicialmente false, dependerá de calibração
    )
    print(f"  Ripples detectados: {len(eventos_ripple)}")
    
    # 2. Resolução espectral (Welch) para as outras bandas
    print("Varrendo janelas espectrais...")
    
    resultados = []
    
    # Vamos pre-calcular as potências para normalização (opcional)
    # Por hora usaremos um limiar simples de teta_ok e energia relativa da janela
    for n_janela, ini in enumerate(inicios):
        fim = ini + win
        trecho = sinal[ini:fim]
        ini_s = ini / fs
        fim_s = fim / fs
        
        # Teta (usamos a função já existente que checa SNR em relação aos flancos)
        teta_ok = int(teta_ok_por_janela(trecho, fs))
        
        # Potência simples para Gama e HG para ter uma métrica contínua
        nperseg = min(len(trecho), int(2 * fs))
        nfft = max(nperseg, int(4 * fs))
        f_w, psd_w = welch(trecho, fs=fs, nperseg=nperseg, nfft=nfft, window='hann')
        
        def pot_banda(lo, hi):
            mask = (f_w >= lo) & (f_w <= hi)
            return float(np.mean(psd_w[mask])) if mask.sum() > 0 else 1e-12
            
        p_gamma = pot_banda(30, 80)
        p_hg = pot_banda(80, 150)
        p_total = np.mean(psd_w) + 1e-12
        
        # Gama e HG ativos se representarem uma proporção considerável da potência?
        # Por enquanto vamos salvar os valores brutos para poder definir limiar depois,
        # ou definir flags simples com base na média. 
        # A triagem de coocorrência pode ser mais branda com gama/HG se focar em HFO.
        
        # HFO / Ripple
        tem_ripple, tem_hfo_cru = resume_eventos_por_janela(
            eventos_ripple, hfo_cru_mask, fs, ini_s, fim_s
        )
        
        resultados.append({
            "janela_ini_s": round(ini_s, 2),
            "janela_fim_s": round(fim_s, 2),
            "teta_ok": teta_ok,
            "gamma_pot_rel": p_gamma / p_total,
            "hg_pot_rel": p_hg / p_total,
            "hfo_cru": int(tem_hfo_cru),
            "ripple": int(tem_ripple)
        })
        
    df = pd.DataFrame(resultados)
    
    # Determinar limiar adaptativo simples para Gamma/HG
    limiar_gamma = df["gamma_pot_rel"].median() * 1.5
    limiar_hg = df["hg_pot_rel"].median() * 1.5
    
    df["gamma_ok"] = (df["gamma_pot_rel"] > limiar_gamma).astype(int)
    df["hg_ok"] = (df["hg_pot_rel"] > limiar_hg).astype(int)
    
    return df

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--origem", required=True, help="Arquivo .mat ou .ns2")
    ap.add_argument("--canal", required=True, help="Canal a carregar")
    ap.add_argument("--saida", default="coocorrencia.csv")
    ap.add_argument("--limiar_dp", type=float, default=4.0, help="Limiar em DP sobre a mediana para ripple")
    ap.add_argument("--duracao_ms", type=float, default=25.0, help="Duracao minima do ripple em ms")
    args = ap.parse_args()

    sinal, fs, nome_canal = carrega_sinal(args.origem, args.canal)
    print(f"Sinal carregado: {len(sinal)} amostras, fs={fs}Hz")
    
    df = triagem_coocorrencia(sinal, fs, args.limiar_dp, args.duracao_ms)
    df.insert(0, "arquivo", os.path.basename(args.origem))
    df.insert(1, "canal", nome_canal)
    
    df.to_csv(args.saida, index=False)
    
    # Resumo agregado
    print("\n=== RESUMO DE COOCORRÊNCIA ===")
    tot = len(df)
    
    teta_gamma = len(df[(df["teta_ok"] == 1) & (df["gamma_ok"] == 1)])
    teta_hg = len(df[(df["teta_ok"] == 1) & (df["hg_ok"] == 1)])
    teta_hfo_cru = len(df[(df["teta_ok"] == 1) & (df["hfo_cru"] == 1)])
    teta_ripple = len(df[(df["teta_ok"] == 1) & (df["ripple"] == 1)])
    
    hfo_sem_ripple = len(df[(df["hfo_cru"] == 1) & (df["ripple"] == 0)])
    hfo_total = len(df[df["hfo_cru"] == 1])
    
    print(f"Total de janelas: {tot}")
    print(f"  Teta + Gamma   : {teta_gamma} ({teta_gamma/tot*100:.1f}%)")
    print(f"  Teta + HG      : {teta_hg} ({teta_hg/tot*100:.1f}%)")
    print(f"  Teta + HFO(cru): {teta_hfo_cru} ({teta_hfo_cru/tot*100:.1f}%)")
    print(f"  Teta + Ripple  : {teta_ripple} ({teta_ripple/tot*100:.1f}%)")
    print("\n--- Diagnóstico do Paradoxo HFO ---")
    print(f"  HFO cru detectado em: {hfo_total} janelas")
    print(f"  Dessas, rebaixadas (HFO sem Ripple): {hfo_sem_ripple}")
    if hfo_total > 0:
        taxa_rejeicao = hfo_sem_ripple / hfo_total * 100
        print(f"  Taxa de rejeição do HFO cru: {taxa_rejeicao:.1f}%")
        if taxa_rejeicao > 95:
            print("  [ALERTA] Detector de ripple rejeitou quase tudo. Calibrar duracao_min_ms/limiar_dp!")
        elif taxa_rejeicao < 5:
            print("  [ALERTA] Detector de ripple aprovou tudo. Pode estar muito permissivo.")
            
    print(f"\nSalvo em {args.saida}")

if __name__ == "__main__":
    main()
