import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Configurações do matplotlib
plt.style.use('default')
sns.set_theme(style="whitegrid", context="talk")

def plotar_resultados_basais(csv_path="C:/acoplamento_theta-gamma/dataset_mestre_final.csv", saida_dir="C:/acoplamento_theta-gamma/analise_basal"):
    os.makedirs(saida_dir, exist_ok=True)
    
    # Carrega dataset
    df = pd.read_csv(csv_path)
    if len(df) == 0:
        print("Dataset vazio!")
        return

    # Garante que o canal seja inteiro para ordenação
    df['canal_num'] = df['canal'].astype(int)
    df = df.sort_values(by='canal_num')

    # 1. Contagem de Janelas Significativas por Canal (Distribuição Espacial)
    plt.figure(figsize=(14, 6))
    ax = sns.countplot(data=df, x='canal_num', hue='par', palette='Set2')
    plt.title("Quantidade de Eventos Significativos por Canal (Estado Basal)", pad=20)
    plt.xlabel("Canal")
    plt.ylabel("Número de Janelas de 10s")
    plt.legend(title='Banda de Acoplamento (PAC)')
    plt.tight_layout()
    plt.savefig(os.path.join(saida_dir, "01_contagem_por_canal.png"), dpi=300)
    plt.close()

    # 2. Força do Acoplamento (Z-score) por Canal
    plt.figure(figsize=(14, 6))
    sns.boxplot(data=df, x='canal_num', y='z_score_refinado', hue='par', palette='Set2', showfliers=False)
    plt.title("Força do Acoplamento (Z-Score Refinado) por Canal (Estado Basal)", pad=20)
    plt.xlabel("Canal")
    plt.ylabel("Z-Score (Força do PAC)")
    plt.axhline(3.0, color='red', linestyle='--', label='Limiar Z=3.0')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(saida_dir, "02_forca_zscore_por_canal.png"), dpi=300)
    plt.close()

    # 3. Força do Acoplamento vs Rigidez de Fase (MVL)
    plt.figure(figsize=(10, 8))
    sns.scatterplot(data=df, x='mvl', y='z_score_refinado', hue='par', alpha=0.6, palette='Set2')
    plt.title("Relação entre Intensidade (Z-Score) e Travamento de Fase (MVL)", pad=20)
    plt.xlabel("Mean Vector Length (MVL) - Rigidez de Fase")
    plt.ylabel("Z-Score (Força do PAC)")
    plt.tight_layout()
    plt.savefig(os.path.join(saida_dir, "03_zscore_vs_mvl.png"), dpi=300)
    plt.close()

    # 4. Distribuição da Assimetria (Skewness / Kurtose)
    plt.figure(figsize=(10, 6))
    sns.histplot(data=df, x='kurtose_banda', hue='par', kde=True, palette='Set2')
    plt.title("Distribuição da Assimetria do Teta (Kurtose da Banda Larga)", pad=20)
    plt.xlabel("Kurtose (Acima de 3 = Caudas Pesadas / Possível Onda Aguda)")
    plt.ylabel("Frequência")
    plt.tight_layout()
    plt.savefig(os.path.join(saida_dir, "04_distribuicao_kurtose.png"), dpi=300)
    plt.close()

    print(f"Gráficos gerados com sucesso na pasta: {saida_dir}")

if __name__ == '__main__':
    plotar_resultados_basais()
