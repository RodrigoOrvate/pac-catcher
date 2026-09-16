"""
plot_basal_results.py -- 4 gráficos descritivos do dataset mestre (contagem e
Z por canal, Z vs MVL, kurtose), filtrados por condição.

    python pipeline/utilitarios/plot_basal_results.py
    python pipeline/utilitarios/plot_basal_results.py --condicao 0h_pos --saida_dir figuras/pos0h
"""
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from pac_core import workspace

# Configurações do matplotlib
plt.style.use('default')
sns.set_theme(style="whitegrid", context="talk")

def plotar_resultados_basais(csv_path=None, saida_dir=None, condicao="basal"):
    csv_path = csv_path or os.path.join(workspace.BASE_RESULTADOS, "dataset_mestre_final.csv")
    saida_dir = saida_dir or workspace.figuras("dataset_mestre", condicao or "todas")
    os.makedirs(saida_dir, exist_ok=True)

    df = pd.read_csv(csv_path, encoding="utf-8-sig")
    if condicao and "condicao" in df.columns:
        df = df[df["condicao"] == condicao]
    if len(df) == 0:
        print(f"Nenhuma linha para condicao={condicao!r} em {csv_path}")
        return
    rotulo = f"Condição: {condicao}" if condicao else "Todas as condições"

    # Garante que o canal seja inteiro para ordenação
    df['canal_num'] = df['canal'].astype(int)
    df = df.sort_values(by='canal_num')

    # 1. Contagem de Janelas Significativas por Canal (Distribuição Espacial)
    plt.figure(figsize=(14, 6))
    ax = sns.countplot(data=df, x='canal_num', hue='par', palette='Set2')
    plt.title(f"Quantidade de Eventos Significativos por Canal ({rotulo})", pad=20)
    plt.xlabel("Canal")
    plt.ylabel("Número de Janelas de 10s")
    plt.legend(title='Banda de Acoplamento (PAC)')
    plt.tight_layout()
    plt.savefig(os.path.join(saida_dir, "01_contagem_por_canal.png"), dpi=300)
    plt.close()

    # 2. Força do Acoplamento (Z-score) por Canal
    plt.figure(figsize=(14, 6))
    sns.boxplot(data=df, x='canal_num', y='z_score_refinado', hue='par', palette='Set2', showfliers=False)
    plt.title(f"Força do Acoplamento (Z-Score Refinado) por Canal ({rotulo})", pad=20)
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
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", default=None,
                    help="dataset mestre (default: resultados/dataset_mestre_final.csv)")
    ap.add_argument("--saida_dir", default=None,
                    help="pasta dos PNGs (default: figuras/dataset_mestre/<condicao>)")
    ap.add_argument("--condicao", default="basal",
                    help="filtra a coluna `condicao` (basal, 0h_pos, 1h_pos, 2h_pos); '' = todas")
    a = ap.parse_args()
    plotar_resultados_basais(a.csv, a.saida_dir, a.condicao)
