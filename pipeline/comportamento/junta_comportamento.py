import pandas as pd
import argparse
import os

def junta_comportamento(dataset_mestre_path, template_path, saida_path):
    print(f"Lendo Dataset Mestre: {dataset_mestre_path}")
    df_mestre = pd.read_csv(dataset_mestre_path)
    
    print(f"Lendo Comportamentos: {template_path}")
    df_comportamento = pd.read_csv(template_path)
    
    # Remove colunas vazias ou irrelevantes do template se existirem (ex: notas vazias)
    # Mantém apenas as chaves e os comportamentos
    colunas_chave = ["sessao", "arquivo", "janela_ini_s", "janela_fim_s"]
    colunas_comportamento = [c for c in df_comportamento.columns if c not in colunas_chave and "video" not in c]
    
    df_comp_limpo = df_comportamento[colunas_chave + colunas_comportamento].copy()
    
    # Como pode haver múltiplas linhas para a mesma janela no mestre (uma por canal e par),
    # usamos um left join. Cada canal/par daquela janela vai receber o mesmo comportamento anotado.
    print("Cruzando dados (Merge)...")
    df_final = pd.merge(df_mestre, df_comp_limpo, on=colunas_chave, how="left")
    
    df_final.to_csv(saida_path, index=False)
    print(f"Sucesso! Dataset final com comportamentos salvo em: {saida_path}")
    
    # Relatório rápido
    if 'comportamento_1' in df_final.columns:
        n_anotados = df_final['comportamento_1'].notna().sum()
        print(f"Total de linhas no dataset: {len(df_final)}")
        print(f"Linhas que receberam anotação de comportamento: {n_anotados} ({(n_anotados/len(df_final))*100:.1f}%)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Junta o dataset mestre PAC com a planilha de comportamentos.")
    parser.add_argument("--mestre", default="../../dataset_mestre_final.csv", help="Caminho do dataset mestre")
    parser.add_argument("--comportamento", default="../../template_comportamento.csv", help="Caminho do CSV de comportamentos")
    parser.add_argument("--saida", default="../../dataset_mestre_COM_COMPORTAMENTO.csv", help="Caminho de saída")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.mestre):
        print(f"Erro: Arquivo mestre não encontrado em {args.mestre}")
    elif not os.path.exists(args.comportamento):
        print(f"Erro: Arquivo de comportamento não encontrado em {args.comportamento}")
    else:
        junta_comportamento(args.mestre, args.comportamento, args.saida)
