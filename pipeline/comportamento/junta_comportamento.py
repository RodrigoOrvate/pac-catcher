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
    
    df_final.to_csv(saida_path, index=False, encoding="utf-8-sig")
    print(f"Sucesso! Dataset final com comportamentos salvo em: {saida_path}")
    
    # Relatório rápido
    if 'comportamento_1' in df_final.columns:
        n_anotados = df_final['comportamento_1'].notna().sum()
        print(f"Total de linhas no dataset: {len(df_final)}")
        print(f"Linhas que receberam anotação de comportamento: {n_anotados} ({(n_anotados/len(df_final))*100:.1f}%)")

if __name__ == "__main__":
    modulo_dir = os.path.dirname(os.path.abspath(__file__))
    script_dir = os.path.abspath(os.path.join(modulo_dir, "..", ".."))
    
    padrao_comp = os.path.join(modulo_dir, "template_comportamento.csv")
    if not os.path.exists(padrao_comp):
        padrao_comp = os.path.join(script_dir, "template_comportamento.csv")

    candidatos_mestre = [
        os.path.join(script_dir, "resultados", "dataset_mestre_final.csv"),
        os.path.join(script_dir, "dataset_mestre_final.csv"),
        os.path.join(script_dir, "dataset_mestre_final_v2.csv"),
    ]
    padrao_mestre = next((c for c in candidatos_mestre if os.path.exists(c)), candidatos_mestre[0])
    padrao_saida = os.path.join(script_dir, "resultados", "dataset_mestre_COM_COMPORTAMENTO.csv")

    parser = argparse.ArgumentParser(description="Junta o dataset mestre PAC com a planilha de comportamentos.")
    parser.add_argument("--mestre", default=padrao_mestre, help="Caminho do dataset mestre")
    parser.add_argument("--comportamento", default=padrao_comp, help="Caminho do CSV de comportamentos")
    parser.add_argument("--saida", default=padrao_saida, help="Caminho de saída")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.mestre):
        print(f"Erro: Arquivo mestre não encontrado em {args.mestre}")
    elif not os.path.exists(args.comportamento):
        print(f"Erro: Arquivo de comportamento não encontrado em {args.comportamento}")
    else:
        junta_comportamento(args.mestre, args.comportamento, args.saida)
