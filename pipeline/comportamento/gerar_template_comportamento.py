import pandas as pd
import os

import argparse

def gerar_template_comportamento(csv_path, saida):
    if not os.path.exists(csv_path):
        print(f"Dataset não encontrado: {csv_path}")
        return
        
    df = pd.read_csv(csv_path)
    if len(df) == 0:
        print("Dataset vazio!")
        return

    # Queremos anotar apenas as janelas de tempo únicas onde HOUVE detecção de PAC
    # Se uma mesma janela de 10s teve PAC no canal 1 e no canal 2, o rato fez o mesmo comportamento.
    # Então agrupamos pelas janelas temporais exclusivas.
    
    janelas_unicas = df[['sessao', 'arquivo', 'janela_ini_s', 'janela_fim_s']].drop_duplicates()
    
    # Ordena cronologicamente
    janelas_unicas = janelas_unicas.sort_values(by=['sessao', 'arquivo', 'janela_ini_s'])
    
    # Adiciona a coluna vazia para o usuário preencher
    janelas_unicas['comportamento'] = ""
    janelas_unicas['observacao_movimento_cabo'] = ""
    
    # Salva o template
    janelas_unicas.to_csv(saida, index=False)
    print(f"Template de comportamento gerado com {len(janelas_unicas)} janelas exclusivas para anotação.")
    print(f"Salvo em: {saida}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Gera o template CSV vazio para anotação de comportamento.")
    parser.add_argument("--csv_mestre", default="dataset_mestre_final.csv", help="Caminho do CSV mestre (entrada)")
    parser.add_argument("--saida", default="template_comportamento.csv", help="Caminho do template (saída)")
    args = parser.parse_args()
    
    gerar_template_comportamento(args.csv_mestre, args.saida)
