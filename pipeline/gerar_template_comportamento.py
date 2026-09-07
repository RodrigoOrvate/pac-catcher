import pandas as pd
import os

def gerar_template_comportamento(csv_path="C:/acoplamento_theta-gamma/dataset_mestre_final.csv", saida="C:/acoplamento_theta-gamma/template_comportamento.csv"):
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
    gerar_template_comportamento()
