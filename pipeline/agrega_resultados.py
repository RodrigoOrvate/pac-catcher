import os
import argparse
import pandas as pd

def agrega_resultados(pasta_saida_base, saida_mestre="dataset_mestre.csv"):
    """
    Varre a árvore RESULTADOS/<sessao>/chan<N>/refinados.csv
    (e seus respectivos harmonico.csv, skewness.csv etc.) para concatenar num dataset final.
    Extrai a condição experimental do nome da sessão.
    """
    print(f"[AGREGADOR] Lendo {pasta_saida_base}...")
    
    if not os.path.exists(pasta_saida_base):
        print("Pasta base não encontrada.")
        return
        
    sessoes = [d for d in os.listdir(pasta_saida_base) if os.path.isdir(os.path.join(pasta_saida_base, d))]
    
    dfs = []
    
    for sessao in sessoes:
        condicao = "basal"
        # Parse rudimentar da condicao
        if "LAC" in sessao.upper(): condicao = "LAC"
        elif "NOCI" in sessao.upper(): condicao = "NOCI"
        
        caminho_sessao = os.path.join(pasta_saida_base, sessao)
        canais = [d for d in os.listdir(caminho_sessao) if d.startswith("chan")]
        
        for canal in canais:
            refinados_path = os.path.join(caminho_sessao, canal, "refinados.csv")
            if not os.path.exists(refinados_path):
                continue
                
            try:
                df = pd.read_csv(refinados_path)
                if len(df) == 0:
                    continue
                    
                df['sessao'] = sessao
                df['condicao'] = condicao
                df['canal'] = canal.replace("chan", "")
                
                dfs.append(df)
            except Exception as e:
                print(f"Erro lendo {refinados_path}: {e}")
                
    if not dfs:
        print("Nenhum dado agregado encontrado.")
        return
        
    df_mestre = pd.concat(dfs, ignore_index=True)
    df_mestre.to_csv(saida_mestre, index=False)
    print(f"[AGREGADOR] Tabela mestre gerada com sucesso: {saida_mestre} ({len(df_mestre)} linhas agregadas).")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Agregador de resultados PAC para Machine Learning.")
    parser.add_argument("--resultados", required=True, help="Pasta contendo as sessões processadas")
    parser.add_argument("--saida", default="dataset_mestre.csv", help="Caminho do CSV de saída")
    args = parser.parse_args()
    
    agrega_resultados(args.resultados, args.saida)
