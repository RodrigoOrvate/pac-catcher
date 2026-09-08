"""
aplica_portao_banda_larga_mestre.py
===================================
Aplica o portao de suspeito_banda_larga no dataset_mestre_final_v2.csv:
Para qualquer candidato com suspeito_banda_larga == True que tenha sido rotulado
como 'Candidato robusto', rebaixa o veredito_refino para
'Revisar: suspeito banda larga / envelope identico ao ruido'.
"""

import sys
import pandas as pd

def main(caminho_csv=r"C:\acoplamento_theta-gamma\dataset_mestre_final_v2.csv"):
    print(f"Lendo dataset: {caminho_csv}")
    df = pd.read_csv(caminho_csv)
    print(f"Total de linhas: {len(df)}")

    if "suspeito_banda_larga" not in df.columns:
        print("[ERRO] Coluna 'suspeito_banda_larga' nao encontrada.")
        return

    # Contar antes
    mask_suspeito = df["suspeito_banda_larga"] == True
    mask_rebaixar = mask_suspeito & (df["veredito_refino"] == "Candidato robusto")
    n_rebaixar = int(mask_rebaixar.sum())

    print(f"Total de candidatos com suspeito_banda_larga=True: {int(mask_suspeito.sum())}")
    print(f"Total de candidatos que serao rebaixados de 'Candidato robusto': {n_rebaixar}")

    # Detalhamento por par
    print("\nRebaixamentos por par:")
    print(df[mask_rebaixar]["par"].value_counts().to_string())

    # Aplicar rebaixamento
    novo_veredito = "Revisar: suspeito banda larga / envelope identico ao ruido"
    df.loc[mask_rebaixar, "veredito_refino"] = novo_veredito

    df.to_csv(caminho_csv, index=False)
    print(f"\nSucesso! {n_rebaixar} candidatos rebaixados e salvos em {caminho_csv}.")

    print("\nDistribuicao final de veredito_refino:")
    print(df["veredito_refino"].value_counts().to_string())

if __name__ == "__main__":
    caminho = sys.argv[1] if len(sys.argv) > 1 else r"C:\acoplamento_theta-gamma\dataset_mestre_final_v2.csv"
    main(caminho)
