"""
aplica_portao_banda_larga_mestre.py
===================================
Shim de compatibilidade: a lógica foi movida para
`enriquece_dataset_mestre.py` (etapa 'portao'). Este script continua
existindo para preservar a CLI (posicional, sobrescreve o próprio
arquivo de entrada) que já existia, delegando para o módulo novo.

Aplica o portao de suspeito_banda_larga no dataset_mestre_final_v2.csv:
Para qualquer candidato com suspeito_banda_larga == True que tenha sido rotulado
como 'Candidato robusto', rebaixa o veredito_refino para
'Revisar: suspeito banda larga / envelope identico ao ruido'.
"""

import sys
import pandas as pd

from enriquece_dataset_mestre import etapa_portao


def main(caminho_csv=r"C:\acoplamento_theta-gamma\dataset_mestre_final_v2.csv"):
    print(f"Lendo dataset: {caminho_csv}")
    df = pd.read_csv(caminho_csv)

    df = etapa_portao(df)

    df.to_csv(caminho_csv, index=False)
    print(f"\nSalvo em {caminho_csv}.")


if __name__ == "__main__":
    caminho = sys.argv[1] if len(sys.argv) > 1 else r"C:\acoplamento_theta-gamma\dataset_mestre_final_v2.csv"
    main(caminho)
