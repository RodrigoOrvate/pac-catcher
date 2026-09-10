"""
atualiza_fooof_mestre.py
========================
Shim de compatibilidade: a lógica foi movida para
`enriquece_dataset_mestre.py` (etapa 'fooof'). Este script continua
existindo para preservar a CLI e os defaults de caminho que já existiam,
delegando para o módulo novo.

Atualiza o dataset_mestre_final.csv com as melhorias do FOOOF (Kuhn et al. 2026):
1. aperiodic_mode='knee'
2. Ajuste particionado (baixa frequencia 2-45 Hz para teta; alta frequencia 35-250 Hz para gama)
3. Extracao do expoente aperiodico (balanco E/I), parametro knee e R^2.

Nao recalcula comodulogramas nem surrogates (preserva as 13h de processamento).
Carrega cada arquivo .ns2 apenas uma vez e processa em paralelo.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from pipeline.dataset_mestre.enriquece_dataset_mestre import etapa_fooof

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Atualiza FOOOF no dataset mestre.")
    parser.add_argument("--mestre_csv", default=r"C:\acoplamento_theta-gamma\dataset_mestre_final.csv",
                        help="CSV mestre atual")
    parser.add_argument("--saida_csv", default=r"C:\acoplamento_theta-gamma\dataset_mestre_final_v2.csv",
                        help="CSV mestre atualizado")
    parser.add_argument("--pasta_dados", default=r"C:\acoplamento_theta-gamma\LAC_NOCI",
                        help="Pasta raiz dos .ns2")
    parser.add_argument("--janela_contexto_s", type=float, default=45.0)
    parser.add_argument("--n_workers", type=int, default=4)
    parser.add_argument("--forca", action="store_true",
                        help="Recalcula mesmo se as colunas FOOOF v2 ja existirem no mestre_csv")
    args = parser.parse_args()

    print(f"Lendo dataset: {args.mestre_csv}")
    df_mestre = pd.read_csv(args.mestre_csv)

    df_mestre_v2 = etapa_fooof(df_mestre, args.pasta_dados, args.janela_contexto_s,
                                args.n_workers, args.forca)

    print(f"Salvando novo dataset mestre em: {args.saida_csv}")
    df_mestre_v2.to_csv(args.saida_csv, index=False)
    print("Salvo com sucesso!")


if __name__ == "__main__":
    main()
