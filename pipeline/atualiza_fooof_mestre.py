"""
atualiza_fooof_mestre.py
========================
Atualiza o dataset_mestre_final.csv com as melhorias do FOOOF (Kuhn et al. 2026):
1. aperiodic_mode='knee'
2. Ajuste particionado (baixa frequencia 2-45 Hz para teta; alta frequencia 35-250 Hz para gama)
3. Extracao do expoente aperiodico (balanco E/I), parametro knee e R^2.

Nao recalcula comodulogramas nem surrogates (preserva as 13h de processamento).
Carrega cada arquivo .ns2 apenas uma vez e processa em paralelo.
"""

import os
import sys
import glob
import time
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

# Ajustar paths
DIR_PIPELINE = os.path.dirname(os.path.abspath(__file__))
DIR_AUDITORIAS = os.path.join(DIR_PIPELINE, "auditorias")
sys.path.insert(0, DIR_PIPELINE)
sys.path.insert(0, DIR_AUDITORIAS)

from ns2_utils import carrega_dados, fatia_janela
from utils_harmonico import extrai_cf_teta_fooof, extrai_cf_gamma_fooof


def processa_arquivo_tarefas(file_path, tarefas_arquivo, janela_contexto_s=45.0, f_linha=60.0):
    """
    Carrega o arquivo .ns2 uma unica vez e roda o FOOOF para todas as janelas/canais solicitados.
    """
    resultados = []
    try:
        dados, fs, canal_ids = carrega_dados(file_path)
    except Exception as e:
        print(f"[ERRO] Falha ao abrir {file_path}: {e}")
        for t in tarefas_arquivo:
            resultados.append({**t, "sucesso": False, "erro": str(e)})
        return resultados

    canal_ids_str = [str(c) for c in canal_ids]

    for t in tarefas_arquivo:
        canal = t["canal"]
        ini = t["janela_ini_s"]
        fim = t["janela_fim_s"]

        # Determinar index do canal
        chan_idx = None
        c_str = f"chan{int(canal)}"
        if c_str in canal_ids_str:
            chan_idx = canal_ids_str.index(c_str)
        elif str(canal) in canal_ids_str:
            chan_idx = canal_ids_str.index(str(canal))
        else:
            try:
                chan_idx = int(canal) - 1
            except Exception:
                pass

        if chan_idx is None or chan_idx < 0 or chan_idx >= dados.shape[1]:
            resultados.append({**t, "sucesso": False, "erro": f"Canal {canal} invalido"})
            continue

        centro = (ini + fim) / 2.0
        ctx_ini = max(0.0, centro - janela_contexto_s / 2.0)
        ctx_fim = centro + janela_contexto_s / 2.0

        try:
            sinal_ctx = fatia_janela(dados, fs, ctx_ini, ctx_fim)[:, chan_idx]
            res_teta = extrai_cf_teta_fooof(sinal_ctx, fs, f_linha=f_linha)
            res_gamma = extrai_cf_gamma_fooof(sinal_ctx, fs, f_linha=f_linha)

            resultados.append({
                "arquivo": t["arquivo"],
                "canal": canal,
                "janela_ini_s": ini,
                "janela_fim_s": fim,
                "sucesso": True,
                "cf_teta_fooof_v2": res_teta["cf_teta"],
                "erro_teta_fooof_v2": res_teta["erro_ajuste"],
                "expoente_teta_fooof": res_teta["expoente_teta"],
                "knee_teta_fooof": res_teta["knee_teta"],
                "offset_teta_fooof": res_teta["offset_teta"],
                "r2_teta_fooof": res_teta["r2_teta"],
                "knee_valido_teta": res_teta["knee_valido_teta"],
                "cf_gamma_fooof_v2": res_gamma["cf_gamma"],
                "erro_gamma_fooof_v2": res_gamma["erro_ajuste"],
                "expoente_gamma_fooof": res_gamma["expoente_gamma"],
                "knee_gamma_fooof": res_gamma["knee_gamma"],
                "offset_gamma_fooof": res_gamma["offset_gamma"],
                "r2_gamma_fooof": res_gamma["r2_gamma"],
                "knee_valido_gamma": res_gamma["knee_valido_gamma"],
            })
        except Exception as e:
            resultados.append({**t, "sucesso": False, "erro": str(e)})

    return resultados


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
    args = parser.parse_args()

    print(f"=== ATUALIZADOR FOOOF V2 (Knee Mode & Piecewise Fit) ===")
    print(f"Lendo dataset: {args.mestre_csv}")
    df_mestre = pd.read_csv(args.mestre_csv)
    print(f"Total de linhas no dataset: {len(df_mestre)}")

    # Indexar arquivos .ns2
    print(f"Mapeando arquivos .ns2 em: {args.pasta_dados}")
    todos_ns2 = glob.glob(os.path.join(args.pasta_dados, "**", "*.ns2"), recursive=True)
    mapa_ns2 = {os.path.basename(p): p for p in todos_ns2}

    # Janelas unicas a processar
    chaves = ["arquivo", "canal", "janela_ini_s", "janela_fim_s"]
    df_unicos = df_mestre[chaves].drop_duplicates()
    print(f"Janelas unicas a analisar: {len(df_unicos)}")

    # Agrupar por arquivo para ler cada arquivo apenas 1x
    tarefas_por_arquivo = {}
    for _, row in df_unicos.iterrows():
        arq = row["arquivo"]
        if arq not in mapa_ns2:
            print(f"[AVISO] Arquivo nao encontrado em disco: {arq}")
            continue
        if arq not in tarefas_por_arquivo:
            tarefas_por_arquivo[arq] = []
        tarefas_por_arquivo[arq].append({
            "arquivo": arq,
            "canal": row["canal"],
            "janela_ini_s": float(row["janela_ini_s"]),
            "janela_fim_s": float(row["janela_fim_s"]),
        })

    print(f"Arquivos a processar: {len(tarefas_por_arquivo)} usando {args.n_workers} workers...")
    t_inicio = time.time()

    resultados_totais = []
    with ProcessPoolExecutor(max_workers=args.n_workers) as executor:
        futures = {}
        for arq, lista_tarefas in tarefas_por_arquivo.items():
            f_path = mapa_ns2[arq]
            f = executor.submit(processa_arquivo_tarefas, f_path, lista_tarefas,
                                args.janela_contexto_s)
            futures[f] = arq

        concluidos = 0
        for f in as_completed(futures):
            arq = futures[f]
            concluidos += 1
            res = f.result()
            resultados_totais.extend(res)
            print(f"[{concluidos}/{len(tarefas_por_arquivo)}] {arq} concluido ({len(res)} janelas).")

    t_fim = time.time()
    print(f"\nProcessamento concluido em {t_fim - t_inicio:.1f} segundos!")

    # Converter resultados em DataFrame
    df_res = pd.DataFrame(resultados_totais)
    print(f"Total de analises FOOOF geradas: {len(df_res)}")

    # Filtrar apenas colunas uteis para o merge
    cols_novas = [
        "arquivo", "canal", "janela_ini_s", "janela_fim_s",
        "cf_teta_fooof_v2", "erro_teta_fooof_v2", "expoente_teta_fooof",
        "knee_teta_fooof", "offset_teta_fooof", "r2_teta_fooof", "knee_valido_teta",
        "cf_gamma_fooof_v2", "erro_gamma_fooof_v2", "expoente_gamma_fooof",
        "knee_gamma_fooof", "offset_gamma_fooof", "r2_gamma_fooof", "knee_valido_gamma"
    ]
    df_res_merge = df_res[cols_novas].drop_duplicates(subset=chaves)

    # Merge no dataset mestre
    print("Mesclando novos dados com o dataset mestre...")
    df_mestre_v2 = df_mestre.merge(df_res_merge, on=chaves, how="left", validate="many_to_one")

    print(f"Salvando novo dataset mestre em: {args.saida_csv}")
    df_mestre_v2.to_csv(args.saida_csv, index=False)
    print("Salvo com sucesso!")

    # Estatisticas
    print("\n--- RESUMO DAS METRICAS FOOOF V2 ---")
    if "expoente_teta_fooof" in df_mestre_v2.columns:
        valid_exp = df_mestre_v2["expoente_teta_fooof"].dropna()
        valid_r2 = df_mestre_v2["r2_teta_fooof"].dropna()
        print(f"Expoente Teta (E/I ratio) - Media: {valid_exp.mean():.3f} +- {valid_exp.std():.3f} [Min: {valid_exp.min():.2f}, Max: {valid_exp.max():.2f}]")
        print(f"R^2 do ajuste Teta (2-45 Hz) - Media: {valid_r2.mean():.4f}")
        if "knee_valido_teta" in df_mestre_v2.columns:
            kv = df_mestre_v2["knee_valido_teta"]
            print(f"Knee valido (f_joelho dentro da faixa ajustada) - Teta: {kv.mean():.1%} "
                  f"({int(kv.sum())}/{kv.notna().sum()}) — USE SO ESSES pra comparar expoente_teta_fooof entre condicoes.")
    if "expoente_gamma_fooof" in df_mestre_v2.columns:
        valid_exp_g = df_mestre_v2["expoente_gamma_fooof"].dropna()
        valid_r2_g = df_mestre_v2["r2_gamma_fooof"].dropna()
        print(f"Expoente Gama (High-freq)  - Media: {valid_exp_g.mean():.3f} +- {valid_exp_g.std():.3f} [Min: {valid_exp_g.min():.2f}, Max: {valid_exp_g.max():.2f}]")
        print(f"R^2 do ajuste Gama (35-250 Hz) - Media: {valid_r2_g.mean():.4f}")
        if "knee_valido_gamma" in df_mestre_v2.columns:
            kvg = df_mestre_v2["knee_valido_gamma"]
            print(f"Knee valido (f_joelho dentro da faixa ajustada) - Gama: {kvg.mean():.1%} "
                  f"({int(kvg.sum())}/{kvg.notna().sum()}) — USE SO ESSES pra comparar expoente_gamma_fooof entre condicoes.")
    print("==================================================")


if __name__ == "__main__":
    main()
