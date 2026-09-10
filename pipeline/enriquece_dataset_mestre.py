"""
enriquece_dataset_mestre.py
==========================================
Rotina única de enriquecimento do dataset mestre, fundindo o que antes
exigia rodar dois scripts em cadeia manualmente:

  1. Etapa 'fooof'  -- adiciona as 14 colunas de decomposição aperiódica
     FOOOF v2 (aperiodic_mode='knee', ajuste particionado teta/gama,
     Kuhn et al. 2026). Não recalcula comodulogramas nem surrogates
     (preserva o processamento caro já feito). Equivalente ao antigo
     atualiza_fooof_mestre.py.
  2. Etapa 'portao'  -- rebaixa veredito_refino para candidatos com
     suspeito_banda_larga=True que ainda estavam como "Candidato
     robusto". Equivalente ao antigo aplica_portao_banda_larga_mestre.py.

As duas etapas são independentes: 'portao' usa só suspeito_banda_larga e
veredito_refino, colunas que já vêm de agrega_resultados.py -- nenhuma
delas é criada pela etapa 'fooof'. A ordem fooof->portao é fixa aqui só
por seguir a convenção histórica dos nomes de arquivo, não por
dependência real de dados.

Semântica de I/O: NÃO-DESTRUTIVA por padrão (grava em --saida, ou num
caminho derivado de --entrada se --saida não for passado). Use --in_place
para sobrescrever a própria entrada (reproduz o comportamento do antigo
aplica_portao_banda_larga_mestre.py).

Uso:
    # As duas etapas, saída em novo arquivo
    python enriquece_dataset_mestre.py --entrada dataset_mestre.csv

    # Só a etapa do portão, sobrescrevendo a entrada
    python enriquece_dataset_mestre.py --entrada dataset_mestre_v2.csv \\
        --etapas portao --in_place

    # Só recalcula FOOOF v2 (aborta se as colunas já existirem, a menos
    # que --forca seja passado)
    python enriquece_dataset_mestre.py --entrada dataset_mestre.csv \\
        --etapas fooof --pasta_dados C:\\acoplamento_theta-gamma\\LAC_NOCI \\
        --saida dataset_mestre_v2.csv
"""

import os
import sys
import glob
import time
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd

DIR_PIPELINE = os.path.dirname(os.path.abspath(__file__))
DIR_AUDITORIAS = os.path.join(DIR_PIPELINE, "auditorias")
sys.path.insert(0, DIR_PIPELINE)
sys.path.insert(0, DIR_AUDITORIAS)

from ns2_utils import carrega_dados, fatia_janela
from utils_harmonico import extrai_cf_teta_fooof, extrai_cf_gamma_fooof

CHAVE_JANELA = ["arquivo", "canal", "janela_ini_s", "janela_fim_s"]

COLS_FOOOF_V2 = [
    "cf_teta_fooof_v2", "erro_teta_fooof_v2", "expoente_teta_fooof",
    "knee_teta_fooof", "offset_teta_fooof", "r2_teta_fooof", "knee_valido_teta",
    "cf_gamma_fooof_v2", "erro_gamma_fooof_v2", "expoente_gamma_fooof",
    "knee_gamma_fooof", "offset_gamma_fooof", "r2_gamma_fooof", "knee_valido_gamma",
]


# ==========================================
# ETAPA FOOOF (ex-atualiza_fooof_mestre.py)
# ==========================================

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


def etapa_fooof(df_mestre, pasta_dados, janela_contexto_s=45.0, n_workers=4, forca=False):
    """
    Adiciona as 14 colunas FOOOF v2 ao df_mestre (merge left por janela
    única). Levanta SystemExit se as colunas já existirem e forca=False
    -- proteção contra o merge do pandas criar colunas _x/_y duplicadas
    silenciosamente ao reenriquecer um dataset que já passou por esta
    etapa.
    """
    existentes = [c for c in COLS_FOOOF_V2 if c in df_mestre.columns]
    if existentes and not forca:
        raise SystemExit(
            f"[ERRO] O dataset ja tem colunas FOOOF v2 ({existentes}). "
            f"Rodar a etapa 'fooof' de novo criaria colunas duplicadas "
            f"(_x/_y) no merge. Use --forca para recalcular do zero "
            f"(as colunas existentes serao substituidas), ou "
            f"--etapas portao para pular esta etapa."
        )
    if existentes and forca:
        print(f"[--forca] Removendo colunas FOOOF v2 existentes antes de recalcular: {existentes}")
        df_mestre = df_mestre.drop(columns=existentes)

    print("=== ETAPA FOOOF (Knee Mode & Piecewise Fit) ===")
    print(f"Total de linhas no dataset: {len(df_mestre)}")

    print(f"Mapeando arquivos .ns2 em: {pasta_dados}")
    todos_ns2 = glob.glob(os.path.join(pasta_dados, "**", "*.ns2"), recursive=True)
    mapa_ns2 = {os.path.basename(p): p for p in todos_ns2}

    df_unicos = df_mestre[CHAVE_JANELA].drop_duplicates()
    print(f"Janelas unicas a analisar: {len(df_unicos)}")

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

    print(f"Arquivos a processar: {len(tarefas_por_arquivo)} usando {n_workers} workers...")
    t_inicio = time.time()

    resultados_totais = []
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {}
        for arq, lista_tarefas in tarefas_por_arquivo.items():
            f_path = mapa_ns2[arq]
            f = executor.submit(processa_arquivo_tarefas, f_path, lista_tarefas, janela_contexto_s)
            futures[f] = arq

        concluidos = 0
        for f in as_completed(futures):
            arq = futures[f]
            concluidos += 1
            res = f.result()
            resultados_totais.extend(res)
            print(f"[{concluidos}/{len(tarefas_por_arquivo)}] {arq} concluido ({len(res)} janelas).")

    t_fim = time.time()
    print(f"\nProcessamento FOOOF concluido em {t_fim - t_inicio:.1f} segundos!")

    df_res = pd.DataFrame(resultados_totais)
    print(f"Total de analises FOOOF geradas: {len(df_res)}")

    cols_novas = CHAVE_JANELA + COLS_FOOOF_V2
    df_res_merge = df_res[cols_novas].drop_duplicates(subset=CHAVE_JANELA)

    print("Mesclando novos dados FOOOF com o dataset mestre...")
    df_mestre = df_mestre.merge(df_res_merge, on=CHAVE_JANELA, how="left", validate="many_to_one")

    print("\n--- RESUMO DAS METRICAS FOOOF V2 ---")
    if "expoente_teta_fooof" in df_mestre.columns:
        valid_exp = df_mestre["expoente_teta_fooof"].dropna()
        valid_r2 = df_mestre["r2_teta_fooof"].dropna()
        print(f"Expoente Teta (E/I ratio) - Media: {valid_exp.mean():.3f} +- {valid_exp.std():.3f} [Min: {valid_exp.min():.2f}, Max: {valid_exp.max():.2f}]")
        print(f"R^2 do ajuste Teta (2-45 Hz) - Media: {valid_r2.mean():.4f}")
        if "knee_valido_teta" in df_mestre.columns:
            kv = df_mestre["knee_valido_teta"]
            print(f"Knee valido (f_joelho dentro da faixa ajustada) - Teta: {kv.mean():.1%} "
                  f"({int(kv.sum())}/{kv.notna().sum()}) — USE SO ESSES pra comparar expoente_teta_fooof entre condicoes.")
    if "expoente_gamma_fooof" in df_mestre.columns:
        valid_exp_g = df_mestre["expoente_gamma_fooof"].dropna()
        valid_r2_g = df_mestre["r2_gamma_fooof"].dropna()
        print(f"Expoente Gama (High-freq)  - Media: {valid_exp_g.mean():.3f} +- {valid_exp_g.std():.3f} [Min: {valid_exp_g.min():.2f}, Max: {valid_exp_g.max():.2f}]")
        print(f"R^2 do ajuste Gama (35-250 Hz) - Media: {valid_r2_g.mean():.4f}")
        if "knee_valido_gamma" in df_mestre.columns:
            kvg = df_mestre["knee_valido_gamma"]
            print(f"Knee valido (f_joelho dentro da faixa ajustada) - Gama: {kvg.mean():.1%} "
                  f"({int(kvg.sum())}/{kvg.notna().sum()}) — USE SO ESSES pra comparar expoente_gamma_fooof entre condicoes.")
    print("==================================================")

    return df_mestre


# ==========================================
# ETAPA PORTAO (ex-aplica_portao_banda_larga_mestre.py)
# ==========================================

def etapa_portao(df):
    """
    Para qualquer candidato com suspeito_banda_larga == True que tenha
    sido rotulado como 'Candidato robusto', rebaixa o veredito_refino
    para 'Revisar: suspeito banda larga / envelope identico ao ruido'.
    """
    print("=== ETAPA PORTAO (suspeito_banda_larga) ===")
    print(f"Total de linhas: {len(df)}")

    if "suspeito_banda_larga" not in df.columns:
        raise SystemExit("[ERRO] Coluna 'suspeito_banda_larga' nao encontrada.")
    if "veredito_refino" not in df.columns:
        raise SystemExit("[ERRO] Coluna 'veredito_refino' nao encontrada.")

    mask_suspeito = df["suspeito_banda_larga"] == True
    mask_rebaixar = mask_suspeito & (df["veredito_refino"] == "Candidato robusto")
    n_rebaixar = int(mask_rebaixar.sum())

    print(f"Total de candidatos com suspeito_banda_larga=True: {int(mask_suspeito.sum())}")
    print(f"Total de candidatos que serao rebaixados de 'Candidato robusto': {n_rebaixar}")

    print("\nRebaixamentos por par:")
    print(df[mask_rebaixar]["par"].value_counts().to_string())

    novo_veredito = "Revisar: suspeito banda larga / envelope identico ao ruido"
    df.loc[mask_rebaixar, "veredito_refino"] = novo_veredito

    print(f"\n{n_rebaixar} candidatos rebaixados.")
    print("\nDistribuicao final de veredito_refino:")
    print(df["veredito_refino"].value_counts().to_string())

    return df


# ==========================================
# MAIN
# ==========================================

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--entrada", required=True, help="CSV mestre de entrada")
    ap.add_argument("--saida", default=None,
                    help="CSV de saida (default: <entrada sem ext>_enriquecido.csv)")
    ap.add_argument("--in_place", action="store_true",
                    help="Sobrescreve a propria entrada (equivalente ao comportamento "
                         "antigo de aplica_portao_banda_larga_mestre.py)")
    ap.add_argument("--etapas", nargs="+", choices=["fooof", "portao"],
                    default=["fooof", "portao"],
                    help="Quais etapas rodar. Ordem de execucao e sempre fooof->portao, "
                         "independente da ordem passada aqui (nao ha dependencia real "
                         "entre as duas).")
    ap.add_argument("--pasta_dados", default=r"C:\acoplamento_theta-gamma\LAC_NOCI",
                    help="Pasta raiz dos .ns2 (usado só se a etapa 'fooof' estiver ativa)")
    ap.add_argument("--janela_contexto_s", type=float, default=45.0)
    ap.add_argument("--n_workers", type=int, default=4)
    ap.add_argument("--forca", action="store_true",
                    help="Permite recalcular a etapa fooof mesmo se as colunas v2 ja existirem")
    ap.add_argument("--dry_run", action="store_true",
                    help="Só relata o que seria feito, não grava nada")
    args = ap.parse_args()

    if args.in_place and args.saida:
        raise SystemExit("[ERRO] --in_place e --saida são mutuamente exclusivos.")

    if args.in_place:
        saida = args.entrada
    elif args.saida:
        saida = args.saida
    else:
        base, ext = os.path.splitext(args.entrada)
        saida = f"{base}_enriquecido{ext}"

    print(f"Lendo dataset: {args.entrada}")
    df = pd.read_csv(args.entrada)
    print(f"Total de linhas no dataset: {len(df)}")

    if args.dry_run:
        print(f"\n[--dry_run] Etapas que seriam executadas (nesta ordem): "
              f"{[e for e in ['fooof', 'portao'] if e in args.etapas]}")
        print(f"[--dry_run] Saida seria: {saida}")
        if "fooof" in args.etapas:
            existentes = [c for c in COLS_FOOOF_V2 if c in df.columns]
            if existentes:
                print(f"[--dry_run] Colunas FOOOF v2 ja presentes: {existentes} "
                      f"-- {'seriam recalculadas (--forca)' if args.forca else 'ABORTARIA sem --forca'}")
        if "portao" in args.etapas:
            if "suspeito_banda_larga" in df.columns and "veredito_refino" in df.columns:
                mask = (df["suspeito_banda_larga"] == True) & (df["veredito_refino"] == "Candidato robusto")
                print(f"[--dry_run] Candidatos que seriam rebaixados pelo portao: {int(mask.sum())}")
        return

    if "fooof" in args.etapas:
        df = etapa_fooof(df, args.pasta_dados, args.janela_contexto_s, args.n_workers, args.forca)
    if "portao" in args.etapas:
        df = etapa_portao(df)

    print(f"\nSalvando dataset enriquecido em: {saida}")
    df.to_csv(saida, index=False)
    print("Salvo com sucesso!")


if __name__ == "__main__":
    main()
