import os
import argparse
import pandas as pd

# Chave que identifica uma janela única (independe do par).
CHAVE_JANELA = ["arquivo", "canal", "janela_ini_s", "janela_fim_s"]
# Chave que identifica uma janela-por-par (a granularidade de refinados.csv).
CHAVE_JANELA_PAR = CHAVE_JANELA + ["par"]


def _le_csv_seguro(caminho):
    if not os.path.exists(caminho):
        return None
    try:
        df = pd.read_csv(caminho)
        return df if len(df) > 0 else None
    except Exception as e:
        print(f"  [AVISO] Erro lendo {caminho}: {e}")
        return None


def agrega_canal(pasta_chan):
    """
    Lê refinados.csv (base, long-format: 1 linha por janela-por-par) e faz
    left-merge com skewness.csv, resumo_comodulogramas.csv, harmonico.csv e
    harmonico_hfo.csv da mesma pasta de canal.

    Retorna (df_mestre_canal, n_base) para permitir checagem de sanidade
    (nenhuma linha deveria ser criada ou perdida nos merges — todos são
    left-join contra a base 'refinados').
    """
    refinados_path = os.path.join(pasta_chan, "refinados.csv")
    df = _le_csv_seguro(refinados_path)
    if df is None:
        return None, 0

    n_base = len(df)

    if "veredito" in df.columns:
        df = df.rename(columns={"veredito": "veredito_refino"})

    # --- skewness.csv: nao depende de 'par', so pode ter ate 3x mais linhas
    # que janelas unicas (uma por par em refinados.csv). Deduplicar antes do
    # merge é obrigatório, senão o merge explode (cada linha de df bate com
    # N linhas duplicadas de skew para a mesma janela).
    df_skew = _le_csv_seguro(os.path.join(pasta_chan, "skewness.csv"))
    if df_skew is not None:
        if "janela_tipo" in df_skew.columns:
            df_skew = df_skew[df_skew["janela_tipo"] == "full"]
        faltando = [c for c in CHAVE_JANELA if c not in df_skew.columns]
        if faltando:
            print(f"  [AVISO] skewness.csv sem colunas {faltando} — pulando merge (schema antigo?).")
        else:
            df_skew = df_skew.drop_duplicates(subset=CHAVE_JANELA)
            cols = CHAVE_JANELA + [c for c in ["skewness", "veredito_skew", "n"] if c in df_skew.columns]
            df = df.merge(df_skew[cols], on=CHAVE_JANELA, how="left", validate="many_to_one")

    # --- resumo_comodulogramas.csv: 1 linha por janela-por-par (mesma
    # granularidade da base), chave precisa incluir 'par'.
    df_comod = _le_csv_seguro(os.path.join(pasta_chan, "resumo_comodulogramas.csv"))
    if df_comod is not None:
        faltando = [c for c in CHAVE_JANELA_PAR if c not in df_comod.columns]
        if faltando:
            print(f"  [AVISO] resumo_comodulogramas.csv sem colunas {faltando} — pulando merge.")
        else:
            dup = df_comod.duplicated(subset=CHAVE_JANELA_PAR).sum()
            if dup:
                print(f"  [AVISO] resumo_comodulogramas.csv tem {dup} linhas duplicadas na chave "
                      f"janela+par — mantendo a primeira ocorrência.")
                df_comod = df_comod.drop_duplicates(subset=CHAVE_JANELA_PAR)
            # 'z_score_refinado' e 'rotulo' sao ecoados de refinados.csv por
            # comodulogram.py — ja existem na base (df), entao precisam ser
            # excluidos aqui, senao colidem e o pandas cria z_score_refinado_x/_y
            # silenciosamente (bug real encontrado em plot_basal_results.py,
            # que assumia 'z_score_refinado' como coluna unica).
            cols_ja_na_base = set(df.columns) - set(CHAVE_JANELA_PAR)
            cols_extra = [c for c in df_comod.columns
                         if c not in CHAVE_JANELA_PAR and c not in cols_ja_na_base]
            df = df.merge(df_comod[CHAVE_JANELA_PAR + cols_extra], on=CHAVE_JANELA_PAR,
                          how="left", validate="many_to_one")

    # --- harmonico.csv (audita_harmonico.py): idem, 1 linha por janela-por-par.
    df_harm = _le_csv_seguro(os.path.join(pasta_chan, "harmonico.csv"))
    if df_harm is not None:
        faltando = [c for c in CHAVE_JANELA_PAR if c not in df_harm.columns]
        if faltando:
            print(f"  [AVISO] harmonico.csv sem colunas {faltando} — pulando merge.")
        else:
            df_harm = df_harm.drop_duplicates(subset=CHAVE_JANELA_PAR)
            cols_ja_na_base = set(df.columns) - set(CHAVE_JANELA_PAR)
            cols_extra = [c for c in df_harm.columns
                         if c not in CHAVE_JANELA_PAR and c not in ("rotulo", "janela")
                         and c not in cols_ja_na_base]
            df = df.merge(df_harm[CHAVE_JANELA_PAR + cols_extra], on=CHAVE_JANELA_PAR,
                          how="left", validate="many_to_one")

    # --- harmonico_hfo.csv (audita_harmonico_hfo.py): so existe/faz sentido
    # para linhas com par == 'theta_hfo'; a chave com 'par' garante que so
    # se anexa as linhas certas (nao contamina theta_gamma/theta_hg).
    df_hfo = _le_csv_seguro(os.path.join(pasta_chan, "harmonico_hfo.csv"))
    if df_hfo is not None:
        faltando = [c for c in CHAVE_JANELA_PAR if c not in df_hfo.columns]
        if faltando:
            print(f"  [AVISO] harmonico_hfo.csv sem colunas {faltando} — pulando merge.")
        else:
            df_hfo = df_hfo.drop_duplicates(subset=CHAVE_JANELA_PAR)
            cols_ja_na_base = set(df.columns) - set(CHAVE_JANELA_PAR)
            cols_extra = [c for c in df_hfo.columns
                         if c not in CHAVE_JANELA_PAR and c not in ("rotulo", "janela")
                         and c not in cols_ja_na_base]
            df = df.merge(df_hfo[CHAVE_JANELA_PAR + cols_extra], on=CHAVE_JANELA_PAR,
                          how="left", validate="many_to_one")

    return df, n_base


def agrega_resultados(pasta_saida_base, saida_mestre="dataset_mestre.csv"):
    """
    Varre a árvore RESULTADOS/<sessao>/chan<N>/refinados.csv, faz o merge
    completo (skewness + comodulograma + harmônico + harmônico_hfo) por
    canal, extrai condição experimental do nome da pasta, e concatena tudo
    numa tabela mestre final.
    """
    print(f"[AGREGADOR] Lendo {pasta_saida_base}...")

    if not os.path.exists(pasta_saida_base):
        print("Pasta base não encontrada.")
        return

    dfs = []
    total_base = 0
    total_final = 0

    for root, dirs, files in os.walk(pasta_saida_base):
        if "refinados.csv" not in files:
            continue

        df_canal, n_base = agrega_canal(root)
        if df_canal is None:
            continue

        # Checagem de sanidade: um merge left-join correto contra a base
        # 'refinados' NUNCA deveria mudar a contagem de linhas. Se mudou,
        # alguma chave de merge está incompleta (duplicação) ou há um bug.
        if len(df_canal) != n_base:
            print(f"  [ALERTA] {root}: {n_base} linhas em refinados.csv -> "
                  f"{len(df_canal)} linhas após merges. Contagem deveria ser "
                  f"IDÊNTICA (merges são left-join 1:1/many:1). Investigar "
                  f"antes de confiar nesse canal.")

        total_base += n_base
        total_final += len(df_canal)

        partes = root.replace("\\", "/").split("/")
        canal = "unknown"
        sessao = "unknown"
        for p in partes:
            if p.startswith("chan"):
                canal = p.replace("chan", "")
            if "MTESC" in p or "Rodada" in p:
                sessao = p

        condicao = "basal"
        if "LAC" in sessao.upper():
            condicao = "LAC"
        elif "NOCI" in sessao.upper():
            condicao = "NOCI"

        df_canal["sessao"] = sessao
        df_canal["condicao"] = condicao
        df_canal["canal"] = canal.replace("chan", "")

        dfs.append(df_canal)

    if not dfs:
        print("Nenhum dado agregado encontrado.")
        return

    df_mestre = pd.concat(dfs, ignore_index=True)
    df_mestre.to_csv(saida_mestre, index=False)

    print(f"\n[AGREGADOR] Tabela mestre gerada: {saida_mestre} "
          f"({len(df_mestre)} linhas agregadas).")
    if total_final != total_base:
        print(f"[AGREGADOR] [ALERTA GLOBAL] Total em refinados.csv somado: {total_base} | "
              f"total após merge: {total_final}. Divergência de {total_final - total_base} "
              f"linhas — NÃO confiar no dataset até resolver.")
    else:
        print(f"[AGREGADOR] Checagem de sanidade OK: {total_base} linhas em refinados.csv "
              f"== {total_final} linhas na tabela mestre (nenhuma linha criada/perdida no merge).")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Agregador de resultados PAC para Machine Learning.")
    parser.add_argument("--resultados", required=True, help="Pasta contendo as sessões processadas")
    parser.add_argument("--saida", default="dataset_mestre.csv", help="Caminho do CSV de saída")
    args = parser.parse_args()

    agrega_resultados(args.resultados, args.saida)
