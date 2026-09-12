"""
consolida_vencedores.py
==========================================
Formaliza o filtro final de "candidatos vencedores" — antes deste script,
essa filtragem só existia como código ad-hoc rodado manualmente (pandas
em conversa). Aplica de uma vez os portões estatístico, FOOOF, harmônico
e comportamental sobre o dataset mestre já enriquecido e anotado
(`dataset_mestre_COM_COMPORTAMENTO.csv`, saída de `junta_comportamento.py`).

Critérios (mesmos já documentados no README, Passo 3.8):
  - Estatístico : veredito_refino == 'Candidato robusto'
                  (já engloba FDR + ajuste Gama + não-suspeito de banda larga)
  - FOOOF       : erro_teta/gamma_fooof_v2 < erro_fooof_max (default 0.15),
                  knee_valido_teta/gamma, cf_teta/gamma_fooof_v2 não-nulo
  - Harmônico   : veredito_harmonico == 'CLEAN'
  - Comportamento: anotado (não-nulo) e diferente de 'Artefato / Cabo'

A saída mantém granularidade canal×par (uma linha por candidato) — a
síntese impressa no console colapsa pseudoreplicação espacial (um evento
por janela, ver CLAUDE.md #3) só para o relatório.

Uso:
    python pipeline/dataset_mestre/consolida_vencedores.py
    # sem argumentos usa os defaults (resultados/dataset_mestre_COM_COMPORTAMENTO.csv
    # -> resultados/candidatos_vencedores_consolidados.csv)
"""
import argparse
import glob
import os
import re

import pandas as pd

CHAVE_JANELA = ["sessao", "arquivo", "janela_ini_s", "janela_fim_s"]


def resolve_rato(sessao, arquivo, base_lac_noci):
    """Extrai o rato (ex.: 'MTESC04') de uma linha do dataset mestre.

    Sessões NOCI trazem o nome do rato literalmente na string de sessao
    (ex.: "MTESC04 -- 2 - infusao..."). Sessões do grupo LAC usam nomes de
    pasta "Rodada-N-DD-MM-2024" que se REPETEM entre MTESC03_LAC e
    MTESC05_LAC (mesmo protocolo, ratos diferentes) -- para essas,
    desambigua verificando em qual pasta *_LAC o arquivo .ns2 realmente
    existe (o nome do arquivo, timestamp da gravação, é único por rato).
    """
    m = re.search(r"MTESC\d+", sessao)
    if m:
        return m.group(0)

    nome_pasta = sessao
    sufixo = "_Basal antes da infusao"
    if nome_pasta.endswith(sufixo):
        nome_pasta = nome_pasta[: -len(sufixo)]

    # Estrutura de pastas do grupo LAC é inconsistente entre datas: a
    # maioria tem uma subpasta "Basal antes da infusao", mas em algumas
    # (ex.: MTESC05_LAC/Rodada-1-04-05-2024, Rodada-2-06-05-2024,
    # Rodada-2-09-05-2024) os .ns2 ficam direto na pasta da sessão --
    # tenta as duas formas.
    candidatos = glob.glob(os.path.join(base_lac_noci, "*_LAC", nome_pasta,
                                        "Basal antes da infusao", str(arquivo)))
    candidatos += glob.glob(os.path.join(base_lac_noci, "*_LAC", nome_pasta, str(arquivo)))
    for c in candidatos:
        m2 = re.search(r"(MTESC\d+)_LAC", c.replace("/", os.sep))
        if m2:
            return m2.group(1)
    return None


def consolida_vencedores(df_mestre, base_lac_noci, erro_fooof_max=0.15):
    robusto = df_mestre["veredito_refino"] == "Candidato robusto"

    fooof_ok = (
        (df_mestre["erro_teta_fooof_v2"] < erro_fooof_max)
        & (df_mestre["erro_gamma_fooof_v2"] < erro_fooof_max)
        & (df_mestre["knee_valido_teta"] == True)
        & (df_mestre["knee_valido_gamma"] == True)
        & df_mestre["cf_teta_fooof_v2"].notna()
        & df_mestre["cf_gamma_fooof_v2"].notna()
    )

    harm_clean = df_mestre["veredito_harmonico"] == "CLEAN"

    comport_puro = df_mestre["comportamento"].notna() & (
        df_mestre["comportamento"] != "Artefato / Cabo"
    )

    vencedores = df_mestre[robusto & fooof_ok & harm_clean & comport_puro].copy()
    vencedores["rato"] = vencedores.apply(
        lambda r: resolve_rato(r["sessao"], r["arquivo"], base_lac_noci), axis=1
    )
    return vencedores


def imprime_sintese(vencedores):
    print(f"\nLinhas vencedoras (canal x par): {len(vencedores)}")
    eventos = (
        vencedores.sort_values("z_score_refinado", ascending=False)
        .drop_duplicates(CHAVE_JANELA)
    )
    print(f"Eventos únicos (colapsando pseudoreplicação espacial): {len(eventos)}")
    print("\n--- Eventos por rato ---")
    print(eventos["rato"].value_counts(dropna=False).to_string())
    print("\n--- Eventos por comportamento ---")
    print(eventos["comportamento"].value_counts().to_string())


def main():
    modulo_dir = os.path.dirname(os.path.abspath(__file__))
    script_dir = os.path.abspath(os.path.join(modulo_dir, "..", ".."))
    workspace_root = os.path.abspath(os.path.join(script_dir, ".."))

    padrao_entrada = os.path.join(script_dir, "resultados", "dataset_mestre_COM_COMPORTAMENTO.csv")
    padrao_saida = os.path.join(script_dir, "resultados", "candidatos_vencedores_consolidados.csv")
    padrao_lac_noci = os.path.join(workspace_root, "LAC_NOCI")

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--entrada", default=padrao_entrada,
                    help="dataset_mestre_COM_COMPORTAMENTO.csv de entrada")
    ap.add_argument("--saida", default=padrao_saida,
                    help="CSV de saída (candidatos_vencedores_consolidados.csv)")
    ap.add_argument("--base_lac_noci", default=padrao_lac_noci,
                    help="Raiz das pastas de dados brutos (p/ desambiguar sessões LAC)")
    ap.add_argument("--erro_fooof_max", type=float, default=0.15,
                    help="Limiar de erro_ajuste do FOOOF (teta e gama), default 0.15")
    args = ap.parse_args()

    if not os.path.exists(args.entrada):
        print(f"Erro: entrada não encontrada em {args.entrada}")
        return

    print(f"Lendo: {args.entrada}")
    df_mestre = pd.read_csv(args.entrada)

    vencedores = consolida_vencedores(df_mestre, args.base_lac_noci, args.erro_fooof_max)

    os.makedirs(os.path.dirname(args.saida), exist_ok=True)
    vencedores.to_csv(args.saida, index=False, encoding="utf-8-sig")
    print(f"Salvo: {args.saida}")

    imprime_sintese(vencedores)


if __name__ == "__main__":
    main()
