"""
gerar_template_comportamento.py -- gera/atualiza o template de anotação a partir
do dataset mestre, filtrado pelos candidatos robustos (veredito_refino ==
"Candidato robusto" -- NÃO a cadeia completa de consolida_vencedores.py, que
já exige comportamento anotado e criaria dependência circular).

Cada janela única (sessao, condicao, arquivo, janela_ini_s, janela_fim_s) vira
uma linha, com os pares/canais que a fizeram entrar (pares_detectados,
n_canais_pac) -- útil pra saber o quão robusto é aquele instante sem abrir o
CSV mestre. Se o arquivo de saída já existir, PRESERVA as anotações já feitas
(comportamento/observacoes/video_tempo_ini/video_tempo_fim), casando pela
chave da janela -- rodar de novo depois de processar mais sessões (ex.: a
infusão, depois de já ter anotado o basal) NUNCA apaga trabalho feito.

    python pipeline/comportamento/gerar_template_comportamento.py
"""
import argparse
import os
import shutil
import datetime

import pandas as pd

CHAVE_JANELA = ["sessao", "arquivo", "janela_ini_s", "janela_fim_s"]
COLUNAS_PRESERVADAS = ["comportamento", "observacoes", "video_tempo_ini", "video_tempo_fim"]


def monta_template(df_mestre):
    robustos = df_mestre[df_mestre["veredito_refino"] == "Candidato robusto"].copy()

    agg = (robustos.groupby(CHAVE_JANELA + ["condicao"])
           .agg(pares_detectados=("par", lambda s: "+".join(sorted(set(s)))),
                n_canais_pac=("canal", "nunique"))
           .reset_index())
    agg = agg.sort_values(["sessao", "arquivo", "janela_ini_s"])

    for col in COLUNAS_PRESERVADAS:
        agg[col] = ""

    colunas = ["sessao", "condicao", "arquivo", "janela_ini_s", "janela_fim_s",
               "video_tempo_ini", "video_tempo_fim", "pares_detectados", "n_canais_pac",
               "comportamento", "observacoes"]
    return agg[colunas]


def preserva_anotacoes(template_novo, saida):
    """Casa pela chave da janela e traz de volta o que já foi anotado no arquivo
    existente em `saida` (se houver). Linhas do template antigo cuja janela
    não sobreviveu no dataset mestre atual são descartadas silenciosamente
    (candidato deixou de ser robusto) -- não há como reintroduzi-las sem
    reprocessar a sessão."""
    if not os.path.exists(saida):
        return template_novo, 0

    antigo = pd.read_csv(saida, encoding="utf-8-sig")
    if not set(CHAVE_JANELA).issubset(antigo.columns):
        return template_novo, 0

    antigo_chave = antigo.set_index(CHAVE_JANELA)
    n_preservadas = 0
    for col in COLUNAS_PRESERVADAS:
        if col not in antigo.columns:
            continue
        valores = template_novo.set_index(CHAVE_JANELA).index.map(
            lambda k: antigo_chave[col].get(k, "") if k in antigo_chave.index else "")
        template_novo[col] = list(valores)
        if col == "comportamento":
            # pd.notna, nao str(v).strip() -- NaN vira a STRING "nan" ao passar por
            # str(), que e' truthy, e inflaria essa contagem com todas as pendentes.
            n_preservadas = template_novo[col].apply(lambda v: pd.notna(v) and str(v).strip() != "").sum()

    return template_novo, n_preservadas


def main():
    modulo_dir = os.path.dirname(os.path.abspath(__file__))
    script_dir = os.path.abspath(os.path.join(modulo_dir, "..", ".."))

    candidatos_mestre = [
        os.path.join(script_dir, "resultados", "dataset_mestre_final.csv"),
        os.path.join(script_dir, "dataset_mestre_final.csv"),
    ]
    padrao_mestre = next((c for c in candidatos_mestre if os.path.exists(c)), candidatos_mestre[0])
    padrao_saida = os.path.join(modulo_dir, "template_comportamento.csv")

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv_mestre", default=padrao_mestre, help="dataset_mestre_final.csv de entrada")
    ap.add_argument("--saida", default=padrao_saida, help="template de anotação (entrada E saída)")
    args = ap.parse_args()

    if not os.path.exists(args.csv_mestre):
        print(f"Dataset não encontrado: {args.csv_mestre}")
        return
    df_mestre = pd.read_csv(args.csv_mestre, encoding="utf-8-sig", low_memory=False)
    if len(df_mestre) == 0:
        print("Dataset vazio!")
        return

    template = monta_template(df_mestre)

    if os.path.exists(args.saida):
        backup_dir = os.path.join(modulo_dir, "backups_comportamento")
        os.makedirs(backup_dir, exist_ok=True)
        ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        shutil.copy2(args.saida, os.path.join(backup_dir, f"template_pre_regeneracao_{ts}.csv"))

    template, n_preservadas = preserva_anotacoes(template, args.saida)

    template.to_csv(args.saida, index=False, encoding="utf-8-sig")
    n_pendentes = template["comportamento"].isna().sum() + (
        template["comportamento"].notna() & (template["comportamento"].astype(str).str.strip() == "")).sum()
    print(f"Template com {len(template)} janelas ({n_preservadas} anotações preservadas, "
          f"{n_pendentes} pendentes).")
    print(f"Salvo em: {args.saida}")


if __name__ == "__main__":
    main()
