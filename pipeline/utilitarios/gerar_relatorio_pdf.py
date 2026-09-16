"""Gera PDF consolidado dos vencedores PAC teta-gama (CSV-driven).

Regra de ouro: NENHUM dado de sessão entra no código. Título, data, ratos,
condições e contagens saem do próprio CSV; as figuras (coluna `png`, nome do
comodulograma gerado pelo pipeline) são localizadas pelo nome sob
RESULTADOS_MESTRADO (ou --pasta_figuras).

Formato do CSV: o dos vencedores consolidados atuais
(`candidatos_vencedores_OURO_PURIFICADO_v2.csv` / `candidatos_vencedores_consolidados.csv`),
uma linha por canal x par. Colunas usadas: sessao, canal, par, arquivo,
janela_ini_s, janela_fim_s, z_score_refinado, fase_pico_hz, amp_pico_hz,
comportamento, png e (opcionais) rato, condicao_x/condicao, observacoes.

    python pipeline/utilitarios/gerar_relatorio_pdf.py
    python pipeline/utilitarios/gerar_relatorio_pdf.py --csv resultados/X.csv --saida relatorio.pdf
"""
import argparse
import datetime
import os
import sys

import pandas as pd
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import (Image, PageBreak, Paragraph, SimpleDocTemplate,
                                Spacer, Table, TableStyle)

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from pac_core.workspace import BASE_RESULTADOS, BASE_RESULTADOS_MESTRADO

COLUNAS_OBRIGATORIAS = ["sessao", "canal", "par", "arquivo", "janela_ini_s", "janela_fim_s",
                        "z_score_refinado", "fase_pico_hz", "amp_pico_hz", "comportamento", "png"]


def indexa_figuras(nomes, pasta):
    """nome do PNG -> caminho completo (primeira ocorrência sob `pasta`)."""
    alvo, idx = set(nomes), {}
    for raiz, _, arquivos in os.walk(pasta):
        for f in arquivos:
            if f in alvo and f not in idx:
                idx[f] = os.path.join(raiz, f)
    return idx


def _txt(v):
    return "" if pd.isna(v) else str(v)


def _tabela(dados, larguras, cor_cab):
    t = Table(dados, colWidths=larguras, repeatRows=1)
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(cor_cab)),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.whitesmoke),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 8),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.HexColor("#ecf0f1"), colors.white]),
        ("GRID", (0, 0), (-1, -1), 0.4, colors.grey),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
    ]))
    return t


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", default=os.path.join(BASE_RESULTADOS, "candidatos_vencedores_OURO_PURIFICADO_v2.csv"),
                    help="CSV de vencedores (default: resultados/candidatos_vencedores_OURO_PURIFICADO_v2.csv, "
                         "o mesmo que o preditor usa)")
    ap.add_argument("--saida", default=None,
                    help="PDF de saída (default: relatorio_<nome do csv>.pdf ao lado do CSV)")
    ap.add_argument("--pasta_figuras", default=BASE_RESULTADOS_MESTRADO,
                    help="Raiz onde procurar os PNGs da coluna `png` (default: RESULTADOS_MESTRADO)")
    args = ap.parse_args()

    df = pd.read_csv(args.csv, encoding="utf-8-sig")
    faltando = [c for c in COLUNAS_OBRIGATORIAS if c not in df.columns]
    if faltando:
        sys.exit(f"CSV sem as colunas {faltando} -- esperado o schema dos vencedores consolidados.")
    df = df.sort_values("z_score_refinado", ascending=False).reset_index(drop=True)

    saida = args.saida or os.path.join(
        os.path.dirname(os.path.abspath(args.csv)),
        f"relatorio_{os.path.splitext(os.path.basename(args.csv))[0]}.pdf")
    figuras = indexa_figuras(df["png"].dropna(), args.pasta_figuras)

    col_cond = next((c for c in ("condicao", "condicao_x") if c in df.columns), None)
    ratos = sorted(df["rato"].dropna().unique()) if "rato" in df.columns else []
    condicoes = sorted(df[col_cond].dropna().unique()) if col_cond else []
    n_eventos = len(df.drop_duplicates(["sessao", "arquivo", "janela_ini_s", "janela_fim_s"]))

    styles = getSampleStyleSheet()
    title = ParagraphStyle("T", parent=styles["Title"], fontSize=18, spaceAfter=15)
    h2 = ParagraphStyle("H2", parent=styles["Heading2"], fontSize=13,
                        textColor=colors.HexColor("#2c3e50"), spaceAfter=8, spaceBefore=12)
    body = ParagraphStyle("B", parent=styles["Normal"], fontSize=10, spaceAfter=5)
    story = [
        Paragraph("Relatório — Vencedores PAC teta-gama", title),
        Paragraph(f"Fonte: <i>{os.path.basename(args.csv)}</i> | "
                  f"Gerado em {datetime.date.today():%d/%m/%Y}", body),
        Paragraph(f"{len(df)} linha(s) canal x par, {n_eventos} janela(s) única(s)"
                  + (f" | Ratos: {', '.join(ratos)}" if ratos else "")
                  + (f" | Condições: {', '.join(map(str, condicoes))}" if condicoes else ""), body),
        Spacer(1, 0.3 * cm),
        Paragraph("Resumo", h2),
        Paragraph("Cada linha passou pelos portões do pipeline: estatístico (Candidato robusto, "
                  "FDR + ajuste Gama), FOOOF, harmônico e conferência comportamental por vídeo.", body),
        Paragraph("Pseudoreplicação espacial: canais com o mesmo pico na mesma janela são UM evento "
                  "biológico -- a contagem de janelas únicas acima já colapsa isso; a tabela abaixo "
                  "mantém a granularidade canal x par.", body),
        Paragraph("Por comportamento (janelas únicas)", h2),
    ]
    por_comp = (df.drop_duplicates(["sessao", "arquivo", "janela_ini_s", "janela_fim_s"])
                ["comportamento"].value_counts())
    story.append(_tabela([["Comportamento", "Janelas"]] + [[k, str(v)] for k, v in por_comp.items()],
                         [8 * cm, 3 * cm], "#2c3e50"))

    story.append(Paragraph("Vencedores (ordenados por Z refinado)", h2))
    linhas = [["Rato", "Canal", "Par", "Janela (s)", "Z", "Pico (Hz)", "Comportamento"]]
    for _, r in df.iterrows():
        linhas.append([_txt(r.get("rato")), _txt(r["canal"]), r["par"],
                       f"{r['janela_ini_s']:.0f}-{r['janela_fim_s']:.0f}",
                       f"{r['z_score_refinado']:.2f}",
                       f"{r['fase_pico_hz']:.0f} x {r['amp_pico_hz']:.0f}", _txt(r["comportamento"])])
    story.append(_tabela(linhas, [1.8 * cm, 1.3 * cm, 2.3 * cm, 2.2 * cm, 1.4 * cm, 2.2 * cm, 5 * cm],
                         "#2c3e50"))

    for _, r in df.iterrows():
        story.append(PageBreak())
        story.append(Paragraph(f"{_txt(r.get('rato'))} — canal {r['canal']} — {r['par']}", h2))
        story.append(Paragraph(_txt(r["sessao"]), body))
        img = figuras.get(r["png"])
        if img:
            story.append(Image(img, width=16 * cm, height=9 * cm, kind="proportional"))
        else:
            story.append(Paragraph(f"[Figura não encontrada: {_txt(r['png'])}]", body))
        story.append(Spacer(1, 0.2 * cm))
        story.append(_tabela([
            ["Métrica", "Valor"],
            ["Arquivo", r["arquivo"]],
            ["Janela", f"{r['janela_ini_s']:.1f}-{r['janela_fim_s']:.1f} s"],
            ["Condição", _txt(r.get(col_cond)) if col_cond else ""],
            ["Z refinado", f"{r['z_score_refinado']:.2f}"],
            ["Pico (fase x amp)", f"{r['fase_pico_hz']:.1f} x {r['amp_pico_hz']:.1f} Hz"],
            ["Comportamento", _txt(r["comportamento"])],
            ["Observações", _txt(r.get("observacoes"))],
        ], [3.5 * cm, 11.5 * cm], "#3498db"))

    SimpleDocTemplate(saida, pagesize=A4, rightMargin=2 * cm, leftMargin=2 * cm,
                      topMargin=2 * cm, bottomMargin=2 * cm).build(story)
    print(f"PDF gerado: {saida} ({len(df)} vencedores, {len(figuras)} figuras encontradas)")


if __name__ == "__main__":
    main()
