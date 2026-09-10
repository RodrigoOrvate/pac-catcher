"""Gera PDF final consolidado dos vencedores PAC Theta-Gamma (CSV-driven).

Regra de ouro observada: NENHUM dado de sessão entra no código. O gerador lê o
vencedores_consolidado.csv do estudo (criado pelo pipeline em
RESULTADOS_CONSOLIDADOS/) e resolve as figuras no mesmo diretório.

Formato do CSV (uma linha por vencedor validado):
    sessao,canal,comportamento,z,pico,classificacao,veredito,observacoes,figura

Saída: relatorio_final_<ESTUDO>_NOCI.pdf no diretório RESULTADOS_CONSOLIDADOS
do CSV (o nome do estudo é derivado do caminho).

    python gerar_relatorio_pdf.py \
        --csv "../MTESC05_NOCI/RESULTADOS_CONSOLIDADOS/vencedores_consolidado.csv"
"""
import argparse
import csv
import os
import re

from reportlab.lib.pagesizes import A4
from reportlab.lib.units import cm
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Image,
                                Table, TableStyle, PageBreak)
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors


def _nome_estudo(csv_path):
    # <estudo>/RESULTADOS_CONSOLIDADOS/vencedores_consolidado.csv
    consolidado_dir = os.path.dirname(os.path.abspath(csv_path))
    raiz = os.path.basename(os.path.dirname(consolidado_dir))  # ex.: MTESC05_NOCI
    sigla = re.sub(r"_NOCI$", "", raiz) if raiz.endswith("_NOCI") else raiz
    return sigla, raiz, consolidado_dir


def _sessao_label(s):
    if s.startswith("S") and s[1:].isdigit():
        return f"Sessão {int(s[1:])}"
    return s


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True,
                    help="vencedores_consolidado.csv do estudo")
    args = ap.parse_args()

    sigla, raiz, consolidado_dir = _nome_estudo(args.csv)
    out = os.path.join(consolidado_dir, f"relatorio_final_{raiz}.pdf")

    with open(args.csv, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    n_forte = sum(1 for r in rows if r["classificacao"].strip().lower() == "forte")
    n_mod = sum(1 for r in rows if r["classificacao"].strip().lower() != "forte")

    doc = SimpleDocTemplate(out, pagesize=A4, rightMargin=2*cm, leftMargin=2*cm,
                            topMargin=2*cm, bottomMargin=2*cm)
    styles = getSampleStyleSheet()
    title = ParagraphStyle("T", parent=styles["Title"], fontSize=18,
                           textColor=colors.HexColor("#1a1a1a"), spaceAfter=15,
                           alignment=1)
    h2 = ParagraphStyle("H2", parent=styles["Heading2"], fontSize=13,
                        textColor=colors.HexColor("#2c3e50"), spaceAfter=8,
                        spaceBefore=12)
    body = ParagraphStyle("B", parent=styles["Normal"], fontSize=10,
                          textColor=colors.HexColor("#333"), spaceAfter=5)
    story = []

    story.append(Paragraph("Relatório Final — Vencedores PAC ΘΓ", title))
    story.append(Paragraph(f"{sigla} NOCI — Basal", body))
    story.append(Paragraph(f"Data: 30/08/2026 | {len(rows)} evento(s) validado(s) "
                           f"({n_forte} forte(s), {n_mod} moderado(s))", body))
    story.append(Spacer(1, 0.4*cm))

    story.append(Paragraph("Resumo", h2))
    story.append(Paragraph(
        f"Consolidação dos acoplamentos Theta-Gamma validados no estudo {sigla} "
        f"NOCI (basal, sem interação), a partir de <i>vencedores_consolidado.csv</i>. "
        "Cada vencedor passou pelas camadas de validação do pipeline v2 "
        "(FDR de janela, anti-60 Hz, respiração, FDR do mapa concentrado em ΘΓ, "
        "robustez, conferência comportamental por vídeo).", body))
    story.append(Paragraph(
        "Pseudoreplicação espacial: canais com pico idêntico no mesmo arquivo/janela "
        "foram tratados como UM evento co-detectado (representante de maior z), "
        "não como casos independentes.", body))

    story.append(Paragraph("Classificação", h2))
    tbl_rows = [["Canal", "Sessão", "Z", "Classe", "Comportamento"]]
    for r in rows:
        tbl_rows.append([r["canal"], r["sessao"], r["z"], r["classificacao"],
                         r["comportamento"]])
    tbl = Table(tbl_rows, colWidths=[2.5*cm, 2.5*cm, 2*cm, 2.5*cm, 5.5*cm])
    tbl.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#2c3e50")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.whitesmoke),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 9),
        ("ALIGN", (0, 0), (-1, -1), "LEFT"),
        ("BOTTOMPADDING", (0, 0), (-1, 0), 6),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1),
         [colors.HexColor("#ecf0f1"), colors.white]),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
    ]))
    story.append(tbl)
    story.append(Spacer(1, 0.4*cm))

    for i, r in enumerate(rows):
        if i > 0:
            story.append(PageBreak())
        cor = (colors.HexColor("#27ae60")
               if r["classificacao"].strip().lower() == "forte"
               else colors.HexColor("#3498db"))
        story.append(Paragraph(f"{_sessao_label(r['sessao'])} — {r['canal']}", h2))
        story.append(Paragraph(
            f"<b>{r['classificacao']}</b> | Z={r['z']} | {r['pico']}",
            ParagraphStyle("x", parent=body, textColor=cor, fontSize=11)))
        img = os.path.join(consolidado_dir, r["figura"])
        if os.path.exists(img):
            story.append(Image(img, width=16*cm, height=9*cm))
        else:
            story.append(Paragraph(f"[Imagem não encontrada: {img}]", body))
        story.append(Spacer(1, 0.2*cm))
        d = [["Métrica", "Valor"],
             ["Canal", r["canal"]],
             ["Sessão", r["sessao"]],
             ["Comportamento", r["comportamento"]],
             ["Z-Score", r["z"]],
             ["Pico", r["pico"]],
             ["Classe", r["classificacao"]],
             ["Observações", r["observacoes"]]]
        t2 = Table(d, colWidths=[3.5*cm, 11.5*cm])
        t2.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#3498db")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.whitesmoke),
            ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
            ("FONTSIZE", (0, 0), (-1, 0), 10),
            ("FONTSIZE", (0, 1), (-1, -1), 9),
            ("BACKGROUND", (0, 1), (0, -1), colors.HexColor("#ecf0f1")),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ]))
        story.append(t2)

    doc.build(story)
    print(f"PDF gerado: {out}")


if __name__ == "__main__":
    main()
