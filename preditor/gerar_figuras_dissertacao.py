"""
gerar_figuras_dissertacao.py
==================================================
Figuras da dissertação geradas A PARTIR DOS DADOS -- nenhum número digitado à
mão, nenhuma curva simulada. (Até 2026-09-14 as figuras 1 e 3 eram curvas
sintéticas de rng.normal, que nem batiam com as AUCs das próprias legendas, e
as figuras 2 e 4 usavam valores copiados à mão.)

1. figura1_comparacao_curvas_roc.png -- ROC fora-da-amostra de três preditores:
   LFP bruto 10 s antes do evento (validar_preditor.py, LOOCV), LFP ultracurto
   (analisar_janelas_ultracurtas.py, 5-fold) e estado comportamental + teta
   (preditor_estado_comportamental.py, 5-fold).
2. figura2_taxas_acoplamento_por_comportamento.png -- % de janelas N com PAC
   robusto por comportamento em N-1, com o qui-quadrado (Hipótese H2).
3. figura3_calibracao_limiar_precisao_recall.png -- curva precisão-recall do
   modelo de estado e o ponto de operação no limiar ótimo por F1.
4. figura4_trajetoria_potencia_teta_pre_evento.png -- potência de teta rápida
   nos 3 s antes do evento vs. controle pareado (Wilcoxon nas inclinações).

Fontes de dados (resultados/):
- figuras 1 (curva de estado), 2 e 3: dataset_mestre_COM_COMPORTAMENTO.csv e
  _preditor_estado_teta.csv (features salvas pelo preditor -- não relê .ns2).
- figura 1 (curvas de LFP): _oof_lfp_10s.csv (rode preditor/validar_preditor.py)
  e _oof_lfp_ultracurto.csv (rode preditor/analisar_janelas_ultracurtas.py).
- figura 4: _trajetoria_teta_pre_evento.csv (analisar_janelas_ultracurtas.py).
Figura cuja fonte não existe é PULADA, com aviso -- nunca preenchida com
dado inventado.

Uso:
    python preditor/gerar_figuras_dissertacao.py [--saida_dir docs/figuras]
        [--cenario_ultracurto "1. Pré 3s puro"] [--permitir_parcial]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency, wilcoxon
from sklearn.metrics import precision_recall_curve, roc_auc_score, roc_curve

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pac_core.workspace import BASE_RESULTADOS, BASE_SCRIPT
from pac_studio import transicao_estado
from preditor.preditor_estado_comportamental import avalia, features_modelo_final

CSV_MESTRE = os.path.join(BASE_RESULTADOS, "dataset_mestre_COM_COMPORTAMENTO.csv")
CSV_FEATURES_ESTADO = os.path.join(BASE_RESULTADOS, "_preditor_estado_teta.csv")
CSV_OOF_10S = os.path.join(BASE_RESULTADOS, "_oof_lfp_10s.csv")
CSV_OOF_ULTRA = os.path.join(BASE_RESULTADOS, "_oof_lfp_ultracurto.csv")
CSV_TRAJETORIA = os.path.join(BASE_RESULTADOS, "_trajetoria_teta_pre_evento.csv")

plt.rcParams["font.sans-serif"] = "Arial"
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["axes.edgecolor"] = "#333333"
plt.rcParams["axes.linewidth"] = 1.0


def _br(valor, fmt):
    """Número com vírgula decimal (padrão das publicações em português)."""
    return format(valor, fmt).replace(".", ",")


def _falta(*caminhos):
    return [c for c in caminhos if not os.path.exists(c)]


def _salva(fig, saida_dir, nome):
    caminho = os.path.join(saida_dir, nome)
    fig.tight_layout()
    fig.savefig(caminho)
    plt.close(fig)
    print(f"Salvo: {caminho}")


def oof_modelo_estado():
    """Probabilidades fora-da-amostra do modelo final de estado, refeitas das
    features salvas (mesma configuração e mesmas sementes do preditor --
    reproduz o limiar gravado em modelo_estado_comportamental.pkl bit a bit)."""
    consec = pd.read_csv(CSV_FEATURES_ESTADO, encoding="utf-8-sig")
    y = consec["vencedor"].astype(int).values
    _, proba, _ = avalia(features_modelo_final(consec).values, y, "modelo de estado (fora-da-amostra)")
    return y, proba


def figura1_roc(saida_dir, estado, cenario_ultra, permitir_parcial):
    faltando = _falta(CSV_OOF_10S, CSV_OOF_ULTRA)
    if faltando and not permitir_parcial:
        print(f"[pulada] figura 1: falta {', '.join(faltando)} -- rode preditor/validar_preditor.py e "
              "preditor/analisar_janelas_ultracurtas.py (ou use --permitir_parcial).")
        return

    curvas = []
    if os.path.exists(CSV_OOF_10S):
        d = pd.read_csv(CSV_OOF_10S, encoding="utf-8-sig")
        curvas.append(("LFP bruto 10 s antes (LOOCV)", d["y_true"], d["y_proba"], "#d9534f", "-"))
    if os.path.exists(CSV_OOF_ULTRA):
        d = pd.read_csv(CSV_OOF_ULTRA, encoding="utf-8-sig")
        print("AUC por cenário ultracurto:",
              {c: round(roc_auc_score(g["y_true"], g["y_proba"]), 3) for c, g in d.groupby("cenario")})
        g = d[d["cenario"] == cenario_ultra]
        if g.empty:
            print(f"[aviso] cenário {cenario_ultra!r} não está em {CSV_OOF_ULTRA}; "
                  f"opções: {sorted(d['cenario'].unique())}")
        else:
            curvas.append((f"LFP ultracurto: {cenario_ultra.split('. ', 1)[-1]}", g["y_true"], g["y_proba"],
                           "#f0ad4e", "-."))
    y_est, p_est = estado
    curvas.append(("Estado comportamental + teta (5-fold)", y_est, p_est, "#2e6da4", "-"))

    fig, ax = plt.subplots(figsize=(6.5, 5.5), dpi=300)
    ax.plot([0, 1], [0, 1], linestyle="--", color="#888888", lw=1.5, label="Chance (AUC = 0,50)")
    for rotulo, y, p, cor, estilo in curvas:
        fpr, tpr, _ = roc_curve(y, p)
        auc = roc_auc_score(y, p)
        ax.plot(fpr, tpr, color=cor, lw=2.2, linestyle=estilo,
                label=f"{rotulo} (AUC = {_br(auc, '.3f')}, n = {len(y)})")
    ax.set_xlim([-0.02, 1.02])
    ax.set_ylim([-0.02, 1.02])
    ax.set_xlabel("Taxa de falso-positivo (1 − especificidade)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Taxa de verdadeiro-positivo (sensibilidade)", fontsize=11, fontweight="bold")
    ax.set_title("Desempenho preditivo fora-da-amostra (ROC)", fontsize=12, fontweight="bold", pad=12)
    ax.legend(loc="lower right", frameon=True, facecolor="#fbfbfb", edgecolor="#cccccc", fontsize=8.5)
    ax.grid(True, linestyle=":", alpha=0.6)
    _salva(fig, saida_dir, "figura1_comparacao_curvas_roc.png")


def figura2_comportamento(saida_dir):
    if _falta(CSV_MESTRE):
        print(f"[pulada] figura 2: falta {CSV_MESTRE}")
        return
    consec = transicao_estado(pd.read_csv(CSV_MESTRE, encoding="utf-8-sig"))
    tabela = pd.crosstab(consec["comport_anterior"], consec["vencedor"])
    chi2, p, dof, _ = chi2_contingency(tabela)
    taxas = (tabela[True] / tabela.sum(axis=1) * 100).sort_values()

    fig, ax = plt.subplots(figsize=(8.5, 5.0), dpi=300)
    cores = ["#c00000" if v == taxas.min() else "#2b5c8f" for v in taxas.values]
    barras = ax.barh(taxas.index, taxas.values, color=cores, edgecolor="#222222", height=0.65)
    for barra in barras:
        w = barra.get_width()
        ax.text(w + 1.0, barra.get_y() + barra.get_height() / 2, f"{_br(w, '.1f')}%",
                va="center", ha="left", fontsize=10, fontweight="bold", color="#222222")
    ax.set_xlim(0, max(60, taxas.max() + 10))
    ax.set_xlabel("Janelas N com acoplamento robusto (%)", fontsize=11, fontweight="bold")
    ax.set_title("Estado comportamental em N-1 prediz acoplamento em N\n"
                 f"(χ² = {_br(chi2, '.1f')}, gl = {dof}, p = {_br(p, '.1e')}, n = {len(consec)} pares)",
                 fontsize=12, fontweight="bold", pad=12)
    ax.grid(axis="x", linestyle=":", alpha=0.6)
    _salva(fig, saida_dir, "figura2_taxas_acoplamento_por_comportamento.png")


def figura3_precisao_recall(saida_dir, estado):
    y, proba = estado
    prec, rec, thr = precision_recall_curve(y, proba)
    f1 = 2 * prec[:-1] * rec[:-1] / (prec[:-1] + rec[:-1] + 1e-9)
    i = int(np.argmax(f1))
    prevalencia = y.mean()

    fig, ax = plt.subplots(figsize=(6.5, 5.5), dpi=300)
    ax.plot(rec, prec, color="#1f77b4", lw=2.5, label="Modelo de estado (fora-da-amostra)")
    ax.scatter([rec[i]], [prec[i]], color="#d9534f", s=120, zorder=5,
               label=(f"Limiar ótimo por F1 = {_br(thr[i], '.3f')}\n"
                      f"Recall = {_br(100 * rec[i], '.1f')}% | Precisão = {_br(100 * prec[i], '.1f')}%"))
    ax.annotate(f"Recall alto (evita perder o evento)\ncusto: {100 * (1 - prec[i]):.0f}% dos disparos "
                "seriam falso-positivos",
                xy=(rec[i], prec[i]), xytext=(max(0.05, rec[i] - 0.6), min(0.95, prec[i] + 0.15)),
                arrowprops=dict(facecolor="#d9534f", shrink=0.08, width=1.5, headwidth=7),
                fontsize=9, bbox=dict(boxstyle="round,pad=0.4", fc="#fff2f2", ec="#d9534f", lw=1))
    ax.axhline(prevalencia, color="#888888", linestyle=":",
               label=f"Prevalência (acaso) = {_br(100 * prevalencia, '.1f')}%")
    ax.set_xlim([0, 1.05])
    ax.set_ylim([0, 1.0])
    ax.set_xlabel("Recall / sensibilidade", fontsize=11, fontweight="bold")
    ax.set_ylabel("Precisão", fontsize=11, fontweight="bold")
    ax.set_title("Compromisso precisão-recall do modelo de estado", fontsize=12, fontweight="bold", pad=12)
    ax.legend(loc="lower left", frameon=True, facecolor="#fbfbfb", edgecolor="#cccccc", fontsize=9)
    ax.grid(True, linestyle=":", alpha=0.6)
    _salva(fig, saida_dir, "figura3_calibracao_limiar_precisao_recall.png")


def figura4_trajetoria_teta(saida_dir):
    if _falta(CSV_TRAJETORIA):
        print(f"[pulada] figura 4: falta {CSV_TRAJETORIA} -- rode preditor/analisar_janelas_ultracurtas.py.")
        return
    d = pd.read_csv(CSV_TRAJETORIA, encoding="utf-8-sig")
    tempos = [-3, -2, -1]
    med_ev = [d[f"ev_teta_{t}s"].median() for t in tempos]
    med_ct = [d[f"ct_teta_{t}s"].median() for t in tempos]
    _, p = wilcoxon(d["ev_slope_teta"], d["ct_slope_teta"])

    fig, ax = plt.subplots(figsize=(6.0, 5.0), dpi=300)
    ax.plot(tempos, med_ev, marker="o", lw=2.5, color="#2b5c8f",
            label=f"Antes do PAC (inclinação mediana = {_br(d['ev_slope_teta'].median(), '+.1f')})")
    ax.plot(tempos, med_ct, marker="s", lw=2.5, color="#888888", linestyle="--",
            label=f"Controle pareado (inclinação mediana = {_br(d['ct_slope_teta'].median(), '+.1f')})")
    ax.set_xticks(tempos)
    ax.set_xticklabels(["[-3, -2] s", "[-2, -1] s", "[-1, 0] s"])
    ax.set_xlabel("Segmento de 1 s antes do início da janela", fontsize=11, fontweight="bold")
    # le_ns2 lê o sinal sem conversão para µV -- unidades brutas do .ns2
    ax.set_ylabel("Potência teta tipo 1, 7–10 Hz (PSD mediana,\nunidades brutas do .ns2)",
                  fontsize=10, fontweight="bold")
    ax.set_title(f"Potência de teta nos 3 s antes do evento\n(Wilcoxon pareado nas inclinações: "
                 f"p = {_br(p, '.3f')}, n = {len(d)} pares)",
                 fontsize=11.5, fontweight="bold", pad=12)
    ax.legend(loc="best", frameon=True, facecolor="#fbfbfb", edgecolor="#cccccc", fontsize=9.5)
    ax.grid(True, linestyle=":", alpha=0.6)
    _salva(fig, saida_dir, "figura4_trajetoria_potencia_teta_pre_evento.png")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--saida_dir", default=os.path.join(BASE_SCRIPT, "docs", "figuras"))
    ap.add_argument("--cenario_ultracurto", default="1. Pré 3s puro",
                    help="cenário de analisar_janelas_ultracurtas.py usado na figura 1")
    ap.add_argument("--permitir_parcial", action="store_true",
                    help="gera a figura 1 mesmo faltando alguma das curvas de LFP")
    args = ap.parse_args()
    os.makedirs(args.saida_dir, exist_ok=True)

    estado = None
    if _falta(CSV_FEATURES_ESTADO):
        print(f"[pulada] figuras 1 e 3: falta {CSV_FEATURES_ESTADO} -- rode "
              "preditor/preditor_estado_comportamental.py.")
    else:
        estado = oof_modelo_estado()

    if estado is not None:
        figura1_roc(args.saida_dir, estado, args.cenario_ultracurto, args.permitir_parcial)
    figura2_comportamento(args.saida_dir)
    if estado is not None:
        figura3_precisao_recall(args.saida_dir, estado)
    figura4_trajetoria_teta(args.saida_dir)


if __name__ == "__main__":
    main()
