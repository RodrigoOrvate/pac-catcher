"""
analisar_janelas_ultracurtas.py
==================================================
Investigação sistemática de precursores ultracurtos (1s, 2s, 3s)
e detecção de onset (1s, 2s) no LFP bruto para os 190 vencedores canônicos.

Objetivo:
1. Testar se o PAC teta-gama possui precursores eletrofisiológicos
   detectáveis em janelas ultracurtas (1 a 3 segundos antes do evento).
2. Testar se a dinâmica temporal imediata (inclinação/aceleração de teta)
   diferencia pré-evento de controle.
3. Testar se os primeiros 1-2 segundos do evento (onset) já são
   suficientes para detecção reativa de circuito fechado (closed-loop).
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
from scipy.signal import welch
from scipy.stats import wilcoxon, mannwhitneyu
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, accuracy_score, precision_recall_fscore_support

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(SCRIPT_DIR, '..'))
from pac_core.io import carrega_dados
from pac_core.filtering import aplica_notch

BASE_LAC_NOCI = r"C:\acoplamento_theta-gamma\LAC_NOCI"
CSV_VENCEDORES = os.path.join(SCRIPT_DIR, '..', 'resultados', 'candidatos_vencedores_OURO_PURIFICADO_v2.csv')
NOTCH_HZ = [60.0, 120.0, 180.0, 240.0]
CONTROL_OFFSET_S = 60.0

BANDAS = {
    "theta_rapida": (7.0, 10.0),   # Vanderwolf tipo 1 (locomoção)
    "theta_lenta": (4.0, 7.0),     # Vanderwolf tipo 2 (sniffing / atenção)
    "gamma_lenta": (30.0, 55.0),
    "gamma_rapida": (60.0, 90.0),
}


def resolve_pasta_basal(sessao_str, arquivo):
    nome_pasta = sessao_str
    sufixo = "_Basal antes da infusao"
    if nome_pasta.endswith(sufixo):
        nome_pasta = nome_pasta[:-len(sufixo)]
    candidatos = glob.glob(os.path.join(BASE_LAC_NOCI, "*", nome_pasta, "Basal antes da infusao"))
    candidatos += [c for c in glob.glob(os.path.join(BASE_LAC_NOCI, "*", nome_pasta)) if c not in candidatos]
    for c in candidatos:
        if os.path.isfile(os.path.join(c, arquivo)):
            return c
    return None


def extrair_potencias(sig, fs):
    """Extrai potências espectrais para janelas curtas."""
    if len(sig) < 50:
        return {b: np.nan for b in BANDAS}
    nperseg = min(500, len(sig))
    f, psd = welch(sig, fs, nperseg=nperseg)
    res = {}
    for nome, (f0, f1) in BANDAS.items():
        m = (f >= f0) & (f <= f1)
        res[nome] = float(np.mean(psd[m])) if np.any(m) else 0.0
    res["razao_tg"] = (res["theta_rapida"] + res["theta_lenta"]) / (res["gamma_lenta"] + res["gamma_rapida"] + 1e-9)
    res["energia_total"] = float(np.var(sig))
    return res


def processar_eventos():
    df = pd.read_csv(CSV_VENCEDORES)
    print(f"Carregando {len(df)} eventos de {CSV_VENCEDORES}...")
    _cache = {}

    registros = []
    for idx, row in df.iterrows():
        pasta = resolve_pasta_basal(str(row["sessao"]), str(row["arquivo"]))
        if pasta is None:
            continue
        ns2_path = os.path.join(pasta, row["arquivo"])
        if ns2_path not in _cache:
            _cache[ns2_path] = carrega_dados(ns2_path)
        dados, fs, canal_ids = _cache[ns2_path]
        canal_idx = int(row["canal"]) - 1
        if not (0 <= canal_idx < dados.shape[1]):
            continue

        t_ev = float(row["janela_ini_s"])
        # Controle pareado
        t_ctrl = t_ev + CONTROL_OFFSET_S
        if (t_ctrl + 5.0) * fs > dados.shape[0]:
            t_ctrl = max(5.0, t_ev - CONTROL_OFFSET_S)

        # Checa limites temporais para janela de 3s pré e 2s pós
        if (t_ev - 3.0) < 0 or (t_ev + 2.0) * fs > dados.shape[0]:
            continue
        if (t_ctrl - 3.0) < 0 or (t_ctrl + 2.0) * fs > dados.shape[0]:
            continue

        # Extrai os LFPs e aplica notch
        sig_ev_full = aplica_notch(dados[int((t_ev - 3.0) * fs):int((t_ev + 2.0) * fs), canal_idx].astype(float), fs, NOTCH_HZ)
        sig_ct_full = aplica_notch(dados[int((t_ctrl - 3.0) * fs):int((t_ctrl + 2.0) * fs), canal_idx].astype(float), fs, NOTCH_HZ)

        # Segmentos pré-evento relativos ao corte de 3s (índice 3*fs é t_ev)
        ponto_zero = int(3.0 * fs)
        pre_3s_ev = sig_ev_full[0:ponto_zero]
        pre_2s_ev = sig_ev_full[int(1.0 * fs):ponto_zero]
        pre_1s_ev = sig_ev_full[int(2.0 * fs):ponto_zero]

        # 3 fatias de 1s para tendência nos 3s antes: [-3, -2], [-2, -1], [-1, 0]
        s1_ev = sig_ev_full[0:int(1.0 * fs)]
        s2_ev = sig_ev_full[int(1.0 * fs):int(2.0 * fs)]
        s3_ev = sig_ev_full[int(2.0 * fs):ponto_zero]

        # Onset (início do PAC): [0, +1s] e [0, +2s]
        onset_1s_ev = sig_ev_full[ponto_zero:ponto_zero + int(1.0 * fs)]
        onset_2s_ev = sig_ev_full[ponto_zero:ponto_zero + int(2.0 * fs)]

        # Mesmos segmentos para o controle
        pre_3s_ct = sig_ct_full[0:ponto_zero]
        pre_2s_ct = sig_ct_full[int(1.0 * fs):ponto_zero]
        pre_1s_ct = sig_ct_full[int(2.0 * fs):ponto_zero]

        s1_ct = sig_ct_full[0:int(1.0 * fs)]
        s2_ct = sig_ct_full[int(1.0 * fs):int(2.0 * fs)]
        s3_ct = sig_ct_full[int(2.0 * fs):ponto_zero]

        onset_1s_ct = sig_ct_full[ponto_zero:ponto_zero + int(1.0 * fs)]
        onset_2s_ct = sig_ct_full[ponto_zero:ponto_zero + int(2.0 * fs)]

        # Features
        feat_pre3_ev = extrair_potencias(pre_3s_ev, fs)
        feat_pre2_ev = extrair_potencias(pre_2s_ev, fs)
        feat_pre1_ev = extrair_potencias(pre_1s_ev, fs)

        feat_pre3_ct = extrair_potencias(pre_3s_ct, fs)
        feat_pre2_ct = extrair_potencias(pre_2s_ct, fs)
        feat_pre1_ct = extrair_potencias(pre_1s_ct, fs)

        feat_on1_ev = extrair_potencias(onset_1s_ev, fs)
        feat_on2_ev = extrair_potencias(onset_2s_ev, fs)
        feat_on1_ct = extrair_potencias(onset_1s_ct, fs)
        feat_on2_ct = extrair_potencias(onset_2s_ct, fs)

        # Tendência de teta rápida nos 3 segundos pré
        teta_traj_ev = [extrair_potencias(s1_ev, fs)["theta_rapida"],
                         extrair_potencias(s2_ev, fs)["theta_rapida"],
                         extrair_potencias(s3_ev, fs)["theta_rapida"]]
        slope_teta_ev = np.polyfit([0, 1, 2], teta_traj_ev, 1)[0]

        teta_traj_ct = [extrair_potencias(s1_ct, fs)["theta_rapida"],
                         extrair_potencias(s2_ct, fs)["theta_rapida"],
                         extrair_potencias(s3_ct, fs)["theta_rapida"]]
        slope_teta_ct = np.polyfit([0, 1, 2], teta_traj_ct, 1)[0]

        registros.append({
            "idx": idx,
            "arquivo": row["arquivo"],
            "canal": row["canal"],
            "comportamento": row.get("comportamento", "N/D"),
            # Pré-evento
            "ev_pre3": feat_pre3_ev,
            "ev_pre2": feat_pre2_ev,
            "ev_pre1": feat_pre1_ev,
            "ev_slope_teta": slope_teta_ev,
            "ev_on1": feat_on1_ev,
            "ev_on2": feat_on2_ev,
            # Controle
            "ct_pre3": feat_pre3_ct,
            "ct_pre2": feat_pre2_ct,
            "ct_pre1": feat_pre1_ct,
            "ct_slope_teta": slope_teta_ct,
            "ct_on1": feat_on1_ct,
            "ct_on2": feat_on2_ct,
        })

    print(f"Total de pares (evento, controle) extraídos com sucesso: {len(registros)} / {len(df)}")
    return registros


def rodar_testes_estatisticos(registros):
    print("\n" + "=" * 75)
    print("TESTE UNIVARIADO: Wilcoxon Pareado e Mann-Whitney (Pré-evento vs Controle)")
    print("=" * 75)

    janelas = [
        ("Pré 3s [-3s, 0s]", "ev_pre3", "ct_pre3"),
        ("Pré 2s [-2s, 0s]", "ev_pre2", "ct_pre2"),
        ("Pré 1s [-1s, 0s]", "ev_pre1", "ct_pre1"),
        ("Onset 1s [0s, +1s]", "ev_on1", "ct_on1"),
        ("Onset 2s [0s, +2s]", "ev_on2", "ct_on2"),
    ]

    metricas = ["theta_rapida", "theta_lenta", "gamma_lenta", "gamma_rapida", "razao_tg", "energia_total"]

    for nome_jan, chave_ev, chave_ct in janelas:
        print(f"\n--- Janela: {nome_jan} (N={len(registros)} pares) ---")
        print(f"{'Métrica':<16} | {'Med. Pré':<10} | {'Med. Ctrl':<10} | {'Diff (%)':<8} | {'p-Wilcoxon':<12} | {'p-MannWhitney':<13} | {'Signif?'}")
        print("-" * 88)
        for m in metricas:
            vals_ev = np.array([r[chave_ev][m] for r in registros])
            vals_ct = np.array([r[chave_ct][m] for r in registros])

            med_ev = np.median(vals_ev)
            med_ct = np.median(vals_ct)
            diff_pct = ((med_ev - med_ct) / (med_ct + 1e-12)) * 100

            # Wilcoxon pareado
            try:
                _, p_wilc = wilcoxon(vals_ev, vals_ct)
            except Exception:
                p_wilc = np.nan

            # Mann-Whitney
            try:
                _, p_mw = mannwhitneyu(vals_ev, vals_ct, alternative='two-sided')
            except Exception:
                p_mw = np.nan

            sig = "***" if p_wilc < 0.001 else ("**" if p_wilc < 0.01 else ("*" if p_wilc < 0.05 else "ns"))
            print(f"{m:<16} | {med_ev:<10.4g} | {med_ct:<10.4g} | {diff_pct:+7.1f}% | {p_wilc:<12.4e} | {p_mw:<13.4e} | {sig}")

    # Teste para a tendência da inclinação
    print(f"\n--- Tendência de Inclinação da Potência Teta [-3s -> -1s] ---")
    slopes_ev = np.array([r["ev_slope_teta"] for r in registros])
    slopes_ct = np.array([r["ct_slope_teta"] for r in registros])
    _, p_w_slope = wilcoxon(slopes_ev, slopes_ct)
    _, p_mw_slope = mannwhitneyu(slopes_ev, slopes_ct)
    print(f"Mediana slope pré-evento: {np.median(slopes_ev):.4g}")
    print(f"Mediana slope controle:   {np.median(slopes_ct):.4g}")
    print(f"Wilcoxon p={p_w_slope:.4e}, Mann-Whitney p={p_mw_slope:.4e}")


def avaliar_classificadores(registros):
    print("\n" + "=" * 75)
    print("CLASSIFICAÇÃO MULTIVARIADA HONESTA (LOOCV / 5-Fold Stratified)")
    print("=" * 75)

    metricas = ["theta_rapida", "theta_lenta", "gamma_lenta", "gamma_rapida", "razao_tg", "energia_total"]

    cenarios = [
        ("1. Pré 3s puro", "ev_pre3", "ct_pre3", False),
        ("2. Pré 2s puro", "ev_pre2", "ct_pre2", False),
        ("3. Pré 1s puro", "ev_pre1", "ct_pre1", False),
        ("4. Pré 1s + Tendência (Slope)", "ev_pre1", "ct_pre1", True),
        ("5. Onset 1s [0s, +1s] (Closed-loop On-line)", "ev_on1", "ct_on1", False),
        ("6. Onset 2s [0s, +2s] (Closed-loop On-line)", "ev_on2", "ct_on2", False),
    ]

    for titulo, chave_ev, chave_ct, inclui_slope in cenarios:
        X_list = []
        y_list = []
        for r in registros:
            # Positivo (Pré-PAC ou Onset)
            f_ev = [np.log(r[chave_ev][m] + 1e-9) for m in metricas]
            if inclui_slope:
                f_ev.append(r["ev_slope_teta"])
            X_list.append(f_ev)
            y_list.append(1)

            # Negativo (Controle Real)
            f_ct = [np.log(r[chave_ct][m] + 1e-9) for m in metricas]
            if inclui_slope:
                f_ct.append(r["ct_slope_teta"])
            X_list.append(f_ct)
            y_list.append(0)

        X = np.array(X_list)
        y = np.array(y_list)

        # Random Forest com Stratified 5-Fold
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        clf = RandomForestClassifier(n_estimators=200, max_depth=4, random_state=42, class_weight="balanced")

        y_proba = cross_val_predict(clf, X, y, cv=cv, method="predict_proba")[:, 1]
        y_pred = (y_proba >= 0.5).astype(int)

        acc = accuracy_score(y, y_pred)
        auc = roc_auc_score(y, y_proba)
        prec, rec, f1, _ = precision_recall_fscore_support(y, y_pred, average='binary', zero_division=0)

        print(f"\n[{titulo}] (N_amostras={len(y)}, N_features={X.shape[1]}):")
        print(f"  Acurácia: {acc:.3f} | AUC-ROC: {auc:.3f} | Precisão (1): {prec:.3f} | Recall (1): {rec:.3f} | F1: {f1:.3f}")

        # Se for Onset, checa calibragem de limiar
        if "Onset" in titulo:
            from sklearn.metrics import precision_recall_curve
            precs, recs, thrs = precision_recall_curve(y, y_proba)
            f1s = 2 * precs[:-1] * recs[:-1] / (precs[:-1] + recs[:-1] + 1e-9)
            best_idx = np.argmax(f1s)
            print(f"  --> Melhor limiar por F1 ({thrs[best_idx]:.3f}): Precisão={precs[best_idx]:.3f}, Recall={recs[best_idx]:.3f}, F1={f1s[best_idx]:.3f}")


def main():
    registros = processar_eventos()
    if not registros:
        print("Nenhum registro extraído. Verifique os caminhos dos arquivos .ns2.")
        return
    rodar_testes_estatisticos(registros)
    avaliar_classificadores(registros)


if __name__ == "__main__":
    main()
