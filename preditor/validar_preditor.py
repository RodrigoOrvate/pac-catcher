"""
validar_preditor.py
==========================================
Validação honesta do preditor de PAC θ-γ.

Diferente do treino (que usa surrogates phase-scrambled como negativos, o que é
otimista), aqui extraímos janelas de CONTROLE REAIS dos mesmos arquivos .ns2
(10 s em instantes longe do evento conhecido) e testamos se o modelo separa
"pré-PAC" de "atividade normal" de verdade.

Métricas reportadas fora do conjunto de treino (hold-out via LOOCV):
  - Acurácia, precisão, recall, F1 por classe
  - Matriz de confusão
  - Threshold: para TTL real, o recall da classe pré-PAC importa mais que a
    acurácia (falso-negativo = perder um PAC que queríamos prevenir).

Uso:
    python validar_preditor.py
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
import joblib
from scipy.signal import welch

BASE = r"C:\acoplamento_theta-gamma"
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(SCRIPT_DIR, '..'))
from pac_core.io import carrega_dados

FS = 1000
PRE_WINDOW = 10
CONTROL_OFFSET_S = 60   # janela de controle começa event_time + offset
FEATURES = ['Energia', 'Theta_Energy', 'Gamma_Energy', 'Ratio_TG', 'Theta_Peak']


def extrair_features(window_data):
    """Janela (n_amostras, n_canais) -> vetor de features (mesma ordem do treino)."""
    sig = np.mean(window_data, axis=1)
    energia = np.var(sig)
    f, psd = welch(sig, FS, nperseg=min(1000, len(sig)))

    tm = (f >= 4) & (f <= 8)
    te = np.mean(psd[tm]) if np.any(tm) else 0.0
    tp = f[tm][np.argmax(psd[tm])] if np.any(tm) else 0.0

    gm = (f >= 30) & (f <= 80)
    ge = np.mean(psd[gm]) if np.any(gm) else 0.0

    return np.array([energia, te, ge, te / (ge + 1e-6), tp])


def localizar_ns2(sessao_dir, filename):
    for sub in ("Basal antes da infusao", ""):
        p = os.path.join(sessao_dir, sub, filename) if sub else os.path.join(sessao_dir, filename)
        if os.path.exists(p):
            return p
    return None


def resolver_canal(canal_name, canal_ids):
    if canal_name in canal_ids:
        return canal_ids.index(canal_name)
    try:
        num = str(canal_name).replace('chan', '')
        for i, cid in enumerate(canal_ids):
            if num in str(cid):
                return i
    except ValueError:
        pass
    return 0


def carregar_positivos_e_controles():
    """Roda sobre todos os vencedores.csv: pré-evento = positivo, e uma janela
    'controle' real no MESMO arquivo (event_time + 60 s) = negativo."""
    X, y, meta = [], [], []
    venc_files = glob.glob(os.path.join(BASE, "**", "RESULTADOS", "vencedores.csv"), recursive=True)

    for v_csv in venc_files:
        sessao_dir = os.path.dirname(os.path.dirname(v_csv))
        try:
            df = pd.read_csv(v_csv)
        except Exception:
            continue
        if 'inicio_s' not in df.columns or 'arquivo' not in df.columns:
            continue

        for _, row in df.iterrows():
            canal_name = row['canal']
            t_evento = float(row['inicio_s'])
            ns2_path = localizar_ns2(sessao_dir, row['arquivo'])
            if not ns2_path:
                print(f"  [skip] sem arquivo p/ {os.path.basename(v_csv)} chan{canal_name}")
                continue

            dados, fs, canal_ids = carrega_dados(ns2_path)
            ch = resolver_canal(canal_name, canal_ids)

            # POSITIVO: 10 s antes do evento
            i0 = int((t_evento - PRE_WINDOW) * fs)
            i1 = i0 + int(PRE_WINDOW * fs)
            if i0 >= 0 and i1 <= dados.shape[0]:
                X.append(extrair_features(dados[i0:i1]))
                y.append(1)
                meta.append((os.path.basename(ns2_path), canal_name, 'pre'))

            # CONTROLE (negativo real): 10 s em outro instante do mesmo arquivo
            c0 = int((t_evento + CONTROL_OFFSET_S) * fs)
            c1 = c0 + int(PRE_WINDOW * fs)
            if c1 <= dados.shape[0]:
                X.append(extrair_features(dados[c0:c1]))
                y.append(0)
                meta.append((os.path.basename(ns2_path), canal_name, 'controle'))
            else:
                # se o offset estoura o arquivo, pega mais cedo
                c0 = max(0, int((t_evento - CONTROL_OFFSET_S - PRE_WINDOW) * fs))
                c1 = c0 + int(PRE_WINDOW * fs)
                if c1 <= dados.shape[0]:
                    X.append(extrair_features(dados[c0:c1]))
                    y.append(0)
                    meta.append((os.path.basename(ns2_path), canal_name, 'controle'))

    return np.array(X), np.array(y), meta


def main():
    X, y, meta = carregar_positivos_e_controles()
    print(f"Amostras reais: {len(X)} (Pré-PAC: {sum(y)}, Controle real: {len(y)-sum(y)})")

    if len(X) < 10:
        print("Pouquíssimas amostras. Não dá para validar com confiança.")
        return

    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import LeaveOneOut
    from sklearn.metrics import confusion_matrix, precision_recall_fscore_support

    clf = RandomForestClassifier(n_estimators=100, random_state=42)

    # LOOCV com eventos reais fora do treino
    loo = LeaveOneOut()
    y_true, y_pred = [], []
    for tr_idx, te_idx in loo.split(X):
        clf.fit(X[tr_idx], y[tr_idx])
        y_true.append(y[te_idx[0]])
        y_pred.append(clf.predict(X[te_idx].reshape(1, -1))[0])
    y_true, y_pred = np.array(y_true), np.array(y_pred)

    print("\n=== Validação LOOCV (controles REAIS, não surrogates) ===")
    acc = np.mean(y_true == y_pred)
    print(f"Acurácia: {acc:.3f}")
    cm = confusion_matrix(y_true, y_pred)
    print("Matriz de confusão (linhas=real, colunas=previsto):")
    print("            prev_Controle  prev_PrePAC")
    print(f"Controle       {cm[0,0]:<9} {cm[0,1]:<7}")
    print(f"Pré-PAC        {cm[1,0]:<9} {cm[1,1]:<7}")

    prec, rec, f1, _ = precision_recall_fscore_support(y_true, y_pred,
                                                       labels=[0, 1],
                                                       zero_division=0)
    for label, p, r, f in zip(['Controle (0)', 'Pré-PAC (1)'], prec, rec, f1):
        print(f"  {label}: precisão={p:.3f} recall={r:.3f} F1={f:.3f}")

    # A classe que importa para o TTL: recall de pré-PAC (não perder eventos)
    print("\nPara TTL, o recall de pré-PAC (evitar falso-negativo) é o que mais importa.")
    print(f"  Recall pré-PAC = {rec[1]:.3f} (fração dos verdadeiros PAC que o modelo pegaria)")

    # Quantos falsos-positivos? (disparo indevido no controle)
    fp_rate = cm[0, 1] / max(1, cm[0, 0] + cm[0, 1])
    print(f"  Taxa de falso-positivo (disparar em atividade normal) = {fp_rate:.3f}")

    # Modelo final treinado em tudo
    clf.fit(X, y)
    out = os.path.join(SCRIPT_DIR, 'modelo_pac.pkl')
    joblib.dump(clf, out)
    print(f"\nModelo final treinado nos dados reais salvo em: {out}")


if __name__ == "__main__":
    main()
