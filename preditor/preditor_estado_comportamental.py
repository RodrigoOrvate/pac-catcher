"""
preditor_estado_comportamental.py
==========================================
Preditor de PAC teta-gama baseado em ESTADO (comportamento + potencia teta),
nao em features espectrais locais do LFP bruto.

Por que existe: validar_preditor.py testou features espectrais/PAC (energia,
MI/MVL no par oficial, footprint espacial) nos 10s antes de cada evento e nao
achou nenhum precursor -- LOOCV no nivel de chance (~0.43), e teste de
Mann-Whitney confirmou nenhuma das 8 features separa pre-evento de controle
real (todas p>0.1). A hipotese testada aqui e diferente: em vez de procurar
uma "rampa" escondida no proprio sinal, testa se o ESTADO COMPORTAMENTAL da
janela anterior (N-1) + a potencia teta dessa mesma janela anterior preveem
se a janela seguinte (N, 10s depois) vira um "Candidato robusto" (passou no
portao estatistico de refino) -- formulacao SEM VAZAMENTO: so usa informacao
que estaria disponivel no momento da decisao, para um disparo de TTL real.

Evidencia que motivou esta abordagem (ver conversa/commits que a precedem):
  - Comportamento da janela atual associa-se ao veredito (qui-quadrado
    p=1.5e-10) -- consistente com Vanderwolf (1969), Tort et al. (2009,
    PNAS), Colgin (2016): PAC teta-gama e dependente de estado, nao um
    evento espontaneo do LFP.
  - O efeito sobrevive restringir a TRANSICOES reais de comportamento
    (p=0.018) -- nao e so autocorrelacao de um estado que persiste.
  - Comportamento anterior adiciona informacao alem da propria potencia
    teta (teste de razao de verossimilhanca, p=1.3e-5) -- nao e redutivel
    a "estava se movendo, logo tinha mais teta".

Resultado (validacao cruzada 5-fold honesta, 2715 amostras reais de
dataset_mestre_COM_COMPORTAMENTO.csv): acuracia=0.584, AUC-ROC=0.597.
Modesto, mas e a PRIMEIRA abordagem desta pasta que bate chance de forma
real (compare com validar_preditor.py, LOOCV ~0.43-0.44, abaixo de chance).
Ainda LONGE de confiavel para disparar hardware (ver README_preditor.md).

Uso:
    python preditor_estado_comportamental.py
"""
import os
import sys
import glob
import numpy as np
import pandas as pd
from scipy.signal import welch

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(SCRIPT_DIR, '..'))
from pac_core.io import carrega_dados

BASE_LAC_NOCI = r"C:\acoplamento_theta-gamma\LAC_NOCI"
CSV_DATASET_MESTRE = os.path.join(SCRIPT_DIR, '..', 'resultados', 'dataset_mestre_COM_COMPORTAMENTO.csv')
CSV_SAIDA = os.path.join(SCRIPT_DIR, '..', 'resultados', '_preditor_estado_teta.csv')
GAP_MAX_S = 15  # so pares onde a janela anterior e' de fato a anterior no tempo (sem buraco)


def resolve_pasta_basal(sessao_str, arquivo):
    """Mesma logica usada em validar_preditor.py / gera_galeria_top5.py."""
    nome_pasta = sessao_str
    sufixo = "_Basal antes da infusao"
    if nome_pasta.endswith(sufixo):
        nome_pasta = nome_pasta[: -len(sufixo)]
    candidatos = glob.glob(os.path.join(BASE_LAC_NOCI, "*", nome_pasta, "Basal antes da infusao"))
    candidatos += [c for c in glob.glob(os.path.join(BASE_LAC_NOCI, "*", nome_pasta)) if c not in candidatos]
    for c in candidatos:
        if os.path.isfile(os.path.join(c, arquivo)):
            return c
    return None


def theta_power(dados, fs, canal_idx, i0, i1):
    if i0 < 0 or i1 > dados.shape[0] or i1 <= i0:
        return np.nan
    sig = dados[i0:i1, canal_idx].astype(float)
    f, psd = welch(sig, fs, nperseg=min(1000, len(sig)))
    tm = (f >= 4) & (f <= 8)
    return float(np.mean(psd[tm])) if np.any(tm) else np.nan


def monta_dataset():
    """Para cada janela N com uma janela N-1 imediatamente anterior anotada,
    calcula a potencia teta de N-1 e junta com o comportamento de N-1.
    Alvo: N virou 'Candidato robusto'?"""
    df = pd.read_csv(CSV_DATASET_MESTRE)
    df["vencedor"] = (df["veredito_refino"] == "Candidato robusto")

    chave = ["sessao", "arquivo", "canal", "par"]
    df = df.sort_values(chave + ["janela_ini_s"])
    df["dur_janela"] = df["janela_fim_s"] - df["janela_ini_s"]
    df["comport_anterior"] = df.groupby(chave)["comportamento"].shift(1)
    df["janela_ini_anterior"] = df.groupby(chave)["janela_ini_s"].shift(1)
    df["gap_s"] = df["janela_ini_s"] - df["janela_ini_anterior"]

    consec = df[(df["gap_s"] <= GAP_MAX_S) & df["comport_anterior"].notna() &
                (df["comport_anterior"] != "Artefato / Cabo")].copy()

    _cache = {}
    thetas = []
    for _, row in consec.iterrows():
        pasta = resolve_pasta_basal(str(row["sessao"]), str(row["arquivo"]))
        if pasta is None:
            thetas.append(np.nan)
            continue
        ns2_path = os.path.join(pasta, row["arquivo"])
        if ns2_path not in _cache:
            _cache[ns2_path] = carrega_dados(ns2_path)
        dados, fs, canal_ids = _cache[ns2_path]
        canal_idx = int(row["canal"]) - 1
        if not (0 <= canal_idx < dados.shape[1]):
            thetas.append(np.nan)
            continue
        i0 = int(float(row["janela_ini_anterior"]) * fs)
        i1 = i0 + int(float(row["dur_janela"]) * fs)
        thetas.append(theta_power(dados, fs, canal_idx, i0, i1))

    consec["theta_power_anterior"] = thetas
    consec = consec.dropna(subset=["theta_power_anterior"])
    consec["log_theta_power_anterior"] = np.log(consec["theta_power_anterior"] + 1e-6)
    return consec


def main():
    consec = monta_dataset()
    print(f"Amostras com feature valida: {len(consec)}")
    consec.to_csv(CSV_SAIDA, index=False, encoding="utf-8-sig")

    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import StratifiedKFold, cross_val_predict
    from sklearn.metrics import confusion_matrix, precision_recall_fscore_support, roc_auc_score
    import joblib

    dummies = pd.get_dummies(consec["comport_anterior"], prefix="comport")
    X = pd.concat([consec[["log_theta_power_anterior"]].reset_index(drop=True),
                   dummies.reset_index(drop=True)], axis=1)
    y = consec["vencedor"].astype(int).values
    feat_names = X.columns.tolist()
    X = X.values

    print(f"Vencedor: {y.sum()}  Nao-vencedor: {len(y) - y.sum()}")
    print(f"Features ({len(feat_names)}): {feat_names}")

    clf = RandomForestClassifier(n_estimators=300, max_depth=5, min_samples_leaf=10,
                                  random_state=42, class_weight="balanced")
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    y_pred = cross_val_predict(clf, X, y, cv=cv)
    y_proba = cross_val_predict(clf, X, y, cv=cv, method="predict_proba")[:, 1]

    print("\n=== Validacao cruzada 5-fold (held-out) ===")
    print(f"Acuracia: {np.mean(y_pred == y):.3f}")
    cm = confusion_matrix(y, y_pred)
    print("Matriz de confusao (linhas=real, colunas=previsto):")
    print("              prev_Nao   prev_Vencedor")
    print(f"Nao-vencedor  {cm[0, 0]:<10} {cm[0, 1]:<10}")
    print(f"Vencedor      {cm[1, 0]:<10} {cm[1, 1]:<10}")
    prec, rec, f1, _ = precision_recall_fscore_support(y, y_pred, labels=[0, 1], zero_division=0)
    for label, p, r, f in zip(["Nao-vencedor", "Vencedor"], prec, rec, f1):
        print(f"  {label}: precisao={p:.3f} recall={r:.3f} F1={f:.3f}")
    print(f"AUC-ROC: {roc_auc_score(y, y_proba):.3f}  (0.5 = chance, 1.0 = perfeito)")

    clf.fit(X, y)
    out = os.path.join(SCRIPT_DIR, 'modelo_estado_comportamental.pkl')
    joblib.dump({"modelo": clf, "features": feat_names}, out)
    print(f"\nModelo salvo em: {out}")

    print("\nImportancia das features (modelo final, todos os dados):")
    for feat, imp in sorted(zip(feat_names, clf.feature_importances_), key=lambda t: -t[1]):
        print(f"  {feat}: {imp:.4f}")


if __name__ == "__main__":
    main()
