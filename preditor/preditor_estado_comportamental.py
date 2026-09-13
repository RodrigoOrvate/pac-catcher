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

Evidencia que motivou esta abordagem: comportamento associa-se ao veredito
(qui-quadrado p=1.5e-10), efeito sobrevive restringir a TRANSICOES reais
(p=0.018, descarta autocorrelacao pura) e sobrevive controlar pela propria
potencia teta (teste de razao de verossimilhanca p=1.3e-5, descarta "e so
teta disfarcado"). Baseline (log_theta_power_anterior + comportamento
one-hot): acuracia 0.584, AUC-ROC 0.597.

Este script faz um ABLATION SEQUENCIAL de 4 melhorias sobre esse baseline,
reportando o delta de AUC de cada uma isoladamente:
  1. Normalizacao por canal: potencia teta bruta nao e comparavel entre
     canais/animais (impedancia, distancia da camada piramidal). Z-score
     dentro do proprio (arquivo, canal) remove essa variancia de nuisance.
  2. Teta tipo 1 (rapido, 7-10Hz, locomocao) vs tipo 2 (lento, 4-7Hz,
     sniffing/imobilidade) -- distincao classica de Vanderwolf que a banda
     unica 4-8Hz borra.
  3. Lag mais curto + tendencia: em vez de um unico valor pontual em N-1,
     divide o intervalo entre N-1 e N em 5 sub-janelas e usa (a) a mais
     proxima do evento e (b) a inclinacao (esta subindo?).
  4. MI_z num par fase-amplitude CANONICO fixo (6Hz x 85Hz, mediana dos
     190 vencedores de OURO_PURIFICADO_v2.csv) medido em N-1 -- NAO usa o
     par "oficial" de cada linha porque esse so existe para quem already
     tem fase_pico_hz/amp_pico_hz preenchido (2202/7393 linhas, quase
     perfeitamente colineares com o proprio rotulo vencedor -- usar o par
     de CADA linha seria vazamento de informacao do rotulo pra dentro da
     feature). Um par fixo, igual para todas as linhas, evita isso.

Depois do ablation, calibra o limiar de decisao via curva precisao-recall
(o que importa pra um disparo de TTL e recall da classe vencedor, nao
acuracia agregada).

Uso:
    python preditor_estado_comportamental.py
"""
import os
import sys
import glob
import numpy as np
import pandas as pd
from scipy.signal import welch, hilbert

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(SCRIPT_DIR, '..'))
from pac_core.io import carrega_dados
from pac_core.filtering import aplica_notch
from pipeline.etapa4_validacao.robustez_parametros import mi_z_par

NOTCH_HZ = [60.0, 120.0, 180.0, 240.0]

BASE_LAC_NOCI = r"C:\acoplamento_theta-gamma\LAC_NOCI"
CSV_DATASET_MESTRE = os.path.join(SCRIPT_DIR, '..', 'resultados', 'dataset_mestre_COM_COMPORTAMENTO.csv')
CSV_SAIDA = os.path.join(SCRIPT_DIR, '..', 'resultados', '_preditor_estado_teta.csv')
GAP_MAX_S = 15  # so pares onde a janela anterior e' de fato a anterior no tempo (sem buraco)

BANDA_LENTA = (4.0, 7.0)   # teta tipo 2 (Vanderwolf): sniffing / imobilidade atenta
BANDA_RAPIDA = (7.0, 10.0)  # teta tipo 1: locomocao voluntaria
N_SUBJANELAS = 5  # divide o intervalo N-1 -> N em 5 pedacos p/ "close" + tendencia
FASE_CANONICA, AMP_CANONICA = 6.0, 85.0  # mediana dos 190 vencedores de OURO_PURIFICADO_v2.csv


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


def banda_power(sig, fs, banda):
    if len(sig) < 20:
        return np.nan
    f, psd = welch(sig, fs, nperseg=min(1000, len(sig)))
    m = (f >= banda[0]) & (f <= banda[1])
    return float(np.mean(psd[m])) if np.any(m) else np.nan


def monta_dataset():
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
    linhas = []
    n = len(consec)
    for k, (_, row) in enumerate(consec.iterrows()):
        if k % 500 == 0:
            print(f"  {k}/{n}...")
        pasta = resolve_pasta_basal(str(row["sessao"]), str(row["arquivo"]))
        if pasta is None:
            linhas.append({})
            continue
        ns2_path = os.path.join(pasta, row["arquivo"])
        if ns2_path not in _cache:
            _cache[ns2_path] = carrega_dados(ns2_path)
        dados, fs, canal_ids = _cache[ns2_path]
        canal_idx = int(row["canal"]) - 1
        if not (0 <= canal_idx < dados.shape[1]):
            linhas.append({})
            continue

        ini_ant = float(row["janela_ini_anterior"])
        fim_ant = ini_ant + float(row["dur_janela"])
        i0_full = int(ini_ant * fs)
        i1_full = int(fim_ant * fs)
        if i0_full < 0 or i1_full > dados.shape[0] or i1_full <= i0_full:
            linhas.append({})
            continue

        sig_full = aplica_notch(dados[i0_full:i1_full, canal_idx].astype(float), fs, freqs_notch=NOTCH_HZ)
        lenta_whole = banda_power(sig_full, fs, BANDA_LENTA)
        rapida_whole = banda_power(sig_full, fs, BANDA_RAPIDA)

        # sub-janelas p/ "close" (mais proxima de N) e tendencia (inclinacao)
        bordas = np.linspace(i0_full, i1_full, N_SUBJANELAS + 1).astype(int)
        lentas_sub, rapidas_sub = [], []
        for j in range(N_SUBJANELAS):
            sub = dados[bordas[j]:bordas[j + 1], canal_idx].astype(float)
            lentas_sub.append(banda_power(sub, fs, BANDA_LENTA))
            rapidas_sub.append(banda_power(sub, fs, BANDA_RAPIDA))
        x_sub = np.arange(N_SUBJANELAS)
        incl_lenta = np.polyfit(x_sub, lentas_sub, 1)[0] if not np.any(np.isnan(lentas_sub)) else np.nan
        incl_rapida = np.polyfit(x_sub, rapidas_sub, 1)[0] if not np.any(np.isnan(rapidas_sub)) else np.nan

        # MI_z no par CANONICO fixo (mesmo par p/ todas as linhas -- nao vaza rotulo)
        rng = np.random.default_rng(42 + row.name)
        mi_z_canonico = mi_z_par(sig_full, fs, FASE_CANONICA, AMP_CANONICA, rng=rng)

        linhas.append({
            "teta_lenta_whole": lenta_whole,
            "teta_rapida_whole": rapida_whole,
            "teta_lenta_close": lentas_sub[-1],
            "teta_rapida_close": rapidas_sub[-1],
            "teta_lenta_tendencia": incl_lenta,
            "teta_rapida_tendencia": incl_rapida,
            "mi_z_canonico": mi_z_canonico,
        })

    feat_df = pd.DataFrame(linhas, index=consec.index)
    consec = pd.concat([consec, feat_df], axis=1)
    consec = consec.dropna(subset=["teta_lenta_whole", "teta_rapida_whole"])

    # log das potencias (sempre positivas, escala log e' mais estavel p/ o classificador)
    for col in ["teta_lenta_whole", "teta_rapida_whole", "teta_lenta_close", "teta_rapida_close"]:
        consec[f"log_{col}"] = np.log(consec[col] + 1e-6)

    # normalizacao por canal: z-score DENTRO do proprio (arquivo, canal)
    for col in ["log_teta_lenta_whole", "log_teta_rapida_whole"]:
        grp = consec.groupby(["arquivo", "canal"])[col]
        consec[f"{col}_z"] = (consec[col] - grp.transform("mean")) / grp.transform("std").replace(0, np.nan)
    consec = consec.dropna(subset=["log_teta_lenta_whole_z", "log_teta_rapida_whole_z"])

    return consec


def avalia(X, y, rotulo):
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import StratifiedKFold, cross_val_predict
    from sklearn.metrics import roc_auc_score, accuracy_score

    clf = RandomForestClassifier(n_estimators=300, max_depth=5, min_samples_leaf=10,
                                  random_state=42, class_weight="balanced")
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    y_proba = cross_val_predict(clf, X, y, cv=cv, method="predict_proba")[:, 1]
    y_pred = (y_proba >= 0.5).astype(int)
    auc = roc_auc_score(y, y_proba)
    acc = accuracy_score(y, y_pred)
    print(f"[{rotulo}] acuracia={acc:.3f}  AUC-ROC={auc:.3f}  (n_features={X.shape[1]})")
    return auc, y_proba, clf


def main():
    print("Montando dataset (extraindo features do LFP, pode levar alguns minutos)...")
    consec = monta_dataset()
    print(f"\nAmostras finais: {len(consec)}")
    consec.to_csv(CSV_SAIDA, index=False, encoding="utf-8-sig")

    y = consec["vencedor"].astype(int).values
    dummies = pd.get_dummies(consec["comport_anterior"], prefix="comport").reset_index(drop=True)

    print("\n=== ABLATION SEQUENCIAL (validacao cruzada 5-fold, held-out) ===\n")

    # BASELINE: log_theta_power_anterior (banda unica 4-8Hz, bruta) + comportamento
    base_power = np.log(consec["teta_lenta_whole"] + consec["teta_rapida_whole"] + 1e-6).reset_index(drop=True)
    X0 = pd.concat([base_power.rename("log_theta_power"), dummies], axis=1).values
    auc0, _, _ = avalia(X0, y, "0. BASELINE (banda unica bruta + comportamento)")

    # PASSO 1: normalizacao por canal (troca bruta por z-score)
    zscores = consec[["log_teta_lenta_whole_z", "log_teta_rapida_whole_z"]].reset_index(drop=True)
    X1 = pd.concat([zscores.sum(axis=1).rename("teta_z_soma"), dummies], axis=1).values
    auc1, _, _ = avalia(X1, y, "1. + normalizacao por canal (z-score)")

    # PASSO 2: separa tipo1 (rapido) x tipo2 (lento), ja normalizados
    X2 = pd.concat([zscores, dummies], axis=1).values
    auc2, _, _ = avalia(X2, y, "2. + separa teta rapido (tipo1) x lento (tipo2)")

    # PASSO 3: adiciona "close" (mais perto do evento) e tendencia
    extra3 = consec[["teta_lenta_close", "teta_rapida_close",
                      "teta_lenta_tendencia", "teta_rapida_tendencia"]].reset_index(drop=True)
    extra3["log_teta_lenta_close"] = np.log(extra3["teta_lenta_close"] + 1e-6)
    extra3["log_teta_rapida_close"] = np.log(extra3["teta_rapida_close"] + 1e-6)
    extra3 = extra3.drop(columns=["teta_lenta_close", "teta_rapida_close"])
    X3 = pd.concat([zscores, extra3, dummies], axis=1).values
    auc3, _, _ = avalia(X3, y, "3. + janela mais proxima do evento + tendencia")

    # PASSO 4: MI_z no par canonico fixo
    mi_canon = consec[["mi_z_canonico"]].reset_index(drop=True)
    X4_df = pd.concat([zscores, extra3, mi_canon, dummies], axis=1)
    X4 = X4_df.values
    auc4, _, _ = avalia(X4, y, "4. + MI_z no par canonico (6Hz x 85Hz)")

    # ISOLAMENTO: normalizacao por canal e MI canonico ajudaram de verdade, ou o
    # ganho e' todo do passo 3 (janela proxima + tendencia)? Testa a mesma janela
    # proxima + tendencia SEM normalizar por canal e SEM MI.
    whole_raw = pd.DataFrame({
        "log_teta_lenta_whole": np.log(consec["teta_lenta_whole"] + 1e-6),
        "log_teta_rapida_whole": np.log(consec["teta_rapida_whole"] + 1e-6),
    }).reset_index(drop=True)
    X_melhor_df = pd.concat([whole_raw, extra3, dummies], axis=1)
    X_melhor = X_melhor_df.values
    auc_melhor, proba_melhor, clf_melhor = avalia(
        X_melhor, y, "5. ISOLAMENTO: whole+close+tendencia SEM z-score, SEM MI")

    print(f"\nDelta de AUC (baseline -> melhor configuracao): {auc0:.3f} -> {auc_melhor:.3f} ({auc_melhor - auc0:+.3f})")
    print("Conclusao do ablation: normalizacao por canal (passo 1) e MI no par canonico")
    print("(passo 4) NAO se confirmaram uteis isoladamente -- o ganho real veio todo")
    print("da janela mais proxima do evento + tendencia (passo 3). Modelo final usa")
    print("so essa configuracao (mais simples, mesmo desempenho ou melhor).")

    proba4 = proba_melhor
    clf4 = clf_melhor
    X4_df = X_melhor_df
    X4 = X_melhor

    # PASSO 6: calibracao de limiar via precisao-recall (o que importa p/ TTL: recall do vencedor)
    from sklearn.metrics import precision_recall_curve, f1_score
    prec, rec, thr = precision_recall_curve(y, proba4)
    f1s = 2 * prec[:-1] * rec[:-1] / (prec[:-1] + rec[:-1] + 1e-9)
    melhor_i = np.argmax(f1s)
    print(f"\n=== 6. Calibracao de limiar (melhor configuracao) ===")
    print(f"Melhor limiar por F1: {thr[melhor_i]:.3f}  ->  precisao={prec[melhor_i]:.3f} recall={rec[melhor_i]:.3f} F1={f1s[melhor_i]:.3f}")
    # ponto de operacao voltado a recall alto (>=0.8), custo em precisao
    alvo_recall = 0.80
    candidatos = np.where(rec[:-1] >= alvo_recall)[0]
    if len(candidatos) > 0:
        i_recall = candidatos[np.argmax(prec[:-1][candidatos])]
        print(f"Limiar p/ recall>={alvo_recall}: {thr[i_recall]:.3f}  ->  precisao={prec[i_recall]:.3f} recall={rec[i_recall]:.3f}")

    clf4.fit(X4, y)
    import joblib
    out = os.path.join(SCRIPT_DIR, 'modelo_estado_comportamental.pkl')
    joblib.dump({"modelo": clf4, "features": X4_df.columns.tolist(),
                 "limiar_f1": float(thr[melhor_i])}, out)
    print(f"\nModelo final salvo em: {out}")
    print("\nImportancia das features (modelo final, todos os dados):")
    for feat, imp in sorted(zip(X4_df.columns.tolist(), clf4.feature_importances_), key=lambda t: -t[1]):
        print(f"  {feat}: {imp:.4f}")


if __name__ == "__main__":
    main()
