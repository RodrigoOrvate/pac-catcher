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

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(SCRIPT_DIR, '..'))
from pac_core.io import carrega_dados
from pac_core.filtering import aplica_notch
from pac_core.workspace import BASE_LAC_NOCI
from pipeline.etapa4_validacao.robustez_parametros import mi_z_par, mvl_z_par

NOTCH_HZ = [60.0, 120.0, 180.0, 240.0]
FOOTPRINT_N_SURR = 50   # varredura grosseira nos outros 31 canais (custo x32)

FS = 1000
PRE_WINDOW = 10
CONTROL_OFFSET_S = 60   # janela de controle começa event_time + offset
FEATURES = ['Energia', 'Theta_Energy', 'Gamma_Energy', 'Ratio_TG', 'Theta_Peak',
            'MI_z_oficial', 'MVL_z_oficial', 'Footprint_n_z3']
CSV_VENCEDORES = os.path.join(SCRIPT_DIR, '..', 'resultados', 'candidatos_vencedores_OURO_PURIFICADO_v2.csv')


def extrair_features_espectrais(window_data):
    """Janela (n_amostras, n_canais) -> 5 features espectrais simples (protótipo original)."""
    sig = np.mean(window_data, axis=1)
    energia = np.var(sig)
    f, psd = welch(sig, FS, nperseg=min(1000, len(sig)))

    tm = (f >= 4) & (f <= 8)
    te = np.mean(psd[tm]) if np.any(tm) else 0.0
    tp = f[tm][np.argmax(psd[tm])] if np.any(tm) else 0.0

    gm = (f >= 30) & (f <= 80)
    ge = np.mean(psd[gm]) if np.any(gm) else 0.0

    return [energia, te, ge, te / (ge + 1e-6), tp]


def extrair_features_pac(dados_janela, canal_idx, fs, fase_pico_hz, amp_pico_hz, rng):
    """3 features PAC-especificas na janela (n_amostras, n_canais_totais):
    - MI_z_oficial / MVL_z_oficial: ja existe acoplamento (no par que a
      sessao eventualmente mostra) se formando NESTE canal, nesta janela?
    - Footprint_n_z3: quantos OUTROS canais ja mostram z>=3 no mesmo par,
      na MESMA janela -- precursor espacial (onda se formando na rede)?
    """
    lfp_alvo = aplica_notch(dados_janela[:, canal_idx].astype(float), fs, freqs_notch=NOTCH_HZ)
    mi_z = mi_z_par(lfp_alvo, fs, fase_pico_hz, amp_pico_hz, rng=rng)
    mvl_z = mvl_z_par(lfp_alvo, fs, fase_pico_hz, amp_pico_hz, rng=rng)

    n_footprint = 0
    for c in range(dados_janela.shape[1]):
        if c == canal_idx:
            continue
        lfp_c = aplica_notch(dados_janela[:, c].astype(float), fs, freqs_notch=NOTCH_HZ)
        z_c = mi_z_par(lfp_c, fs, fase_pico_hz, amp_pico_hz, n_surr=FOOTPRINT_N_SURR, rng=rng)
        if z_c >= 3.0:
            n_footprint += 1

    return [mi_z, mvl_z, n_footprint]


def extrair_features(dados_janela, canal_idx, fs, fase_pico_hz, amp_pico_hz, rng):
    return np.array(extrair_features_espectrais(dados_janela) +
                     extrair_features_pac(dados_janela, canal_idx, fs, fase_pico_hz, amp_pico_hz, rng))


def resolve_pasta_basal(sessao_str, arquivo):
    """Mesma lógica usada em gera_galeria_top5.py / checa_desalinhamento_190.py."""
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


def carregar_positivos_e_controles():
    """Le o resultado ATUAL do pipeline (candidatos_vencedores_OURO_PURIFICADO_v2.csv,
    190 eventos pos-purificacao por notch) em vez do legado RESULTADOS/vencedores.csv
    (schema antigo, ~10 eventos). Pre-evento = positivo, janela 'controle' real no
    MESMO arquivo (event_time + 60 s) = negativo."""
    X, y, meta = [], [], []
    df = pd.read_csv(CSV_VENCEDORES)
    _cache = {}

    for _, row in df.iterrows():
        canal_idx = int(row['canal']) - 1  # convencao 1-based do dataset mestre
        t_evento = float(row['janela_ini_s'])
        fp = float(row['fase_pico_hz'])
        fa = float(row['amp_pico_hz'])
        rng = np.random.default_rng(42 + row.name)
        pasta = resolve_pasta_basal(str(row['sessao']), str(row['arquivo']))
        if pasta is None:
            print(f"  [skip] sessao={row['sessao']!r} arquivo={row['arquivo']!r} nao resolvido")
            continue
        ns2_path = os.path.join(pasta, row['arquivo'])

        if ns2_path not in _cache:
            _cache[ns2_path] = carrega_dados(ns2_path)
        dados, fs, canal_ids = _cache[ns2_path]
        if not (0 <= canal_idx < dados.shape[1]):
            print(f"  [skip] canal {row['canal']} fora do range em {ns2_path}")
            continue

        # POSITIVO: 10 s antes do evento
        i0 = int((t_evento - PRE_WINDOW) * fs)
        i1 = i0 + int(PRE_WINDOW * fs)
        if i0 >= 0 and i1 <= dados.shape[0]:
            X.append(extrair_features(dados[i0:i1, :], canal_idx, fs, fp, fa, rng))
            y.append(1)
            meta.append((os.path.basename(ns2_path), row['canal'], 'pre'))

        # CONTROLE (negativo real): 10 s em outro instante do mesmo arquivo
        c0 = int((t_evento + CONTROL_OFFSET_S) * fs)
        c1 = c0 + int(PRE_WINDOW * fs)
        if c1 <= dados.shape[0]:
            X.append(extrair_features(dados[c0:c1, :], canal_idx, fs, fp, fa, rng))
            y.append(0)
            meta.append((os.path.basename(ns2_path), row['canal'], 'controle'))
        else:
            # se o offset estoura o arquivo, pega mais cedo
            c0 = max(0, int((t_evento - CONTROL_OFFSET_S - PRE_WINDOW) * fs))
            c1 = c0 + int(PRE_WINDOW * fs)
            if c1 <= dados.shape[0]:
                X.append(extrair_features(dados[c0:c1, :], canal_idx, fs, fp, fa, rng))
                y.append(0)
                meta.append((os.path.basename(ns2_path), row['canal'], 'controle'))

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

    print("\nImportancia das features (modelo final, todos os dados):")
    for feat, imp in sorted(zip(FEATURES, clf.feature_importances_), key=lambda t: -t[1]):
        print(f"  {feat}: {imp:.4f}")


if __name__ == "__main__":
    main()
