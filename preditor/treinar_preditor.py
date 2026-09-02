"""
treinar_preditor.py
Treina um modelo de Machine Learning para prever acoplamento theta-gamma
com base em features espectrais simples (Razão Theta/Gamma, Energia, etc.).
"""
import pandas as pd
import numpy as np
import glob
import os
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
import joblib

FS = 1000
PRE_WINDOW = 10

def extract_features(window_data):
    """Extrai features de uma janela de sinal (amostras x canais)."""
    # Média dos canais para simplificar (ou usar o canal principal)
    # window_data vem como (n_amostras, n_canais) -> média ao longo das colunas
    sig = np.mean(window_data, axis=1)

    # 1. Energia Total (Variância)
    energy = np.var(sig)

    # 2. Espectro de Potência (Welch)
    from scipy.signal import welch
    # Aumentar nperseg para melhor resolução em frequência (1s = 1000 amostras)
    f, psd = welch(sig, FS, nperseg=min(1000, len(sig)))

    # 3. Energia Theta (4-8 Hz)
    theta_mask = (f >= 4) & (f <= 8)
    if np.any(theta_mask):
        theta_energy = np.mean(psd[theta_mask])
        theta_peak = f[theta_mask][np.argmax(psd[theta_mask])]
    else:
        theta_energy = 0.0
        theta_peak = 0.0

    # 4. Energia Gamma (30-80 Hz)
    gamma_mask = (f >= 30) & (f <= 80)
    if np.any(gamma_mask):
        gamma_energy = np.mean(psd[gamma_mask])
    else:
        gamma_energy = 0.0

    # 5. Razão Theta/Gamma (Feature Principal)
    ratio_tg = theta_energy / (gamma_energy + 1e-6)

    return [energy, theta_energy, gamma_energy, ratio_tg, theta_peak]

def _surrogate_phase_scramble(data, seed_offset=0):
    """Surrogate: embaralha a fase de cada canal via FFT, destruindo qualquer
    acoplamento temporal mas preservando o espectro de potência. Serve de
    baseline realista (não é ruído branco puro)."""
    rng = np.random.default_rng(42 + seed_offset)
    out = np.empty_like(data)
    for c in range(data.shape[1]):
        spec = np.fft.rfft(data[:, c])
        phase = rng.uniform(0, 2 * np.pi, spec.shape)
        spec_permut = spec * np.exp(1j * phase)
        out[:, c] = np.fft.irfft(spec_permut, n=data.shape[0])
    return out

def load_training_data():
    """Carrega dados positivos (pré-evento) e gera negativos (baseline surrogado)."""
    base = r"C:\acoplamento_theta-gamma\ANALISE_PRE_EVENTO"
    raw_files = glob.glob(os.path.join(base, "raw_*.csv"))

    X, y = [], []

    for i, f in enumerate(raw_files):
        data = pd.read_csv(f).values

        # 1. Dado POSITIVO (os 10s antes do evento)
        X.append(extract_features(data))
        y.append(1) # 1 = PAC iminente

        # 2. Dado NEGATIVO (surrogate phase-scrambled: mesmo espectro, sem acoplamento)
        surrogate = _surrogate_phase_scramble(data, seed_offset=i)
        X.append(extract_features(surrogate))
        y.append(0) # 0 = Sem PAC

    return np.array(X), np.array(y)

def train_model():
    from sklearn.model_selection import LeaveOneOut, cross_val_score, StratifiedKFold

    print("Carregando e processando dados...")
    X, y = load_training_data()

    print(f"Amostras: {len(X)} (Positivos: {sum(y)}, Negativos: {len(y)-sum(y)})")

    if len(X) < 20:
        print("Dados insuficientes para train/test split (n < 20). Usando Leave-One-Out Cross-Validation (LOOCV).")
        cv = LeaveOneOut()
    else:
        print("Usando Stratified 5-Fold Cross-Validation.")
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    print("Treinando Random Forest com Cross-Validation...")
    clf = RandomForestClassifier(n_estimators=100, random_state=42)

    scores = cross_val_score(clf, X, y, cv=cv, scoring='accuracy')
    print(f"\nResultados da Validação Cruzada:")
    print(f"Acurácia por fold: {scores}")
    print(f"Acurácia Média: {np.mean(scores):.3f} (+/- {np.std(scores):.3f})")

    # Treinar o modelo final com TODOS os dados para salvar
    print("\nTreinando modelo final com todos os dados...")
    clf.fit(X, y)

    # Salvar modelo
    output_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "modelo_pac.pkl")
    joblib.dump(clf, output_path)
    print(f"Modelo final (treinado em todos os dados) salvo em: {output_path}")

    # Importância das features no modelo final
    features = ['Energia', 'Theta_Energy', 'Gamma_Energy', 'Ratio_TG', 'Theta_Peak']
    importances = clf.feature_importances_
    print("\nImportância das features (modelo final):")
    for f, imp in zip(features, importances):
        print(f"  {f}: {imp:.4f}")

if __name__ == "__main__":
    train_model()
