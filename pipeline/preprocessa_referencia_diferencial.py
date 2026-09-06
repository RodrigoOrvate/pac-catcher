import numpy as np
from scipy.signal import butter, filtfilt, hilbert

def seleciona_pool_referencia(dados, fs, flags_coocorrencia_por_canal, canal_alvo_idx,
                              banda_ripple=(150, 250), n_canais_pool=8, 
                              limiar_amplitude_minima=1e-6): # um chute inicial
    """
    Seleciona um subconjunto de canais "silenciosos" para servir de
    referência diferencial, sem exigir mapeamento anatômico manual.
    canal_alvo_idx eh EXCLUIDO do pool.
    """
    n_canais = dados.shape[1]
    candidatos = []

    for c in range(n_canais):
        # 0. Excluir o proprio canal alvo para evitar auto-subtracao
        if c == canal_alvo_idx:
            continue

        # 1. Excluir canais já sinalizados (pode ser dict vazio se não rodou triagem completa)
        flags = flags_coocorrencia_por_canal.get(f'chan{c+1}', {})
        if flags.get("teta_hg", False) or flags.get("teta_gama", False):
            continue  # provavelmente dentro do hipocampo

        sinal_c = dados[:, c]
        amplitude_rms_total = np.sqrt(np.mean(sinal_c**2))
        
        # 2. Excluir canal morto
        if limiar_amplitude_minima and amplitude_rms_total < limiar_amplitude_minima:
            continue

        # 3. Potência na banda de ripple
        b, a = butter(4, [banda_ripple[0]/(fs/2), banda_ripple[1]/(fs/2)], btype='band')
        sinal_filtrado = filtfilt(b, a, sinal_c)
        env_ripple = np.abs(hilbert(sinal_filtrado))
        potencia_ripple = np.mean(env_ripple**2)

        candidatos.append((c, potencia_ripple, amplitude_rms_total))

    # Ordenar por menor potência de ripple primeiro
    candidatos.sort(key=lambda x: x[1])
    
    # Previne erro se sobraram menos canais que o n_canais_pool
    num_selecionar = min(n_canais_pool, len(candidatos))
    
    if num_selecionar == 0:
        return []

    pool_idx = [c for c, _, _ in candidatos[:num_selecionar]]
    return pool_idx

def constroi_referencia(dados, pool_idx):
    """
    Media dos canais do pool.
    """
    if not pool_idx:
        return np.zeros(dados.shape[0])
    return dados[:, pool_idx].mean(axis=1)

def aplica_referencia_diferencial(sinal_alvo, sinal_referencia):
    """
    Subtrai o sinal de referência do sinal alvo.
    """
    if sinal_referencia is None:
        return sinal_alvo
    return sinal_alvo - sinal_referencia
