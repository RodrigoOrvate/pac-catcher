"""
pac_core/io.py
==========================================
Módulo compartilhado de leitura de dados eletrofisiológicos e salvamento
de CSV, usado por todo o pipeline (triagem_pac, refina_candidatos,
comodulogram, auditorias/*, preditor/*).

Lê .ns2 (Blackrock) diretamente via `neo`, sem precisar do extrator.exe.
Também aceita .bin legado (int16 puro, já extraído), para quem ainda
tiver arquivos processados pelo pipeline antigo, e .mat (MATLAB/tetrodo),
para o fluxo de `triagem_pac_mat.py`.

Requer: neo, numpy, scipy
    pip install neo numpy scipy
"""

import os
import numpy as np


def le_ns2(caminho_arquivo):
    """
    Lê um .ns2 usando neo.rawio.BlackrockRawIO (leitura preguiçosa).
    load_nev=False: não carrega o .nev (eventos TTL) junto. Isso evita o erro
    "Inconsistent ns2 and nev file" quando a gravação foi pausada/retomada e
    os dois arquivos ficam com número de segmentos diferente -- para a
    triagem de PAC só precisamos do sinal contínuo, não dos eventos.
    Retorna: dados (n_amostras x n_canais), fs (Hz), nomes_canais (list)
    """
    import neo.rawio as neorawio

    reader = neorawio.BlackrockRawIO(filename=caminho_arquivo, load_nev=False)
    reader.parse_header()

    fs = reader.get_signal_sampling_rate(stream_index=0)
    n_amostras = reader.get_signal_size(block_index=0, seg_index=0, stream_index=0)

    dados = reader.get_analogsignal_chunk(
        block_index=0, seg_index=0,
        i_start=0, i_stop=n_amostras,
        stream_index=0,
    )
    canal_ids = list(reader.header["signal_channels"]["name"])
    return dados, float(fs), canal_ids


def le_bin_legado(caminho_arquivo, n_canais, fs=1000.0):
    """
    Lê o formato .bin antigo (int16 puro, extraído pelo extrator.exe).
    Retorna: dados (n_amostras x n_canais), fs (Hz), nomes_canais (list)
    """
    dados = np.fromfile(caminho_arquivo, dtype=np.int16)
    if n_canais > 1:
        dados = dados.reshape(-1, n_canais)
    else:
        dados = dados.reshape(-1, 1)
    canal_ids = [f"Ch{i+1}" for i in range(dados.shape[1])]
    return dados, float(fs), canal_ids


def fs_do_mat(m):
    """
    Detecta a taxa de amostragem em um dicionário do scipy.io.loadmat.
    Tenta variáveis comuns: fs, Fs, srate, samplingRate.
    Retorna 1000.0 se não encontrar.
    """
    for chave in ["fs", "Fs", "srate", "samplingRate", "SR"]:
        if chave in m:
            val = m[chave]
            try:
                return float(np.squeeze(val))
            except Exception:
                pass
    return 1000.0


def le_mat(path, chave=None):
    """
    Carrega um LFP de um arquivo .mat (MATLAB/tetrodo).

    Se 'chave' é None, tenta na ordem:
      lfp, LFP, lfpBruto, lfpHG, lfpHFO, signal, data

    Retorna (array_float64, fs_float, chave_usada).

    NOTA: não garante 1D -- se a chave encontrada apontar para um array 2D
    genuinamente multi-canal (ex.: FOOOF/CA1_example.mat, shape (N,4)),
    o array retornado também é 2D (np.squeeze não achata dimensões >1).
    Esse comportamento é herdado de adapta_lfp_mat.py sem alteração --
    hoje nenhum consumidor real aponta para arquivos desse formato por
    este caminho.
    """
    import scipy.io

    m = scipy.io.loadmat(path)
    chaves_candidatas = (
        [chave] if chave else
        ["lfp", "LFP", "lfpBruto", "lfpHG", "lfpHFO", "signal", "data"]
    )
    for c in chaves_candidatas:
        if c in m:
            arr = np.squeeze(m[c]).astype(np.float64)
            fs = fs_do_mat(m)
            return arr, fs, c
    raise KeyError(
        f"Nenhuma das chaves {chaves_candidatas} encontrada em {path}.\n"
        f"Chaves disponíveis: {[k for k in m if not k.startswith('_')]}"
    )


def normaliza_sinal(arr):
    """
    Normalização robusta: (x - median) / MAD.
    Mantém a forma do sinal, remove outliers de offset, escala para ~unidades.
    Preferida ao z-score mean/std quando há outliers (artefatos).
    """
    med = np.median(arr)
    mad = np.median(np.abs(arr - med))
    if mad < 1e-12:
        # Fallback para std se MAD for zero (sinal constante)
        std = arr.std()
        if std < 1e-12:
            return arr - med  # sinal constante
        return (arr - med) / std
    return (arr - med) / mad


def info_sinal(arr, nome="sinal", fs=None):
    """
    Imprime diagnóstico de escala do sinal. Essencial para detectar
    problemas de normalização antes de passar pelo pipeline.
    """
    duracao = len(arr) / fs if fs else None
    print(f"  [{nome}]")
    print(f"    shape: {arr.shape}, dtype: {arr.dtype}")
    if duracao:
        print(f"    duração: {duracao:.1f} s @ {fs:.0f} Hz")
    print(f"    min={arr.min():.4g}  max={arr.max():.4g}")
    print(f"    mean={arr.mean():.4g}  std={arr.std():.4g}")
    med = np.median(arr)
    mad = np.median(np.abs(arr - med))
    print(f"    median={med:.4g}  MAD={mad:.4g}")
    # Alerta de escala
    if arr.std() > 5000:
        print(f"    *** AVISO: std={arr.std():.1f} muito alto — provável escala em nV ou "
              f"contagem ADC. Considere --normaliza ou converta para µV.")
    elif arr.std() < 0.001:
        print(f"    *** AVISO: std={arr.std():.6f} muito baixo — provável escala em V.")
    else:
        print(f"    Escala aparenta ser µV (std razoável para LFP).")


def carrega_dados(caminho_arquivo, n_canais_bin=16, fs_bin=1000.0):
    """
    Dispatcher automático: escolhe o leitor certo pela extensão do arquivo.
    - .ns2         -> le_ns2 (leitura direta, sem extrator.exe)
    - .bin/.dat    -> le_bin_legado (formato antigo já extraído)
    - .mat         -> le_mat (MATLAB/tetrodo)

    n_canais_bin e fs_bin só são usados no caminho .bin legado, onde essa
    informação não está no cabeçalho do arquivo (precisa ser informada
    manualmente, como nos scripts originais).
    """
    ext = os.path.splitext(caminho_arquivo)[1].lower()
    if ext == ".ns2":
        return le_ns2(caminho_arquivo)
    elif ext in (".bin", ".dat"):
        return le_bin_legado(caminho_arquivo, n_canais_bin, fs_bin)
    elif ext == ".mat":
        arr, fs, chave = le_mat(caminho_arquivo)
        if arr.ndim == 1:
            dados = arr.reshape(-1, 1)
            canal_ids = [chave]
        else:
            dados = arr
            canal_ids = [f"{chave}{i}" for i in range(dados.shape[1])]
        return dados, fs, canal_ids
    else:
        raise ValueError(
            f"Extensão '{ext}' não reconhecida. Use .ns2 (direto), "
            f".bin/.dat (formato legado extraído pelo extrator.exe) ou "
            f".mat (MATLAB/tetrodo)."
        )


def fatia_janela(dados, fs, t_inicio_seg, t_fim_seg):
    """Recorta uma janela de tempo [t_inicio_seg, t_fim_seg) de `dados`."""
    idx_inicio = int(t_inicio_seg * fs)
    idx_fim = int(t_fim_seg * fs)
    return dados[idx_inicio:idx_fim]


def concatena_sessao(pasta_basal, fs_esperado=1000.0):
    """
    Lê todos os arquivos .ns2 (ou .bin) da pasta `pasta_basal` e concatena
    em uma matriz única. Usado pelo passo 0 (exploração visual) e pelo
    `comodulogram_interativo.py` para ter o registro completo de uma sessão
    em uma única estrutura.

    Parâmetros
    ----------
    pasta_basal : str
        Caminho da pasta com os arquivos .ns2 (ou .bin) consecutivos
        (ex.: a pasta "Basal antes da infusao" de uma sessão MTESC04_NOCI).
    fs_esperado : float
        Frequência de amostragem esperada (default 1000 Hz). Usado só para
        sanity check (assertion) — se os arquivos tiverem fs diferente, a
        função avisa mas não falha.

    Retorna
    -------
    dados : np.ndarray, shape (n_amostras_total, n_canais)
    fs : float
    ids_canais : list of str
    offsets : list of int
        Índice (em amostras) de início de cada arquivo dentro de `dados`.
        Ex.: offsets = [0, 300000, 600000] significa que o segundo arquivo
        começa em 300 000 amostras e o terceiro em 600 000.
    """
    arquivos = sorted([
        os.path.join(pasta_basal, f)
        for f in os.listdir(pasta_basal)
        if f.lower().endswith(('.ns2', '.bin', '.dat'))
    ])
    if not arquivos:
        raise FileNotFoundError(
            f"Nenhum .ns2/.bin/.dat em {pasta_basal}"
        )

    dados_lista = []
    ids_canais = None
    fs = None
    for arq in arquivos:
        d, f, ids = carrega_dados(arq)
        if ids_canais is None:
            ids_canais = ids
            fs = f
        else:
            # sanity check: mesma lista de canais e mesmo fs
            if f != fs:
                raise ValueError(
                    f"Frequência de amostragem inconsistente: {f} Hz em "
                    f"{arq} vs {fs} Hz nos anteriores"
                )
        dados_lista.append(d)

    dados = np.vstack(dados_lista)
    offsets = np.cumsum([0] + [d.shape[0] for d in dados_lista[:-1]]).tolist()

    if fs is not None and abs(fs - fs_esperado) > 1.0:
        print(f"  AVISO: fs={fs} Hz (esperado {fs_esperado} Hz)")

    return dados, float(fs), ids_canais, offsets


def salva_csv(df, caminho, **kwargs):
    """
    Salva um DataFrame em CSV com encoding utf-8-sig (BOM) por padrão, para
    que o Excel no Windows abra corretamente caracteres acentuados em
    português (Exploração, Locomoção, Imóvel, Atenção...).

    Helper opcional: nenhum `to_csv()` existente no pipeline foi migrado
    para usar isto automaticamente. Use em scripts novos ou ao migrar um
    call site existente deliberadamente.
    """
    kwargs.setdefault("encoding", "utf-8-sig")
    df.to_csv(caminho, **kwargs)
