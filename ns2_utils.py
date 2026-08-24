"""
ns2_utils.py
==========================================
Módulo compartilhado de leitura de dados eletrofisiológicos, usado por
leitor_bin.py, comodulogram.py e triagem_pac.py.

Lê .ns2 (Blackrock) diretamente via `neo`, sem precisar do extrator.exe.
Também aceita .bin legado (int16 puro, já extraído), para quem ainda
tiver arquivos processados pelo pipeline antigo.

Requer: neo, numpy
    pip install neo numpy
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


def carrega_dados(caminho_arquivo, n_canais_bin=16, fs_bin=1000.0):
    """
    Dispatcher automático: escolhe o leitor certo pela extensão do arquivo.
    - .ns2         -> le_ns2 (leitura direta, sem extrator.exe)
    - .bin/.dat    -> le_bin_legado (formato antigo já extraído)

    n_canais_bin e fs_bin só são usados no caminho .bin legado, onde essa
    informação não está no cabeçalho do arquivo (precisa ser informada
    manualmente, como nos scripts originais).
    """
    ext = os.path.splitext(caminho_arquivo)[1].lower()
    if ext == ".ns2":
        return le_ns2(caminho_arquivo)
    elif ext in (".bin", ".dat"):
        return le_bin_legado(caminho_arquivo, n_canais_bin, fs_bin)
    else:
        raise ValueError(
            f"Extensão '{ext}' não reconhecida. Use .ns2 (direto) ou "
            f".bin/.dat (formato legado extraído pelo extrator.exe)."
        )


def fatia_janela(dados, fs, t_inicio_seg, t_fim_seg):
    """Recorta uma janela de tempo [t_inicio_seg, t_fim_seg) de `dados`."""
    idx_inicio = int(t_inicio_seg * fs)
    idx_fim = int(t_fim_seg * fs)
    return dados[idx_inicio:idx_fim]