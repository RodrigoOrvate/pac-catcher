"""
triagem_pac.py
==========================================
Varredura automática de arquivos .ns2 (Blackrock) para flagar janelas
de tempo com suspeita de Acoplamento Fase-Amplitude (PAC) Theta-Gamma,
ANTES de checar o vídeo comportamental.

DIFERENÇA CRÍTICA em relação aos scripts anteriores (simu_proc_PAC.py,
comodulogram.py, leitor_bin.py): aqueles calculam o MI "bruto". Isso é
correto para validar o algoritmo em sinal simulado, mas é insuficiente
para decidir "esse arquivo/janela tem acoplamento real" em dados in vivo,
porque:

  1. MI bruto não tem escala interpretável — não dá pra saber se 0.03 é
     muito ou pouco sem saber o que ruído puro geraria na mesma janela.
  2. Theta hipocampal é assimétrico (sawtooth), não senoidal, e isso
     por si só gera MI espúrio (Kramer et al. 2008; Cole & Voytek 2017).
  3. Ruído muscular/movimento contamina exatamente a banda Gamma (30-80Hz)
     durante os comportamentos que você quer flagar (sniffing, rearing,
     locomoção) — o que cria risco de circularidade: flag alto -> vídeo
     mostra movimento -> "confirma" acoplamento que na verdade é artefato.

Este script resolve (1) e (2) comparando o MI observado contra uma
distribuição NULA gerada por surrogates (deslocamento circular do
envelope de Gamma), reportando um z-score e p-valor empírico em vez do
MI bruto. Para (3), como você não tem canal de EMG/acelerômetro, o script
calcula um PROXY heurístico de artefato motor (potência de banda larga em
alta frequência, ~150-450 Hz) por janela — isso NÃO substitui um canal de
EMG real, mas permite descartar candidatos onde o z-MI alto coincide com
alta energia de banda larga (assinatura clássica de artefato muscular).

Requer: neo, numpy, scipy, pandas
    pip install neo numpy scipy pandas

Uso:
    # Varrer uma pasta inteira de .ns2, todos os canais
    python triagem_pac.py --pasta /caminho/para/ns2 --saida resultados.csv

    # Só um canal específico, janelas de 15s com passo de 5s
    python triagem_pac.py --pasta /caminho/para/ns2 --canal 3 \
        --janela 15 --passo 5 --saida resultados.csv

    # Rodar o autoteste com dados sintéticos (sem precisar de arquivo real)
    python triagem_pac.py --demo
"""

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
import scipy.signal as signal

# ==========================================
# FUNÇÕES DE PROCESSAMENTO (mesma lógica dos scripts originais)
# ==========================================

def filtra_sinal(sinal, lowcut, highcut, fs, order=3):
    nyq = 0.5 * fs
    low = max(lowcut / nyq, 1e-6)
    high = min(highcut / nyq, 0.999)
    b, a = signal.butter(order, [low, high], btype="bandpass")
    return signal.filtfilt(b, a, sinal)


def calcula_mi(fase, envelope, n_bins=18):
    bins = np.linspace(-np.pi, np.pi, n_bins + 1)
    bin_idx = np.clip(np.digitize(fase, bins) - 1, 0, n_bins - 1)
    return _mi_de_bin_idx(bin_idx, envelope, n_bins)


def _mi_de_bin_idx(bin_idx, envelope, n_bins):
    """Núcleo vetorizado: soma/conta por bin via bincount em vez de loop Python."""
    soma_bins = np.bincount(bin_idx, weights=envelope, minlength=n_bins)
    cont_bins = np.bincount(bin_idx, minlength=n_bins)
    media_bins = np.divide(soma_bins, cont_bins, out=np.zeros(n_bins), where=cont_bins > 0)

    soma = np.sum(media_bins)
    if soma <= 0:
        return 0.0
    P = media_bins / soma
    H = -np.sum(P * np.log(P + 1e-10))
    return (np.log(n_bins) - H) / np.log(n_bins)


def mi_com_surrogates(fase, envelope, fs, n_surr=200, n_bins=18,
                       shift_min_s=1.0, rng=None):
    """
    Calcula o MI observado e o compara contra n_surr surrogates gerados
    por deslocamento circular do envelope de Gamma (quebra a relação
    temporal fase-amplitude mas preserva o espectro de cada sinal).

    Versão vetorizada: o índice de bin da fase é calculado uma única vez
    (não muda entre surrogates, só o envelope é deslocado), e a soma por
    bin usa np.bincount em vez de loop Python -- ~10-20x mais rápido que
    a versão ingênua com np.where em loop.

    Retorna: mi_obs, z_score, p_empirico, mi_surr_media, mi_surr_dp
    """
    if rng is None:
        rng = np.random.default_rng()

    bins = np.linspace(-np.pi, np.pi, n_bins + 1)
    bin_idx = np.clip(np.digitize(fase, bins) - 1, 0, n_bins - 1)

    mi_obs = _mi_de_bin_idx(bin_idx, envelope, n_bins)

    n = len(envelope)
    shift_min = int(shift_min_s * fs)
    if n <= 2 * shift_min:
        shift_min = max(1, n // 10)

    mi_surr = np.empty(n_surr)
    deslocamentos = rng.integers(shift_min, n - shift_min, size=n_surr)
    for i, desloc in enumerate(deslocamentos):
        env_shift = np.roll(envelope, desloc)
        mi_surr[i] = _mi_de_bin_idx(bin_idx, env_shift, n_bins)

    media = np.mean(mi_surr)
    dp = np.std(mi_surr)
    z = (mi_obs - media) / dp if dp > 0 else 0.0
    p_emp = np.mean(mi_surr >= mi_obs)

    return mi_obs, z, p_emp, media, dp


def proxy_artefato_motor(sinal, fs, banda=(150, 450)):
    """
    Proxy heurístico de artefato muscular/movimento: potência relativa
    em banda larga de alta frequência (acima de Gamma), onde EMG
    tipicamente vaza no LFP. NÃO substitui EMG real -- serve só para
    marcar candidatos suspeitos de contaminação motora.
    """
    nyq = 0.5 * fs
    alta = filtra_sinal(sinal, banda[0], min(banda[1], nyq * 0.98), fs, order=3)
    pot_alta = np.mean(alta ** 2)
    pot_total = np.mean(sinal.astype(float) ** 2) + 1e-12
    return pot_alta / pot_total


# ==========================================
# VARREDURA DE UM CANAL / JANELA
# ==========================================

def varre_canal(sinal, fs, window_s=10.0, step_s=5.0, n_surr=200,
                 theta_band=(4, 8), gamma_band=(30, 80), rng=None,
                 rotulo_progresso=None):
    win = int(window_s * fs)
    step = int(step_s * fs)
    resultados = []

    inicios = list(range(0, len(sinal) - win + 1, step))
    total = len(inicios)

    for n_janela, ini in enumerate(inicios):
        fim = ini + win
        trecho = sinal[ini:fim].astype(float)

        lfp_theta = filtra_sinal(trecho, *theta_band, fs)
        lfp_gamma = filtra_sinal(trecho, *gamma_band, fs)

        fase = np.angle(signal.hilbert(lfp_theta))
        env = np.abs(signal.hilbert(lfp_gamma))

        mi_obs, z, p_emp, mi_surr_m, mi_surr_dp = mi_com_surrogates(
            fase, env, fs, n_surr=n_surr, rng=rng
        )
        artefato = proxy_artefato_motor(trecho, fs)

        resultados.append({
            "janela_ini_s": ini / fs,
            "janela_fim_s": fim / fs,
            "mi_observado": mi_obs,
            "mi_surrogate_media": mi_surr_m,
            "mi_surrogate_dp": mi_surr_dp,
            "z_score": z,
            "p_empirico": p_emp,
            "proxy_artefato_motor": artefato,
        })

        if rotulo_progresso and (n_janela % 10 == 0 or n_janela == total - 1):
            print(f"  {rotulo_progresso}: janela {n_janela + 1}/{total}", end="\r")

    if rotulo_progresso:
        print()  # nova linha ao terminar o canal

    return pd.DataFrame(resultados)


# ==========================================
# LEITURA DE .ns2 VIA NEO (Blackrock)
# ==========================================

from ns2_utils import le_ns2  # leitura compartilhada com leitor_bin.py / comodulogram.py


def varre_arquivo(caminho, canais=None, window_s=10.0, step_s=5.0,
                   n_surr=200, seed=None):
    rng = np.random.default_rng(seed)
    dados, fs, nomes_canais = le_ns2(caminho)
    n_canais_total = dados.shape[1]

    if canais is None:
        canais_idx = range(n_canais_total)
    else:
        canais_idx = canais

    todos = []
    for ch in canais_idx:
        sinal = dados[:, ch]
        nome_canal = nomes_canais[ch] if ch < len(nomes_canais) else ch
        rotulo = f"{os.path.basename(caminho)} | canal {nome_canal}"
        df = varre_canal(sinal, fs, window_s, step_s, n_surr, rng=rng,
                          rotulo_progresso=rotulo)
        df.insert(0, "arquivo", os.path.basename(caminho))
        df.insert(1, "canal", nome_canal)
        todos.append(df)

    return pd.concat(todos, ignore_index=True) if todos else pd.DataFrame()


# ==========================================
# MODO DEMO (validação em dados sintéticos, sem precisar de .ns2)
# ==========================================

def roda_demo():
    print("=== DEMO: validando o critério de flagging em dados sintéticos ===\n")
    print("(usando Theta com jitter de frequência/fase, como no LFP real -- um Theta\n"
          " perfeitamente periódico quebra o surrogate por deslocamento circular,\n"
          " porque o envelope de Gamma vira periódico com o mesmo período do Theta\n"
          " e o shift deixa de destruir o acoplamento. Isso não é hipotético: foi o\n"
          " que aconteceu na primeira versão deste demo, sem jitter.)\n")
    fs = 1000.0
    dur = 60.0
    t = np.arange(0, dur, 1 / fs)
    rng = np.random.default_rng(0)

    # Theta com jitter biológico: frequência instantânea variando (random walk lento)
    # entre ~5 e ~7 Hz, em vez de uma senoide cravada em 6.0 Hz.
    freq_inst = 6.0 + np.cumsum(rng.normal(0, 0.002, len(t)))
    freq_inst = np.clip(freq_inst, 5.0, 7.0)
    fase_theta_acum = 2 * np.pi * np.cumsum(freq_inst) / fs
    theta = np.cos(fase_theta_acum)
    modulador = (-theta + 1) / 2.0
    gamma_acoplado = modulador * np.sin(2 * np.pi * 60.0 * t)

    meio = len(t) // 2
    sinal = np.zeros_like(t)
    sinal[:meio] = theta[:meio] + rng.normal(0, 0.4, meio)  # sem gamma acoplado -> só ruído
    sinal[meio:] = theta[meio:] + gamma_acoplado[meio:] + rng.normal(0, 0.4, len(t) - meio)

    df = varre_canal(sinal, fs, window_s=10.0, step_s=5.0, n_surr=200, rng=rng)
    df["esperado"] = np.where(df["janela_ini_s"] < dur / 2, "ruído (sem acoplamento)",
                               "acoplamento real")

    pd.set_option("display.width", 120)
    print(df[["janela_ini_s", "janela_fim_s", "mi_observado", "z_score",
              "p_empirico", "esperado"]].to_string(index=False))

    print("\nCheck esperado: janelas 'ruído' devem ter z baixo (~0-2) e p alto;")
    print("janelas 'acoplamento real' devem ter z alto (>>3) e p baixo (~0).")
    print("Se isso bater, o critério de corte por z-score está funcionando")
    print("antes de aplicar nos arquivos .ns2 reais.")


# ==========================================
# MAIN
# ==========================================

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pasta", help="Pasta contendo os arquivos .ns2")
    ap.add_argument("--canal", type=int, default=None,
                     help="Índice de um único canal (default: todos)")
    ap.add_argument("--janela", type=float, default=10.0, help="Tamanho da janela (s)")
    ap.add_argument("--passo", type=float, default=5.0, help="Passo entre janelas (s)")
    ap.add_argument("--n_surr", type=int, default=200, help="Número de surrogates")
    ap.add_argument("--z_corte", type=float, default=3.0,
                     help="z-score mínimo para considerar candidato")
    ap.add_argument("--saida", default="resultados_triagem.csv")
    ap.add_argument("--demo", action="store_true", help="Roda autoteste sintético")
    args = ap.parse_args()

    if args.demo:
        roda_demo()
        return

    if not args.pasta:
        print("Erro: forneça --pasta com os arquivos .ns2, ou use --demo.")
        sys.exit(1)

    arquivos = sorted(glob.glob(os.path.join(args.pasta, "*.ns2")))
    if not arquivos:
        print(f"Nenhum .ns2 encontrado em {args.pasta}")
        sys.exit(1)

    canais = [args.canal] if args.canal is not None else None

    todos_resultados = []
    for i, arq in enumerate(arquivos):
        print(f"[{i+1}/{len(arquivos)}] Processando {os.path.basename(arq)} ...")
        try:
            df = varre_arquivo(arq, canais=canais, window_s=args.janela,
                                step_s=args.passo, n_surr=args.n_surr, seed=i)
            todos_resultados.append(df)
        except Exception as e:
            print(f"  -> ERRO ao processar {arq}: {e}")

    if not todos_resultados:
        print("Nenhum resultado gerado.")
        sys.exit(1)

    resultado_final = pd.concat(todos_resultados, ignore_index=True)
    resultado_final = resultado_final.sort_values("z_score", ascending=False)
    resultado_final.to_csv(args.saida, index=False)

    candidatos = resultado_final[resultado_final["z_score"] >= args.z_corte]
    print(f"\nSalvo: {args.saida}")
    print(f"Total de janelas analisadas: {len(resultado_final)}")
    print(f"Candidatos com z >= {args.z_corte}: {len(candidatos)}")

    if len(candidatos) > 0:
        print("\nTop candidatos (revisar proxy_artefato_motor antes de confiar):")
        print(candidatos.head(15)[[
            "arquivo", "canal", "janela_ini_s", "janela_fim_s",
            "z_score", "p_empirico", "proxy_artefato_motor"
        ]].to_string(index=False))
        print("\nATENÇÃO: se proxy_artefato_motor estiver alto nos mesmos candidatos"
              " de z-score alto, desconfie de artefato muscular antes de ir ao vídeo.")


if __name__ == "__main__":
    main()