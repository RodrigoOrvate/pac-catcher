"""
triagem_pac_mat.py
==========================================
Wrapper para rodar o pipeline de triagem PAC multi-banda diretamente
sobre arquivos .mat (sem precisar do adaptador CSV intermediário).

Carrega o LFP do .mat, aplica verificação de escala automática,
e roda varre_canal() com os três pares Theta-Gamma / Theta-HG / Theta-HFO
em uma única passagem sobre o sinal.

DIFERENÇAS do fluxo antigo (adapta_lfp_mat → tmp.csv → triagem_pac.py):
  1. Sem CSV intermediário — o sinal vai direto para o pipeline.
  2. Verificação de escala antes do processamento (evita overflow de z-score).
  3. teta_ok calculado por janela (não cacheado).
  4. Todos os três pares numa única passagem — o rato pode estar em exploração
     (teta_gamma alto) ou repouso (theta_hfo alto); o pipeline registra os três.
  5. ratio_hfo_gamma incluído — proxy de harmônico Gamma→HFO por janela.

Saída: FLAGS_<nome>_TRIPLO.csv com colunas:
  janela_ini_s, janela_fim_s,
  teta_ok,
  mi_theta_gamma, z_theta_gamma, p_theta_gamma,
  mi_theta_hg,    z_theta_hg,    p_theta_hg,
  mi_theta_hfo,   z_theta_hfo,   p_theta_hfo,
  ratio_hfo_gamma,
  transiente_detectado, frac_transiente, max_diff_z, max_amp_z,
  mvl, correlacao_ruido, suspeito_banda_larga, proxy_artefato_motor

Uso:
    # Todos os três pares (default)
    python triagem_pac_mat.py --mat DADOS_EXEMPLO/LFP_HG_HFO.mat

    # Só theta_gamma e theta_hfo
    python triagem_pac_mat.py --mat DADOS_EXEMPLO/LFP_HG_HFO.mat \\
        --pares theta_gamma theta_hfo

    # Forçar normalização da escala antes de processar
    python triagem_pac_mat.py --mat DADOS_EXEMPLO/LFP_HG_HFO.mat --normaliza

    # Especificar qual variável do .mat usar
    python triagem_pac_mat.py --mat DADOS_EXEMPLO/LFP_HG_HFO.mat \\
        --canal_mat lfpHG --pares theta_hg

Requer: numpy, scipy, pandas + fooof (para audita_harmonico_hfo.py)
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

# Adiciona o pipeline e a raiz do SCRIPT ao path (triagem_pac + pac_core)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from triagem_pac import BAND_PAIRS, varre_canal

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


# ==========================================
# LEITURA E DIAGNÓSTICO DO .mat
# ==========================================

def carrega_mat_com_diagnostico(caminho_mat, canal_mat=None, normaliza=False,
                                  normaliza_limiar_std_alto=5000,
                                  normaliza_limiar_std_baixo=0.001):
    """
    Carrega LFP de um .mat e imprime diagnóstico de escala.

    Se normaliza=True: aplica normalização robusta median/MAD.
    Se normaliza=False mas escala fora do range: imprime AVISO e sugere --normaliza.

    Retorna: (sinal_1d_float64, fs, nome_chave_usada)
    """
    from pac_core.io import le_mat as carrega_lfp_mat, info_sinal, normaliza_sinal, fs_do_mat

    print(f"Carregando: {caminho_mat}")
    arr, fs, chave = carrega_lfp_mat(caminho_mat, canal_mat)
    print(f"  Variável: '{chave}'  |  fs: {fs:.0f} Hz  |  "
          f"duração: {len(arr)/fs:.1f} s  |  dtype: {arr.dtype}")

    info_sinal(arr, nome=chave, fs=fs)

    std_sinal = float(arr.std())
    if normaliza:
        print("\n  [--normaliza] Aplicando normalização robusta (median/MAD)...")
        arr = normaliza_sinal(arr)
        print(f"  Após normalização: std={arr.std():.4f}")
    elif std_sinal > normaliza_limiar_std_alto:
        print(f"\n  [AVISO] std={std_sinal:.1f} muito alto. "
              f"Use --normaliza para evitar overflow de z-score.")
        print("  Continuando SEM normalização (z-scores podem ser absurdos).")
    elif std_sinal < normaliza_limiar_std_baixo:
        print(f"\n  [AVISO] std={std_sinal:.6f} muito baixo. "
              f"Use --normaliza para escala compatível com o pipeline.")

    return arr, fs, chave


# ==========================================
# RELATÓRIO DE FLAGS COMPACTO
# ==========================================

def _imprime_resumo_flags(df, pares_usados, z_corte):
    """Imprime resumo dos flags por par e detecta padrões anômalos."""
    print(f"\n{'='*60}")
    print(f"RESUMO — {len(df)} janelas")
    print(f"{'='*60}")

    for par in pares_usados:
        z_col = f"z_{par}"
        if z_col not in df.columns:
            continue
        z_vals = df[z_col].dropna()
        n_cand = int((z_vals >= z_corte).sum())
        # Alerta de variância zero
        if z_vals.std() < 0.01:
            aviso = "⚠ VARIÂNCIA ZERO — provável cache/bug de escala"
        elif n_cand == len(df):
            aviso = "⚠ 100% candidatos — limiar provavelmente permissivo demais"
        elif n_cand == 0:
            aviso = "⚠ 0 candidatos — z sempre baixo, checar escala do sinal"
        else:
            aviso = "OK"
        print(f"  {par:<15}  z med={z_vals.median():.2f}  "
              f"std={z_vals.std():.2f}  "
              f"candidatos(z>={z_corte})={n_cand}/{len(df)}  {aviso}")

    if "teta_ok" in df.columns:
        n_teta = int(df["teta_ok"].sum())
        pct = 100 * n_teta / len(df) if len(df) > 0 else 0
        if df["teta_ok"].std() < 0.01:
            aviso_teta = "⚠ teta_ok CONSTANTE — verifique cálculo por janela"
        else:
            aviso_teta = "OK (variação esperada)"
        print(f"  teta_ok         {n_teta}/{len(df)} janelas ({pct:.0f}%)  {aviso_teta}")

    if "ratio_hfo_gamma" in df.columns:
        r = df["ratio_hfo_gamma"]
        n_suspeito = int((r > 0.3).sum())
        print(f"  ratio_hfo_gamma med={r.median():.3f}  "
              f"suspeito_harmonico(>0.3): {n_suspeito} janelas")
        if n_suspeito > 0:
            print("    → Execute audita_harmonico_hfo.py para confirmar harmônico Gamma→HFO")

    print(f"{'='*60}\n")


# ==========================================
# MAIN
# ==========================================

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mat", required=True,
                     help="Arquivo .mat com o LFP")
    ap.add_argument("--canal_mat", default=None,
                     help="Variável do .mat a usar (default: auto-detecta na ordem "
                          "lfp, LFP, lfpBruto, lfpHG, lfpHFO)")
    ap.add_argument("--fs", type=float, default=None,
                     help="Taxa de amostragem (Hz). Sobrescreve a detectada no .mat")
    ap.add_argument("--pares", nargs="+",
                     default=["theta_gamma", "theta_hg", "theta_hfo"],
                     choices=list(BAND_PAIRS.keys()),
                     help="Pares a calcular (default: todos os três)")
    ap.add_argument("--janela", type=float, default=10.0, help="Janela (s, default 10)")
    ap.add_argument("--passo",  type=float, default=5.0,  help="Passo  (s, default 5)")
    ap.add_argument("--n_surr", type=int,   default=200,  help="Surrogates (default 200)")
    ap.add_argument("--z_corte", type=float, default=3.0, help="z mínimo para candidato")
    ap.add_argument("--normaliza", action="store_true",
                     help="Normaliza sinal (median/MAD) antes de processar")
    ap.add_argument("--saida", default=None,
                     help="Arquivo CSV de saída (default: FLAGS_<nome>_TRIPLO.csv)")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    if not os.path.isfile(args.mat):
        print(f"[ERRO] Arquivo não encontrado: {args.mat}")
        sys.exit(1)

    # Carrega e diagnostica
    sinal, fs, chave = carrega_mat_com_diagnostico(
        args.mat, canal_mat=args.canal_mat, normaliza=args.normaliza)

    if args.fs is not None:
        print(f"  [--fs] Sobrescrevendo fs para {args.fs} Hz")
        fs = args.fs

    # Verifica suporte de banda vs. fs
    nyq = fs * 0.5
    pares_sel = {}
    for par in args.pares:
        cfg = BAND_PAIRS[par]
        amp_hi = cfg["amp"][1]
        if cfg["amp"][0] >= nyq * 0.98:
            print(f"  [AVISO] Par '{par}' incompatível com fs={fs:.0f} Hz "
                  f"(banda {cfg['amp']} Hz > Nyquist). Ignorado.")
        else:
            pares_sel[par] = cfg

    if not pares_sel:
        print("[ERRO] Nenhum par compatível com a fs do arquivo.")
        sys.exit(1)

    print(f"\nPares a calcular: {list(pares_sel.keys())}")
    print(f"Janela {args.janela}s / passo {args.passo}s / {args.n_surr} surrogates")

    # Estima nº de janelas
    win  = int(args.janela * fs)
    step = int(args.passo  * fs)
    n_janelas_est = max(0, (len(sinal) - win) // step + 1)
    print(f"Janelas estimadas: {n_janelas_est}")
    print()

    rng = np.random.default_rng(args.seed)
    df = varre_canal(
        sinal, fs,
        window_s=args.janela,
        step_s=args.passo,
        n_surr=args.n_surr,
        pares=pares_sel,
        rng=rng,
        rotulo_progresso=os.path.basename(args.mat),
    )

    # Adiciona metadados
    nome_base = os.path.splitext(os.path.basename(args.mat))[0]
    df.insert(0, "arquivo", os.path.basename(args.mat))
    df.insert(1, "canal",   chave)

    # Relatório de flags
    _imprime_resumo_flags(df, list(pares_sel.keys()), args.z_corte)

    # Gera FLAGS resumido (formato FLAGS_*.csv existente + novas colunas)
    saida = args.saida or f"FLAGS_{nome_base}_TRIPLO.csv"
    df.to_csv(saida, index=False)
    print(f"Salvo: {saida}  ({len(df)} janelas × {len(df.columns)} colunas)")

    # Dica de próximo passo
    z_col_hfo = "z_theta_hfo"
    if z_col_hfo in df.columns:
        n_hfo = int((df[z_col_hfo] >= args.z_corte).sum())
        if n_hfo > 0:
            print(f"\nPróximo passo para {n_hfo} janelas com z_theta_hfo >= {args.z_corte}:")
            print(f"  python pipeline/auditorias/audita_harmonico_hfo.py \\")
            print(f"      --csv {saida} \\")
            print(f"      --mat {args.mat} \\")
            print(f"      --z_corte {args.z_corte}")


if __name__ == "__main__":
    main()
