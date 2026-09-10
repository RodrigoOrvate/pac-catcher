"""
refina_candidatos.py
==========================================
ETAPA 2 do pipeline: refinamento estatístico dos candidatos gerados pelo
triagem_pac.py para UM ou MAIS pares fase-amplitude.

SUPORTA OS TRÊS PARES (detectados automaticamente pelo CSV de entrada):
  theta_gamma  : Theta 4-8 Hz (fase) x Gamma 30-80 Hz (amplitude)
  theta_hg     : Theta 4-8 Hz (fase) x High-Gamma 80-150 Hz (amplitude)
  theta_hfo    : Theta 4-8 Hz (fase) x HFO 150-250 Hz (amplitude)

O QUE FAZ:
  - p-valor PARAMÉTRICO: ajusta distribuição Gama aos surrogates do par
    correto (KL-MI é não-negativo e assimétrico à direita — Gama apropriada).
  - Correção FDR (Benjamini-Hochberg) sobre TODAS as janelas originais.
  - Co-ocorrência entre canais por par: indica ruído de modo comum vs.
    gerador focal.
  - Kurtose filtrada NA BANDA DE AMPLITUDE do par (não Gamma fixo 30-80 Hz):
    detecta transientes musculares na banda que realmente importa.
  - ratio_hfo_gamma (quando disponível): herdado da triagem, indica suspeita
    de harmônico Gamma→HFO antes mesmo da auditoria FOOOF.

Detecção automática de pares:
  Se o CSV contiver colunas z_theta_gamma, z_theta_hg, z_theta_hfo,
  processa cada par que tiver candidatos com z >= z_pre_filtro.
  Se o CSV for legado (coluna z_score), trata como theta_gamma.

Uso:
    python refina_candidatos.py --csv resultados_triplo.csv \
        --pasta_ns2 /caminho/ns2 \
        --z_pre_filtro 2.0 --n_surr 1000 --fdr_q 0.05 \
        --saida resultados_refinados.csv

Requer: neo, numpy, scipy, pandas
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import scipy.signal as signal
import scipy.stats as stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ns2_utils import le_ns2, fatia_janela
from triagem_pac import BAND_PAIRS, detecta_transiente, correlacao_gama_ruido
from pac_core.filtering import filtra_sinal

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


# ==========================================
# FUNÇÕES DE PROCESSAMENTO
# filtra_sinal agora vem de pac_core.filtering (import no topo)
# ==========================================


def _mi_de_bin_idx(bin_idx, envelope, n_bins):
    soma_bins = np.bincount(bin_idx, weights=envelope, minlength=n_bins)
    cont_bins = np.bincount(bin_idx, minlength=n_bins)
    media_bins = np.divide(soma_bins, cont_bins, out=np.zeros(n_bins), where=cont_bins > 0)
    soma = np.sum(media_bins)
    if soma <= 0:
        return 0.0
    P = media_bins / soma
    H = -np.sum(P * np.log(P + 1e-10))
    return (np.log(n_bins) - H) / np.log(n_bins)


def mi_com_surrogates_parametrico(fase, envelope, fs, n_surr=1000, n_bins=18,
                                   shift_min_s=1.0, rng=None):
    """
    Gera n_surr surrogates por deslocamento circular, ajusta distribuição
    Gama à distribuição nula e devolve p-valor ANALÍTICO (maior resolução
    que empírico 1/n_surr) e z-score para comparabilidade.
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

    mi_surr_pos = mi_surr[mi_surr > 0]
    if len(mi_surr_pos) < 20 or mi_obs <= 0:
        p_analitico = float(np.mean(mi_surr >= mi_obs))
    else:
        try:
            forma, loc, escala = stats.gamma.fit(mi_surr_pos, floc=0)
            p_analitico = float(np.clip(
                stats.gamma.sf(mi_obs, forma, loc=loc, scale=escala), 0.0, 1.0))
        except Exception:
            p_analitico = float(np.mean(mi_surr >= mi_obs))

    return mi_obs, z, p_analitico, media, dp


def kurtose_banda(sinal_in, fs, banda):
    """
    Kurtose (excesso) do sinal filtrado NA BANDA DE AMPLITUDE do par.
    Rajadas musculares têm kurtose alta mesmo após o filtro; oscilação
    genuína tem kurtose perto de 0. Agora genérica para qualquer banda.
    """
    nyq = fs * 0.5
    lo, hi = banda
    hi = min(hi, nyq * 0.98)
    if lo >= hi:
        return 0.0
    lfp_filtrado = filtra_sinal(sinal_in, lo, hi, fs)
    return float(stats.kurtosis(lfp_filtrado, fisher=True))


def proxy_saturacao(sinal_bruto, limiar_fracao=0.98):
    """Fração de amostras perto do limite do int16 (indicativo de clipping)."""
    limite = 32767 * limiar_fracao
    return float(np.mean(np.abs(sinal_bruto.astype(float)) >= limite))


def bh_fdr(pvalues, m_total, alpha=0.05):
    """
    Correção Benjamini-Hochberg. m_total = família completa (todas as
    janelas da triagem original) — escolha conservadora e correta.
    """
    p = np.asarray(pvalues)
    n = len(p)
    ordem = np.argsort(p)
    p_ordenado = p[ordem]
    limiar = (np.arange(1, n + 1) / m_total) * alpha
    passou = p_ordenado <= limiar
    significativo = np.zeros(n, dtype=bool)
    if np.any(passou):
        corte = np.max(np.where(passou))
        significativo[ordem[:corte + 1]] = True
    return significativo


# ==========================================
# CO-OCORRÊNCIA ENTRE CANAIS
# ==========================================

def adiciona_cooccorrencia(df, z_col, limiar_z=2.0):
    """
    Para cada (arquivo, janela), conta quantos canais também passaram do
    limiar no mesmo par e janela. Alto = suspeita de ruído de modo comum.
    """
    if z_col not in df.columns:
        df["n_canais_simultaneos"] = 0
        return df
    contagem = (
        df[df[z_col] >= limiar_z]
        .groupby(["arquivo", "janela_ini_s"])
        .size()
        .rename("n_canais_simultaneos")
    )
    return df.merge(contagem, on=["arquivo", "janela_ini_s"], how="left").fillna(
        {"n_canais_simultaneos": 0}
    )


# ==========================================
# DETECÇÃO AUTOMÁTICA DE PARES
# ==========================================

def detecta_pares_no_csv(df):
    """
    Identifica quais pares estão presentes no CSV de triagem.
    Retorna lista de tuplas (nome_par, z_col, banda_fase, banda_amp).
    Suporta CSV novo (z_theta_gamma etc.) e CSV legado (z_score = theta_gamma).
    """
    pares_encontrados = []
    for nome_par, cfg in BAND_PAIRS.items():
        z_col = f"z_{nome_par}"
        if z_col in df.columns:
            pares_encontrados.append((
                nome_par, z_col,
                cfg["fase"], cfg["amp"]
            ))
    # Suporte legado: CSV antigo com coluna "z_score" (= theta_gamma)
    if not pares_encontrados and "z_score" in df.columns:
        cfg = BAND_PAIRS["theta_gamma"]
        pares_encontrados.append((
            "theta_gamma", "z_score",
            cfg["fase"], cfg["amp"]
        ))
        print("  [info] CSV legado detectado (z_score) — tratando como theta_gamma.")
    return pares_encontrados


# ==========================================
# VEREDITO
# ==========================================

def gera_veredito(row, banda_nome):
    """Veredito interpretável combinando significância + alertas de artefato."""
    if not row["significativo_fdr"]:
        return "Nao significativo apos FDR"
    alertas = []
    if row.get("kurtose_banda", 0) > 3:
        alertas.append(f"transiente na banda {banda_nome} (kurtose alta)")
    if row.get("proxy_saturacao", 0) > 0.001:
        alertas.append("sinal saturado")
    if row.get("n_canais_simultaneos", 0) >= 8:
        alertas.append(f"{int(row['n_canais_simultaneos'])} canais simultaneos (possivel ruido comum)")
    if row.get("ratio_hfo_gamma", 0) > 0.3 and "hfo" in banda_nome:
        alertas.append("ratio_hfo_gamma alto (possivel harmonico Gamma->HFO)")
    if row.get("suspeito_banda_larga", False):
        alertas.append("suspeito banda larga / envelope identico ao ruido")
    if alertas:
        return "Revisar: " + "; ".join(alertas)
    return "Candidato robusto"


# ==========================================
# MAIN
# ==========================================

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True,
                    help="CSV gerado pelo triagem_pac.py (suporta multi-par e legado)")
    ap.add_argument("--pasta_ns2", required=True, help="Pasta com os .ns2 originais")
    ap.add_argument("--z_pre_filtro", type=float, default=2.0,
                    help="Corte leniente da triagem para entrar no refinamento (default 2.0)")
    ap.add_argument("--n_surr", type=int, default=1000,
                    help="Surrogates no refinamento (default 1000)")
    ap.add_argument("--fdr_q", type=float, default=0.05,
                    help="Nivel FDR Benjamini-Hochberg (default 0.05)")
    ap.add_argument("--janela", type=float, default=10.0,
                    help="Tamanho da janela usada na triagem (s)")
    ap.add_argument("--saida", default="resultados_refinados.csv")
    args = ap.parse_args()

    print(f"Lendo {args.csv} ...")
    df_orig = pd.read_csv(args.csv)
    m_total = len(df_orig)
    print(f"  -> {m_total} janelas na triagem original")

    pares = detecta_pares_no_csv(df_orig)
    if not pares:
        print("[ERRO] Nenhum par detectado no CSV.")
        print(f"  Colunas encontradas: {list(df_orig.columns)}")
        sys.exit(1)

    print(f"  Pares detectados: {[p[0] for p in pares]}")

    rng = np.random.default_rng(42)
    todos_resultados = []

    for nome_par, z_col, fase_band, amp_band in pares:
        print(f"\n{'='*55}")
        print(f"PAR: {nome_par}  [{fase_band[0]}-{fase_band[1]} Hz x {amp_band[0]}-{amp_band[1]} Hz]")
        print(f"{'='*55}")

        df_par = df_orig.copy()
        df_par = adiciona_cooccorrencia(df_par, z_col, limiar_z=args.z_pre_filtro)
        candidatos = df_par[df_par[z_col] >= args.z_pre_filtro].copy()

        if "arquivo" not in candidatos.columns:
            print(f"  [AVISO] Coluna 'arquivo' nao encontrada — sem .ns2 para refinar {nome_par}.")
            continue

        print(f"  {len(candidatos)} candidatos com {z_col} >= {args.z_pre_filtro}")
        if len(candidatos) == 0:
            continue

        resultados_par = []

        for arquivo, grupo in candidatos.groupby("arquivo"):
            caminho = os.path.join(args.pasta_ns2, arquivo)
            if not os.path.isfile(caminho):
                print(f"  [AVISO] {caminho} nao encontrado — pulando.")
                continue

            print(f"\n  Carregando {arquivo} ({len(grupo)} candidatos) ...")
            try:
                dados, fs, nomes_canais = le_ns2(caminho)
            except Exception as e:
                print(f"  [ERRO] {e} — pulando {arquivo}")
                continue

            mapa_canal = {str(nome): i for i, nome in enumerate(nomes_canais)}

            # Verifica compatibilidade da banda amp com Nyquist
            nyq = fs * 0.5
            amp_hi = min(amp_band[1], nyq * 0.98)
            if amp_band[0] >= amp_hi:
                print(f"  [AVISO] Banda {amp_band} Hz incompativel com fs={fs:.0f} Hz — pulando par.")
                continue
            banda_efetiva = (amp_band[0], amp_hi)

            for n_i, row in enumerate(grupo.itertuples(), 1):
                canal_idx = mapa_canal.get(str(row.canal))
                if canal_idx is None:
                    print(f"  aviso: canal {row.canal} nao encontrado em {arquivo}, pulando")
                    continue

                trecho_bruto = fatia_janela(
                    dados[:, canal_idx], fs, row.janela_ini_s, row.janela_fim_s
                )
                trecho = trecho_bruto.astype(float)

                # Fase: sempre theta 4-8 Hz
                lfp_fase = filtra_sinal(trecho, fase_band[0], fase_band[1], fs)
                fase = np.angle(signal.hilbert(lfp_fase))

                # Amplitude: banda específica do par
                lfp_amp = filtra_sinal(trecho, banda_efetiva[0], banda_efetiva[1], fs)
                env = np.abs(signal.hilbert(lfp_amp))

                mi_obs, z_novo, p_analitico, mi_surr_m, mi_surr_dp = \
                    mi_com_surrogates_parametrico(
                        fase, env, fs, n_surr=args.n_surr, rng=rng
                    )

                # Kurtose na banda de amplitude do par (não Gamma fixo)
                kurt = kurtose_banda(trecho, fs, banda_efetiva)
                sat = proxy_saturacao(trecho_bruto)

                # Camadas anti-falso-positivo
                trans_info = detecta_transiente(trecho, fs)
                banda_info = correlacao_gama_ruido(trecho, fs, theta_band=fase_band, gamma_band=banda_efetiva)

                linha = {
                    "par": nome_par,
                    "arquivo": arquivo,
                    "canal": row.canal,
                    "janela_ini_s": row.janela_ini_s,
                    "janela_fim_s": row.janela_fim_s,
                    "z_triagem": float(getattr(row, z_col, float('nan'))),
                    "mi_observado": round(mi_obs, 6),
                    "z_score_refinado": round(z_novo, 4),
                    "p_analitico": p_analitico,
                    "n_canais_simultaneos": int(row.n_canais_simultaneos)
                        if hasattr(row, 'n_canais_simultaneos') else 0,
                    "kurtose_banda": round(kurt, 3),
                    "proxy_saturacao": round(sat, 6),
                    # Camada 1: Transiente
                    "transiente_detectado": trans_info["transiente_encontrado"],
                    "frac_transiente": round(trans_info["frac_transiente"], 4),
                    # Camada 2: MVL
                    "mvl": round(banda_info["mvl"], 4),
                    # Camada 4: Correlacao gamma<->ruido
                    "correlacao_ruido": round(banda_info["correlacao_ruido"], 3),
                    "suspeito_banda_larga": banda_info["suspeito_banda_larga"],
                }

                # ratio_hfo_gamma se disponivel (herda da triagem)
                if hasattr(row, 'ratio_hfo_gamma'):
                    linha["ratio_hfo_gamma"] = float(row.ratio_hfo_gamma)

                # teta_ok se disponivel
                if hasattr(row, 'teta_ok'):
                    linha["teta_ok"] = int(row.teta_ok)

                resultados_par.append(linha)

                if n_i % 20 == 0 or n_i == len(grupo):
                    print(f"  refinado {n_i}/{len(grupo)}", end="\r")
            print()

        if not resultados_par:
            print(f"  Nenhum resultado para {nome_par}.")
            continue

        df_par_result = pd.DataFrame(resultados_par)

        # FDR sobre a família completa m_total
        df_par_result["significativo_fdr"] = bh_fdr(
            df_par_result["p_analitico"].values,
            m_total=m_total,
            alpha=args.fdr_q
        )

        df_par_result["veredito"] = df_par_result.apply(
            lambda r: gera_veredito(r, nome_par), axis=1
        )

        n_sig = int(df_par_result["significativo_fdr"].sum())
        n_robusto = int((df_par_result["veredito"] == "Candidato robusto").sum())
        print(f"  {nome_par}: {len(df_par_result)} refinados | {n_sig} significativos FDR "
              f"| {n_robusto} candidatos robustos")

        todos_resultados.append(df_par_result)

    if not todos_resultados:
        print("\nNenhum resultado gerado.")
        sys.exit(0)

    resultado_final = pd.concat(todos_resultados, ignore_index=True)
    resultado_final = resultado_final.sort_values(
        ["par", "significativo_fdr", "z_score_refinado"],
        ascending=[True, False, False]
    )
    resultado_final.to_csv(args.saida, index=False)

    print(f"\n{'='*55}")
    print(f"Salvo: {args.saida}")
    print(f"Total refinados: {len(resultado_final)} ({resultado_final['par'].nunique()} pares)")
    print(f"Significativos FDR (q={args.fdr_q}): {int(resultado_final['significativo_fdr'].sum())}")
    n_rob = int((resultado_final['veredito'] == 'Candidato robusto').sum())
    print(f"Candidatos robustos (sig + sem alertas): {n_rob}")
    print("\nDistribuicao por par:")
    resumo = resultado_final.groupby('par').agg(
        n_total=('z_score_refinado', 'count'),
        n_sig=('significativo_fdr', 'sum'),
        n_robusto=('veredito', lambda x: (x == 'Candidato robusto').sum()),
        z_mediana=('z_score_refinado', 'median'),
    ).round(2)
    print(resumo.to_string())


if __name__ == "__main__":
    main()