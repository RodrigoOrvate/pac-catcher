"""
refina_candidatos.py
==========================================
ETAPA 2 do pipeline: refinamento estatístico dos candidatos gerados pelo
triagem_pac.py. A triagem inicial (z-score com 200 surrogates) é rápida
mas insuficiente para decidir "é acoplamento real ou não", por 4 motivos
que apareceram na análise do resultados.csv:

  1. 5472 janelas testadas sem correção de múltiplas comparações -> boa
     parte dos "candidatos" com z>=3 é esperada por puro acaso.
  2. Resolução do p-valor empírico travada em 1/200=0.005 (muitos
     candidatos empatados em p=0.0) -- fino demais pra ranquear ou
     aplicar correção FDR de forma confiável.
  3. O proxy de artefato motor original (150-450 Hz) não sobrepõe a
     banda de Gamma (30-80 Hz) analisada -- não detecta contaminação
     na banda que realmente importa.
  4. Não havia nenhuma checagem de CO-OCORRÊNCIA entre canais: um
     gerador focal real aparece em poucos canais vizinhos; se muitos
     canais disparam juntos na mesma janela, é assinatura clássica de
     ruído de modo comum (movimento, cabo, referência), não de dipolo
     neural -- mas também pode ser Theta genuinamente difuso (ex. REM).
     O script agora reporta isso em vez de esconder.

O QUE MUDA AQUI:
  - p-valor PARAMÉTRICO: ajusta uma distribuição Gama aos surrogates
    (KL-MI é não-negativo e assimétrico à direita -- Gama é apropriada,
    diferente da suposição implícita de normalidade do z-score) e tira
    o p-valor da cauda dessa distribuição ajustada. Isso dá resolução
    praticamente arbitrária sem precisar de milhões de permutações.
  - Correção de FDR (Benjamini-Hochberg) sobre TODAS as 5472 janelas
    originais (m = total de testes da triagem, não só os refinados) --
    é a escolha conservadora e correta para controlar falsos positivos
    na família inteira de testes.
  - Coluna de co-ocorrência: quantos canais do mesmo arquivo/janela
    também passaram do corte -- para você julgar "dipolo focal" vs.
    "ruído difuso" com um número na mão, não no olho.
  - Kurtose do sinal filtrado em Gamma (30-80Hz): rajadas musculares
    são transientes de alta kurtose mesmo depois do filtro; oscilação
    genuína tem kurtose baixa. Detecta contaminação NA banda certa,
    ao contrário do proxy antigo.
  - Checagem de saturação/clipping do sinal bruto.

Uso:
    python refina_candidatos.py --csv resultados.csv --pasta_ns2 /caminho/ns2 \
        --z_pre_filtro 2.0 --n_surr 1000 --fdr_q 0.05 \
        --saida resultados_refinados.csv

Requer: neo, numpy, scipy, pandas
"""

import argparse
import os
import numpy as np
import pandas as pd
import scipy.signal as signal
import scipy.stats as stats

from ns2_utils import le_ns2, fatia_janela


# ==========================================
# FUNÇÕES DE PROCESSAMENTO (mesmo núcleo do triagem_pac.py)
# ==========================================

def filtra_sinal(sinal_in, lowcut, highcut, fs, order=3):
    nyq = 0.5 * fs
    low = max(lowcut / nyq, 1e-6)
    high = min(highcut / nyq, 0.999)
    b, a = signal.butter(order, [low, high], btype="bandpass")
    return signal.filtfilt(b, a, sinal_in)


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


def calcula_mi(fase, envelope, n_bins=18):
    bins = np.linspace(-np.pi, np.pi, n_bins + 1)
    bin_idx = np.clip(np.digitize(fase, bins) - 1, 0, n_bins - 1)
    return _mi_de_bin_idx(bin_idx, envelope, n_bins)


def mi_com_surrogates_parametrico(fase, envelope, fs, n_surr=1000, n_bins=18,
                                   shift_min_s=1.0, rng=None):
    """
    Gera n_surr surrogates por deslocamento circular, ajusta uma distribuição
    Gama à distribuição nula resultante, e devolve um p-valor ANALÍTICO
    (cauda da Gama ajustada) em vez do p-valor empírico bruto (que trava em
    1/n_surr de resolução). Também devolve o z-score, por comparabilidade.
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

    # Ajuste paramétrico: Gama sobre os valores de MI dos surrogates.
    # Pequeno deslocamento (loc) fixo em 0 -- MI é estritamente >= 0.
    mi_surr_pos = mi_surr[mi_surr > 0]
    if len(mi_surr_pos) < 20 or mi_obs <= 0:
        p_analitico = np.mean(mi_surr >= mi_obs)  # cai pro empírico se não der pra ajustar
    else:
        try:
            forma, loc, escala = stats.gamma.fit(mi_surr_pos, floc=0)
            p_analitico = stats.gamma.sf(mi_obs, forma, loc=loc, scale=escala)
            p_analitico = float(np.clip(p_analitico, 0.0, 1.0))
        except Exception:
            p_analitico = np.mean(mi_surr >= mi_obs)

    return mi_obs, z, p_analitico, media, dp


def kurtose_gamma(sinal_in, fs, gamma_band=(30, 80)):
    """
    Kurtose (excesso) do sinal já filtrado na banda de Gamma. Rajadas
    musculares/movimento geram transientes que mantêm alta kurtose mesmo
    após o filtro; oscilação genuína tem kurtose baixa (perto de 0).
    Diferente do proxy antigo, mede DENTRO da banda que importa.
    """
    lfp_gamma = filtra_sinal(sinal_in, *gamma_band, fs)
    return float(stats.kurtosis(lfp_gamma, fisher=True))


def proxy_saturacao(sinal_bruto_int16, limiar_fracao=0.98):
    """Fração de amostras perto do limite do int16 (indicativo de clipping/saturação)."""
    limite = 32767 * limiar_fracao
    return float(np.mean(np.abs(sinal_bruto_int16.astype(float)) >= limite))


def bh_fdr(pvalues, m_total, alpha=0.05):
    """
    Correção de Benjamini-Hochberg. m_total é o número de testes da
    família COMPLETA (todas as 5472 janelas da triagem original), não
    só os que foram refinados aqui -- essa é a escolha conservadora:
    tratamos os que não entraram no refinamento como automaticamente
    não-significativos.
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

def adiciona_cooccorrencia(df, limiar_z=2.0):
    """
    Para cada (arquivo, janela), conta quantos canais também passaram do
    limiar na mesma janela de tempo. Números altos = suspeita de ruído
    de modo comum (ou Theta genuinamente difuso -- o número não decide
    sozinho, mas te dá o dado pra julgar em vez de descobrir só no vídeo).
    """
    contagem = (
        df[df["z_score"] >= limiar_z]
        .groupby(["arquivo", "janela_ini_s"])
        .size()
        .rename("n_canais_simultaneos")
    )
    return df.merge(contagem, on=["arquivo", "janela_ini_s"], how="left").fillna(
        {"n_canais_simultaneos": 0}
    )


# ==========================================
# MAIN
# ==========================================

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True, help="CSV gerado pelo triagem_pac.py")
    ap.add_argument("--pasta_ns2", required=True, help="Pasta com os .ns2 originais")
    ap.add_argument("--z_pre_filtro", type=float, default=2.0,
                     help="Corte (leniente) do z-score da triagem para entrar no refinamento")
    ap.add_argument("--n_surr", type=int, default=1000,
                     help="Nº de surrogates no refinamento (mais que a triagem)")
    ap.add_argument("--fdr_q", type=float, default=0.05, help="Nível de FDR (Benjamini-Hochberg)")
    ap.add_argument("--janela", type=float, default=10.0, help="Tamanho da janela usada na triagem (s)")
    ap.add_argument("--saida", default="resultados_refinados.csv")
    args = ap.parse_args()

    print(f"Lendo {args.csv} ...")
    df = pd.read_csv(args.csv)
    m_total = len(df)
    print(f"  -> {m_total} janelas na triagem original")

    df = adiciona_cooccorrencia(df, limiar_z=args.z_pre_filtro)

    candidatos = df[df["z_score"] >= args.z_pre_filtro].copy()
    print(f"  -> {len(candidatos)} candidatos com z >= {args.z_pre_filtro} vão para refinamento")

    if len(candidatos) == 0:
        print("Nenhum candidato para refinar. Encerrando.")
        return

    rng = np.random.default_rng(42)
    resultados = []

    for arquivo, grupo in candidatos.groupby("arquivo"):
        caminho = os.path.join(args.pasta_ns2, arquivo)
        print(f"\nCarregando {arquivo} para refinamento ({len(grupo)} candidatos) ...")
        dados, fs, nomes_canais = le_ns2(caminho)
        mapa_canal = {nome: i for i, nome in enumerate(nomes_canais)}

        for n_i, row in enumerate(grupo.itertuples(), 1):
            canal_idx = mapa_canal.get(row.canal)
            if canal_idx is None:
                print(f"  aviso: canal {row.canal} não encontrado em {arquivo}, pulando")
                continue

            trecho_bruto = fatia_janela(
                dados[:, canal_idx], fs, row.janela_ini_s, row.janela_fim_s
            )
            trecho = trecho_bruto.astype(float)

            lfp_theta = filtra_sinal(trecho, 4.0, 8.0, fs)
            lfp_gamma = filtra_sinal(trecho, 30.0, 80.0, fs)
            fase = np.angle(signal.hilbert(lfp_theta))
            env = np.abs(signal.hilbert(lfp_gamma))

            mi_obs, z_novo, p_analitico, mi_surr_m, mi_surr_dp = mi_com_surrogates_parametrico(
                fase, env, fs, n_surr=args.n_surr, rng=rng
            )
            kurt = kurtose_gamma(trecho, fs)
            sat = proxy_saturacao(trecho_bruto)

            resultados.append({
                "arquivo": arquivo,
                "canal": row.canal,
                "janela_ini_s": row.janela_ini_s,
                "janela_fim_s": row.janela_fim_s,
                "mi_observado": mi_obs,
                "z_score_refinado": z_novo,
                "p_analitico": p_analitico,
                "n_canais_simultaneos": row.n_canais_simultaneos,
                "kurtose_gamma": kurt,
                "proxy_saturacao": sat,
                "proxy_artefato_motor_150_450hz": row.proxy_artefato_motor,
            })

            if n_i % 20 == 0 or n_i == len(grupo):
                print(f"  refinado {n_i}/{len(grupo)}", end="\r")
        print()

    resultado_final = pd.DataFrame(resultados)

    # Correção FDR sobre a família completa de testes originais (m_total)
    resultado_final["significativo_fdr"] = bh_fdr(
        resultado_final["p_analitico"].values, m_total=m_total, alpha=args.fdr_q
    )

    # Veredito interpretável, combinando significância + sinais de artefato/co-ocorrência
    def veredito(row):
        if not row["significativo_fdr"]:
            return "Não significativo após correção FDR"
        alertas = []
        if row["kurtose_gamma"] > 3:  # oscilação senoidal limpa tem kurtose perto de 0
            alertas.append("possível transiente muscular (kurtose alta em Gamma)")
        if row["proxy_saturacao"] > 0.001:
            alertas.append("sinal saturado/clipado")
        if row["n_canais_simultaneos"] >= 8:
            alertas.append(f"{int(row['n_canais_simultaneos'])} canais simultâneos (possível ruído comum)")
        if alertas:
            return "Revisar: " + "; ".join(alertas)
        return "Candidato robusto"

    resultado_final["veredito"] = resultado_final.apply(veredito, axis=1)
    resultado_final = resultado_final.sort_values(
        ["significativo_fdr", "z_score_refinado"], ascending=[False, False]
    )
    resultado_final.to_csv(args.saida, index=False)

    n_sig = resultado_final["significativo_fdr"].sum()
    n_robusto = (resultado_final["veredito"] == "Candidato robusto").sum()
    print(f"\nSalvo: {args.saida}")
    print(f"Candidatos refinados: {len(resultado_final)}")
    print(f"Significativos após FDR (q={args.fdr_q}): {n_sig}")
    print(f"'Candidato robusto' (significativo + sem alertas): {n_robusto}")
    print("\nEsses são os que valem o tempo de conferir no vídeo primeiro.")


if __name__ == "__main__":
    main()