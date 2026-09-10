"""
audita_janela.py
==========================================
Orquestrador forense: para UM canal/janela, roda numa só passada as 4
auditorias que compartilham o mesmo formato de entrada (canal + tempos
de início/fim) -- skewness, footprint espacial, transientes/despike e
razão harmônica teta-gama -- e devolve/grava um relatório consolidado.

Papel diferente de `processa_sessao.py`: aquele é o orquestrador EM LOTE
(roda por canal, sobre um CSV de candidatos, produz os CSVs que
alimentam agrega_resultados.py). Este é o orquestrador POR JANELA:
"tenho um caso suspeito -- me diga tudo sobre ele numa chamada só".
Os dois coexistem sem sobreposição.

3 auditorias ficam DE FORA por terem pré-requisitos estruturais
incompatíveis com "canal + janela genérica" (não integradas aqui):
  - audita_held_out.py: exige janela reancorada (ilha vs. resto).
  - audita_segmentos.py: exige sub-segmentos explícitos definidos à mão,
    e recomputa o comodulograma completo (~275 células x 200 surrogates
    POR segmento) -- ordens de magnitude mais caro que os 4 daqui.
  - audita_harmonico_hfo.py: CSV/fonte de dados distintos (triagem HFO,
    aceita .mat), já integrado à rota de audita_harmonico.py separado.

Ganho de I/O: cada arquivo .ns2 é lido UMA vez (cache por arquivo), e
as 3 variantes de sinal que as 4 auditorias precisam (bruto sem notch,
com notch 60Hz, janela de contexto de 45s p/ FOOOF) são derivadas em
memória -- hoje, rodar os 4 scripts individualmente relê o mesmo
arquivo pelo menos 3-4 vezes.

Custo (em unidades de mi_z_par, n_surr=200): skewness ~0, harmônico ~0
em MI (1 fit FOOOF + PLV opcional), transientes 8x, footprint 32x (o
mais caro, ~4x tudo o resto somado) -- use --pula footprint para
descartar o mais caro.

Uso:
    # Janela única
    python audita_janela.py --pasta_ns2 "<sessao>/<basal>" \\
        --arquivo 20240708-123605-003.ns2 --canal chan22 \\
        --inicio 20 --fim 30 --fp 5 --fa 35

    # Modo lote, reaproveitando um vencedores.csv já existente
    python audita_janela.py --pasta_ns2 "<sessao>/<basal>" \\
        --csv "<sessao>/RESULTADOS/vencedores.csv"

    # Pulando a auditoria mais cara (footprint, 32x)
    python audita_janela.py --pasta_ns2 ... --csv ... --pula footprint
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ns2_utils import le_ns2, fatia_janela

from audita_skewness import skewness_de_sinal, classifica_skew
from audita_transientes import audita_transientes_de_sinal, contexto_amplitude
from audita_footprint import footprint_de_janela
from audita_harmonico import avalia_harmonico

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


# ==========================================
# NÚCLEO: um caso (canal/janela) já com dados em memória
# ==========================================

def audita_janela(dados, fs, nomes, canal, ini, fim, fp, fa,
                  janela_contexto_s=45.0, notch_hz=60.0,
                  n_surr=200, n_bins=18,
                  tol_rel=0.1, limiar_plv=0.8, limiar_skew=0.5,
                  modo_preprocesso="hibrido", f_linha=60.0,
                  limiar_skewness=0.5, pula_footprint=False):
    """
    Roda as 4 auditorias sobre um canal/janela já em memória (dados de
    TODOS os canais do arquivo, já carregados). Devolve um dict
    estruturado com o resultado bruto de cada auditoria + um veredito
    consolidado conservador.

    Nota sobre `rng`: propositalmente NÃO thread um Generator compartilhado
    entre as chamadas de audita_transientes_de_sinal/footprint_de_janela --
    ambas dependem, por dentro (via mi_z_par), do default embutido
    `default_rng(42)` recriado a cada chamada, que é o que garante os
    MESMOS deslocamentos de surrogate entre baseline/despike/sub-janelas
    (documentado em audita_transientes.py). Passar um rng compartilhado
    aqui quebraria essa garantia.
    """
    idx = nomes.index(canal) if canal in nomes else int(canal.replace("chan", "")) - 1

    # Para skewness/transientes/contexto: janela já em float (mesma convenção
    # que audita_transientes.py/audita_footprint.py já usam internamente).
    lfp_bruto = fatia_janela(dados[:, idx], fs, ini, fim).astype(float)

    centro = (ini + fim) / 2.0
    ctx_ini = max(0.0, centro - janela_contexto_s / 2.0)
    ctx_fim = centro + janela_contexto_s / 2.0

    # Para o harmônico: SEM astype(float) -- audita_harmonico.py original
    # passa sinal_ctx/sinal_cand com o dtype BRUTO de le_ns2 (int16) direto
    # para extrai_cf_teta_fooof/compute_plv_harmonico. Confirmado empiricamente
    # que o fit do FOOOF diverge na 9a casa decimal se o array já vier
    # convertido para float64 antes -- não é diferença de precisão trivial,
    # é uma dependência real de algum passo interno (Welch/detrend) ao dtype
    # de entrada. Manter o dtype bruto aqui é o que garante paridade bit-exata
    # com harmonico.csv gerado por audita_harmonico.py.
    sinal_cand_raw = fatia_janela(dados[:, idx], fs, ini, fim)
    sinal_ctx_raw = fatia_janela(dados[:, idx], fs, ctx_ini, ctx_fim)

    resultado = {"canal": canal, "janela_ini_s": ini, "janela_fim_s": fim,
                 "fase_pico_hz": fp, "amp_pico_hz": fa}

    # 1) skewness (sinal bruto, sem notch -- mesmo pré-processo do script individual)
    skew, n_skew = skewness_de_sinal(lfp_bruto, fs)
    veredito_skew = classifica_skew(skew, n_skew, limiar_skewness)
    resultado["skewness"] = {"skew": skew, "n": n_skew, "veredito": veredito_skew}

    # 2) transientes (aplica notch internamente)
    res_trans = audita_transientes_de_sinal(lfp_bruto, fs, fp, fa, notch_hz=notch_hz,
                                            n_bins=n_bins, n_surr=n_surr)
    resultado["transientes"] = res_trans

    # 2b) contexto de amplitude (usa o arquivo inteiro do canal + janela de todos os canais)
    coluna_canal_inteira = dados[:, idx]
    dados_janela_todos = fatia_janela(dados, fs, ini, fim)
    resultado["contexto_amplitude"] = contexto_amplitude(
        lfp_bruto, coluna_canal_inteira, fs, dados_janela_todos, idx, fa)

    # 3) footprint espacial (opcional -- o mais caro, 32x)
    if not pula_footprint:
        resultado["footprint"] = footprint_de_janela(
            dados, fs, nomes, ini, fim, fp, fa, notch_hz=notch_hz)
    else:
        resultado["footprint"] = None

    # 4) harmônico (sinal_cand = mesma janela, dtype bruto -- ver nota acima)
    resultado["harmonico"] = avalia_harmonico(
        sinal_ctx_raw, sinal_cand_raw, fs, fa, skew=skew,
        tol_rel=tol_rel, limiar_plv=limiar_plv, limiar_skew=limiar_skew,
        modo_preprocesso=modo_preprocesso, f_linha=f_linha)

    resultado["veredito_consolidado"], resultado["motivos"] = _consolida_veredito(resultado)
    return resultado


def _consolida_veredito(resultado):
    """
    Regra conservadora e transparente: começa em LIMPO, acumula um
    motivo por sinal de alerta. Reusa os limiares que cada auditoria
    individual já usa (--limiar do skewness, z>=3 do footprint,
    --limiar_plv/--tol_rel do harmônico). Os 2 limiares de transientes
    abaixo (razão de queda no despike, nº de sub-janelas com z>=3) são
    heurísticos de TRIAGEM, não critério científico -- mesmo espírito do
    aviso que audita_harmonico.py já traz no cabeçalho.
    """
    motivos = []

    if resultado["skewness"]["veredito"].startswith("SUSPECT"):
        motivos.append("skew_suspeito")

    t = resultado["transientes"]
    z_base = t["z_baseline"]
    z_k6 = t["despike"][6.0]["z"]
    razao_k6 = (z_k6 / z_base) if z_base not in (0, None) and not np.isnan(z_base) and z_base != 0 else np.nan
    if not np.isnan(razao_k6) and razao_k6 < 0.5:
        motivos.append(f"queda_forte_no_despike_k6 (razao={razao_k6:.2f})")
    n_sub_ge_3 = int(np.sum(np.array(t["z_sub"]) >= 3.0))
    if n_sub_ge_3 <= 1:
        motivos.append(f"acoplamento_nao_persiste_nas_subjanelas ({n_sub_ge_3}/5 com z>=3)")

    if resultado["footprint"] is not None:
        foot = resultado["footprint"]
        n_canais = len(foot["nomes"])
        if foot["n_z_ge_3"] > n_canais / 3:
            motivos.append(f"footprint_difuso ({foot['n_z_ge_3']}/{n_canais} canais com z>=3)")

    if resultado["harmonico"]["veredito_harmonico"] not in ("CLEAN", "SEM_REFERENCIA_TETA"):
        motivos.append(f"harmonico:{resultado['harmonico']['veredito_harmonico']}")

    if not motivos:
        return "LIMPO", []
    return "REVISAR: " + "; ".join(motivos), motivos


# ==========================================
# CLI: carrega arquivo(s), monta relatório(s)
# ==========================================

def _linha_csv(rotulo, arquivo, res):
    t = res["transientes"]
    return {
        "rotulo": rotulo, "arquivo": arquivo, "canal": res["canal"],
        "janela_ini_s": res["janela_ini_s"], "janela_fim_s": res["janela_fim_s"],
        "fase_pico_hz": res["fase_pico_hz"], "amp_pico_hz": res["amp_pico_hz"],
        "skew_skewness": res["skewness"]["skew"], "skew_n": res["skewness"]["n"],
        "skew_veredito": res["skewness"]["veredito"],
        "trans_z_baseline": t["z_baseline"], "trans_mi_baseline": t["mi_baseline"],
        "trans_z_despike_k6": t["despike"][6.0]["z"], "trans_frac_k6": t["despike"][6.0]["frac"],
        "trans_z_despike_k5": t["despike"][5.0]["z"], "trans_frac_k5": t["despike"][5.0]["frac"],
        **{f"trans_z_sub{i+1}": z for i, z in enumerate(t["z_sub"])},
        "ctx_rms_janela": res["contexto_amplitude"]["rms_janela"],
        "ctx_pct_rms_vs_arquivo": res["contexto_amplitude"]["pct_rms_vs_arquivo"],
        "ctx_rank_rms_entre_32_canais": res["contexto_amplitude"]["rank_rms_entre_32_canais"],
        "ctx_kurtose_bruto": res["contexto_amplitude"]["kurtose_bruto"],
        "ctx_kurtose_gama": res["contexto_amplitude"]["kurtose_gama"],
        "foot_n_z_ge_3": res["footprint"]["n_z_ge_3"] if res["footprint"] else None,
        "foot_n_canais": len(res["footprint"]["nomes"]) if res["footprint"] else None,
        "harm_cf_teta_fooof": res["harmonico"]["cf_teta_fooof"],
        "harm_erro_ajuste_fooof": res["harmonico"]["erro_ajuste_fooof"],
        "harm_plv": res["harmonico"]["plv_harmonico"],
        "harm_veredito": res["harmonico"]["veredito_harmonico"],
        "veredito_consolidado": res["veredito_consolidado"],
        "motivos": "; ".join(res["motivos"]),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pasta_ns2", required=True)
    ap.add_argument("--csv", default=None,
                    help="modo lote: vencedores.csv (colunas canal, arquivo, "
                         "janela_ini_s/inicio_s, janela_fim_s/fim_s, fase_pico_hz, amp_pico_hz)")
    ap.add_argument("--arquivo", default=None, help="modo direto: .ns2 dentro de --pasta_ns2")
    ap.add_argument("--canal", default=None)
    ap.add_argument("--inicio", type=float, default=None)
    ap.add_argument("--fim", type=float, default=None)
    ap.add_argument("--fp", type=float, default=None, help="fase de pico (Hz)")
    ap.add_argument("--fa", type=float, default=None, help="amplitude de pico (Hz)")
    ap.add_argument("--rotulo", default=None)

    ap.add_argument("--pula", nargs="+", default=[], choices=["footprint"],
                    help="pula auditorias caras (hoje só 'footprint', 32x o custo)")
    ap.add_argument("--n_surr", type=int, default=200)
    ap.add_argument("--n_bins", type=int, default=18)
    ap.add_argument("--janela_contexto_s", type=float, default=45.0)
    ap.add_argument("--tol_rel", type=float, default=0.1)
    ap.add_argument("--limiar_plv", type=float, default=0.8)
    ap.add_argument("--limiar_skew", type=float, default=0.5)
    ap.add_argument("--modo_preprocesso", default="hibrido",
                    choices=["gaussiana", "cirurgica", "hibrido"])
    ap.add_argument("--f_linha", type=float, default=60.0)
    ap.add_argument("--saida_dir", default="auditoria_janela")
    args = ap.parse_args()

    if not args.csv and not (args.arquivo and args.canal and args.inicio is not None
                              and args.fim is not None and args.fp is not None and args.fa is not None):
        raise SystemExit("[ERRO] Use --csv (modo lote) OU --arquivo/--canal/--inicio/--fim/--fp/--fa (modo direto).")

    if args.csv:
        df = pd.read_csv(args.csv)
        casos = []
        for _, r in df.iterrows():
            casos.append(dict(
                rotulo=r.get("rotulo"), arquivo=str(r["arquivo"]), canal=str(r["canal"]),
                ini=float(r.get("janela_ini_s", r.get("inicio_s", 0))),
                fim=float(r.get("janela_fim_s", r.get("fim_s", 0))),
                fp=float(r["fase_pico_hz"]), fa=float(r["amp_pico_hz"]),
            ))
    else:
        casos = [dict(rotulo=args.rotulo or f"{args.canal}@{args.inicio:.0f}-{args.fim:.0f}",
                      arquivo=args.arquivo, canal=args.canal, ini=args.inicio, fim=args.fim,
                      fp=args.fp, fa=args.fa)]

    pula_footprint = "footprint" in args.pula

    os.makedirs(args.saida_dir, exist_ok=True)
    cache = {}
    linhas = []
    linhas_footprint = []

    for caso in casos:
        path = os.path.join(args.pasta_ns2, caso["arquivo"])
        if path not in cache:
            dados, fs, nomes = le_ns2(path)
            cache[path] = (dados, fs, [str(n) for n in nomes])
        dados, fs, nomes = cache[path]

        print(f"\n=== {caso['rotulo']} — {caso['arquivo']} @ {caso['canal']} "
              f"{caso['ini']:.0f}-{caso['fim']:.0f}s — pico {caso['fp']:g}x{caso['fa']:g} Hz ===")

        res = audita_janela(
            dados, fs, nomes, caso["canal"], caso["ini"], caso["fim"], caso["fp"], caso["fa"],
            janela_contexto_s=args.janela_contexto_s, n_surr=args.n_surr, n_bins=args.n_bins,
            tol_rel=args.tol_rel, limiar_plv=args.limiar_plv,
            limiar_skew=args.limiar_skew, modo_preprocesso=args.modo_preprocesso,
            f_linha=args.f_linha, limiar_skewness=args.limiar_skew,
            pula_footprint=pula_footprint,
        )

        print(f"  skew={res['skewness']['skew']:+.3f} ({res['skewness']['veredito']})")
        t = res["transientes"]
        print(f"  transientes: baseline z={t['z_baseline']:.2f} | despike_k6 z={t['despike'][6.0]['z']:.2f} "
              f"| sub-janelas z=" + ", ".join(f"{z:.1f}" for z in t["z_sub"]))
        if res["footprint"] is not None:
            print(f"  footprint: {res['footprint']['n_z_ge_3']}/{len(res['footprint']['nomes'])} canais com z>=3")
            for nome, z in res["footprint"]["top6"]:
                linhas_footprint.append({"caso": caso["rotulo"], "canal": nome, "z": round(float(z), 2)})
        h = res["harmonico"]
        print(f"  harmonico: cf_teta={h['cf_teta_fooof']} veredito={h['veredito_harmonico']}")
        print(f"  >>> VEREDITO CONSOLIDADO: {res['veredito_consolidado']}")

        linhas.append(_linha_csv(caso["rotulo"], caso["arquivo"], res))

    saida_csv = os.path.join(args.saida_dir, "auditoria_janela.csv")
    pd.DataFrame(linhas).to_csv(saida_csv, index=False)
    print(f"\nSalvo: {saida_csv}")

    if linhas_footprint:
        saida_foot = os.path.join(args.saida_dir, "auditoria_janela_footprint.csv")
        pd.DataFrame(linhas_footprint).to_csv(saida_foot, index=False)
        print(f"Salvo: {saida_foot}")


if __name__ == "__main__":
    main()
