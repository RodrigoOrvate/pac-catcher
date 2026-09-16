"""Fonte unica de verdade para os caminhos do workspace (a pasta
acoplamento_theta-gamma, que contem SCRIPT/, LAC_NOCI/, RESULTADOS_MESTRADO/
como irmas -- ver .claude/CLAUDE.md). Nunca faca hardcode de caminho de disco
num script -- importe daqui.

O workspace padrao e' a pasta-mae de SCRIPT/ (deduzida da localizacao deste
arquivo), entao mover o workspace inteiro de disco -- como em 13/09/2026, de
C: para D: -- nao exige editar nada. Para usar dados em outro lugar, defina a
variavel de ambiente ACOPLAMENTO_BASE."""
import os
import re

BASE_SCRIPT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE_WORKSPACE = os.environ.get("ACOPLAMENTO_BASE", os.path.dirname(BASE_SCRIPT))
BASE_LAC_NOCI = os.path.join(BASE_WORKSPACE, "LAC_NOCI")
BASE_RESULTADOS_MESTRADO = os.path.join(BASE_WORKSPACE, "RESULTADOS_MESTRADO")
BASE_RESULTADOS = os.path.join(BASE_SCRIPT, "resultados")
BASE_FIGURAS = os.path.join(BASE_SCRIPT, "figuras")  # figuras de diagnóstico (fora do git)

_indices_ns2 = {}


def _indexa_ns2(base):
    if base not in _indices_ns2:
        indice = {}
        for raiz, _, arquivos in os.walk(base):
            for nome in arquivos:
                if nome.lower().endswith(".ns2"):
                    indice.setdefault(nome, []).append(os.path.join(raiz, nome))
        _indices_ns2[base] = indice
    return _indices_ns2[base]


def localiza_ns2(arquivo, dica=None, base=None):
    """Caminho completo do .ns2 bruto `arquivo` dentro de LAC_NOCI (ou `base`).

    Acha pelo nome do arquivo (carimbo de data/hora da gravação), sem depender
    de como a pasta da sessão se chama -- serve pra basal ("Basal antes da
    infusao") e pras condições de infusão ("0h pos infusao - 30min", ...).

    O nome NÃO é único no disco: em 2026-09-14 a pós-infusão do LAC tinha as
    mesmas gravações copiadas em MTESC03_LAC e MTESC05_LAC. Quando há mais de
    uma cópia, desempata pelo rato (`MTESC\\d+`) encontrado em `dica` (use a
    coluna `sessao` da linha). Se ainda assim for ambíguo, levanta ValueError
    com os candidatos -- nunca escolhe um em silêncio.

    Retorna None se o arquivo não existir.
    """
    candidatos = _indexa_ns2(base or BASE_LAC_NOCI).get(os.path.basename(str(arquivo)), [])
    if len(candidatos) <= 1:
        return candidatos[0] if candidatos else None
    rato = re.search(r"MTESC\d+", dica or "")
    if rato:
        do_rato = [c for c in candidatos if rato.group(0) in c]
        if len(do_rato) == 1:
            return do_rato[0]
        if do_rato:
            candidatos = do_rato
    raise ValueError(
        f"{arquivo} existe em {len(candidatos)} lugares e a dica {dica!r} não desempata: "
        + "; ".join(candidatos))


def figuras(*subpastas):
    """Pasta de saída de figuras de diagnóstico (SCRIPT/figuras/<subpastas>),
    criada se não existir."""
    caminho = os.path.join(BASE_FIGURAS, *subpastas)
    os.makedirs(caminho, exist_ok=True)
    return caminho
