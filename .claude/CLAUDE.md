# CLAUDE.md — SCRIPT central do projeto (acoplamento theta-gamma)

## O que é esta pasta

**Código único do grupo de estudos** `C:\acoplamento_theta-gamma\`,
compartilhado por TODOS os estudos e sessões (MTESC04/MTESC05 × NOCI/LAC).
Nada de específico de sessão mora aqui — nem dados, nem saídas, nem listas de
vencedores. Rodar os comandos DE DENTRO desta pasta, apontando para a sessão
por caminho relativo (ver `README_pipeline.md`).

```
C:\acoplamento_theta-gamma\
├── SCRIPT\                          ← ESTA pasta (código + README + CLAUDE.md)
├── MTESC04_NOCI\
│   ├── MTESC04 -- 1 - infusao - 08-07-2024\
│   │   ├── Basal antes da infusao\  ← 3 .ns2 (gravação 20240708-123605)
│   │   └── RESULTADOS\              ← SAÍDAS: resultados*.csv, comodulogramas*/,
│   │                                    figuras/, auditoria/, vencedores.csv,
│   │                                    registro_resultados.md, .claude/
│   │                                    (as pastas 1–2 chamavam-se "SCRIPT",
│   │                                     renomeadas para RESULTADOS em 24/08/2026)
│   ├── MTESC04 -- 2 - infusao - 09-07-2024\   (idem; gravação 20240709-141215)
│   └── MTESC04 -- 3..6 - infusao - 11/15/22/23-07-2024\   (pendentes)
├── MTESC04_LAC\
├── MTESC05_NOCI\
└── MTESC05_LAC\
```

Regra de ouro: **sessões têm dados e saídas; o SCRIPT central tem código.**
Se um script precisar de informação específica da sessão (vencedores,
caminhos), ela entra por CLI/CSV — nunca editar listas no código.

## Estrutura canônica (padronizada em 30/08/2026)

### SCRIPT/ (código — único, nada de sessão)
```
SCRIPT/
├── triagem_pac.py        → results: resultados.csv
├── refina_candidatos.py  → resultados_refinados.csv
├── comodulogram.py       → comodulogramas_fdr/
├── figura_apresentacao.py→ figuras/<rotulo>_<canal>.png
├── robustez_parametros.py→ robustez_parametros.csv
├── extrair_picos.py      → extracao_picos_v1_vs_v2.csv (re-verifica pico com banda restrita + FDR)
├── ns2_utils.py
├── gerar_relatorio_pdf.py← CSV-driven (lê vencedores_consolidado.csv)
└── auditorias/           ← validação (um por função):
    ├── audita_transientes.py   audita_segmentos.py   audita_footprint.py
    ├── diagnostico_janela.py   audita_skewness.py    audita_held_out.py
```

### RESULTADOS/ por sessão (canônico)
```
RESULTADOS/  (antigas pastas "SCRIPT" das sessões 1–2, renomeadas 24/08/2026)
├── figuras/                  ← só figuras dos vencedores VALIDADOS
├── comodulogramas_fdr/  diagnosticos/   auditoria/
├── resultados.csv            (triagem)
├── resultados_refinados.csv  (refinamento)
├── robustez_parametros.csv
├── vencedores.csv            (corrente, SEM sufixo — sempre a versão final validada)
├── rejeitados.csv            (rejeitados com motivo — pseudoreplicação, FDR, fora da banda...)
└── registro_resultados.md
└── auditoria/                ← artefatos de validação (não confundir com resultados):
    ├── skewness.csv  held_out.csv  footprint.csv  extracao_picos_v1_vs_v2.csv
    ├── candidatos_auditoria.csv     (lista de candidatos testados na re-auditoria)
    ├── *.log                        (logs de execução das etapas)
    └── figuras/mapas de candidatos REJEITADOS   (evidência da auditoria)
```
Regra: a **raiz de RESULTADOS só tem saídas de pipeline** (resultados*, robustez,
vencedores, rejeitados, registro). Tudo que é validação/auditoria/log/figura de
candidato rejeitado vai para `auditoria/`. O **consolidado por estudo também só
carrega figuras de vencedores validados** — as de rejeitados ficam na auditoria
da sessão (ex.: MTESC05 consolidado = só `S5_chan32_sniffing.png`, 1 evento).

### RESULTADOS_CONSOLIDADOS/ (por estudo — MTESC04_NOCI/, MTESC05_NOCI/)
```
RESULTADOS_CONSOLIDADOS/
├── relatorio_final_<ESTUDO>_NOCI.pdf
├── vencedores_consolidado.csv      ← alimenta gerar_relatorio_pdf.py
├── vencedores_consolidado_NOCI.md  ← revisão crítica + linguagem de manuscrito
├── limitacoes_PAC_NOCI.md
├── registro_resultados_S{1..6}.md  (um por sessão)
└── S{1..6}_chanXX_comportamento.png   (figura canônica: <S#>_<canal>_<comportamento>)
```

### Regras de nomes (não gerar lixo)
- **Nunca** salvar variantes de pipeline com sufixo: `_atualizado`, `_picos_reais`,
  `_v2`, `resultados_triagem`, `candidatos_refinados`. O nome canônico é o da
  estrutura acima; versões intermediárias só se forem evidência nomeada
  (`extracao_picos_v1_vs_v2.csv`), nunca lixo de pipeline.
- `vencedores.csv` da sessão = versão **final validada** (sem sufixo). O
  histórico/cadeia de versões é preservado em `BACKUP_PRE_PADRONIZACAO_*` e no
  git, não em arquivos `_vN` soltos.
- Consolidado: um `vencedores_consolidado.csv` por estudo é a fonte da verdade
  para o PDF.

## Convenções do estudo

- **Pipeline**: 8 passos (+2.5 skewness e 7.5 held-out) documentados no
  `README_pipeline.md` (triagem → refinamento → comodulogramas → conferência
  comportamental manual → respiração → FDR mapa → robustez+figuras → registro).
- **Validação = 5 etapas**: FDR de janela; não-60 Hz (notch sempre);
  respiração; FDR do mapa "concentrado em ΘΓ"; robustez (n_bins/larguras/
  MVL) + coerência comportamental no vídeo.
- **`vencedores.csv`** na pasta RESULTADOS de cada sessão (formato no README,
  passo 7): alimenta `robustez_parametros.py` e `figura_apresentacao.py`.
  Colunas-base: `rotulo,arquivo,canal,inicio_s,fim_s,fase_pico_hz,amp_pico_hz,
  comportamento` (MTESC05 acrescenta `pico,z_score`).
- **`vencedores_consolidado.csv`** no RESULTADOS_CONSOLIDADOS (por estudo) =
  lista final validada (sessao, canal, comportamento, z, pico, classificacao,
  veredito, observacoes, figura); alimenta `gerar_relatorio_pdf.py` (CSV-driven).
- Nome de PNG na sessão = `<rotulo>_<canal>.png`; no consolidado =
  `<S#>_chanXX_comportamento.png`.
- **`registro_resultados.md`** por sessão = registro formal consolidável
  (identificação + offset vídeo↔ns2 + números das etapas + checklist +
  rejeitados com motivo).
- **`.claude/CLAUDE.md`** por sessão = contexto vivo detalhado daquela
  sessão (lições, auditorias, decisões).
- Tempo: usuário fala em MM:SS GLOBAL da gravação (cada arquivo ≈ 5 min);
  converter para janela LOCAL do arquivo antes de rodar. Offset vídeo↔ns2 é
  POR SESSÃO (08/07: +321 s; 09/07: ≈0; 11/07: +100 s).
- Notch 60 Hz em tudo (rede brasileira). Nula de surrogates: deslocamento
  circular ≥1 s, 200 repetições, semente 42 nos scripts de validação.

## Estado das sessões (24/08/2026)

| # | Data | Gravação | Estado | Vencedores |
|---|---|---|---|---|
| 1 | 08/07 | 123605 | **VALIDADA** (+ auditoria de transientes, chan16 z=18,3, re-ancoragem do grooming) | 4 canais / 3 episódios |
| 2 | 09/07 | 141215 | Etapas automáticas completas; faltam sub-segmentação fina e pontuação de sniffing | 4 episódios (V1 destaque z=11,17) |
| 3 | 11/07 | 121046 | **VALIDADA** (+ auditoria transientes+pegada espacial; ilhas 2–4 s; offset +100 s; dados chegaram duplicados da sessão 2 — substituídos) | 2 canais / 15 episódios |
| 4 | 15/07 | — | **EXCLUÍDA** (erro de sincronia vídeo↔ns2) | — |
| 5–6 | 22/23-07 | — | pendentes | — |

Detalhes e lições por sessão: ler o `.claude/CLAUDE.md` e
`registro_resultados.md` DA SESSÃO (não duplicar aqui).

## Lições já aprendidas (válidas para todas as sessões)

1. 60 Hz contamina gamma — comparar sempre com/sem notch.
2. Janela mista (dois estados comportamentais) infla/esconde z — ancorar a
   janela num estado puro antes de publicar número.
3. Proxy respiratório (0,5–3 Hz) é cego para sniffing 4–8 Hz — desempate por
   vídeo + teste de segmento.
4. Transiente espalha no mapa FDR; acoplamento genuíno concentra em ΘΓ.
5. A varredura própria pode perder acoplamentos fortes mascarados pelo 60 Hz
   (chan16 da sessão 1) — pegada espacial de 32 canais é rede de segurança.
6. MVL fraco com MI alto acontece (distribuição vs 1º momento) — registrar ⚠,
   não rejeitar sozinho.
7. Conferir o ID da gravação contra a data da sessão ANTES de rodar em
   "dados novos" (sessão 3 chegou com cópia exata da sessão 2 — MD5 idêntico).
8. Acoplamentos podem viver em ILHAS de 2–4 s dentro da janela de 10 s:
   âncora fina com `audita_segmentos.py`; despike derruba até PAC genuíno
   bursty — a PEGADA ESPACIAL (gradiente entre vizinhos) é o árbitro final.
9. **Pseudoreplicação espacial**: canais vizinhos com pico idêntico no mesmo
   arquivo/janela = UM evento co-detectado (condução de volume/referência
   comum), não N casos independentes. Pegada difusa (`audita_footprint.py`,
   muitos canais z≥3) = difuso/co-detecção; focal (poucos vizinhos) = fonte
   local. Ex.: "7 vencedores" do MTESC05 vira **1 evento** (chan32 8×50).
10. **Pico deve ser validado na banda + FDR, não argmax cru**: pico fora da
    banda canônica (4–8 × 30–80 Hz) ou não-significativo sob FDR na célula
    vencedora é desqualificado (ex.: chan16 8×115 e estático 5×20, MTESC05).
    Rodar `extrair_picos.py` antes de publicar o par.
