# SCRIPT TABLE OF CONTENTS — Acoplamento theta-gamma

Este é o repositório central de **código** do grupo de estudos
`C:\acoplamento_theta-gamma\`. Compartilhado por todos os estudos e
sessões (MTESC04/05 × NOCI/LAC). Nada de específico de sessão mora aqui —
nem dados, nem saídas, nem listas de vencedores (entram por CLI/CSV).

```
acoplamento_theta-gamma/
├── SCRIPT/          ← ESTA pasta: código único + documentação (você está aqui)
├── MTESC04_NOCI/
├── MTESC05_LAC/
├── MTESC05_NOCI/
└── EXPLORACAO_OBJETOS/
```

## Leitura rápida — qual documento ler?

| Situação | Leia |
|---|---|
| **Primeira vez aqui / orientar-se** | Este `README.md` |
| Rodar o pipeline PAC numa sessão nova | [`pipeline/README_pipeline.md`](pipeline/README_pipeline.md) |
| Entender cada script em detalhe | [`scripts_explicados.md`](scripts_explicados.md) |
| Projeto de prever PAC p/ optogenética (TTL) | [`preditor/README_preditor.md`](preditor/README_preditor.md) |
| Estado/lições por sessão | `.claude/CLAUDE.md` e `RESULTADOS/registro_resultados.md` **da sessão** |

## Estrutura de pastas

```
SCRIPT/
├── README.md                    ← este guia (índice geral)
├── scripts_explicados.md        ← explicação script a script do pipeline
├── requirements.txt
├── pipeline/                    ← código do pipeline PAC
│   ├── README_pipeline.md       ← como rodar uma sessão nova (LER PRIMEIRO)
│   ├── triagem_pac.py            Etapa 1: varre janelas, 3 pares Θ×Γ/HG/HFO
│   ├── triagem_pac_mat.py        Wrapper para .mat (sem .ns2) — multi-banda
│   ├── refina_candidatos.py      Etapa 2: FDR de janela + filtros de artefato
│   ├── comodulogram.py           Etapa 3: mapas z-scoredos ±notch + FDR do mapa
│   ├── figura_apresentacao.py    Etapa 7: figuras dos vencedores
│   ├── robustez_parametros.py    Etapa 7: sweep de parâmetros + MVL
│   ├── extrair_picos.py          re-verifica pico (banda restrita + FDR)
│   ├── diagnostico_janela.py     Etapa 5: 5 painéis + MI resp vs Θ×Γ
│   ├── exploracao_minuto.py      Passo 0: painel visual por minuto
│   ├── ns2_utils.py              leitura compartilhada de .ns2
│   ├── adapta_lfp_mat.py         carrega .mat com diagnóstico de escala
│   ├── gerar_relatorio_pdf.py    PDF CSV-driven
│   └── auditorias/               validação pós-hoc:
│       ├── audita_transientes.py    despike + sub-janelas
│       ├── audita_segmentos.py      localização temporal do acoplamento
│       ├── audita_footprint.py      pegada espacial (32 canais)
│       ├── audita_skewness.py       assimetria do theta
│       ├── audita_held_out.py       validação hold-out
│       ├── audita_harmonico.py      razão harmônica Θ→Γ (FOOOF Kuhn)
│       ├── audita_harmonico_hfo.py  razão harmônica Γ→HFO (FOOOF)  ← NOVO
│       ├── diagnostico_janela.py
│       ├── linha_noise_kuhn.py      limpeza de linha 60 Hz
│       └── compara_preprocesso_linha.py
├── FOOOF/                       ← referência Kuhn et al. 2026 (não é produção)
├── preditor/                    ← PROJETO: prever PAC → disparar TTL
│   ├── README_preditor.md       ← contexto, estado e comandos (LER!)
│   ├── analisar_pre_evento.py
│   ├── treinar_preditor.py
│   ├── validar_preditor.py
│   ├── prever_pac_tempo_real.py
│   └── modelo_pac.pkl
└── DADOS_EXEMPLO_LFP_HG_HFO/   ← dados de teste do pipeline multi-banda
    ├── LFP_HG_HFO.mat           LFP real 300s @1kHz (lfpHG + lfpHFO)
    └── FLAGS_TRIPLO_NOVO.csv    resultado da triagem: 3 pares, 59 janelas
```

## Três acoplamentos suportados

O pipeline detecta os três acoplamentos fase-amplitude relevantes em CA1:

| Par | Banda de fase | Banda de amplitude | Estado comportamental |
|---|---|---|---|
| `theta_gamma` | Theta 4–8 Hz | Gamma 30–80 Hz | Exploração locomotora |
| `theta_hg` | Theta 4–8 Hz | High-Gamma 80–150 Hz | Misto |
| `theta_hfo` | Theta 4–8 Hz | HFO 150–250 Hz | Repouso / SWR-associado |

```bash
# Varredura completa em .mat (todos os 3 pares, ~8s para 300s de sinal):
python pipeline/triagem_pac_mat.py --mat DADOS_EXEMPLO_LFP_HG_HFO/LFP_HG_HFO.mat

# Varredura em .ns2 com todos os 3 pares:
python pipeline/triagem_pac.py --pasta <sessao>/BASAL \
    --pares theta_gamma theta_hg theta_hfo \
    --saida <sessao>/RESULTADOS/resultados_triplo.csv

# Auditoria harmônico Gamma→HFO (pós-triagem):
python pipeline/auditorias/audita_harmonico_hfo.py \
    --csv FLAGS_TRIPLO_NOVO.csv --mat LFP_HG_HFO.mat --z_corte 3.0
```

## Os dois arquivos "núcleo" de documentação

- **`pipeline/README_pipeline.md`** — guia operacional: como rodar a triagem em
  uma sessão nova, as etapas, os critérios de validação (5 etapas) e os
  parâmetros que mudam por sessão (caminhos, offset vídeo↔ns2, rótulos).
- **`scripts_explicados.md`** — referência por script: o que cada um faz, por
  que daquela forma (z-score e não MI bruto, 200 surrogates, notch 60 Hz,
  double-dipping, integer ratio etc.).

## Regra de ouro (vale para tudo)

**Sessões têm dados e saídas; o `SCRIPT` central tem código.**
Se um script precisar de informação específica da sessão (vencedores,
caminhos, offsets), ela entra por CLI/CSV — **nunca editar listas no código**.

## Pastas de referência (não são código de pipeline)

- **`FOOOF/`** — dados e código de referência de Kuhn et al. 2026
  (`fit.py`, `funcs.py`, `CA1_example.mat`, `DG_example.mat`).
  Usado por `audita_harmonico.py` e `compara_preprocesso_linha.py`.
  **Não é código de produção** — não editar.

- **`preditor/`** — protótipo para prever PAC ~10 s antes e disparar TTL.
  **Ainda não validado**; ver aviso em `preditor/README_preditor.md`.

- **`DADOS_EXEMPLO_LFP_HG_HFO/`** — dados de teste do pipeline multi-banda.
  Contém `LFP_HG_HFO.mat` (LFP real, 300s, 1 kHz) e o CSV de resultado
  `FLAGS_TRIPLO_NOVO.csv` (59 janelas × 3 pares). Os arquivos anteriores
  (`FLAGS_LFP_HG_HFO.csv`, `FLAGS_LFP_HG_HFO_CORRIGIDO.csv`,
  `RESULTADOS_PIPELINE_LFP_HG_HFO.csv`) foram mantidos aqui como histórico —
  representam resultados com bugs corrigidos em 2026-09-05.
