# SCRIPT TABLE OF CONTENTS — Acoplamento theta-gamma

Este é o repositório central de **código** do grupo de estudos
`C:\acoplamento_theta-gamma\`. Ele é compartilhado por todos os estudos e
sessões (MTESC04/05 × NOCI/LAC). Nada de específico de sessão mora aqui —
nem dados, nem saídas, nem listas de vencedores (entram por CLI/CSV).

```
acoplamento_theta-gamma/
├── SCRIPT/          ← ESTA pasta: código único + documentação (você está aqui)
├── MTESC04_NOCI/
├── MTESC04_LAC/
├── MTESC05_NOCI/
└── MTESC05_LAC/     ← cada estudo tem sessões com dados, vídeo e RESULTADOS/
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
├── README.md                  ← este guia (índice geral)
├── scripts_explicados.md      ← explicação script a script do pipeline
├── requirements.txt
├── pipeline/                  ← código do pipeline PAC (deteção/validação)
│   ├── README_pipeline.md     ← como rodar uma sessão nova
│   ├── triagem_pac.py          triagem → resultados.csv
│   ├── refina_candidatos.py    refinamento → resultados_refinados.csv
│   ├── comodulogram.py         mapas FDR
│   ├── figura_apresentacao.py  figuras dos vencedores
│   ├── robustez_parametros.py  robustez + MVL
│   ├── extrair_picos.py        re-verifica pico (banda restrita + FDR)
│   ├── ns2_utils.py            leitura compartilhada de .ns2
│   ├── gerar_relatorio_pdf.py  PDF CSV-driven
│   └── auditorias/             validação pós-hoc:
│       ├── audita_transientes.py  audita_segmentos.py
│       ├── audita_footprint.py    audita_skewness.py
│       ├── audita_held_out.py     diagnostico_janela.py
│       ├── audita_harmonico.py    teste de razão harmônica (FOOOF Kuhn)
│       ├── linha_noise_kuhn.py    3 configs de limpeza de linha (60Hz)
│       └── compara_preprocesso_linha.py comparador gaussiana/cirurgica/hibrido
└── preditor/                  ← PROJETO NOVO: prever PAC → disparar TTL
    ├── README_preditor.md     ← contexto, estado e comandos (LER!)
    ├── analisar_pre_evento.py  extrai os 10 s antes de cada vencedor
    ├── treinar_preditor.py     RandomForest → modelo_pac.pkl
    ├── validar_preditor.py     valida com controles reais (LOOCV)
    ├── prever_pac_tempo_real.py janela deslizante → dispara_ttl()
    └── modelo_pac.pkl
```

## Os dois arquivos "núcleo" de documentação

- **`pipeline/README_pipeline.md`** — guia operacional: como rodar a triagem de
  PAC numa sessão nova, as etapas, os critérios de validação (5 etapas) e os
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
  (LFP_FOOOF): `.mat` de exemplo (CA1/DG, tetrode 4 ch, 1000 Hz),
  `rem_noise.m` (implementação MATLAB da limpeza de linha),
  `Figure_1c.m`, `Example_fitting.ipynb`, `fit.py`, `funcs.py`.
  Usado por `compara_preprocesso_linha.py` para comparação com
  dados reais. **Não é código de produção** — os `.mat` carregam o
  campo `lfp` (não `mlfp`, que tem NaN). Manter para referência
  e reprodução; não editar.

## Projetos dentro do SCRIPT

1. **Pipeline PAC** (maduro, validado) — ideal para seguir o guia do
   `pipeline/README_pipeline.md`.
2. **Preditor PAC + optogenética** (protótipo, `preditor/`) — prever o
   acoplamento ~10 s antes e disparar um pulso TTL. **Ainda não validado**;
   ver o aviso no `preditor/README_preditor.md`.
