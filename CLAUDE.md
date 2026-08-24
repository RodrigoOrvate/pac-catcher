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
│   │   └── SCRIPT\                  ← SAÍDAS: resultados*.csv, comodulogramas*/,
│   │                                    figuras/, auditoria/, vencedores.csv,
│   │                                    registro_resultados.md, .claude/
│   ├── MTESC04 -- 2 - infusao - 09-07-2024\   (idem; gravação 20240709-141215)
│   └── MTESC04 -- 3..6 - infusao - 11/15/22/23-07-2024\   (pendentes)
├── MTESC04_LAC\
├── MTESC05_NOCI\
└── MTESC05_LAC\
```

Regra de ouro: **sessões têm dados e saídas; o SCRIPT central tem código.**
Se um script precisar de informação específica da sessão (vencedores,
caminhos), ela entra por CLI/CSV — nunca editar listas no código.

## Convenções do estudo

- **Pipeline**: 8 passos documentados no `README_pipeline.md` (triagem →
  refinamento → comodulogramas → conferência comportamental manual →
  respiração → FDR mapa → robustez+figuras → registro).
- **Validação = 5 etapas**: FDR de janela; não-60 Hz (notch sempre);
  respiração; FDR do mapa "concentrado em ΘΓ"; robustez (n_bins/larguras/
  MVL) + coerência comportamental no vídeo.
- **`vencedores.csv`** na pasta SCRIPT de cada sessão (formato no README,
  passo 7): alimenta `robustez_parametros.py` e `figura_apresentacao.py`.
  Nome de PNG gerado = `<rotulo>_<canal>.png`.
- **`registro_resultados.md`** por sessão = registro formal consolidável
  (identificação + offset vídeo↔ns2 + números das etapas + checklist +
  rejeitados com motivo).
- **`.claude/CLAUDE.md`** por sessão = contexto vivo detalhado daquela
  sessão (lições, auditorias, decisões).
- Tempo: usuário fala em MM:SS GLOBAL da gravação (cada arquivo ≈ 5 min);
  converter para janela LOCAL do arquivo antes de rodar. Offset vídeo↔ns2 é
  POR SESSÃO (08/07: +321 s; 09/07: ≈0).
- Notch 60 Hz em tudo (rede brasileira). Nula de surrogates: deslocamento
  circular ≥1 s, 200 repetições, semente 42 nos scripts de validação.

## Estado das sessões (24/08/2026)

| # | Data | Gravação | Estado | Vencedores |
|---|---|---|---|---|
| 1 | 08/07 | 123605 | **VALIDADA** (+ auditoria de transientes, chan16 z=18,3, re-ancoragem do grooming) | 4 canais / 3 episódios |
| 2 | 09/07 | 141215 | Etapas automáticas completas; faltam sub-segmentação fina e pontuação de sniffing | 4 episódios (V1 destaque z=11,17) |
| 3–6 | 11/15/22/23-07 | — | pendentes | — |

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
