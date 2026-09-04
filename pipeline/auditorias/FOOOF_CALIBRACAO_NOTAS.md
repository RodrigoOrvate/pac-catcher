---
name: fooof-calibracao
description: Notas de calibração do FOOOF em audita_harmonico.py — descobertas e pendências em 04/09/2026
---

# Notas de calibração FOOOF (audita_harmonico.py)

## TL;DR (atualizado em 04/09/2026)

- `aperiodic_mode='knee'` é o **default de produção** e deve permanecer.
- A configuração de produção atual **NÃO detecta teta** em sinal sintético
  realista (45s, 1/f knee=28Hz, SNR=20dB). Sintético bem comportado falha.
- O teste sintético usa kwargs relaxados (`nperseg=1.0s`, `pwl=(1,4)`)
  APENAS para isolar a validação da lógica razão+PLV. Esta diferença é
  DELIBERADA, DOCUMENTADA, e **NÃO representa a produção**.
- A integração do FOOOF de produção com o portão `erro_ajuste<0.15`
  ainda não foi validada ponta a ponta.

## Por que `knee` e não `fixed`

LFP real de CA1/DG tem `knee frequency` real:
- ~28 Hz em CA1
- ~70 Hz em DG

(fonte: Kuhn et al. 2026, LFP_FOOOF).

`fixed` assume 1/f puro sem quebra. Em LFP real, `fixed` superestima a
potência de banda larga e mascara teta/gamma. Logo, `fixed` só serve em
sinais sem componente 1/f ou em testes sintéticos sem estrutura.

## O que o teste sintético exige (FOOOF isolado)

Sweep de parâmetros sobre sinal sintético 3x, 1/f realista, 45s, SNR=20dB:

| nperseg_s | freq_res | pwl          | resultado               |
|-----------|----------|--------------|-------------------------|
| 0.5s      | 2.0Hz    | qualquer     | FALHA                   |
| 1.0s      | 1.0Hz    | (1,4)        | cf=7.99, erro=0.30      |
| 1.0s      | 1.0Hz    | (2,5)        | cf=7.99, erro=0.30      |
| **1.2s (produção)** | 0.83Hz | (2,5) | **FALHA — nem detecta** |
| 2.0s      | 0.5Hz    | (1,4)        | cf=7.99, erro=0.20      |

**Conclusão**: A combinação `nperseg=1.2s + pwl=(2,5)` é instável em
4-12 Hz. Com 9-10 bins no intervalo, a Gaussiana do teta (largura 2-5 Hz
= 2.5-6 bins) gera ajuste numérico instável.

## Portão de qualidade (`erro_ajuste < 0.15`)

Chute não-calibrado. Nos testes sintéticos bem comportados:
- `nperseg=1.0s`: erro = 0.30 (portão REJEITA)
- `nperseg=2.0s`: erro = 0.20-0.54 (portão REJEITA)
- `nperseg=1.2s` produção: não detecta, retorna None

Conclusão: o portão de 0.15 **rejeita o sintético bem comportado** em
qualquer configuração. Isso significa que:
- O portão pode estar calibrado em LFP real (que tem estruturas
  diferentes do sintético)
- OU o portão está muito restritivo

**Não dá para concluir sem validação em sessões reais**.

## Como isso apareceu

Sessão de validação em 04/09/2026. Reporte inicial marcou a lógica
razão+PLV como "robusta" — esta formulação foi corrigida depois que o
usuário apontou que:

1. O teste de quase-coincidência γ=25Hz (3×8=24) só validava a borda
   da tolerância, **não o valor do PLV** (porque a razão já rejeitava).
2. O cenário crítico que **prova o valor do PLV** é o oposto: γ dentro
   da tolerância MAS com fase LIVRE. Este caso foi adicionado.
3. O "veredito robusto" misturava validação de uma peça com validação
   do pipeline inteiro — frase corrigida para "lógica validada
   isoladamente; integração com portão FOOOF de produção ainda em
   teste".

## Estado atual (frase honesta para apresentação)

> "Lógica de discriminação (razão+PLV) validada isoladamente em sinal
> sintético realista: PLV rejeita oscilador com fase LIVRE mesmo
> dentro da tolerância de frequência (γ=23.7Hz vs 3×8=24Hz, |Δ|=0.3Hz
> dentro de 0.8Hz tol, PLV=0.27 < 0.5). Integração com o portão de
> qualidade do FOOOF de produção ainda em teste — a configuração
> atual (nperseg=1.2s, pwl=(2,5)) não detecta teta no sintético
> bem comportado. Validação biológica em sessões MTESC04/05 pendente."

## Pendências

- [ ] Calibração empírica do portão `erro_ajuste<0.15` em LFP real
- [ ] Investigar se `nperseg=1.2s + pwl=(2,5)` está correto para LFP
      real (pode estar otimizado para teta mais largo de CA1, mas
      pode falhar em DG com teta mais fino)
- [ ] Validação biológica em sessões MTESC04/05 (não sintética)
- [ ] Alternativa: usar `specparam` em vez de `fooof` (fooof está
      deprecado)
