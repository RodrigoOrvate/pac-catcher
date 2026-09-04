---
name: fooof-calibracao
description: Notas de calibração do FOOOF em audita_harmonico.py — descobertas e pendências em 04/09/2026
---

# Notas de calibração FOOOF (audita_harmonico.py)

## TL;DR (atualizado em 04/09/2026)

- `aperiodic_mode='knee'` é o **default de produção** e deve permanecer.
- **Resolvido em 04/09/2026**: adição de `nfft=4000` (zero-padding) na
  chamada `welch` restaura a detecção de teta em sinal sintético
  realista com a configuração de produção (`nperseg=1.2s`, `pwl=(2,5)`).
  A causa era uma diferença de implementação: Kuhn et al. 2026 usa
  `nfft=4000` explicitamente, e o `welch` do scipy sem `nfft` produzia
  grade de frequência com resolução 0.83 Hz/bin (instável para Gaussiana
  de 2-5 Hz). Com `nfft=4000`, a grade cai para 0.25 Hz/bin e o
  ajuste converge.
- O teste sintético usa kwargs relaxados (`nperseg=1.0s`, `pwl=(1,4)`)
  apenas para isolar a validação da lógica razão+PLV do ajuste FOOOF.
  Após a correção, **produção e teste usam o mesmo `nfft=4000`**.
- A integração do FOOOF de produção com o portão `erro_ajuste<0.15`
  ainda exige calibração empírica — o limiar rejeita o sintético bem
  comportado (erro=0.18-0.21 vs limiar 0.15).

## Por que `knee` e não `fixed`

LFP real de CA1/DG tem `knee frequency` real:
- ~28 Hz em CA1
- ~70 Hz em DG

(fonte: Kuhn et al. 2026, LFP_FOOOF).

`fixed` assume 1/f puro sem quebra. Em LFP real, `fixed` superestima a
potência de banda larga e mascara teta/gamma. Logo, `fixed` só serve em
sinais sem componente 1/f ou em testes sintéticos sem estrutura.

## A correção `nfft=4000`

Sweep de parâmetros sobre sinal sintético 3x, 1/f realista, 45s, SNR=20dB:

| nperseg_s | nfft  | pwl        | cf_teta | erro    | has_model |
|-----------|-------|------------|---------|---------|-----------|
| 1.0s      | None  | (1,4)      | 7.99    | 0.297   | ✓         |
| 1.0s      | 4000  | (1,4)      | 7.99    | **0.222** | ✓       |
| 1.0s      | 4000  | (2,5)      | 7.99    | **0.222** | ✓       |
| **1.2s (produção)** | None | (2,5) | -       | -       | ✗         |
| **1.2s (produção)** | **4000** | (2,5) | **7.99** | **0.180** | **✓** |
| 2.0s      | None  | (1,4)      | 7.99    | 0.204   | ✓         |
| 2.0s      | None  | (2,5)      | 7.96    | 0.540   | ✓         |
| 2.0s      | 4000  | (2,5)      | 7.98    | 0.424   | ✓         |

A linha chave é a **5ª**: com `nfft=4000`, a configuração de produção
agora detecta teta com erro 0.18 (próximo do limiar 0.15). Kuhn et al.
2026 relatam erro < 0.1% para esta configuração; estamos em ~2%, mas o
sintético não tem a mesma estrutura fina de harmônicos biológicos que o
LFP real.

**Implementação em `extrai_cf_teta_fooof`**:

```python
nperseg = int(nperseg_s * fs)
nfft = 4000 if nperseg <= 4000 else nperseg
freqs, psd = welch(sinal, fs=fs, window='hann',
                    nperseg=nperseg, noverlap=nperseg // 2,
                    nfft=nfft)
```

O `nfft=4000` é zero-padding: Welch ainda janelado em 1.2s (a resolução
estatística real do espectro continua limitada pela duração da janela),
mas o FFT é calculado sobre 4000 pontos, interpolando o espectro para
grade fina de ~0.25 Hz/bin. Isso dá ao FOOOF pontos suficientes para
convergir numa Gaussiana de 2-5 Hz.

## Portão de qualidade (`erro_ajuste < 0.15`)

Após a correção `nfft=4000`:
- `nperseg=1.2s` produção: erro = 0.18 (portão ainda REJEITA)
- `nperseg=1.0s` teste: erro = 0.22
- SNR=0dB: erro = 0.14 (portão ACEITA em SNR muito baixa)
- SNR=20dB: erro = 0.18 (portão REJEITA)
- SNR=30dB: erro = 0.32 (portão REJEITA — mais sinal, mais resíduo)

O portão de 0.15 ainda é chute não-calibrado, e há evidência de que em
sintético ele rejeita o caso ideal. Kuhn et al. 2026 relatam erro < 0.1%
em LFP real — possivelmente o limiar correto está mais perto de 0.20
para sintético, ou a estrutura fina de harmônicos no LFP real torna o
ajuste mais fácil do que no sintético.

**Não dá para calibrar o limiar sem validação em sessões reais**. Mas
agora temos a configuração de detecção que produz cf_teta correto;
falta apenas calibrar o portão.

## Validação de ponta a ponta (com `nfft=4000`)

| Cenário                              | cf    | erro  | razão | PLV   | veredito |
|--------------------------------------|-------|-------|-------|-------|----------|
| A2_2x (harmônico travado)            | 7.99  | 0.18  | sim, n=2 | 0.996 | suspeito |
| A3_3x (harmônico travado)            | 7.99  | 0.18  | sim, n=3 | 0.990 | suspeito |
| A4_4x (harmônico travado)            | 8.00  | 0.19  | sim, n=4 | 0.984 | suspeito |
| A5_5x (harmônico travado)            | 8.00  | 0.19  | sim, n=5 | 0.978 | suspeito |
| B (genuíno 35Hz)                     | 7.99  | 0.21  | não  | -     | clean    |
| C (fase livre 23.7Hz)                | 7.99  | 0.21  | sim, n=3 | 0.327 | rejeitado por PLV |

## Estado atual (frase honesta para apresentação)

> "A lógica de discriminação (razão + PLV) foi validada no cenário
> mais difícil: um oscilador com fase livre cuja frequência cai por
> acaso dentro da tolerância harmônica (γ=23.7Hz vs 3×θ=24Hz, |Δ|=0.3Hz
> dentro de 0.8Hz tol) foi corretamente rejeitado, com PLV=0.327
> contra um limiar de 0.5. A integração com o portão de qualidade do
> FOOOF de produção está em ajuste final — identificamos que nossa
> implementação do Welch não replicava o `nfft=4000` usado no artigo
> original de Kuhn et al. 2026, e a correção restaura a detecção de
> teta em sinal sintético realista. Falta calibrar empiricamente o
> limiar do portão de qualidade (erro < 0.15) em sessões reais."

## Como isso apareceu

Sessão de validação em 04/09/2026. Sequência de correções:

1. **Reverte `aperiodic_mode='fixed'` para `'knee'`** — escolha de
   default não pode ser ditada por teste sintético incompleto.
2. **Adiciona teste de quase-coincidência com fase LIVRE** — sem ele,
   a lógica `suspeito AND skew AND plv` não tinha evidência de que o
   PLV era o discriminador.
3. **Adiciona 1/f realista no sintético** — sem ele, FOOOF não
   tinha estrutura para ajustar.
4. **Identifica `nfft=4000` faltando** — causa-raiz da falha de
   detecção do FOOOF de produção, resolvida por fidelidade ao método
   publicado (Kuhn et al. 2026, Methods).

## Pendências

- [x] Justificar por que o teste sintético usa limiar PLV=0.5 e a
      produção usa 0.8 (default CLI)
- [ ] Calibração empírica do portão `erro_ajuste<0.15` em LFP real
- [ ] Validação biológica em sessões MTESC04/05 (não sintética)
- [ ] Alternativa: usar `specparam` em vez de `fooof` (fooof está
      deprecado)

