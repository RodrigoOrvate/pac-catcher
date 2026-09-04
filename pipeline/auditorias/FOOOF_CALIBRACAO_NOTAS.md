---
name: fooof-calibracao
description: Notas sobre calibração do FOOOF em audita_harmonico.py — por que `knee` é o default, e como o sintético difere do LFP real
---

# Notas de calibração FOOOF (audita_harmonico.py)

## TL;DR

- `aperiodic_mode='knee'` é o **default de produção** e deve permanecer.
- Não trocar o default só porque o teste sintético falha — isso é um sinal
  de que o sintético está incompleto, não de que `knee` é ruim.

## Por que `knee` e não `fixed`

LFP real de CA1/DG tem uma `knee frequency` real:
- ~28 Hz em CA1
- ~70 Hz em DG

(fonte: Kuhn et al. 2026, LFP_FOOOF).

O modelo `knee` ajusta essa quebra espectral; `fixed` assume 1/f puro sem
quebra. Em LFP real, `fixed` superestima a potência de banda larga e
mascara teta/gamma. Logo, `fixed` só serve em sinais sem componente 1/f
ou em testes sintéticos.

## O que o teste sintético exige para ser realista

Para o FOOOF detectar teta em sinal sintético, o sinal precisa ter:
- **Fundo 1/f realista** (pink noise com knee ≈ 28 Hz, slope ≈ 1.2)
- **Resolução de frequência compatível** com `peak_width_limits`:
  - `nperseg=1.2s` → freq_res = 0.83 Hz; `peak_width_limits=(2,5)` exige
    2-5 bins. Em 4-12 Hz cabem apenas 10 bins. Ajuste falha.
  - `nperseg=1.0s` → freq_res = 1 Hz; `peak_width_limits=(1,4)` ajusta.

Os parâmetros acima (relaxados) são para o **teste sintético apenas**.
O código de produção usa `nperseg=1.2s` e `peak_width_limits=(2,5)`,
calibrados para LFP real.

## Portão de qualidade (`erro_ajuste < 0.15`)

Chute não-calibrado. No teste sintético com fundo 1/f realista, o erro
de ajuste ficou em 0.28-0.44 (acima do limiar). O portão rejeitaria
todos os sinais sintéticos — confirmando que **o limiar precisa ser
calibrado empiricamente em sessões reais**, não no sintético.

## Como isso apareceu

Sessão de validação em 04/09/2026. Reporte inicial sugeria trocar
`knee` por `fixed` para fazer o sintético passar; revertido depois de
constatar que a estrutura 1/f do sintético estava ausente.
