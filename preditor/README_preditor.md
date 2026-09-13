# Preditor de PAC θ-γ (tempo real / optogenética)

**Contexto:** o grupo quer não só *detectar* acoplamento theta-gamma, mas
**prever** que ele vai acontecer ~10 s antes e, em tempo real, disparar um pulso
**TTL** para ligar um laser de optogenética. Esta pasta contém o protótipo
desse preditor. É código `SCRIPT/` (nada de sessão).

> ⚠️ **Estado 13/09/2026: PROTÓTIPO não validado.** Nenhum modelo aqui é
> confiável para acionar luz ainda. Duas linhas de investigação, resultado
> oposto:
> - **Features espectrais/PAC locais no LFP bruto (`validar_preditor.py`)**:
>   agora com o resultado atual do pipeline (190 eventos de
>   `candidatos_vencedores_OURO_PURIFICADO_v2.csv`, não mais o schema legado
>   de ~10 eventos) — LOOCV com controles reais dá acurácia ≈0.43, ao nível
>   de chance, mesmo incluindo MI/MVL no par oficial e footprint espacial.
>   Mann-Whitney confirma: nenhuma das 8 features separa pré-evento de
>   controle real (p>0.1 em todas). **Não há precursor detectável no LFP
>   bruto nos 10s antes do evento** com essas métricas.
> - **Estado comportamental + potência teta (`preditor_estado_comportamental.py`,
>   novo)**: em vez de procurar uma rampa no sinal, testa se comportamento da
>   janela anterior + potência teta dessa janela preveem a janela seguinte.
>   Baseline (banda 4-8Hz bruta + comportamento): AUC-ROC 0.602. Rodado um
>   ablation de 4 melhorias candidatas — só uma se confirmou: separar a
>   potência numa janela mais próxima do evento + a tendência (subindo ou
>   caindo) ao longo dos últimos ~10s levou a AUC-ROC **0.618**.
>   Normalização por canal (z-score) e MI num par fase-amplitude fixo NÃO
>   ajudaram isoladamente (testado e descartado, ver histórico de commits).
>   Consistente com Vanderwolf (1969)/Tort et al. (2009)/Colgin (2016): PAC
>   é dependente de estado comportamental, não um evento espontâneo do LFP
>   local — e a feature mais forte agora é a potência teta RÁPIDA (7-10Hz,
>   tipo 1, ligada a locomoção), não a banda genérica 4-8Hz.
>
> **Não usar nenhum dos dois para disparar hardware** até validação mais
> robusta (AUC 0.62 está longe do necessário para controlar luz em tempo
> real). Ver seção "Limitação conhecida" abaixo.

## Fluxo dos dados

```
vencedores.csv (todas as sessões)
      │
      ▼
analisar_pre_evento.py ──► D:\acoplamento_theta-gamma\ANALISE_PRE_EVENTO\
      │                       raw_<estudo>_<sessao>_<canal>.csv  (10 s antes do evento)
      │                       pre_*.png  (espectrogramas)
      ▼
treinar_preditor.py  ──► SCRIPT/preditor/modelo_pac.pkl
      │                    (RandomForest; positivos = pré-evento real;
      │                     negativos = surrogate phase-scramble)
      ▼
validar_preditor.py   ──► valida contra janelas de CONTROLE REAIS (mesmo .ns2,
      │                     evento +60 s), LOOCV honesto. Re-salva o modelo.
      ▼
prever_pac_tempo_real.py ─► janela deslizante 10 s / passo 1 s → P(PAC) →
                             dispara_ttl() quando P ≥ limiar (0.7)
```

## Scripts

| Script | O que faz | Entrada → Saída |
|---|---|---|
| `analisar_pre_evento.py` | Extrai os 10 s anteriores a cada vencedor de todos os `vencedores.csv` (schema legado, não usado pelos scripts abaixo) | `vencedores.csv` → `raw_*.csv` + `pre_*.png` |
| `treinar_preditor.py` | Treina RandomForest com surrogates (otimista, schema legado) | `ANALISE_PRE_EVENTO/raw_*.csv` → `modelo_pac.pkl` |
| `validar_preditor.py` | Features espectrais/PAC (energia, MI/MVL no par oficial, footprint) vs. controles reais, LOOCV honesto | `resultados/candidatos_vencedores_OURO_PURIFICADO_v2.csv` + `.ns2` → relatório; re-salva `modelo_pac.pkl` |
| `analisar_janelas_ultracurtas.py` | Testa precursores ultracurtos (1s, 2s, 3s) e onset no LFP bruto vs controles pareados | `candidatos_vencedores_OURO_PURIFICADO_v2.csv` + `.ns2` → relatório estatístico |
| `preditor_estado_comportamental.py` | Estado comportamental + potência teta da janela anterior, validação cruzada 5-fold | `resultados/dataset_mestre_COM_COMPORTAMENTO.csv` + `.ns2` → relatório; salva `modelo_estado_comportamental.pkl` |
| `prever_pac_tempo_real.py` | Gating de estado em tempo real/replay usando o modelo consolidado; `dispara_ttl()` é plugável | `.ns2` / streaming → P(PAC) por janela + disparo TTL |

Para síntese teórica completa formatada para a dissertação de mestrado, consulte `docs/sintese_investigacao_preditor_pac.md`.

## Features usadas

**`validar_preditor.py`** (8 features, por canal/janela, na convenção
1-based do dataset mestre): `Energia (variância)`, `Theta_Energy (4–8 Hz)`,
`Gamma_Energy (30–80 Hz)`, `Ratio_TG`, `Theta_Peak` (Welch, nperseg 1000) +
`MI_z_oficial`, `MVL_z_oficial` (no par fase/amp oficial do evento, via
`mi_z_par`/`mvl_z_par` de `robustez_parametros.py`) + `Footprint_n_z3`
(quantos dos outros 31 canais já mostram z≥3 no mesmo par).

**`preditor_estado_comportamental.py`** (13 features, config. vencedora do
ablation): potência teta em duas sub-bandas — lenta/tipo2 (4–7Hz,
sniffing/imobilidade) e rápida/tipo1 (7–10Hz, locomoção, distinção de
Vanderwolf) — cada uma em 3 versões (janela N-1 inteira, sub-janela mais
próxima de N, e tendência/inclinação ao longo de 5 sub-janelas) + one-hot
do `comportamento` anotado de N-1 (7 categorias). Todas medidas em N-1 para
prever se N vira vencedor — formulação sem vazamento. O script roda um
ablation completo antes de treinar o modelo final e imprime o efeito
isolado de cada melhoria testada (inclusive as que NÃO ajudaram:
normalização por canal e MI num par canônico fixo).

## Dependências

`neo` (leitura .ns2 via `pac_core.io`), `numpy`, `scipy`, `pandas`,
`scikit-learn`, `joblib`, `matplotlib`.

## Comandos

```bash
# Da pasta SCRIPT/preditor/
python validar_preditor.py                  # features espectrais/PAC vs. controle real (honesto)
python preditor_estado_comportamental.py    # estado comportamental + potencia teta (honesto)
python prever_pac_tempo_real.py --ns2 <arquivo.ns2> --canal chan16 --replay
python prever_pac_tempo_real.py --simular   # demo do loop de disparo
```

## Limitação conhecida / lições

- **Surrogates ≠ controle real**: o treino com fase-scrambled (`treinar_preditor.py`)
  reporta 0.83, mas isso é otimista — sempre validar com janelas reais de
  controle (`validar_preditor.py`/`preditor_estado_comportamental.py`).
- **LFP bruto local não basta**: mesmo com MI/MVL no par oficial e footprint
  espacial (não só energia espectral genérica), e com 373 amostras reais
  (18x o protótipo original de ~10), nenhuma feature derivada só do sinal
  elétrico local separa pré-evento de controle. O acoplamento parece mais
  um burst transiente que uma rampa gradual visível no LFP isolado.
- **Estado comportamental funciona, modestamente**: comportamento da janela
  anterior + potência teta chegam a AUC 0.602 — sobrevive a restringir só
  transições reais de comportamento (p=0.018, descarta autocorrelação pura)
  e a controlar pela própria potência teta (teste de razão de
  verossimilhança p=1.3e-5, descarta "é só teta disfarçado"). Ainda é uma
  discriminação fraca-a-moderada, não confiável para hardware.
- **Ablation de melhorias (13/09/2026)**: testadas 4 hipóteses sobre o
  baseline de estado. Só uma se confirmou: **janela mais próxima do evento
  + tendência de subida/descida** ao longo dos últimos ~10s (AUC 0.602 →
  0.618). As outras três NÃO ajudaram isoladamente, e isso foi verificado
  explicitamente (não só assumido): normalização de potência teta por canal
  (z-score) é neutra, não prejudicial — o resultado inicial pior era um
  artefato de testá-la sem as outras features; MI num par fase-amplitude
  canônico fixo (6Hz×85Hz, mediana dos 190 vencedores) piora levemente; e
  separar teta rápida (tipo 1, locomoção)/lenta (tipo 2, sniffing) não deu
  ganho isolado, mas a rápida virou a feature de maior importância no
  modelo final combinado. Limiar calibrado por F1 (0.431): recall 0.896 do
  vencedor com precisão 0.452 — útil se a prioridade é não perder o evento
  (custo: quase metade dos disparos seria falso-positivo).
- **Próximos candidatos não testados**: janelas de pré-evento mais curtas
  que os ~10s atuais (testar 1–3s, ainda mais perto do início); combinar
  estado+teta com as features PAC locais de `validar_preditor.py` (MI/MVL/
  footprint) num único modelo, já que nunca foram testadas juntas; mais
  sub-bandas teta (a divisão 4–7/7–10Hz é uma primeira aproximação grosseira
  da distinção de Vanderwolf, não uma calibração fina).
- **TTL**: `dispara_ttl()` em `prever_pac_tempo_real.py` é um stand-in — plugar
  o backend (NI-DAQ `nidaqmx`, Arduino serial etc.) quando houver a placa.
