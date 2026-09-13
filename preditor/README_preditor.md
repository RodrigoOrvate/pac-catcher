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
>   Validação cruzada 5-fold: acurácia 0.584, **AUC-ROC 0.597** — modesto,
>   mas a primeira abordagem que bate chance de forma real. Consistente com
>   Vanderwolf (1969)/Tort et al. (2009)/Colgin (2016): PAC é dependente de
>   estado comportamental, não um evento espontâneo do LFP local.
>
> **Não usar nenhum dos dois para disparar hardware** até validação mais
> robusta (AUC 0.6 está longe do necessário para controlar luz em tempo
> real). Ver seção "Limitação conhecida" abaixo.

## Fluxo dos dados

```
vencedores.csv (todas as sessões)
      │
      ▼
analisar_pre_evento.py ──► C:\acoplamento_theta-gamma\ANALISE_PRE_EVENTO\
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
| `preditor_estado_comportamental.py` | Estado comportamental + potência teta da janela anterior, validação cruzada 5-fold | `resultados/dataset_mestre_COM_COMPORTAMENTO.csv` + `.ns2` → relatório; salva `modelo_estado_comportamental.pkl` |
| `prever_pac_tempo_real.py` | Preditor em tempo real/replay; `dispara_ttl()` é plugável | `.ns2` → P(PAC) por janela + disparo TTL |

`analisar_pre_evento.py`/`treinar_preditor.py` ficaram no schema antigo
(`RESULTADOS/vencedores.csv`, ~10 eventos) e não foram atualizados — o
caminho validado hoje é `validar_preditor.py` (retreina e re-salva o modelo
final sozinho, sem precisar dos dois primeiros) e o novo
`preditor_estado_comportamental.py`.

## Features usadas

**`validar_preditor.py`** (8 features, por canal/janela, na convenção
1-based do dataset mestre): `Energia (variância)`, `Theta_Energy (4–8 Hz)`,
`Gamma_Energy (30–80 Hz)`, `Ratio_TG`, `Theta_Peak` (Welch, nperseg 1000) +
`MI_z_oficial`, `MVL_z_oficial` (no par fase/amp oficial do evento, via
`mi_z_par`/`mvl_z_par` de `robustez_parametros.py`) + `Footprint_n_z3`
(quantos dos outros 31 canais já mostram z≥3 no mesmo par).

**`preditor_estado_comportamental.py`** (8 features): `log_theta_power_anterior`
(Welch 4–8Hz da janela N-1) + one-hot do `comportamento` anotado da janela
N-1 (7 categorias). Ambos medidos em N-1 para prever se N vira vencedor —
formulação sem vazamento (só usa o que estaria disponível no momento da
decisão).

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
  anterior + potência teta chegam a AUC 0.597 — sobrevive a restringir só
  transições reais de comportamento (p=0.018, descarta autocorrelação pura)
  e a controlar pela própria potência teta (teste de razão de
  verossimilhança p=1.3e-5, descarta "é só teta disfarçado"). Ainda é uma
  discriminação fraca-a-moderada, não confiável para hardware.
- **Próximos candidatos**: um proxy de movimento em tempo real (a própria
  potência teta de banda larga já serve, sem precisar de vídeo/EMG);
  janelas de pré-evento mais curtas que 10s (testar 1–3s, mais perto do
  início) para ver se o precursor comportamental fica mais forte perto do
  evento; combinar as duas abordagens (features PAC locais + estado) num
  único modelo.
- **TTL**: `dispara_ttl()` em `prever_pac_tempo_real.py` é um stand-in — plugar
  o backend (NI-DAQ `nidaqmx`, Arduino serial etc.) quando houver a placa.
