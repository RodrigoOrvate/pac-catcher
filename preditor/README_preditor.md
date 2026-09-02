# Preditor de PAC θ-γ (tempo real / optogenética)

**Contexto:** o grupo quer não só *detectar* acoplamento theta-gamma, mas
**prever** que ele vai acontecer ~10 s antes e, em tempo real, disparar um pulso
**TTL** para ligar um laser de optogenética. Esta pasta contém o protótipo
desse preditor. É código `SCRIPT/` (nada de sessão).

> ⚠️ **Estado 31/08/2026: PROTÓTIPO não validado.** O modelo atual NÃO é
> confiável para acionar luz. Com janelas de controle reais (não surrogates),
> a validação LOOCV dá acurácia ≈0.45 (≈ chance), recall pré-PAC ≈0.40 e taxa
> de falso-positivo ≈0.50. **Não usar para disparar hardware até ter mais
> dados e melhor validação.** Ver `validar_preditor.py`.

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
| `analisar_pre_evento.py` | Extrai os 10 s anteriores a cada vencedor de todos os `vencedores.csv` | `vencedores.csv` → `raw_*.csv` + `pre_*.png` |
| `treinar_preditor.py` | Treina RandomForest (features espectrais simples) | `ANALISE_PRE_EVENTO/raw_*.csv` → `modelo_pac.pkl` |
| `validar_preditor.py` | Valida contra controles reais (LOOCV, matriz de confusão, recall/FDR) | `.ns2` + `vencedores.csv` → relatório; re-salva modelo |
| `prever_pac_tempo_real.py` | Preditor em tempo real/replay; `dispara_ttl()` é plugável | `.ns2` → P(PAC) por janela + disparo TTL |

## Features usadas (ordem fixa)

`Energia (variância)`, `Theta_Energy (4–8 Hz)`, `Gamma_Energy (30–80 Hz)`,
`Ratio_TG`, `Theta_Peak`. Calculadas na **média dos canais** de cada janela
via Welch (nperseg 1000 para resolução θ). Estão hardcode nos 3 scripts
(extração/treino/validação); se mudar uma, mudar nas três.

## Dependências

`neo` (leitura .ns2 via `ns2_utils` em `SCRIPT/pipeline/`), `numpy`,
`scipy`, `pandas`, `scikit-learn`, `joblib`, `matplotlib`.

## Comandos

```bash
# Da pasta SCRIPT/preditor/
python analisar_pre_evento.py         # gerar/regenerar pré-eventos
python treinar_preditor.py            # treinar com surrogates (otimista)
python validar_preditor.py            # validar com controles reais (honesto)
python prever_pac_tempo_real.py --ns2 <arquivo.ns2> --canal chan16 --replay
python prever_pac_tempo_real.py --simular          # demo do loop de disparo
```

## Limitação conhecida / lições

- **Surrogates ≠ controle real**: o treino com fase-scrambled reporta 0.83,
  mas isso é otimista — o modelo não separa pré-PAC de atividade normal real
  (0.45). Sempre validar com janelas reais de controle.
- **Dados insuficientes**: só ~10 eventos. Precisamos de dezenas/hipercentenas
  antes de qualquer modelo funcionar.
- **Features ainda fracas**: média do canal + Welch não capturam o que difere
  o pré-PAC. Próximos candidatos: MVL de PAC em janela curta, espectrograma
  2D fase×amp achatado, e o **footprint espacial** (canais vizinhos) que já
  funciona no pipeline.
- **TTL**: `dispara_ttl()` em `prever_pac_tempo_real.py` é um stand-in — plugar
  o backend (NI-DAQ `nidaqmx`, Arduino serial etc.) quando houver a placa.
