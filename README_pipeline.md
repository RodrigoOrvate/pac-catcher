# README — Pipeline PAC theta-gamma (MTESC04 / MTESC05)

Como encontrar e validar acoplamento fase-amplitude theta-gamma no registro
**Basal** de uma sessão (3 arquivos .ns2 consecutivos, 32 canais, 1000 Hz).

Validado nas sessões 08/07 e 09/07/2024 do MTESC04 NOCI (ver o
`registro_resultados.md` de cada sessão e o `CLAUDE.md` desta pasta). Este
guia descreve como rodar uma **nova sessão**.

## Arquitetura (30/08/2026): um SCRIPT central, dados por sessão

O código mora SOMENTE aqui (`C:\acoplamento_theta-gamma\SCRIPT\`) — na raiz
do grupo, acima das pastas de estudo — e é compartilhado por todas as
sessões de todos os estudos (MTESC04/05 × NOCI/LAC). Cada pasta de sessão
guarda apenas dados, saídas e os arquivos de estado:

```
acoplamento_theta-gamma/
├── SCRIPT/                        ← ESTA pasta: todo o código (única cópia)
├── MTESC04_NOCI/
│   ├── MTESC04 -- 1 - infusao - 08-07-2024/
│   │   ├── Basal antes da infusao/    ← os 3 .ns2
│   │   ├── video .MPG                 ← etapa comportamental (manual)
│   │   └── RESULTADOS/                ← SAÍDAS da sessão: resultados*.csv,
│   │      (era "SCRIPT")                 comodulogramas*/, figuras/, logs,
│   │                                      vencedores.csv, registro_resultados.md, .claude/
│   ├── MTESC04 -- 2 - infusao - 09-07-2024/   (idem)
│   └── ... (sessões 3–6)
├── MTESC04_LAC/
├── MTESC05_NOCI/
└── MTESC05_LAC/
```

Todos os comandos abaixo são dados DENTRO desta pasta `SCRIPT\`, apontando
para a sessão por caminho relativo à raiz do grupo.

## Pré-requisitos

```
pip install -r requirements.txt   # numpy scipy matplotlib pandas neo
```

- Pasta com os .ns2 da sessão.
- Vídeo da sessão (para a etapa comportamental — manual).

## Começando em uma nova sessão (do zero)

NADA é copiado para a sessão — o código é este SCRIPT central. Basta que a
pasta da sessão tenha os dados e crie as saídas dentro do próprio `SCRIPT/`
da sessão:

```
<ESTUDO>/<SESSAO>/              ex.: MTESC04_NOCI/MTESC04 -- 3 - infusao - 11-07-2024
├── Basal <nome>/          ← os 3 .ns2 da nova sessão
├── video da sessao.MPG    ← para a etapa comportamental (manual)
└── RESULTADOS/            ← criar vazia; recebe resultados*.csv,
                              comodulogramas*/, figuras/, vencedores.csv,
                              registro_resultados.md, .claude/
```

Nos comandos do fluxo, troque `<ESTUDO>/<SESSAO>` pelo caminho da sessão a
partir da raiz do grupo (ex.: `MTESC04_NOCI/MTESC04 -- 3 - infusao -
11-07-2024`) e `<BASAL>` pela pasta de dados dentro dela.

## Fluxo completo (7 passos)

### 1. Triagem estatística — varre TODAS as janelas
```bash
python triagem_pac.py --pasta "<SESSAO>/<BASAL>" \
    --saida "<SESSAO>/RESULTADOS/resultados.csv"
```
Janelas de 10 s a cada 5 s, todos os canais; KL-MI com 200 surrogates
(deslocamento circular), z-score, proxy de artefato motor.
Padrões: `--janela 10 --passo 5 --n_surr 200 --z_corte 3.0`.

### 2. Refinamento — FDR de janela + filtros de artefato
```bash
python refina_candidatos.py --csv "<SESSAO>/RESULTADOS/resultados.csv" \
    --pasta_ns2 "<SESSAO>/<BASAL>" \
    --saida "<SESSAO>/RESULTADOS/resultados_refinados.csv"
```
p-valor paramétrico (Gama), FDR-BH sobre todas as janelas (m≈5472),
co-ocorrência entre canais (descarta ruído de modo comum), kurtose em gamma,
saturação. Linhas com veredito **"Candidato robusto"** seguem adiante.

### 2.5. Teste de skewness do theta (anti-harmônicos) — `audita_skewness.py`
```bash
python auditorias/audita_skewness.py \
    --csv "<SESSAO>/RESULTADOS/vencedores.csv" \
    --pasta_ns2 "<SESSAO>/<BASAL>" \
    --saida "<SESSAO>/RESULTADOS/auditoria/skewness.csv"
```
Verifica se a banda de fase é senoide limpa. Quando a razão amp/fase é inteira
(5×20=4,0; 5×35; 5×45), uma theta assimétrica ("dente de serra") gera MI
espúrio por carregar um harmônico na banda de amplitude. Skewness
Fisher-Pearson |g|>0,5 (padrão `--limiar`) = **SUSPECT**; ~0 = **CLEAN**.
Roda na janela cheia e, se houver, na ilha reancorada (`inicio_ilha_s/fim_ilha_s`).
**Atenção:** limpar não prova ausência de problema — só descarta a hipótese
harmônica (ex.: o par 5×20 do MTESC05 era CLEAN mas foi rejeitado por
instabilidade temporal + FDR, não por harmônico).

### 3. Comodulogramas em lote (com notch de 60 Hz)
```bash
python comodulogram.py --csv "<SESSAO>/RESULTADOS/resultados_refinados.csv" \
    --pasta_ns2 "<SESSAO>/<BASAL>" \
    --saida_dir "<SESSAO>/RESULTADOS/comodulogramas_notch" \
    --veredito_prefixo "Candidato robusto" --notch 60
```
MI z-scoredo por célula do mapa (fase 4–14 Hz × amplitude 30–150 Hz).
`resumo_comodulogramas.csv` traz o pico ΘΓ de cada janela (z, frequências e
MI bruto). **Sempre rodar também sem notch** (`--saida_dir comodulogramas`)
e comparar: pico que cai >1,5z com o notch é rede elétrica, não acoplamento.

**⚠️ NOTCH HARMÔNICOS:** Além do notch em 60 Hz, aplicar também em **120 Hz** (e
opcionalmente 180 Hz) para controlar harmônicos da rede elétrica brasileira.

### 4. Conferência comportamental — MANUAL
- Derivar o mapeamento vídeo↔ns2 da NOVA sessão (o offset muda!):
  anotar o instante do vídeo em que o rato é colocado = início do 1º .ns2.
  Na sessão 08/07: `tempo_vídeo = tempo_ns2 + 321 s`.
- Assistir ao vídeo nos instantes dos candidatos do resumo e anotar o
  comportamento (rearing, grooming, walking, imóvel...). Este passo não é
  automatizável — é o rótulo comportamental que valida ou rejeita.

**⚠️ REGRAS DE SUB-JANELA (evitar double dipping):**
- Se a janela inicial de 10s contiver dois estados comportamentais, a
  sub-janela deve ser escolhida por **ENERGIA THETA** (não por MI) para
  evitar circularidade: "Sempre pegar o sub-segmento de 3s de maior energia
  theta, independente do MI."
- O MI/Z-score deve ser recalculado nessa sub-janela, nunca reutilizando
  o valor da janela original.
- **Split-half (recomendado para fechar o viés):** usar primeira metade da
  ilha (ex: 3s) só para decidir o ponto de ancoragem (por energia theta),
  e calcular MI/Z apenas na segunda metade (os 3s restantes), nunca vista
  na etapa de seleção. Isso quebra circularidade de fato.

### 5. Teste de respiração (descarta a fase respiratória)
```bash
python diagnostico_janela.py \
    --arquivo "<SESSAO>/<BASAL>/ARQUIVO.ns2" \
    --canal chan20 --inicio 85 --fim 95 --notch 60 \
    --saida_png "<SESSAO>/RESULTADOS/diagnosticos/diag_chan20_85-95s.png"
```
Para cada candidato: se `MI resp×γ` também for alto, o "theta" é sniffing/
respiração (6–10 Hz), não theta. Sobrevive quem tem θ×γ alto e resp×γ ≈ 0.
Checar também o fator de crista do theta (>2 = transientes, suspeito).
Atenção: a banda respiratória do proxy é 0,5–3 Hz — CEGA para sniffing
4–8 Hz; desempate fino exige vídeo.

### 6. FDR sobre o mapa do comodulograma
```bash
python comodulogram.py --csv "<SESSAO>/RESULTADOS/resultados_refinados.csv" \
    --pasta_ns2 "<SESSAO>/<BASAL>" \
    --saida_dir "<SESSAO>/RESULTADOS/comodulogramas_fdr" \
    --veredito_prefixo "Candidato robusto" \
    --notch 60 --fdr_q 0.05
```
Contorno preto nas células significantes. Classificação: **"concentrado em
ΘΓ"** (≥2 células sig, ≥50% em ΘΓ) = acoplamento genuíno e estreito;
"esparso/fora de ΘΓ" (muitas células espalhadas) = transientes ritmados.

### 7. Robustez de parâmetros + figura dos vencedores
Criar o **`vencedores.csv` na pasta RESULTADOS da sessão** (uma linha por
vencedor; sem vírgulas no texto de comportamento):

```
rotulo,arquivo,canal,inicio_s,fim_s,fase_pico_hz,amp_pico_hz[,comportamento]
vencedor1_rearing,20240708-123605-001.ns2,chan20,85,95,7,75,Rearing — vídeo 6:46-6:56
```

O nome do PNG gerado é `<rotulo>_<canal>.png`. Se `fase_pico_hz`/`amp_pico_hz`
ficarem vazios, o par de pico é lido do resumo FDR (`--resumo_fdr`). Rodar da
pasta central:

```bash
python robustez_parametros.py \
    --pasta "<SESSAO>/<BASAL>" \
    --vencedores "<SESSAO>/RESULTADOS/vencedores.csv" \
    [--resumo_fdr "<SESSAO>/RESULTADOS/comodulogramas_fdr/resumo_comodulogramas.csv"] \
    --saida_csv "<SESSAO>/RESULTADOS/robustez_parametros.csv"

python figura_apresentacao.py \
    --pasta_ns2 "<SESSAO>/<BASAL>" \
    --vencedores "<SESSAO>/RESULTADOS/vencedores.csv" \
    --saida_dir "<SESSAO>/RESULTADOS/figuras"
```
Critério: z estável (≥3) em todo o sweep de n_bins, pico do mapa imóvel,
MVL confirmando.

> **Extrair picos com banda + FDR (recomendado ANTES de publicar o par):**
> `python extrair_picos.py --csv "<SESSAO>/RESULTADOS/vencedores.csv"
> --pasta_ns2 "<SESSAO>/<BASAL>"
> --saida "<SESSAO>/RESULTADOS/auditoria/extracao_picos_v1_vs_v2.csv"`.
> Pico fora da banda 4–8 × 30–80 Hz ou não-significativo sob FDR na célula
> vencedora é desqualificado (nada de argmax cru) — roda sempre na banda
> padrão; adicionar `--busca_ampla` só para auditar fora dela.

### 7.5. Held-out anti-double-dipping — `audita_held_out.py`
```bash
python auditorias/audita_held_out.py \
    --csv "<SESSAO>/RESULTADOS/vencedores.csv" \
    --pasta_ns2 "<SESSAO>/<BASAL>" \
    --saida "<SESSAO>/RESULTADOS/auditoria/held_out.csv"
```
OBRIGATÓRIO quando uma janela foi **reancorada** (`inicio_ilha_s/fim_ilha_s`
no vencedores.csv). Compara o MI da ILHA reportada contra o RESTO-das-janela
(+ janela adjacente), z vs 200 surrogates — confirma que o efeito é focal na
ilha e não inflado pela circularidade da seleção. "Validado" só se z≥3 na
ilha e <3 no resto. Se nenhuma janela foi reancorada, grava `N/A` e não é
necessário.

### 8. Registro de resultados — `registro_resultados.md`
Ao fim da sessão, criar/atualizar `registro_resultados.md` na pasta
RESULTADOS da sessão: identificação (rato/dia/arquivos/vídeo + offset),
números das etapas, vencedores com o checklist das 5 etapas e rejeitados
com motivo.
Um arquivo por sessão — o índice global dos animais é consolidado depois a
partir deles. Ver exemplo preenchido na sessão 09/07/2024.

## Checklist de validação (as 5 etapas)

Um acoplamento só é "VALIDADO" com **todas**:
1. FDR de janela (etapa 2) ✓
2. Não é 60 Hz (etapa 3, com/sem notch) ✓
3. Não é respiração (etapa 5) ✓
4. FDR do mapa concentrado em ΘΓ (etapa 6) ✓
5. Robustez de parâmetros + métrica alternativa (etapa 7) ✓
   — e o comportamento no vídeo é coerente (etapa 4).

## O que muda entre sessões (e o que não muda)

| Item | Muda? | Observação |
|---|---|---|
| Caminhos (`<SESSAO>`, `<BASAL>`) | Sim | todo comando |
| Offset vídeo↔ns2 | **Sim** | re-derivar assistindo ao vídeo (08/07: +321 s; 09/07: ≈0) |
| Rótulos de comportamento | **Sim** | manual, etapa 4 |
| `vencedores.csv` da sessão | Sim | criar após escolher os vencedores |
| Parâmetros (janela 10 s, bandas 4–8 × 30–80 Hz, n_bins 18, 200 surrogates) | **Não** | padrões de rato; revisar só se a sessão tiver estado farmacológico muito diferente (ex.: NOCI reduz as frequências) |
| Notch 60 Hz | **Não** | rede brasileira; usar 50 se dados de outro país |
| Critérios de validação | **Não** | as 5 etapas acima |

## Limitações a lembrar

- N = 1 animal/sessão por análise: validação interna, não generalização.
- Sem EMG: risco de circularidade movimento→EMG→PAC é mitigado, não eliminado.
- A varredura cobre theta 4–8 × gamma 30–80 em janelas de 10 s; acoplamentos
  mais curtos ou fora dessas bandas não seriam detectados.
- **Sniffing 4–8 Hz:** proxy respiratório (0,5–3 Hz) é cego para sniffing na
  banda theta; desempate requer vídeo qualitativo ou termistor.
- **Double dipping:** sub-janelas escolhidas por energia theta, não por MI,
  para evitar circularidade (regra: energia theta máxima, não MI máximo).
- **Integer ratio check:** razão amp/fase inteira (ex: 5×35, 5×45) pode indicar
  theta não-senoidal gerando MI espúrio; checar skewness do theta filtrado.
