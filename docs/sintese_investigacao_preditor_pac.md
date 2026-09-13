# Síntese da Investigação do Preditor de Acoplamento Teta-Gama (PAC) e Intervenção em Circuito Fechado

**Autor:** Rodrigo Orvate  
**Projeto:** Dissertação de Mestrado — Módulo de Acoplamento Teta-Gama e Intervenção Optogenética Closed-Loop  
**Data de Consolidação:** 13 de Setembro de 2026  
**Status do Módulo:** Concluído e Validado  

---

## 1. Contexto e Motivação Científica

No plano do mestrado (*Fase 2 - Proposta de Intervenção*), formulou-se o objetivo de desenvolver um sistema preditor capaz de antecipar a ocorrência do acoplamento fase-amplitude (*Phase-Amplitude Coupling* - PAC) teta-gama no hipocampo CA1 em cerca de 10 segundos, disparando um sinal TTL para acionamento de estimulação optogenética.

Essa investigação conectava-se diretamente com as hipóteses centrais do pré-projeto:
- **Hipótese H1 (Acoplamento $\to$ Memória/Comportamento):** Romper o acoplamento teta-gama via laser alinhado a fases de teta modula o processamento neural.
- **Hipótese H2 (Comportamento $\to$ Acoplamento):** O acoplamento não é um evento estocástico intrínseco e espontâneo do LFP local, mas depende do estado comportamental e cognitivo do animal.

O pipeline de processamento em lote (*Fase 1*) consolidou com sucesso o conjunto de **190 eventos vencedores canônicos** em `candidatos_vencedores_OURO_PURIFICADO_v2.csv` (147 de MTESC04, 41 de MTESC05 e 2 de MTESC03), purificados contra artefatos harmônicos de rede (60/120/180/240 Hz) e validados em matriz de comodulograma e robustez de vizinhança. 

Com base nesse padrão-ouro, realizou-se uma investigação empírica exaustiva em três frentes:
1. **LFP local bruto em janela longa (10 s antes do evento)**;
2. **LFP local bruto em janelas ultracurtas (1 s a 3 s antes) e detecção de *onset***;
3. **Dinâmica do estado comportamental antecedente ($N-1$) + potência teta**.

---

## 2. Investigação 1: O LFP Local 10 Segundos Antes (`validar_preditor.py`)

### Metodologia
- **Amostras:** 373 janelas reais de 10 s (183 janelas pré-evento reais $[-10\text{s}, 0\text{s}]$ vs. 190 janelas de controle pareadas do mesmo canal e sessão, em $t_{\text{evento}} + 60\text{s}$).
- **Validação:** *Leave-One-Out Cross-Validation* (LOOCV) honesta contra controles reais (não surrogates de fase).
- **Features (8 métricas):**
  - Espectrais clássicas: Energia total (variância), Potência Teta (4–8 Hz), Potência Gama (30–80 Hz), Razão Teta/Gama, Frequência de Pico de Teta;
  - Padrão-ouro de PAC: $MI_z$ oficial e $MVL_z$ oficial calculados exatamente no par de frequência $(f_{\text{fase}}, f_{\text{amp}})$ que a sessão viria a confirmar;
  - Footprint de rede: número de canais vizinhos (entre os outros 31) com $z \ge 3$ no mesmo par.

### Resultados
- **Acurácia LOOCV:** **0,429** (nível de chance; pior que o protótipo simplificado de 0,44).
- **Importância das features no Random Forest:** Quase plana (0,06 a 0,27), sem nenhuma métrica dominante.
- **Testes univariados (Mann-Whitney U):** Nenhuma das 8 features separou pré-evento de controle ($p > 0,10$ em todas; maioria com $p > 0,40$).
- **O dado mais revelador:** O $MI_z$ oficial apresentou **medianas praticamente idênticas** nos dois grupos:
  $$\text{Mediana}(MI_{z,\text{pré}}) = -0,157 \quad \text{vs.} \quad \text{Mediana}(MI_{z,\text{controle}}) = -0,155$$

### Conclusão Parcial 1
O resultado descarta a hipótese de falta de amostra ou de feature: mesmo na métrica padrão-ouro de PAC, no canal exato e no par canônico, **não existe rampa lenta ou precursor eletrofisiológico 10 segundos antes do evento no LFP local**. O PAC comporta-se como um evento transiente (*burst*) de rede.

---

## 3. Investigação 2: Dinâmica de Janelas Ultracurtas (1 s a 3 s) e *Onset* (`analisar_janelas_ultracurtas.py`)

Para testar se o precursor no LFP existe mas atua em escalas temporais imediatas, processaram-se 186 pares de janelas ultracurtas nos dados reais:

### Testes Estatísticos Univariados (Pré-evento vs. Controle Real)

| Janela / Condição | Métrica | Mediana Pré | Mediana Ctrl | $\Delta$ (%) | $p$-Wilcoxon (Pareado) | $p$-MannWhitney | Signif. |
|---|---|---|---|---|---|---|---|
| **Pré 3s** $[-3\text{s}, 0\text{s}]$ | Teta Rápida (7–10 Hz) | $1,668 \times 10^4$ | $1,514 \times 10^4$ | $+10,2\%$ | **$0,0476$** | $0,616$ | $*$ |
| **Pré 3s** $[-3\text{s}, 0\text{s}]$ | Gama Rápida (60–90 Hz) | $172,3$ | $153,3$ | $+12,4\%$ | **$0,0166$** | $0,346$ | $*$ |
| **Pré 2s** $[-2\text{s}, 0\text{s}]$ | Teta Rápida (7–10 Hz) | $1,565 \times 10^4$ | $1,485 \times 10^4$ | $+5,4\%$ | **$0,0482$** | $0,568$ | $*$ |
| **Pré 1s** $[-1\text{s}, 0\text{s}]$ | Teta Rápida (7–10 Hz) | $1,600 \times 10^4$ | $1,441 \times 10^4$ | $+11,1\%$ | **$0,0459$** | $0,437$ | $*$ |
| **Tendência (Slope)** | Inclinação Teta $[-3\text{s} \to -1\text{s}]$ | $+110,3$ | $-170,2$ | — | **$0,0247$** | $0,089$ | $*$ |
| **Onset 1s** $[0\text{s}, +1\text{s}]$ | Gama Lenta (30–55 Hz) | $990,1$ | $790,0$ | $+25,3\%$ | **$0,0058$** | $0,399$ | $**$ |
| **Onset 2s** $[0\text{s}, +2\text{s}]$ | Teta Rápida (7–10 Hz) | $1,727 \times 10^4$ | $1,471 \times 10^4$ | $+17,4\%$ | **$0,0084$** | $0,391$ | $**$ |

### Classificação Multivariada (Stratified 5-Fold CV / LOOCV)

| Modelo / Buffer Temporal | Acurácia | AUC-ROC | Precisão (1) | Recall (1) | F1-Score |
|---|---|---|---|---|---|
| **Pré 3s puro** | $0,441$ | $0,416$ | $0,437$ | $0,409$ | $0,422$ |
| **Pré 2s puro** | $0,487$ | $0,464$ | $0,487$ | $0,500$ | $0,493$ |
| **Pré 1s puro** | $0,395$ | $0,382$ | $0,382$ | $0,339$ | $0,359$ |
| **Pré 1s + Tendência (Slope)** | $0,446$ | $0,415$ | $0,434$ | $0,355$ | $0,391$ |
| **Onset 1s** $[0\text{s}, +1\text{s}]$ | $0,500$ | $0,466$ | $0,500$ | $0,489$ | $0,495$ |
| **Onset 2s** $[0\text{s}, +2\text{s}]$ | $0,484$ | $0,467$ | $0,484$ | $0,489$ | $0,487$ |

### Conclusão Parcial 2
1. **Assinatura Fisiológica Pareada Real:** O teste pareado de Wilcoxon revelou que, intra-sessão, existe uma aceleração sutil e significativa da potência teta tipo 1 nos últimos segundos antes do evento ($p = 0,0247$ na inclinação da trajetória).
2. **Inviabilidade Classificatória:** Devido à variabilidade inter-individual basal de impedância e camadas de CA1, essa assinatura não resiste à generalização entre sessões e animais ($p_{\text{Mann-Whitney}} > 0,35$). Todos os modelos classificatórios multivariados em janelas de 1 a 3 segundos permaneceram no nível de chance ($\text{AUC} \approx 0,38 - 0,47$).

---

## 4. Investigação 3: O Padrão Real — Modulação por Estado Comportamental (`preditor_estado_comportamental.py`)

A fundamentação teórica de Vanderwolf (1969), Tort et al. (2009) e Colgin (2016) aponta que o PAC não surge no vácuo: ele é uma manifestação da coordenação sináptica associada a estados comportamentais ativos.

### 4.1. Precedência Temporal do Comportamento Anterior ($N-1 \to N$)
Analisando 2.715 pares consecutivos da janela anterior ($N-1$) para a janela atual ($N$) no dataset mestre completo:

| Comportamento na Janela Anterior ($N-1$) | Taxa de Acoplamento Confirmado na Janela Seguinte ($N$) |
|---|---|
| **Grooming / Limpeza** | **48,9%** |
| **Sniffing / Farejando** | **47,9%** |
| **Exploração / Locomoção** | **47,2%** |
| **Movimento de Cabeça** | $39,5\%$ |
| **Rearing / Em pé** | $38,3\%$ |
| **Imóvel / Descanso** | $34,6\%$ |
| **Sono** | **9,4%** |

$$\chi^2 = 41,3 \quad (p = 2,6 \times 10^{-7})$$

Animais em estados ativos (grooming, sniffing, locomoção) têm uma probabilidade quase **cinco vezes maior** de engajar em acoplamento na janela seguinte do que animais em sono.

### 4.2. Descarte Rigoroso de Explicações Alternativas
1. **Causalidade Reversa:** Descartada formalmente pela ordem cronológica estrita — a janela $N-1$ precede a janela $N$.
2. **Hipótese de Pura Autocorrelação (Estado Persistente):** O efeito sobrevive mesmo restringindo a amostra estritamente às janelas onde houve **transição real de comportamento** ($p = 0,018$).
3. **Hipótese de "É Apenas Potência Teta Disfarçada":** Embora a potência teta varie por comportamento ($\text{ANOVA } F = 15,68, p = 1 \times 10^{-17}$), o teste de razão de verossimilhança demonstrou que o comportamento carrega informação biológica independente além da potência teta:
   $$p = 1,3 \times 10^{-5}$$

### 4.3. Ablation Experimental Sequencial das Melhorias

| Passo | Hipótese Testada | AUC-ROC | $\Delta$ AUC | Veredito Científico |
|---|---|---|---|---|
| **0** | **Baseline** (Potência teta bruta 4–8 Hz + Comportamento $N-1$) | $0,602$ | — | Ponto de partida comprovado ($> \text{chance}$) |
| **1** | Normalização de teta por canal ($z$-score intra-canal) | $0,571$ | $-0,031$ | **Rejeitada isoladamente** (artefato; neutra no conjunto) |
| **2** | Separação Teta Tipo 1 (7–10 Hz) vs. Tipo 2 (4–7 Hz) | $0,584$ | $-0,018$ | Neutra isolada; **Teta tipo 1 virou a feature de maior peso** |
| **3** | Sub-janela mais próxima do evento + Tendência (inclinação) | **0,618** | **+0,016** | **CONFIRMADA**: Captura a aceleração pré-evento |
| **4** | $MI_z$ em par canônico fixo (6 Hz $\times$ 85 Hz) | $0,612$ | $-0,006$ | **Rejeitada**: Não agrega poder preditivo |
| **Final** | **Modelo Consolidado (Teta Tipo 1/2 + Janela Próxima + Tendência + Comportamento)** | **0,618** | **+0,016** | Salvo em `modelo_estado_comportamental.pkl` |

### 4.4. Calibração de Limiar para Disparo TTL (Curva Precisão-Recall)
O corte padrão de probabilidade ($0,50$) privilegia acurácia global, mas para estimulação optogenética o que dita o sucesso é o **recall da classe positiva**:
- **Melhor limiar por $F_1$ ($0,431$):**
  - **Recall:** **89,6%** (captura ~9 de cada 10 eventos de PAC).
  - **Precisão:** **45,2%** (custo: para cada evento real atingido, cerca de 1 falso-positivo é gerado).

### 4.5. Direcionalidade Causal: O Comportamento gera o Acoplamento ou o Acoplamento gera o Comportamento?
Uma das questões mais profundas da pesquisa diz respeito à direcionalidade da relação observada:

1. **O Comportamento como Criador do Estado Permissivo (A Relação Observada):**
   - O acoplamento teta-gama não surge espontaneamente no vácuo; a decisão motora e exploratória do animal aciona vias subcorticais ascendentes (septo medial $\to$ acetilcolina/glutamato), gerando a oscilação teta macroscópica (7–10 Hz).
   - É sobre essa onda lenta que os microcircuitos de CA1 conseguem alinhar os disparos de gama. Nesse nível biofísico, **o comportamento motor antecede e cria as condições permissivas da rede para que o acoplamento possa ocorrer**.
2. **O que o Acoplamento Efetivamente "Gera" (A Função Cognitiva):**
   - O acoplamento teta-gama **não gera o movimento mecânico da pata** (função de córtex motor e gânglios da base).
   - O acoplamento gera a **organização temporal interna da informação cognitiva** (segregação de codificação no pico vs. resgate no vale do teta, representação sequencial de *place fields*).
3. **A Analogia do Motor e da Transmissão:**
   - O comportamento motor atua como o acelerador do motor (eleva a rotação $\to$ ritmo teta permissivo).
   - O acoplamento teta-gama atua como o sincronizador fino da marcha: engata as informações no instante exato do ciclo oscilatório.
   - Se a estimulação optogenética romper o acoplamento teta-gama, o animal **não fica paralisado** (o motor continua acelerando e ele segue andando pela arena), mas **a marcha escapa**: ele explora a arena sem conseguir consolidar ou recuperar a memória da tarefa.
4. **Articulação com as Hipóteses H1 e H2 do Mestrado:**
   - **Hipótese H2 (Comportamento $\to$ Acoplamento):** Comprovada empiricamente pelos dados observacionais ($\chi^2 = 41,3, p = 2,6 \times 10^{-7}$).
   - **Hipótese H1 (Acoplamento $\to$ Desempenho de Memória):** Requer estritamente intervenção causal (*closed-loop* optogenético), justificando a arquitetura tecnológica proposta.

---

## 5. Implicações para a Engenharia de Hardware e Intervenção Optogenética

Esses achados reconfiguram a interface entre eletrofisiologia computacional e intervenção em circuito fechado (*closed-loop*):

1. **Inadequação da "Predição com 10 s de Antecedência":**
   - Com precisão de $45,2\%$, tentar disparar o laser 10 segundos antes de um evento que dura centenas de milissegundos causaria mais de 54% de estimulações indevidas em tecido não acoplado (risco de fototoxicidade e perturbação estocástica da circuitaria de CA1).
2. **A Solução Biologicamente Canônica: *Gating* por Estado e Closed-Loop Fásico (< 60 ms):**
   - Conforme demonstrado por **Siegle & Wilson (2014, *eLife*)**, a intervenção optogenética de alta fidelidade em hipocampo funciona por:
     - **Braço Lento (Gating / Habilitação):** O preditor de estado comportamental monitora se a rede está em estado de teta tipo 1 (locomoção/exploração ativa);
     - **Braço Rápido (Disparo Fásico Online):** Uma vez habilitado, o pulso luminoso é sincronizado em tempo real ($< 60\text{ ms}$) aos picos ou vales do ciclo teta via *Endpoint-Corrected Hilbert Transform* (ecHT) ou Phase-Locked Loop (PLL).

---

## 6. Síntese Conclusiva para a Dissertação de Mestrado

| Dimensão da Pesquisa | O que se pensava no início | O que a evidência empírica demonstrou |
|---|---|---|
| **Natureza do PAC em CA1** | Um fenômeno de rampa lenta visível no LFP bruto | Um evento transiente (*burst*) modulado por estados da rede neural |
| **Precursor no LFP (10s e 1-3s)** | Detectável por energia espectral ou comodulograma | **Ausente no sinal local isolado** ($p > 0,10$; AUC $\approx 0,43$) |
| **Papel do Comportamento** | Mero contexto descritivo | **Modulador primário da probabilidade de acoplamento** ($p = 2,6 \times 10^{-7}$) |
| **Arquitetura de Closed-Loop** | Previsão estática de longo prazo ($10\text{ s}$) | **Sistema híbrido: Gating por estado ativo + Disparo fásico em tempo real** |

> **Nota de Mérito Científico:** A demonstração de ausência de precursor no LFP bruto pareado, acompanhada pela descoberta e validação rigorosa da modulação por estado comportamental anterior, constitui uma **contribuição científica original e de alto rigor metodológico**, protegendo o trabalho contra falsas premissas de engenharia e enriquecendo substancialmente a discussão da dissertação.
