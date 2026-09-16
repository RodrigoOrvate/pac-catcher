# Contexto Estratégico: Da Validação do Preditor à Criação da Ferramenta

**Data:** 13 de Setembro de 2026  
**Documento de Alinhamento e Transição de Fase**  
**Projeto:** Mestrado em Eletrofisiologia / Acoplamento Teta-Gama e Circuito Fechado  

---

## 1. Fechamento da Fase de Investigação do Preditor

Conceitualmente e estatisticamente, o núcleo desta fase está completo, validado e esgotado com o máximo de rigor que os dados permitem:

1. **Ausência de Precursor Eletrofisiológico Local no LFP Bruto:**
   - **Janelas Longas (10 s):** LOOCV honesto contra controles reais resultou em acurácia de ~0,43 (nível de chance), com $MI_z$ oficial idêntico entre pré-evento e controle ($-0,157$ vs. $-0,155$, $p > 0,10$).
   - **Janelas Ultracurtas (1 s a 3 s) e Onset:** Teste em 186 pares de janelas demonstrou elevação fisiológica sutil de teta tipo 1 pareada intra-sessão ($p = 0,0247$ na inclinação da trajetória), mas sem capacidade discriminatória multivariada entre sessões e animais (AUC $\approx 0,42$).
   - **Conclusão Biológica:** O acoplamento teta-gama em CA1 comporta-se como um evento transiente (*burst*) de rede, e não como uma oscilação que "carrega" lentamente ao longo de segundos no sinal elétrico local.

2. **O Precursor Real está no Estado Comportamental:**
   - A precedência temporal do comportamento na janela $N-1$ prediz a ocorrência de PAC na janela $N$ com $\chi^2 = 41,3$ ($p = 2,6 \times 10^{-7}$).
   - Estados ativos (Grooming 48,9%, Sniffing 47,9%, Locomoção 47,2%) quase quintuplicam a chance de acoplamento em relação ao Sono (9,4%).
   - O efeito sobrevive a transições reais de comportamento ($p = 0,018$) e ao controle direto pela potência teta ($p = 1,3 \times 10^{-5}$).
   - No *ablation*, a potência teta rápida tipo 1 (7–10 Hz, Vanderwolf) com dinâmica de janela próxima + tendência levou o modelo a **AUC-ROC 0,618**, com limiar calibrado para recall de ~90% (precisão ~45%).

3. **Arremates Executados:**
   - **Figuras Científicas:** 4 gráficos em alta resolução gerados em `docs/figuras/` (curvas ROC, taxas por comportamento, precisão-recall e dinâmica de teta).
   - **Replay em Sessão Real:** Executado com sucesso no arquivo `20240708-123605-001.ns2` (canal 30, 272 s).
   - **Síntese Teórica:** Consolidada em `docs/sintese_investigacao_preditor_pac.md`.
   - **Versionamento:** Commit `d7b21bc` com árvore de trabalho 100% limpa.

---

## 2. A Pergunta Científica do Aluno: Farmacologia e Escopo *In Vivo*

> **Pergunta levantada pelo autor:**
> *"Se eu pegar os dados da infusão (altera o comportamento, como nociceptina que reduz o walking e vai reduzir os acoplamentos também, é o esperado) vai comprovar mais algo e se eu precisaria fazer algo in vivo para confirmar mais algo? Tentar ainda prever o acoplamento (ou de fato não dá) ou prever o comportamento?"*

### 2.1. O Impacto dos Dados da Infusão (Nociceptina)
- **No Basal (Dado Observacional):** A relação locomoção $\to$ acoplamento é correlacional. O animal escolhe espontaneamente andar ou parar.
- **Na Infusão de Nociceptina (Perturbação Exógena / Causal):** A nociceptina atua nos receptores NOP em CA1/septo, reduzindo a excitabilidade neuronal e suprimindo o *walking*. Se a taxa de acoplamento cair na mesma proporção, obtém-se uma **comprovação causal e reversível da Hipótese H2 (Comportamento $\to$ Acoplamento)**.
- **Validação Cruzada do Preditor:** Aplicar o modelo de estado (treinado no basal) nas sessões de infusão permite atestar sua generalização biológica e farmacológica.

### 2.2. A Necessidade de Experimentos *In Vivo* Adicionais
- **Para o Mestrado: NÃO é necessário realizar novos experimentos cirúrgicos.**
  - O banco de dados já possui eletrofisiologia multicanal *in vivo* (32 canais), vídeo sincronizado com anotação comportamental quadro a quadro, e manipulação farmacológica real (Nociceptina vs. Lactato).
  - Um experimento cirúrgico completo de optogenética em malha fechada (*in vivo closed-loop* com AAV, canulação, implante crônico de fibra óptica + matriz de eletrodos e latência estrita $< 60$ ms) é escopo de **Doutorado inteiro** (2 a 3 anos).
  - No pré-projeto de mestrado, a **Fase 1 era a "Garantida"** e a **Fase 2 era "Em Exploração / Proposta"**. O software e a modelagem teórica desenvolvidos cobrem integralmente o escopo proposto.

### 2.3. Prever Acoplamento vs. Prever Comportamento
- **Prever o acoplamento determinístico 10s antes no LFP bruto:** De fato **não é viável**, pois decorre de uma limitação biofísica da rede neural (PAC é um *burst* discreto).
- **Prever comportamento:** Para fins de intervenção em malha fechada, não é necessário prever o comportamento futuro; basta a detecção *on-line* instantânea do estado ativo (a potência teta rápida de 7–10 Hz atua como *proxy* imediato).

---

## 3. Decisão Estratégica: Delimitação de Escopo e Desenvolvimento da Ferramenta

> **Diretriz definida pelo autor:**
> *"Se isso de fato é apenas coisa de doutorado, não dá tempo de fazer no mestrado, então vamos limitar isso mesmo. O que podemos fazer é iniciar uma confecção de uma ferramenta. Todo esse processo gerado em uma ferramenta. A gente meio que responde uma hipótese e como bônus gera uma ferramenta."*

### O Valor Científico e Acadêmico dessa Abordagem:
1. **Responde formalmente a uma hipótese biológica:** Comprova a **Hipótese H2** (o acoplamento teta-gama é governado pelo estado comportamental e motor da rede neural, e não um oscilador elétrico local independente).
2. **Entrega um produto tangível de engenharia de software biomédico:** Transforma scripts dispersos em uma plataforma integrada (*Neuroinformatics Toolbox*), abrindo caminho para publicação em periódicos de métodos (*Journal of Neuroscience Methods*, *Frontiers in Neuroinformatics*, *Neuroinformatics*, *SoftwareX*).
3. **Deixa um legado tecnológico para o laboratório:** Novos integrantes do grupo poderão executar triagens, purificação de harmônicos e análises de circuito fechado com interface amigável e reprodutível.

---

## 4. Arquitetura Proposta para a Ferramenta (`ThetaGamma-Suite` / `PAC-Studio`)

### Estrutura Modular dos 4 Componentes:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                            THETA-GAMMA SUITE                                │
├─────────────────────────┬─────────────────────────┬─────────────────────────┤
│  Módulo 1: Ingestão     │  Módulo 2: PAC & Refino │  Módulo 3: Comportamento│
│  • Leitor Blackrock .ns2│  • Comodulograma fase×amp│ • Sincronização vídeo  │
│  • Filtros Notch/Banda  │  • Tort MI, PLV, FOOOF  │ • Matriz de transição   │
│  • Inspeção 32 canais   │  • Auditoria harmônicos │ • Teste Qui-quadrado    │
├─────────────────────────┴─────────────────────────┴─────────────────────────┤
│                     Módulo 4: Circuito Fechado (Closed-Loop)                │
│  • Gating por Estado Ativo (Teta Rápido 7-10 Hz)                            │
│  • Replay em Tempo Real de Sessão com Indicador de Disparo TTL              │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Formato Híbrido (O Padrão Ouro de Engenharia):
- **Camada de Base (Backend):** O núcleo dos algoritmos consolidado como biblioteca modular e testável em Python (`pac_core` e `pipeline`).
- **Camada de Apresentação (Frontend):** programa desktop local ThetaGamma-Studio (`dashboard_desktop/`, Tkinter + ttkbootstrap) que roda as etapas do pipeline (triagem, refinamento, comodulograma, auditorias, agregação) com log ao vivo e reúne a exploração de sessões, comodulogramas, cruzamentos comportamentais, anotador de vídeo e simulação do laço fechado com disparo TTL; relatórios gerados por `pipeline/utilitarios/gerar_relatorio_pdf.py`.
