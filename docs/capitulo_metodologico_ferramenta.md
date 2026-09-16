# Capítulo Metodológico — ThetaGamma-Studio

*Rascunho para incorporação ao capítulo de Métodos da dissertação. Descreve
a ferramenta de software desenvolvida como parte da metodologia de
neuroinformática do estudo, não como um resultado científico em si.*

## 1. Motivação

A análise de acoplamento fase-amplitude (PAC) teta-gama conduzida neste
trabalho envolveu mais de 30 scripts independentes, cada um cobrindo uma
etapa do pipeline (triagem, refinamento, comodulograma, auditorias de
qualidade, dataset mestre, predição de estado). Essa fragmentação, comum
em pipelines de neurociência de dados construídos incrementalmente ao
longo de um mestrado, dificulta tanto a reprodução externa dos resultados
quanto a inspeção visual rápida de um candidato específico por terceiros
(orientador, banca, colaboradores). O ThetaGamma-Studio foi desenvolvido
para consolidar esse pipeline numa ferramenta única, sem alterar nenhuma
das métricas científicas já validadas.

## 2. Arquitetura

A ferramenta segue uma arquitetura híbrida em duas camadas:

1. **`pac_core`** (núcleo matemático) e **`pac_studio`** (API de
   consolidação): módulos Python puros, sem interface, que orquestram as
   funções já existentes de leitura de dados (`pac_core.io`), filtragem
   (`pac_core.filtering`) e cálculo de modulação de fase-amplitude
   (`pac_core.pac_metrics`). Nenhuma métrica foi recalculada com lógica
   nova — `pac_studio` apenas expõe, como chamadas de uma linha, sequências
   que antes existiam duplicadas manualmente em cada script do pipeline em
   lote (notch → filtro de fase/amplitude → transformada de Hilbert → MI
   com surrogates).
2. **`dashboard_desktop`** (programa desktop nativo, Tkinter com tema
   ttkbootstrap): camada de interface em três seções. **Pipeline** executa
   as etapas do processamento em lote — triagem, refinamento com correção
   FDR, comodulogramas, auditorias de qualidade e construção do dataset
   mestre — chamando os mesmos scripts usados pela linha de comando, com o
   registro de execução exibido ao vivo. **Análise de Dados** reúne a
   inspeção dos resultados: traçados brutos multicanal, comodulogramas
   interativos, checklist dos portões de qualidade, galeria dos candidatos
   vencedores, cruzamento comportamento-acoplamento e replay do preditor de
   estado. **Anotador de Vídeo** abre a ferramenta de anotação
   comportamental já existente no projeto. O programa roda inteiramente na
   máquina local, sem servidor, navegador, nuvem ou banco de dados externo.

## 3. Princípio de reuso e fidelidade científica

A decisão de arquitetura mais relevante do ponto de vista metodológico foi
a de **não reimplementar** nenhum cálculo já validado no pipeline em lote.
Toda métrica exibida na interface (MI de Tort, z-score contra surrogates,
vereditos de harmônico/skewness/transiente, o teste qui-quadrado de
transição de estado comportamental) é produzida pela mesma função Python
já usada na análise em lote — a interface só decide como e quando chamá-la
e como apresentar o resultado.

Essa decisão foi verificada, não apenas assumida: `tests/test_pac_studio.py`
inclui um caso de **paridade bit-a-bit** entre a API de conveniência e a
chamada manual equivalente via `pac_core` (mesmo padrão de verificação já
estabelecido em `tests/test_pac_metrics_parity.py`, usado para validar a
migração do núcleo matemático para módulo compartilhado). Um segundo caso
de teste reproduz numericamente a Hipótese H2 do estudo (comportamento →
acoplamento: χ²=41,3, p=2,6×10⁻⁷, n=2715 pares consecutivos) inteiramente
a partir da ferramenta, confirmando que a interface não introduz nenhum
desvio numérico em relação à análise original.

## 4. Limitações documentadas

Dois pontos são deliberadamente sinalizados na própria interface, em vez
de omitidos:

- O **preditor de estado comportamental** usado na aba de replay
  (AUC-ROC 0,618) é explicitamente rotulado como protótipo não validado
  para acionamento de hardware real — a aba de closed-loop é um replay
  offline, não um controle em tempo real.
- O **painel farmacológico** (Basal vs. Nociceptina vs. Lactato) degrada
  para uma mensagem informativa quando o dataset mestre ainda não tem as
  colunas de grupo/condição agregadas (situação em curso no momento da
  escrita deste capítulo, com o processamento da fase de infusão ainda em
  andamento) — a ferramenta nunca preenche esse vácuo com dado inventado.

## 5. Critério de usabilidade

O critério de aceite adotado para a interface foi que um membro de banca
sem formação em programação devesse conseguir, sem assistência, carregar
uma sessão, inspecionar os 32 canais e visualizar o comodulograma
correspondente. O `docs/manual_usuario.md` documenta esse fluxo passo a
passo e serve como material de apoio à defesa.
