# PAC Catcher 🎣

Bem-vindo ao **PAC Catcher**! Este repositório centraliza o pipeline computacional para detectar, auditar e validar **Acoplamento Fase-Amplitude (Phase-Amplitude Coupling - PAC)** em registros eletrofisiológicos (LFP).

## 🎯 Objetivo

O foco principal do pipeline é extrair e validar com rigor estatístico os acoplamentos **Teta-Gama (30-80 Hz)**, **Teta-High Gamma (HG, 80-150 Hz)** e **Teta-HFO (150-250 Hz)**. 

Devido à vasta quantidade de falsos positivos reportados na literatura de PAC (gerados por harmônicos de fase, transientes e ruído), este pipeline não se baseia apenas no tradicional *Modulation Index*. Ele cruza estatística de deslocamento circular, paramétrica, comportamento de rede espacial e decomposição espectral avançada para assegurar que o acoplamento detectado é um fenômeno neurofisiológico genuíno.

### 📚 Fundamentos e Citações

O desenvolvimento deste fluxo é fortemente embasado em literatura de referência:
- **Separação Aperiódica/Oscilatória:** Utilizamos o algoritmo **FOOOF (SpecParam)**. A implementação em `pipeline/auditorias/utils_harmonico.py` e `pipeline/auditorias/linha_noise_kuhn.py` é uma adaptação para LFP do método descrito por Kuhn et al. no artigo *"Aperiodicity in Mouse CA1 and DG Power Spectra"* (Kuhn et al., 2026) — o código-fonte original dos autores (MATLAB/Python, incluindo `rem_noise.m`) está disponível no repositório oficial: [LFP_FOOOF](https://github.com/huRashidy/LFP_FOOOF). Essa abordagem é essencial para High-Gamma e HFO, garantindo que a elevação do *spiking* (disparo neuronal) não seja confundida com oscilação real de banda estreita. **Nota:** O pacote `fooof` foi substituído upstream pelo `specparam`, mas nossa implementação adaptou o core original para as necessidades de LFP.
- **Dados de exemplo / validação sintética:** O dataset `LFP_HG_HFO.mat`, usado para validar `triagem_pac_mat.py`, tem origem no repositório [phase-amplitude-coupling](https://github.com/tortlab/phase-amplitude-coupling) (Tort Lab). Esses dados não são versionados neste repositório (`.gitignore` exclui `DADOS_EXEMPLO*/`) — baixe-os do repositório original se precisar reproduzir os testes que os usam. (`testa_metodos_matlab.py` também referencia um `Example_plv_modindex.mat` de terceiros, na mesma pasta — confirmar a origem exata antes de documentá-la aqui.)
- **Ambiguidade Estatística no HFO:** Ao testar harmônicos em ordens muito altas (ex: Teta 8Hz contra HFO 200Hz = $n=25$), o volume de comparações múltiplas gera uma taxa de Falso Positivo puramente estocástica de **~31%** (medida via Monte Carlo). Isso significa que 1 em cada 3 osciladores genuínos e aleatórios do HFO cairia na tolerância harmônica de Teta por puro acaso. Isso torna a validação pelo travamento de fase (**PLV**) absolutamente indispensável nas altas frequências.
  - **Limitação Conceitual do PLV**: O PLV alto descarta osciladores independentes (ruído) que caíram na janela de harmônico por acaso. No entanto, se um acoplamento Teta-HFO biológico for excepcionalmente forte e travado em fase, ele produzirá um PLV idêntico ao de um harmônico matemático. O pipeline sinaliza esses casos como `REVISAR_FASE_TRAVADA`, assumindo implicitamente que PLVs > 0.8 são artefatos matemáticos, mas a inspeção visual humana é necessária para descartar acoplamentos biológicos bizarramente robustos.
  - **Detecção de Ripples (150-250Hz)**: O critério clássico de 25ms pode rejeitar 100% dos eventos em gravações basais sem referenciamento diferencial. Taxas de detecção a 10ms mostraram-se indistinguíveis de ruído de banda larga (idênticas na exploração e imobilidade). Qualquer evento candidato em limiares menores (ex: 15-20ms) deve ser validado com `inspeciona_evento.py` em vez de ser aceito cegamente como SWR fisiológico, já que a taxa de eventos basais pode ser virtualmente nula dependendo da montagem do eletrodo.
    > **Nota Específica (Dataset MTESC04):** A feature de Teta-HFO / Ripple encontra-se arquivada em **standby permanente** para este animal/sessão. Testes exaustivos de referenciamento diferencial (incluindo exclusão *Leave-One-Out* de canais silenciosos) provaram estatisticamente e morfologicamente que os poucos candidatos detectados eram artefatos (sem sharp-wave coincidente e sem o envelope *spindle*). A conclusão isolada é que **a montagem física da sonda neste animal não captou a camada piramidal de CA1** com proximidade suficiente. **Esta limitação é de hardware/gravação, não do pipeline:** o detector de SWR, o teste de harmônico e a correção por PLV estão validados e perfeitamente funcionais. Se outra sessão ou animal tiver um posicionamento ideal na camada, o pipeline atual de Teta-HFO extrairá o PAC legitimamente sem necessidade de nenhuma alteração no código.
- **Modulação Teta-HG e Teta-HFO:** As bases e os métodos de distribuição de graus de fase (Tort's Modulation Index, Polar Plots) têm fortes raízes em trabalhos como *"Theta phase modulates multiple layer-specific oscillations in the CA1 region"* de Adriano Tort.

---

## 🚀 Como Rodar o Pipeline

Todo o código reside aqui (`SCRIPT/`). Há duas formas de rodar as mesmas etapas — por baixo, os scripts são os mesmos:

- **Programa desktop (ThetaGamma-Studio)** — formulários com log ao vivo, sem digitar comandos. Ver a seção *Programa desktop* no fim deste README.
- **Linha de comando** — os Passos 1 a 5 abaixo, apontando os scripts para as pastas com os seus dados `.ns2`.

> [!NOTE]
> Para uma explicação exaustiva, teórica e matemática de **TODOS** os scripts, consulte **[`docs/scripts_explicados.md`](./docs/scripts_explicados.md)**. Ele é o verdadeiro manual técnico interno do pipeline.

### Ordem real de execução

Os Passos estão numerados por tema; a dependência de dados é esta:

1. **Por sessão:** Passos 1 → 2 → 3 + auditorias de skewness/harmônico — ou tudo de uma vez com `processa_sessao.py` (ver *Atalho* abaixo).
2. **Dataset mestre:** Passo 3.6 — agrega todas as sessões e enriquece com FOOOF + portão de banda larga.
3. **Comportamento:** Passo 3.5 — template → anotação no vídeo → junta ao dataset mestre. Depende do dataset do Passo 3.6.
4. **Vencedores:** Passo 3.8 — consolidação final.
5. **Figuras e exploração:** Passos 4 e 5.

## 📂 Estrutura do Repositório (Guia de Scripts)

O código é unificado e vive na pasta raiz (`SCRIPT/`). Abaixo, o mapa de ferramentas:

- **`pac_core/`**: Núcleo matemático compartilhado (extraído em 2026-09, ver `docs/scripts_explicados.md` para o histórico da refatoração). Nenhum outro módulo deve reimplementar o que está aqui.
  - `io.py`: Leitura de `.ns2` (Blackrock), `.bin` legado e `.mat` (MATLAB/tetrodo), dispatcher único por extensão via `carrega_dados()`.
  - `filtering.py`: Filtro Butterworth passa-faixa (`filtra_sinal`) e notch multi-harmônico (`aplica_notch`) canônicos.
  - `pac_metrics.py`: KL-MI de Tort + surrogates por deslocamento circular (`calcula_mi_com_surrogates`, `z_score_mi`, `gera_deslocamentos`).
  - `workspace.py`: fonte única do caminho-base do workspace (`BASE_WORKSPACE`, `BASE_LAC_NOCI`, `BASE_RESULTADOS_MESTRADO`) — nunca faça hardcode de `D:\acoplamento_theta-gamma`; sobrescreva com a variável de ambiente `ACOPLAMENTO_BASE` se a pasta mudar de lugar.
  - `pipeline/ns2_utils.py` (raiz de `pipeline/`) e `pipeline/dataset_mestre/{atualiza_fooof_mestre,aplica_portao_banda_larga_mestre}.py` são **shims finos** que re-exportam daqui/de `enriquece_dataset_mestre.py`, preservando CLIs antigas — novo código deve importar de `pac_core` diretamente.

- **`pac_studio/`**: API de 1 linha sobre `pac_core` (`analisa_janela()`, `transicao_estado()`) — consolida a sequência notch → filtro → Hilbert → MI + surrogates, e a reconstrução dos pares comportamentais N-1 → N da Hipótese H2. Backend do programa desktop.

- **`dashboard_desktop/`**: o programa desktop ThetaGamma-Studio (Tkinter + ttkbootstrap) — ver a seção *Programa desktop* no fim deste README.

- **`pipeline/`**: O núcleo duro do PAC Catcher, organizado em pastas de etapa (reorganização 2026-09 — ver `docs/scripts_explicados.md` para detalhes da refatoração e a convenção de import qualificado usada):
  - **`processa_sessao.py`** e **`ns2_utils.py`** ficam na raiz de `pipeline/` (orquestrador principal e shim de I/O, não pertencem a uma etapa específica).
  - **`etapa1_triagem/`** (Passo 1 + 0.5): `triagem_pac.py` / `triagem_pac_mat.py` (varredura inicial `.ns2`/`.mat` em busca de Teta-Gama, Teta-HG, Teta-HFO), `triagem_coocorrencia.py` (Passo 0.5, resolução amostral HFO/Ripple), `deteccao_ripple.py`, `preprocessa_referencia_diferencial.py`, `inspeciona_evento.py`.
  - **`etapa2_refinamento/`**: `refina_candidatos.py` — p-valor paramétrico (Gama), FDR de Benjamini-Hochberg, filtros de kurtose.
  - **`etapa3_comodulograma/`**: `comodulogram.py` — mapas de calor 2D (fase x amplitude) com filtros notch.
  - **`etapa4_validacao/`** (Passo 4, ex-`etapa7_validacao`): `robustez_parametros.py` & `figura_apresentacao.py` — estabilidade (varredura de n_bins, filtros) e figuras finais (STFT, polar plots, FOOOF).
  - **`etapa5_exploracao/`** (Passo 5, ex-`etapa8_exploracao`): `exploracao_minuto.py`, `comodulogram_interativo.py` — dissecação interativa com FOOOF dos eventos campeões (notebook ou CLI).
  - **`dataset_mestre/`**: `agrega_resultados.py` → `enriquece_dataset_mestre.py` (constrói e depois enriquece o dataset mestre — FOOOF v2 + portão de banda larga, ver Passo 3.6) → `consolida_vencedores.py` (filtro final de vencedores, ver Passo 3.8) + os 2 shims `atualiza_fooof_mestre.py`/`aplica_portao_banda_larga_mestre.py`.
  - **`comportamento/`**: `gerar_template_comportamento.py` (template de janelas exclusivas), `anotador_comportamento.py` (GUI de sincronização com vídeo), `junta_comportamento.py` (mescla anotações ao dataset mestre).
  - **`utilitarios/`**: `extrair_picos.py`, `plot_basal_results.py`, `gerar_relatorio_pdf.py`, `gera_plot_fooof.py`.

- **`pipeline/auditorias/`** (Passo 3.7): Filtros e testes secundários rigorosos para falsos positivos.
  - `audita_janela.py`: **orquestrador forense** — roda skewness + transientes + footprint + harmônico numa só passada para um canal/janela (lê o `.ns2` uma única vez). Recomendado para investigar um caso específico; os 4 scripts abaixo continuam existindo individualmente (e são os que `processa_sessao.py` chama em lote).
  - `audita_harmonico.py` / `audita_harmonico_hfo.py`: Usa o FOOOF para separar 1/f e confirmar se os picos de amplitude não são harmônicos matemáticos da fase.
  - `audita_skewness.py`: Checa se a assimetria (dente-de-serra) da onda lenta forjou o acoplamento.
  - `audita_footprint.py`: Checa se a distribuição do acoplamento pelos 32 canais é focal (verdadeira) ou difusa (condução de volume/artefato).
  - `diagnostico_janela.py`: Cruza o ritmo Teta com faixas respiratórias/olfatórias (0.5-3Hz / 4-8Hz) para descartar artefatos respiratórios.
  - `audita_transientes.py`, `audita_segmentos.py`: Garantem que o acoplamento não é dirigido por *spikes* (espigões) de ruído mecânico.
  - `audita_held_out.py`: Testa double-dipping em janelas reancoradas (ilha vs. resto) — pré-requisito estrutural incompatível com `audita_janela.py`, fica separado.

- **`preditor/`**: Protótipo de previsão de PAC para closed-loop (`preditor_estado_comportamental.py`, `prever_pac_tempo_real.py`, `validar_preditor.py`, entre outros). **Não validado para acionar hardware** (AUC-ROC 0,618) — ver `preditor/README_preditor.md`.
- **`tests/`**: Scripts standalone — rode `python tests/<arquivo>.py` (rodar `pytest` não testa nada). `test_pac_metrics_parity.py` e `test_pac_studio.py` são os gates: paridade bit a bit do núcleo e da API com o pipeline em lote, saem com erro se falharem. Os demais simulam LFPs sintéticos (`test_synthetic_harmonico.py`, etc.) ou rodam em dados reais do MTESC04 (`test_integracao_pipeline.py`, `test_referencia_diferencial.py`, `sweep_ripple.py`) e imprimem o resultado para conferência. Tabela completa em `docs/scripts_explicados.md`, seção 23.
- **`docs/`**: Documentação — `scripts_explicados.md` (manual técnico de todos os scripts), `manual_usuario.md` (programa desktop), sínteses de investigação, decisões e o capítulo metodológico da ferramenta.

---

### Pré-requisitos
Instale as bibliotecas necessárias:
```bash
pip install -r requirements.txt
```
*(Certifique-se de usar Python 3.9+ e de ter suas sessões `.ns2` organizadas nas pastas correspondentes).*

---

### Atalho: sessão inteira de uma vez (`processa_sessao.py`)
É o que o processamento em lote usa: roda os Passos 1 → 3 e as auditorias de skewness/harmônico canal a canal, para uma pasta de sessão (ou de condição de infusão):
```bash
python pipeline/processa_sessao.py \
    --pasta "<LAC_NOCI>/MTESC05_NOCI/MTESC05 -- 2 - infusao - 09-07-2024/0h pos infusao - 30min" \
    --saida "<RESULTADOS_MESTRADO>/MTESC05_NOCI_2_09-07-2024" \
    --min_janelas 3
```
- **Estágio 0:** lê os primeiros 30 s do primeiro `.ns2`, descarta canais com amplitude morta e sinaliza (sem descartar) outliers de RMS.
- **Estágio 1:** `triagem_coocorrencia.py` conta, arquivo a arquivo, as janelas com teta + gama e teta + HG. O canal só segue para o pipeline pesado se uma das duas contagens chegar a `--min_janelas` (padrão 3). A co-ocorrência teta + HFO é registrada, mas **não** entra nesse corte.
- **Estágio 2:** triagem → refinamento → comodulograma (3 pares) → `audita_skewness.py` / `audita_harmonico.py` / `audita_harmonico_hfo.py`.

A saída fica em `<saida>/<nome da pasta>/chanN/`, mais um `resumo_canais.csv` com as contagens de co-ocorrência de todos os canais (útil para entender uma sessão com zero candidatos). Em `RESULTADOS_MESTRADO`, o basal fica em `basal/LAC/` e `basal/NOCI/`, e cada rodada de infusão em `<RATO>_NOCI_<N>_<DATA>/<condição>/`.

---

### Passo 1: Triagem Estatística
Varre todas as janelas do registro e calcula o Z-score do Modulation Index (KL-MI) contra *surrogates* de deslocamento circular, buscando os 3 pares simultaneamente.
```bash
python pipeline/etapa1_triagem/triagem_pac.py \
    --pasta "../SESSAO_EXEMPLO/Basal" \
    --pares theta_gamma theta_hg theta_hfo \
    --saida "../SESSAO_EXEMPLO/RESULTADOS/resultados.csv"
```
*Gera uma tabela bruta com Z-scores, p-valores empíricos e SNR espectral de todas as janelas (ex: janelas de 10s com sobreposição).*

---

### Passo 2: Refinamento e Correção FDR
Aplica o *False Discovery Rate* (FDR - Benjamini-Hochberg) em todas as janelas, descarta falsos positivos por co-ocorrência multicanal e kurtose na banda alta.
```bash
python pipeline/etapa2_refinamento/refina_candidatos.py \
    --csv "../SESSAO_EXEMPLO/RESULTADOS/resultados.csv" \
    --pasta_ns2 "../SESSAO_EXEMPLO/Basal" \
    --saida "../SESSAO_EXEMPLO/RESULTADOS/resultados_refinados.csv"
```
*Quem sobrevive a esse funil estatístico ganha o selo "Candidato robusto".*

---

### Passo 3: Comodulogramas e Controle de Linha
Gera os mapas bidimensionais (fase x amplitude) para os candidatos robustos. **Crucial:** aplica filtros *notch* em 60 Hz e 120 Hz para evitar que harmônicos da rede elétrica forjem acoplamento no High-Gamma.

> [!WARNING]
> **Frequência da Rede Elétrica:** O comando abaixo utiliza `--notch 60`, que é o padrão para a rede elétrica do Brasil e EUA. Se os seus dados foram coletados em um país que utiliza **50 Hz** (ex: Europa, parte da Ásia), você **DEVE** alterar este parâmetro para `--notch 50`, caso contrário, o ruído elétrico não será removido e poderá gerar falsos positivos no High-Gamma.

```bash
python pipeline/etapa3_comodulograma/comodulogram.py \
    --csv "../SESSAO_EXEMPLO/RESULTADOS/resultados_refinados.csv" \
    --pasta_ns2 "../SESSAO_EXEMPLO/Basal" \
    --saida_dir "../SESSAO_EXEMPLO/RESULTADOS/comodulogramas" \
    --notch 60 --fdr_q 0.05
```

---

### Passo 3.5: Anotação Comportamental Sincronizada com Vídeo
Para correlacionar os episódios de acoplamento detectados com o comportamento real do animal (exploração, sniffing, grooming, imobilidade/descanso, rearing):

1. **Gerar o Template de Janelas Exclusivas:**
   Agrupa as janelas onde houve detecção de PAC em qualquer canal, evitando anotações repetidas da mesma janela temporal:
   ```bash
   python pipeline/comportamento/gerar_template_comportamento.py
   # sem argumentos usa os defaults (pasta autocontida, reorg 2026-09):
   #   --csv_mestre resultados/dataset_mestre_final.csv
   #   --saida      pipeline/comportamento/template_comportamento.csv
   ```
   > [!WARNING]
   > Depende do `dataset_mestre_final.csv` do **Passo 3.6** — rode-o antes. O gerador lista **toda** janela do CSV que receber, inclusive as "Não significativo após FDR" (o agregador não filtra por veredito). Passe um CSV já filtrado para `veredito_refino == "Candidato robusto"`, senão o template ganha centenas de janelas de ruído estatístico para anotar à toa. O template atual foi gerado assim.

2. **Anotar via Interface Gráfica:**
   Abre o aplicativo gráfico que sincroniza o vídeo contínuo (`.MPG`, `.mp4`) com os blocos de gravação neural (`.ns2`):
   ```bash
   python pipeline/comportamento/anotador_comportamento.py
   ```
   - **Múltiplos Arquivos por Sessão:** Gerencia automaticamente sessões com múltiplos blocos (ex.: `001.ns2`, `002.ns2`, `003.ns2`), permitindo offsets independentes no `config_anotador.json` e compensando pausas ou gaps de gravação da máquina.
   - **Playback Ágil:** Toca cada janela de 10s em loop contínuo com reprodução automática ao avançar.
   - **Atalhos Rápidos:** Pressione as teclas de atalho numéricas `[1]` a `[8]` ou digite livremente e tecle `Enter` para salvar e pular para o próximo momento.
   - **Compatibilidade Excel:** Salva nativamente com codificação `utf-8-sig` (UTF-8 com BOM), garantindo integridade de acentos no Windows/Excel.

3. **Mesclar Anotações ao Dataset Mestre:**
   ```bash
   python pipeline/comportamento/junta_comportamento.py
   # sem argumentos usa os defaults:
   #   --mestre        resultados/dataset_mestre_final.csv
   #   --comportamento pipeline/comportamento/template_comportamento.csv
   #   --saida         resultados/dataset_mestre_COM_COMPORTAMENTO.csv
   ```

---

### Passo 3.6: Construção e Enriquecimento do Dataset Mestre
Depois que `refina_candidatos.py` gerou `refinados.csv` por canal (e as auditorias do Passo 4 já rodaram, se for o caso), duas etapas consolidam tudo numa tabela única:

```bash
# 1. Agrega refinados.csv + skewness/comodulograma/harmônico(_hfo) de cada canal
python pipeline/dataset_mestre/agrega_resultados.py --resultados "<sessao>/RESULTADOS" --saida "<sessao>/dataset_mestre.csv"

# 2. Enriquece: FOOOF v2 (aperiodic_mode='knee', ajuste particionado) + portão de banda larga
python pipeline/dataset_mestre/enriquece_dataset_mestre.py --entrada "<sessao>/dataset_mestre.csv" \
    --saida "<sessao>/dataset_mestre_v2.csv"
```
Por padrão `enriquece_dataset_mestre.py` roda as duas etapas (`--etapas fooof portao`) e **não sobrescreve a entrada** (`--in_place` reproduz o comportamento antigo, se precisar). Os comandos antigos (`atualiza_fooof_mestre.py`, `aplica_portao_banda_larga_mestre.py`) continuam funcionando idênticos, como atalhos finos para este script.

O agregador varre a pasta **recursivamente** e tira do caminho de cada canal as colunas `sessao` (a parte com o prefixo do rato, ex. `MTESC05_NOCI_2_09-07-2024`), `grupo` (`NOCI`/`LAC`/`VEH` — nas sessões do MTESCnn_LAC vale o que foi infundido, `-lac_hem...` ou `-veh_hem...` no nome da pasta) e `condicao` (`basal`, `0h_pos`, `1h_pos`, `2h_pos`). Pastas que começam com `_` (`_CONTAMINADO_nao_usar/`, `_LAC_basal_antigo_nao_usar/`) são ignoradas. As pastas de sessão do LAC levam o rato no nome desde 2026-09-14 (`MTESC03 -- Rodada-1-02-05-2024-lac_hemdir`; mapeamento em `LAC_NOCI/_renomeacoes_pastas_LAC.csv`). Lote do basal: `RESULTADOS_MESTRADO/roda_lote_mtesc.ps1 -ratos "MTESC0?_LAC" -listar` mostra o plano sem rodar.

---

### Passo 3.7: Auditorias Específicas
O pipeline conta com auditorias separadas para blindar os resultados contra falhas físicas e matemáticas do sinal. Para investigar um canal/janela específico de uma vez, use o orquestrador:
```bash
python pipeline/auditorias/audita_janela.py --pasta_ns2 "<sessao>/<BASAL>" \
    --arquivo <arquivo>.ns2 --canal chan22 --inicio 20 --fim 30 --fp 5 --fa 35
```
Ele roda as 4 auditorias abaixo numa só passada (lendo o `.ns2` uma única vez) e devolve um veredito consolidado. Cada uma também pode ser rodada isoladamente:
- **`audita_skewness.py`**: Avalia a assimetria (formato dente de serra) do Teta para evitar geração de harmônicos de fase.
- **`audita_harmonico.py` / `audita_harmonico_hfo.py`**: Usa o algoritmo **FOOOF** para verificar se o Gama é harmônico do Teta ou se o HFO é harmônico do Gama.
- **`audita_footprint.py`**: Garante que o gerador da oscilação é focal e descarta propagação por condução de volume (co-detecção massiva nos 32 canais).
- **`audita_transientes.py`**: Testa se o acoplamento é dirigido por espigões (*spikes*) de ruído mecânico via despike + sub-janelas.
- **`diagnostico_janela.py`**: Confirma se o que estamos chamando de Teta não é na verdade respiração/sniffing do roedor (ritmo olfatório).

`audita_held_out.py`, `audita_segmentos.py` e `audita_harmonico_hfo.py` ficam de fora do orquestrador (pré-requisitos estruturais próprios: janela reancorada, sub-segmentos manuais, e fonte de dados/CSV distinta, respectivamente) — continuam scripts separados.

---

### Passo 3.8: Consolidação dos Vencedores
Filtra `dataset_mestre_COM_COMPORTAMENTO.csv` (Passo 3.5) pelos portões estatístico, FOOOF, harmônico e comportamental de uma vez, gerando a tabela final que os Passos 4 e 5 consomem:
```bash
python pipeline/dataset_mestre/consolida_vencedores.py
# sem argumentos usa os defaults:
#   --entrada resultados/dataset_mestre_COM_COMPORTAMENTO.csv
#   --saida   resultados/candidatos_vencedores_consolidados.csv
```
Critérios aplicados por padrão:
- **Estatístico:** `veredito_refino == 'Candidato robusto'`.
- **FOOOF:** erro de ajuste < 0,15 e pico periódico real nas duas bandas (teta e gama). Knee válido **não** é exigido — é diagnóstico de identificabilidade do expoente 1/f, não de autenticidade do PAC; `--exigir_knee_valido` reativa o filtro antigo.
- **Harmônico:** `CLEAN` ou `REVISAR_RAZAO_INTEIRA` (coincidência numérica sem travamento de fase). Reprova todo o resto — evidência real de contaminação (PLV > 0,8, ambiguidade de ordem, teta assimétrico) ou sem referência de teta; `--exigir_harmonico_clean` reativa o filtro antigo.
- **Comportamento:** anotado e diferente de `Artefato / Cabo`.

Imprime uma síntese de eventos por rato e por comportamento (já colapsando pseudoreplicação espacial só para o relatório — o CSV de saída mantém a granularidade canal×par).

---

### Passo 4: Geração de Figuras Finais
Ao final, você seleciona os vencedores consolidados (Passo 3.8) e gera as imagens (STFT, distribuição polar, LFP bruto e painéis FOOOF teta/gama) para apresentação, *lab meetings* ou artigo.
```bash
python pipeline/etapa4_validacao/figura_apresentacao.py \
    --vencedores resultados/candidatos_vencedores_consolidados.csv \
    --pasta_ns2 "<sessao>/Basal antes da infusao" \
    --saida_dir resultados/figuras
```
`--canal` na tabela de vencedores segue a convenção 1-based do dataset mestre (mesma do Passo 5); `rotulo` é opcional (gerado automaticamente a partir de arquivo+canal+janela se ausente).

*(Imagens de exemplo serão adicionadas aqui em breve para facilitar a visualização).*

---

### Passo 5: Exploração e Dissecação Interativa dos Vencedores (Opcional)
Alternativa interativa ao Passo 4: em vez de gerar todas as figuras de uma vez, inspeciona visualmente qualquer evento de `candidatos_vencedores_consolidados.csv` com 4 painéis: LFP bruto + teta + gama (10s), FOOOF banda baixa (2-45 Hz, knee + pico de teta), FOOOF banda alta (35 Hz-~0,95×Nyquist, com limpeza de linha Kuhn) e o comodulograma fase×amplitude daquela janela.

Duas formas de usar a mesma lógica (`pipeline/etapa5_exploracao/comodulogram_interativo.py`):

**A) Notebook — dashboard com dropdowns Rato → Comportamento → Janela campeã:**
```bash
jupyter notebook pipeline/etapa5_exploracao/notebooks/exploracao_interativo.ipynb
```

**B) Linha de comando — gera direto uma janela específica, sem abrir GUI:**
```bash
python pipeline/etapa5_exploracao/comodulogram_interativo.py \
    --pasta_ns2 "<sessao>/Basal antes da infusao" \
    --canal 13 --par theta_hg \
    --saida resultados/figuras_campeoes/evento_x \
    --zoom_t_center 50
```

> **Atenção ao canal:** `--canal` usa a MESMA convenção 1-based da coluna `canal` do dataset mestre (linha com `canal=13` no CSV → `--canal 13` no comando); o script converte internamente para o índice 0-based do array de dados. Não confundir com o nome nativo do canal no `.ns2` (ex.: `chan26`), que só aparece no título da figura para conferência cruzada.

---

## 🖥️ Programa desktop (ThetaGamma-Studio)

Interface gráfica para rodar o pipeline e inspecionar os resultados sem digitar comandos. Não reimplementa nenhum cálculo: cada botão chama os mesmos scripts dos Passos acima (a paridade com o pipeline em lote é verificada por `tests/test_pac_studio.py`).

```bash
python dashboard_desktop/app.py
```

A janela abre maximizada, com três seções:

**Pipeline** — cada sub-aba é um formulário com os parâmetros reais do script e um log ao vivo; **Parar** encerra a etapa e os processos que ela disparou.

| Sub-aba | Script | Equivale a |
|---|---|---|
| Sessão Completa | `processa_sessao.py` | *Atalho* (sessão inteira) |
| Triagem | `triagem_pac.py` (botão **Demo** = teste sintético rápido) | Passo 1 |
| Refinamento e Filtros | `refina_candidatos.py` + gráfico da distribuição de vereditos | Passo 2 |
| Comodulograma (lote) | `comodulogram.py --csv` | Passo 3 |
| Auditorias | `audita_skewness.py`, `audita_harmonico.py`, `audita_harmonico_hfo.py` | Passo 3.7 |
| Agregação | `agrega_resultados.py` → `enriquece_dataset_mestre.py` → `junta_comportamento.py` → `consolida_vencedores.py` | Passos 3.6, 3.5 (junção) e 3.8 |

`audita_transientes.py` e `audita_footprint.py` ficam fora de propósito: a linha de comando delas é presa a uma sessão específica (casos passados como texto). As saídas padrão da Agregação gravam em `*_novo.csv` e toda etapa pergunta antes de sobrescrever um arquivo existente.

**Análise de Dados** — inspeção de resultados já processados: traçados LFP multicanal com filtros, comodulograma de uma janela, checklist dos portões de qualidade de um candidato, galeria dos 190 vencedores, comportamento (incluindo o χ² da Hipótese H2) e replay do preditor (offline, não aciona hardware).

**Anotador de Vídeo** — abre `anotador_comportamento.py` (Passo 3.5) numa janela separada e mostra o progresso de anotação por sessão.

Código em `dashboard_desktop/` (`app.py` monta as seções; `aba_pipeline.py`, `aba_analise.py`, `aba_anotador.py`; `runner.py` executa os scripts com log ao vivo; `estilo.py` tem o tema). Guia passo a passo: [`docs/manual_usuario.md`](./docs/manual_usuario.md).
