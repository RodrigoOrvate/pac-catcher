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

Abaixo estão as etapas principais do pipeline. Todo o código reside aqui. Você apenas apontará os comandos para as pastas com os seus dados `.ns2`.

> [!NOTE]
> Para uma explicação exaustiva, teórica e matemática de **TODOS** os scripts listados abaixo, consulte o arquivo **[`scripts_explicados.md`](./scripts_explicados.md)**. Ele é o verdadeiro manual técnico interno do pipeline.

## 📂 Estrutura do Repositório (Guia de Scripts)

O código é unificado e vive na pasta raiz (`SCRIPT/`). Abaixo, o mapa de ferramentas:

- **`pac_core/`**: Núcleo matemático compartilhado (extraído em 2026-09, ver `scripts_explicados.md` para o histórico da refatoração). Nenhum outro módulo deve reimplementar o que está aqui.
  - `io.py`: Leitura de `.ns2` (Blackrock), `.bin` legado e `.mat` (MATLAB/tetrodo), dispatcher único por extensão via `carrega_dados()`.
  - `filtering.py`: Filtro Butterworth passa-faixa (`filtra_sinal`) e notch multi-harmônico (`aplica_notch`) canônicos.
  - `pac_metrics.py`: KL-MI de Tort + surrogates por deslocamento circular (`calcula_mi_com_surrogates`, `z_score_mi`, `gera_deslocamentos`).
  - `pipeline/ns2_utils.py` (raiz de `pipeline/`) e `pipeline/dataset_mestre/{atualiza_fooof_mestre,aplica_portao_banda_larga_mestre}.py` são **shims finos** que re-exportam daqui/de `enriquece_dataset_mestre.py`, preservando CLIs antigas — novo código deve importar de `pac_core` diretamente.

- **`pipeline/`**: O núcleo duro do PAC Catcher, organizado em pastas de etapa (reorganização 2026-09 — ver `scripts_explicados.md` para detalhes da refatoração e a convenção de import qualificado usada):
  - **`processa_sessao.py`** e **`ns2_utils.py`** ficam na raiz de `pipeline/` (orquestrador principal e shim de I/O, não pertencem a uma etapa específica).
  - **`etapa0_exploracao/`**: `exploracao_minuto.py`, `comodulogram_interativo.py` — navegação visual e inspeção prévia dos dados antes da triagem cega.
  - **`etapa1_triagem/`** (Passo 1 + 0.5): `triagem_pac.py` / `triagem_pac_mat.py` (varredura inicial `.ns2`/`.mat` em busca de Teta-Gama, Teta-HG, Teta-HFO), `triagem_coocorrencia.py` (Passo 0.5, resolução amostral HFO/Ripple), `deteccao_ripple.py`, `preprocessa_referencia_diferencial.py`, `inspeciona_evento.py`.
  - **`etapa2_refinamento/`**: `refina_candidatos.py` — p-valor paramétrico (Gama), FDR de Benjamini-Hochberg, filtros de kurtose.
  - **`etapa3_comodulograma/`**: `comodulogram.py` — mapas de calor 2D (fase x amplitude) com filtros notch.
  - **`etapa7_validacao/`**: `robustez_parametros.py` & `figura_apresentacao.py` — estabilidade (varredura de n_bins, filtros) e figuras finais (STFT, polar plots).
  - **`dataset_mestre/`**: `agrega_resultados.py` → `enriquece_dataset_mestre.py` (constrói e depois enriquece o dataset mestre — FOOOF v2 + portão de banda larga, ver Passo 3.6) + os 2 shims `atualiza_fooof_mestre.py`/`aplica_portao_banda_larga_mestre.py`.
  - **`comportamento/`**: `gerar_template_comportamento.py` (template de janelas exclusivas), `anotador_comportamento.py` (GUI de sincronização com vídeo), `junta_comportamento.py` (mescla anotações ao dataset mestre).
  - **`utilitarios/`**: `extrair_picos.py`, `plot_basal_results.py`, `gerar_relatorio_pdf.py`, `gera_plot_fooof.py`.

- **`pipeline/auditorias/`**: Filtros e testes secundários rigorosos para falsos positivos.
  - `audita_janela.py`: **orquestrador forense** — roda skewness + transientes + footprint + harmônico numa só passada para um canal/janela (lê o `.ns2` uma única vez). Recomendado para investigar um caso específico; os 4 scripts abaixo continuam existindo individualmente (e são os que `processa_sessao.py` chama em lote).
  - `audita_harmonico.py` / `audita_harmonico_hfo.py`: Usa o FOOOF para separar 1/f e confirmar se os picos de amplitude não são harmônicos matemáticos da fase.
  - `audita_skewness.py`: Checa se a assimetria (dente-de-serra) da onda lenta forjou o acoplamento.
  - `audita_footprint.py`: Checa se a distribuição do acoplamento pelos 32 canais é focal (verdadeira) ou difusa (condução de volume/artefato).
  - `diagnostico_janela.py`: Cruza o ritmo Teta com faixas respiratórias/olfatórias (0.5-3Hz / 4-8Hz) para descartar artefatos respiratórios.
  - `audita_transientes.py`, `audita_segmentos.py`: Garantem que o acoplamento não é dirigido por *spikes* (espigões) de ruído mecânico.
  - `audita_held_out.py`: Testa double-dipping em janelas reancoradas (ilha vs. resto) — pré-requisito estrutural incompatível com `audita_janela.py`, fica separado.

- **`preditor/`**: Ferramentas experimentais de Machine Learning (`treinar_preditor.py`, `prever_pac_tempo_real.py`, `analisar_pre_evento.py`) para prever ocorrência de PAC em tempo real.
- **`tests/`**: Suite de testes automatizados (`test_synthetic_harmonico.py`, etc.) que simulam LFPs sintéticos ruidosos para garantir que a matemática do pipeline não falha sob *stress*.

---

### Pré-requisitos
Instale as bibliotecas necessárias:
```bash
pip install -r requirements.txt
```
*(Certifique-se de usar Python 3.9+ e de ter suas sessões `.ns2` organizadas nas pastas correspondentes).*

---

### Passo 0: Exploração Interativa (Opcional, mas recomendado)
Antes de rodar a varredura cega, você pode navegar pelo sinal bruto, STFT e PSD concatenado de toda a sessão para identificar visualmente eventos de interesse.
```bash
# Executado via Jupyter Notebook ou interface interativa
python pipeline/etapa0_exploracao/exploracao_minuto.py ...
```

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
   python pipeline/comportamento/gerar_template_comportamento.py --csv_mestre ../dataset_mestre_final_v2.csv --saida ../template_comportamento.csv
   ```

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
   python pipeline/comportamento/junta_comportamento.py --template ../template_comportamento.csv --csv_mestre ../dataset_mestre_final_v2.csv --saida ../dataset_mestre_final_comportamento.csv
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

---

### Passo 4: Auditorias Específicas
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

### Passo 5: Geração de Figuras Finais
Ao final, você seleciona os verdadeiros vencedores e passa para a geração das imagens (STFT, distribuição polar e LFP bruto) para apresentação, *lab meetings* ou artigo.
```bash
python pipeline/etapa7_validacao/figura_apresentacao.py \
    --vencedores "../SESSAO_EXEMPLO/RESULTADOS/vencedores.csv" \
    --pasta_ns2 "../SESSAO_EXEMPLO/Basal" \
    --saida_dir "../SESSAO_EXEMPLO/RESULTADOS/figuras"
```

*(Imagens de exemplo serão adicionadas aqui em breve para facilitar a visualização).*
