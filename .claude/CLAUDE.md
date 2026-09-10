# CLAUDE.md — SCRIPT central do projeto (acoplamento theta-gamma)

## O que é esta pasta

**Código único do grupo de estudos** `C:\acoplamento_theta-gamma\`,
compartilhado por TODOS os estudos e sessões (MTESC04/MTESC05 × NOCI/LAC e EXPLORACAO_OBJETOS).
Nada de específico de sessão mora aqui — nem dados brutos (.ns2, .mpg), nem saídas.
Rodar os comandos preferencialmente de dentro da pasta `SCRIPT/` ou apontando caminhos relativos ao workspace.

```
C:\acoplamento_theta-gamma\
├── SCRIPT\                          ← ESTA pasta (código central + pipeline/ + README + CLAUDE.md)
│   ├── pipeline/                    ← Núcleo do pipeline PAC, comodulogramas e módulo comportamental
│   │   ├── triagem_pac.py
│   │   ├── refina_candidatos.py
│   │   ├── comodulogram.py
│   │   ├── robustez_parametros.py
│   │   ├── figura_apresentacao.py
│   │   ├── anotador_comportamento.py  ← GUI de anotação de vídeo sincronizada com janelas PAC
│   │   ├── gerar_template_comportamento.py
│   │   ├── junta_comportamento.py
│   │   ├── ns2_utils.py              ← shim fino, reexporta de pac_core/io.py
│   │   ├── agrega_resultados.py     ← constrói o dataset mestre (merge por canal)
│   │   ├── enriquece_dataset_mestre.py  ← FOOOF v2 + portão de banda larga
│   │   └── auditorias/              ← Validação e filtros de artefatos
│   │       └── audita_janela.py     ← orquestrador forense (skew+transientes+footprint+harmônico)
│   ├── pac_core/                    ← Núcleo matemático compartilhado (io/filtering/pac_metrics)
│   ├── preditor/                    ← Previsão de PAC em tempo real / ML
│   └── tests/                       ← Testes unitários e dados sintéticos
├── LAC_NOCI\                        ← Dados brutos eletrofisiológicos (.ns2, .nev) e vídeos (.MPG)
│   ├── MTESC04_NOCI\
│   │   ├── MTESC04 -- 1 - infusao - 08-07-2024\
│   │   │   ├── Basal antes da infusao\  ← 3 arquivos .ns2/.nev (gravação 20240708-123605)
│   │   │   └── basal antes da infusão - duracao 15m.MPG (vídeo contínuo de 20 min)
│   │   └── MTESC04 -- 2..6 ...
│   └── MTESC05_NOCI\
├── EXPLORACAO_OBJETOS\              ← Experimentos de exploração de objetos
├── RESULTADOS_MESTRADO\             ← Saídas consolidadas de análises por sessão e canal
├── template_comportamento.csv       ← Planilha de anotação comportamental (utf-8-sig)
└── config_anotador.json             ← Configurações e offsets individuais de vídeo por arquivo
```

Regra de ouro: **sessões têm dados e vídeos; a pasta SCRIPT tem o código.**

---

## Módulo Comportamental e Sincronização Vídeo ↔ LFP

### 1. Ferramentas de Anotação (`SCRIPT/pipeline/`)
- **`gerar_template_comportamento.py`**:
  Varre o dataset mestre e extrai uma lista deduplicada de janelas temporais de 10s onde houve detecção de PAC em qualquer canal. Evita redundância de ter que anotar a mesma janela várias vezes quando múltiplos canais cruzam o limiar. Gera `template_comportamento.csv`.
- **`anotador_comportamento.py`**:
  Aplicativo gráfico (Tkinter + OpenCV + Pillow) para visualização e anotação produtiva.
  - Toca a janela de 10s em loop contínuo com autoplay ao avançar.
  - Atalhos rápidos numéricos `[1]` a `[8]` ou digitação livre com autocompletar e salvamento ao teclar `Enter`.
  - Transição de momentos e atualização automática de status (`⏳` pendente $\rightarrow$ `✔` anotado).
- **`junta_comportamento.py`**:
  Mescla as colunas `comportamento` e `observacoes` anotadas de volta ao dataset mestre de resultados.

### 2. Sincronia Vídeo ↔ Eletrofisiologia & Gaps de Gravação
- **Arquivos `.nev` e Ausência de Hardware TTL:**
  O arquivo `.nev` (*Neural Event*) não contém eventos na porta digital (`digital_input_port: 0 events`). A câmera e o sistema Blackrock foram acionados manualmente de forma independente.
- **Carimbos de Horário Real (`SYSTEMTIME`):**
  Os cabeçalhos binários dos arquivos `.ns2` e `.nev` registram o horário do relógio da máquina com precisão de milissegundos:
  - Exemplo Sessão 1 (`MTESC04 -- 1 - infusao - 08-07-2024`):
    - `001.ns2`: Início às 15:36:05 (durou 272s até interrupção por `Data Loss` registrado no `.nev` em $t=273s$).
    - `002.ns2`: Início às 15:41:07 (iniciou 302,5s após o início do 001 $\rightarrow$ gap de 30,4s sem LFP enquanto o vídeo continuou gravando).
    - `003.ns2`: Início às 15:46:08 (iniciou 301,0s após o início do 002 $\rightarrow$ transição automática contínua sem gap).
- **Offsets por Arquivo Neural (`offsets_arquivos`):**
  Como a filmadora gravou um arquivo contínuo de 20 minutos (`.MPG`), cada bloco neural de 5 minutos tem seu próprio marco temporal:
  - `config_anotador.json` armazena `offsets_arquivos` individualmente:
    ```json
    "offsets_arquivos": {
      "20240708-123605-001.ns2": 325.0,
      "20240708-123605-002.ns2": 627.5,
      "20240708-123605-003.ns2": 928.5
    }
    ```
  - O aplicativo recalcula e sincroniza `video_tempo_ini` e `video_tempo_fim` respeitando o offset do arquivo correspondente.

### 3. Regra Crítica de Codificação CSV (`utf-8-sig`)
- Todos os arquivos CSV de anotação comportamental (`template_comportamento.csv`, etc.) **DEVEM** ser salvos e lidos com codificação **`utf-8-sig`** (UTF-8 com BOM).
- O Excel no Windows falha ao abrir UTF-8 puro sem BOM, corrompendo caracteres portugueses como `Exploração`, `Locomoção`, `Imóvel`, `Atenção`. O uso de `utf-8-sig` garante compatibilidade nativa sem erros.

---

## Estrutura Canônica dos Scripts (`SCRIPT/`)

### `SCRIPT/pac_core/` — núcleo matemático compartilhado (refatoração 2026-09)
Extraído das cópias duplicadas que existiam em `triagem_pac.py`/`comodulogram.py`/`refina_candidatos.py`/`auditorias/*.py`. **Qualquer função de filtro, notch, leitura de arquivo ou KL-MI/surrogates nova deve ir aqui, nunca reimplementada num script individual.**
- **`io.py`**: `le_ns2`, `le_bin_legado`, `le_mat`, dispatcher único `carrega_dados()` por extensão (`.ns2`/`.bin`/`.dat`/`.mat`), `fatia_janela`, `concatena_sessao`, `salva_csv` (utf-8-sig opcional).
- **`filtering.py`**: `filtra_sinal` (Butterworth passa-faixa canônico, `order=3`) e `aplica_notch` (notch multi-harmônico). Variantes com `order`/clamp diferentes (`deteccao_ripple.py`, `audita_held_out.py`, `utils_harmonico.py::narrow_band`, etc.) foram deixadas **de propósito** fora daqui — unificá-las é decisão científica, não refatoração.
- **`pac_metrics.py`**: `_mi_de_bin_idx`, `calcula_mi_com_surrogates`, `gera_deslocamentos`, `z_score_mi`/`z_score_mi_mapa`. `rng` é sempre parâmetro explícito (nunca semeado internamente por padrão) — cada consumidor decide sua própria política de reprodutibilidade. `audita_held_out.py` (métrica de MI diferente) e `mvl_z_par`/`mvl_bruto_e_rayleigh` (MVL de Canolty, não KL-MI) ficam de fora, documentado no próprio módulo.
- Teste de paridade: `tests/test_pac_metrics_parity.py` (cópias congeladas + valores numéricos literais capturados antes da migração — qualquer mudança de comportamento no núcleo deve ser validada contra ele).
- **Padrão de shim**: `pipeline/ns2_utils.py`, `pipeline/atualiza_fooof_mestre.py`, `pipeline/aplica_portao_banda_larga_mestre.py` são wrappers finos que preservam CLI/nome antigo mas delegam para `pac_core`/`enriquece_dataset_mestre.py` — não remover, scripts/hábitos antigos dependem deles continuarem funcionando. (`adapta_lfp_mat.py` foi removido em 2026-09 por não ter mais uso — `triagem_pac_mat.py` já lê `.mat` direto via `pac_core.io.le_mat`.)

### `SCRIPT/pipeline/`
- **`triagem_pac.py`** / **`triagem_pac_mat.py`**: Varredura em janelas deslizantes (10s com passo 5s) calculando KL-MI e Z-score vs 200 surrogates para `theta_gamma` (30-80 Hz), `theta_hg` (80-150 Hz) e `theta_hfo` (150-250 Hz) — a versão `_mat` lê `.mat` de tetrodo/terceiros em vez de pasta `.ns2` (CLIs deliberadamente separadas, não fundidas — ver `pac_core/io.py`).
- **`refina_candidatos.py`**: Refinamento paramétrico com distribuição Gama, correção FDR Benjamini-Hochberg, análise de coocorrência multicanal e kurtose na banda alta.
- **`comodulogram.py`**: Mapas bidimensionais de calor (fase × amplitude) com filtros notch em 60 Hz e harmônicos (120, 180, 240 Hz). Hub histórico — vários módulos ainda importam `z_pico_par`/`calcula_comodulograma_z` daqui.
- **`robustez_parametros.py`**: Varredura de estabilidade de parâmetros (variação de n_bins, filtros e métrica MVL).
- **`figura_apresentacao.py`**: Plota STFT, polar plots e LFP bruto dos vencedores validados.
- **`agrega_resultados.py`**: Constrói o dataset mestre bruto (merge de `refinados.csv` + `skewness.csv`/`resumo_comodulogramas.csv`/`harmonico*.csv` por canal, renomeia `veredito`→`veredito_refino`).
- **`enriquece_dataset_mestre.py`**: Enriquece o dataset mestre com FOOOF v2 (etapa `fooof`) e o portão de banda larga (etapa `portao`) — não-destrutivo por padrão (`--in_place` p/ sobrescrever). Substitui a cadeia manual `atualiza_fooof_mestre.py` + `aplica_portao_banda_larga_mestre.py` (ambos viram shims).
- **`gerar_relatorio_pdf.py`**: Compila relatório consolidado do estudo a partir de `vencedores_consolidado.csv`.

### `SCRIPT/pipeline/auditorias/`
- **`audita_janela.py`**: Orquestrador forense — roda skewness + transientes + footprint + harmônico numa só passada para um canal/janela (cache de arquivo, deriva as 3 variantes de sinal que cada uma precisa). `audita_held_out.py`/`audita_segmentos.py`/`audita_harmonico_hfo.py` ficam de fora (pré-requisitos estruturais incompatíveis: janela reancorada, sub-segmentos manuais, CSV/fonte de dados distinta).
- **`audita_harmonico.py` / `audita_harmonico_hfo.py`**: Usa o algoritmo FOOOF para separar 1/f e confirmar se os picos de amplitude não são harmônicos matemáticos da fase. **CLI/schema de CSV congelados** — `processa_sessao.py` chama os dois em produção via subprocess.
- **`audita_skewness.py`**: Checa se a assimetria (dente-de-serra) da onda lenta forjou o acoplamento. **CLI/schema de CSV congelados** (idem, `processa_sessao.py` chama em produção).
- **`audita_footprint.py`**: Checa se a distribuição do acoplamento pelos 32 canais é focal (verdadeira) ou difusa (condução de volume/artefato de referência comum).
- **`diagnostico_janela.py`**: Cruza o ritmo Teta com faixas respiratórias/olfatórias (0.5–3 Hz / 4–8 Hz) para descartar artefatos de respiração/sniffing.
- **`audita_transientes.py`, `audita_segmentos.py`**: Garantem que o acoplamento não é dirigido por espigões transitórios (*spikes*) de ruído mecânico.
- **`audita_held_out.py`**: Testa double-dipping (viés de reancoragem) comparando a "ilha" reancorada contra o resto da janela — métrica de MI própria (histograma 2D), não usa `pac_core.pac_metrics` de propósito.

---

## Convenções do Estudo e Regras de Validação

1. **Notch 60 Hz em tudo**: Harmônicos da rede em 60, 120, 180, 240 Hz devem ser filtrados para não inflar artificialmente o High-Gamma.
2. **Janelas Mistas vs Puras**: Janelas com múltiplos estados comportamentais diluem ou inflam z-scores; priorizar momentos comportamentais puros na validação final.
3. **Pseudoreplicação Espacial**: Canais vizinhos com o mesmo pico na mesma janela representam **um único evento biológico**, não descobertas independentes.
4. **Sniffing vs Teta**: Sniffing na faixa de 4–8 Hz não é descartado cegamente; a pegada espacial (`audita_footprint.py`) + confirmação no vídeo definem se a oscilação é focal neural ou artefato olfatório/motor.
5. **Nula de Surrogates**: Deslocamento circular $\ge 1\text{ s}$, 200 permutações, semente fixa 42 para reproducibilidade.
6. **Reuso de `pac_core/`**: filtro, notch, leitura de arquivo e KL-MI/surrogates já existem em `pac_core/`. Antes de escrever uma função nova desse tipo, verificar se já existe lá. Migração de um consumidor para `pac_core` **sempre** exige teste de paridade numérica (bit-exata, `np.testing.assert_array_equal`, nunca `assert_allclose`) contra dados reais antes/depois — ver `tests/test_pac_metrics_parity.py` como referência de formato.
