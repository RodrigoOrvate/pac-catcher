# Scripts do SCRIPT central — explicação rápida

> Pasta `CONTEXTO` reúne documentação que não muda entre sessões e
> serve de referência para qualquer sessão futura do estudo
> (MTESC04/MTESC05 × NOCI/LAC). Nada de específico de sessão aqui
> — caminhos, offsets e vencedores entram por CLI/CSV, nunca no código.

## Visão geral do pipeline (passo 0 + 7 passos)

```
PASSO 0: exploracao_interativo.py (NOVO — PRIMEIRO, nao opcional)
  → LFP dos 3 .ns2 juntos (concatena_sessao); navegacao assíncrona
     (VisPy/Neuroglancer/Bqplot no Jupyter) em milhões de pontos
  → visual: espectrograma + PSD + LFP bruto navegavel
  → usuario marca instante interessante → carimbo (t_start, t_end, canal)
     enviado ao pipeline por CLI/CSV — nao hardcoded

triagem_pac.py          → etapa 1: varre TODAS as janelas (10 s / 5 s)
refina_candidatos.py    → etapa 2: FDR de janela + filtros de artefato
comodulogram.py         → etapa 3: mapas z...scoredos ±notch + FDR do mapa
diagnostico_janela.py   → etapa 5: 5 painéis + MI respiração vs θ×γ
robustez_parametros.py  → etapa 7 (validação): sweep de parâmetros + MVL
figura_apresentacao.py  → etapa 7 (opcional): figuras dos vencedores
agrega_resultados.py         → constrói o dataset mestre (merge por canal)
enriquece_dataset_mestre.py  → FOOOF v2 + portão de banda larga no mestre
audita_*.py             → auditorias pós...hoc (transientes, segmentos,
                          pegada espacial, held-out, harmônico...)
audita_janela.py        → orquestrador forense: as 4 auditorias
                          compatíveis (skew+transientes+footprint+
                          harmônico) numa só passada por canal/janela
pac_core/               → núcleo matemático compartilhado (refatoração
                          2026-09): io.py, filtering.py, pac_metrics.py.
                          ns2_utils.py e demais "shims" listados abaixo
                          reexportam daqui — ver seção dedicada.
```

**Refatoração 2026-09 (Fatia 1+2):** o pipeline passou por uma consolidação
para eliminar duplicação de código entre os scripts (filtro, notch,
leitura de arquivo, KL-MI+surrogates viviam copiados em 4-6 lugares cada).
Toda a matemática compartilhada agora mora em `pac_core/` (ver seção
dedicada mais abaixo) e os scripts antigos viraram wrappers finos —
**nenhum resultado numérico mudou** (cada migração foi validada com
paridade bit-exata contra dados reais antes/depois, `tests/
test_pac_metrics_parity.py` é o gate permanente). Se você rodou o
pipeline antes dessa data, não precisa rodar de novo — os CSVs já
gerados continuam válidos.

**Reorganização em pastas de etapa (2026-09, mesma rodada):** os nomes de
script acima são os mesmos, mas os arquivos fisicos moveram para
`pipeline/etapa<N>_<nome>/` (ver mapa completo e exemplos de CLI com
caminho completo na seção "Estrutura Canônica" abaixo, ou em
`README.md`). `pipeline/` e cada pasta de etapa viraram pacotes Python
de verdade (`__init__.py`), e todo import interno é qualificado
(`from pipeline.etapa1_triagem.triagem_pac import X`) em vez de
resolvido por `sys.path` implícito. `auditorias/` não mudou de lugar.

.........

## 1. `triagem_pac.py` — Varredura estatística (etapa 1)

**O que faz:** varre todas as janelas de 10 s (passo 5 s) de todos os
arquivos .ns2 de uma sessão, calcula o **KL...MI** (informação mútua)
entre a *fase* do theta (4–8 Hz) e o *envelope* do gamma (30–80 Hz),
e reporta um **z...score** e p...valor empírico — não o MI bruto.

**Por que z...score e não MI bruto:**
... MI bruto não tem escala (0,03 é muito ou pouco sem contexto);
... Theta hipocampal é assimétrico (sawtooth) — gera MI espúrio
  mesmo sem acoplamento (Kramer et al. 2008; Cole & Voytek 2017);
... Ruído motor/contaminação EMG infla o gamma exatamente nas
  frequências que queremos olhar.

**A nula:** 200 *surrogates* por **deslocamento circular** do envelope
de gamma (≥ 1 s de deslocamento, semente 42). Isso quebra a relação
temporal fase...amplitude mas preserva o espectro de cada sinal,
simulando "nenhum acoplamento". O z...score é (MI_obs − MI_surr_mean)
/ MI_surr_std.

**Proxy de artefato motor:** potência relativa em banda larga de alta
frequência (150–450 Hz), onde o EMG costuma vazar no LFP. Não
substitui EMG real — serve só para *marcar* candidatos suspeitos
de contaminação motora (desconfie se z alto + proxy alto no
mesmo candidato).

**CLI principal:**
```bash
python pipeline/etapa1_triagem/triagem_pac.py ......pasta "<sessao>/<BASAL>" \
    ......saida "<sessao>/RESULTADOS/resultados.csv"
# janela 10 s / passo 5 s / 200 surrogates / z_corte 3.0 (padrão)
```

**Saída:** `resultados.csv` — uma linha por (arquivo, canal, janela):
`janela_ini_s`, `janela_fim_s`, `mi_observado`, `mi_surrogate_media`,
`mi_surrogate_dp`, `z_score`, `p_empirico`, `proxy_artefato_motor`.

**Notas:**
... `......demo` roda um teste com dados sintéticos (theta com jitter de
  frequência + acoplamento genuíno na metade) para validar o corte
  por z antes de usar em dados reais.
... `filtra_sinal` (filtro butterworth bandpass, `filtfilt`) vem de
  `pac_core/filtering.py` — mesma implementação em todos os scripts
  que a usam, não mais cópias locais (ver seção "Núcleo Compartilhado").

.........

## 2. `refina_candidatos.py` — Refinamento FDR (etapa 2)

**O que faz:** pega os candidatos da triagem (z ≥ `......z_pre_filtro`,
default 2,0) e aplica filtros estatísticos mais rigorosos:

1. **p...valor paramétrico (Gama):** ajusta uma distribuição Gama aos
   200 surrogates de *cada janela* e tira o p...valor analítico da
   cauda. Resolve o problema do p...valor empírico travado em
   1/200 = 0,005 (resolução fina sem precisar de milhões de
   permutações). KL...MI é não...negativo e assimétrico → Gama é a
   família adequada.
2. **Correção FDR...Benjamini...Hochberg** sobre **todas** as janelas da
   triagem original (m = 5472 na sessão 1, 2832 na #4). Conservador:
   trata as que não foram refinadas como automaticamente não...
   significantes.
3. **Co...ocorrência entre canais:** quantos canais do mesmo arquivo
   e janela também passaram do corte. Número alto (= "ruído de modo
   comum", ex. 12–15 canais simultâneos) → provavelmente cabo/
   movimento/referência, não dipolo neural.
4. **Kurtose do sinal filtrado em gamma (30–80 Hz):** transientes
   musculares mantêm kurtose alta mesmo após o filtro; oscilação
   genuína tem kurtose baixa (~0). Detecta contaminação *dentro*
   da banda de amplitude — ao contrário do proxy antigo de 150–450 Hz.
5. **Checagem de saturação/clipping:** fração de amostras próximas
   do limite do int16 (32767). Sinal clipado = artefato.

**Veredito por linha:**
... `"Não significativo após correção FDR"` → não sobreviveu ao BH;
... `"Revisar: ..."` → significativo mas com alerta (kurtose alta /
  saturação / muitos canais simultâneos);
... `"Candidato robusto"` → significativo + sem alertas.

**CLI:**
```bash
python pipeline/etapa2_refinamento/refina_candidatos.py ......csv resultados.csv \
    ......pasta_ns2 "<sessao>/<BASAL>" \
    ......saida "<sessao>/RESULTADOS/resultados_refinados.csv"
# ......z_pre_filtro 2.0 ......n_surr 1000 ......fdr_q 0.05 (padrão)
```

**Saída:** `resultados_refinados.csv` com colunas: `arquivo`, `canal`,
`janela_ini_s`, `janela_fim_s`, `mi_observado`, `z_score_refinado`,
`p_analitico`, `n_canais_simultaneos`, `kurtose_gamma`,
`proxy_saturacao`, `proxy_artefato_motor_150_450hz`, `significativo_fdr`,
`veredito`.

**Nota:** `filtra_sinal` e `_mi_de_bin_idx` vêm de `pac_core/filtering.py`
e `pac_core/pac_metrics.py` respectivamente — núcleo único, importado
(não copiado) por praticamente todos os scripts do pipeline. `bh_fdr`
permanece local a cada script que a usa (não migrada para `pac_core`).

.........

## 3. `comodulogram.py` — Mapas z...scoredos (etapa 3)

**O que faz:** gera um **mapa de calor** (heatmap) do z...score do MI
para cada par (frequência de fase × frequência de amplitude),
antes e depois do notch 60 Hz. Opcionalmente aplica FDR por célula
do mapa e classifica o padrão.

**Dois modos:**
... **Lote (`......csv`):** lê `resultados_refinados.csv`, gera um PNG por
  candidato. Usa `......top_n N` para pegar só os N melhores por
  `z_score_refinado` (exploratório). Filtra por `......veredito_prefixo`.
... **Janela única (`......arquivo`):** exploração rápida de um trecho
  específico — `......canal`, `......inicio`, `......fim`.

**O que muda vs. a versão antiga:** MI agora é **z...scoredo** contra
surrogates (a mesma nula do triagem). Colormap `RdBu_r` centrado em
0 — azul = abaixo do acaso, vermelho = acima, branco = nada. Isso
permite distinguir acoplamento genuíno (pico focal no quadrante
θ×γ) de transientes (coluna inteira de 8 Hz, sem banda focal).

**Notch 60 Hz (`......notch 60`):** rejeita a rede elétrica e harmônicos
(iirnotch) antes de qualquer filtragem. Sempre rodar **com e sem**
notch e comparar: pico que cai > 1,5z com o notch é rede elétrica,
não acoplamento.

**FDR do mapa (`......fdr_q 0.05`):** Benjamini...Hochberg sobre as 275
células do mapa (fases 4–14 Hz × amplitudes 30–150 Hz, 5 Hz de passo).
Classifica cada janela em:
... **"concentrado em ΘΓ"** — ≥ 2 células sig, ≥ 50± dentro de ΘΓ
  (≥2 células sig, ≥ 50± dentro do quadrante theta...gamma) =
  acoplamento genuíno e estreito;
... **"esparso/fora de ΘΓ"** — poucas células significativas, fora do
  quadrante = transientes ritmados;
... **"nada sobrevive ao FDR"** — zero células significantes.

**Discriminador validado:** a janela rejeitada tem *mais* células
significantes mas menor ± em ΘΓ; acoplamento genuíno concentra.
Janelas de 60 Hz têm 0± em ΘΓ.

**Saída:** PNGs em `<saida_dir>/<arquivo>_<canal>_<ini>...<fim>s_zcomodo.png`
+ `resumo_comodulogramas.csv` com: `canal`, `janela_ini_s`,
`janela_fim_s`, `z_pico_theta_gamma`, `fase_pico_hz`, `amp_pico_hz`,
`mi_pico`, `png`, e (com FDR) `classe_fdr`, `n_sig_fdr`, etc.

**Funções...chave internas:** `calcula_comodulograma_z` (mapa z por
célula), `p_valores_por_celula` (p por célula via Gama),
`bh_fdr_mapa` (BH sobre o mapa), `resume_cluster_fdr` (classifica
o padrão), `z_pico_theta_gamma` (pico no quadrante θ×γ).

.........

## 4. `diagnostico_janela.py` — Figura diagnóstica (etapa 5)

**O que faz:** para um (arquivo, canal, janela) específico, plota 5
painéis alinhados no tempo + 4 números que discriminam acoplamento
genuíno de artefato de transiente ou respiração.

**Os 5 painéis:**
1. LFP bruto (com notch, se `......notch`);
2. Banda respiratória (0,5–3 Hz) — deflexões lentas grandes =
   artefato/potencial respiratório;
3. Theta (4–12 Hz) — visual;
4. Gamma (30–80 Hz) + envelope;
5. Espectrograma (STFT) com bandas θ e γ marcadas.

**Os 4 números:**
... **MI θ×γ z** — o acoplamento original (deve ser alto e estável);
... **MI resp×γ z** — se alto, a "fase" que organiza o gamma é
  respiratória (~1–2 Hz), não theta;
... **MI resp×θ z** — respiração modula o "theta"? (sobreposição
  sniffing 4–8 Hz ≈ theta);
... **Fator de crista do theta filtrado** — senoide pura ≈ 1,4;
  transientes / ondas agudas >> 2 (discrimina "oscilação" de
  "bursts").

**Limitação declarada:** o proxy respiratório (0,5–3 Hz) é **CEGO
para sniffing a 4–8 Hz**, que cai dentro da banda "theta". Para
desempate fino, o vídeo é o árbitro (etapa 4).

**CLI:**
```bash
python pipeline/auditorias/diagnostico_janela.py ......arquivo "<BASAL>/<arquivo>.ns2" \
    ......canal chan20 ......inicio 85 ......fim 95 ......notch 60
```

**Saída:** PNG diagonal + 4 linhas no console.

.........

## 5. `robustez_parametros.py` — Validação de robustez (etapa 7)

**O que faz:** testa se o acoplamento de cada *vencedor* sobrevive a
mudanças de parâmetro e a uma métrica alternativa — sem bins.

**O que sweepa para cada vencedor (no par de pico fixo):**
... **A. n_bins** ∈ {10, 12, 15, 18, 24, 30} — o parâmetro do KL...MI.
  Acoplamento genuíno não pode desaparecer com 12 ou 24 bins;
  o esperado é um **plató de z alto** em torno de n_bins canônico
  (±5 Hz), não imunidade aos extremos;
... **B. Largura do filtro de fase** ∈ {±0,7, ±1,0, ±1,5, ±2,5} Hz;
... **C. Largura do filtro de amplitude** ∈ {±2,5, ±3,5, ±5,0, ±7,5,
  ±10,0} Hz;
... **D. Mapa completo recomputado com n_bins=12 e n_bins=24** — o
  pico ΘΓ deve continuar no mesmo lugar (±1 Hz fase, ±5 Hz amp);
... **E. MVL (mean vector length, Canolty et al. 2006)** — métrica
  *sem bins*: |média(envelope · e^{i·fase})|. Confirma sem o
  binning do KL...MI. MVL z < z do KL...MI é normal (captura só o 1º
  momento, não a distribuição).

**Critério "ROBUSTO" (honesto):** z ≥ 3 em *todo* o sweep de n_bins +
pico do mapa estável entre n_bins=12 e n_bins=24 + MVL confirmando
(z ≥ 3). As larguras de filtro *não* exigem z ≥ 3 nos extremos —
filtro mais estreito que o evento de gamma (±2,5 Hz) perde o evento
por construção.

**Nota importante:** MVL fraco com MI alto *pode* acontecer (MI captura
a distribuição inteira, MVL só o 1º momento). Registrar ⚠, não
rejeitar sozinho — lição registrada na sessão #1 (chan30, MVL 1,9).

**CLI:**
```bash
python pipeline/etapa7_validacao/robustez_parametros.py ......pasta "<sessao>/<BASAL>" \
    ......vencedores "<sessao>/vencedores.csv" \
    ......resumo_fdr "<sessao>/RESULTADOS/comodulogramas_fdr/resumo_comodulogramas.csv" \
    ......saida_csv "<sessao>/RESULTADOS/robustez_parametros.csv"
```

**Saída:** console + `robustez_parametros.csv` com linhas:
`janela`, `canal`, `par_pico`, `teste`, `parametro`, `z`.

.........

## 6. `figura_apresentacao.py` — Figura dos vencedores (etapa 7, opcional)

**O que faz:** gera uma figura de apresentação para cada vencedor do
`vencedores.csv`. Cada figura tem 4 séries temporais + comodulograma
+ polar fase×amplitude (estilo Tort et al. 2010).

**Layout (4×2):**
... Esquerda (séries alinhadas, janela completa):
  1. LFP bruto (notch 60 Hz);
  2. theta filtrado no par de pico (±1 Hz);
  3. gamma filtrado no par de pico (±5 Hz) + envelope;
  4. espectrograma (STFT) com bandas θ/γ marcadas;
... Direita:
  5. comodulograma z...scoredo (recalculado aqui, com notch) com o
     par de pico marcado;
  6. distribuição polar fase×amplitude (18 bins), amplitude de γ
     normalizada por bin de fase do theta — o MI polar é reportado
     abaixo do círculo.

**Nota:** o comodulograma é recalculado aqui com o par de pico do
vencedor — garante consistência visual com os números da robustez.

**CLI:**
```bash
python pipeline/etapa7_validacao/figura_apresentacao.py ......pasta_ns2 "<sessao>/<BASAL>" \
    ......vencedores "<sessao>/vencedores.csv" \
    ......saida_dir "<sessao>/RESULTADOS/figuras"
```

**Saída:** `<saida_dir>/<rotulo>_<canal>.png` (nome do PNG = rótulo +
canal, conforme convenção).

**Depende de:** `vencedores.csv` (formato: `rotulo,arquivo,canal,
inicio_s,fim_s,fase_pico_hz,amp_pico_hz[,comportamento]`).

.........

## 0 (NOVO). `exploracao_interativo.py` — Passo 0: exploração visual em Jupyter

**Arquitetura:** navegador assíncrono para milhões de pontos,
integrado ao Jupyter. O pipeline canônico (etapas 1–7) varre
**cego** — processa todas as janelas e devolve números. O passo 0
inverte isso: o pesquisador **vê o sinal primeiro**, marca o que
interessa, e só então o pipeline processa. O pipeline não acha
acoplamento? Provavelmente a janela была errada — você viu o
acoplamento no LFP e precisa informar o carimbo.

**Tecnologia:** VisPy (CanvasSci, GPU...accelerated) ou Bqplot (d3.js
no browser) — lidam com milhões de pontos sem travar, zoom livre,
scroll, atualização em tempo real. Alternativa: Neuroglancer (se
houver voxels 3D) ou painel MNE...Python com TimeSeriesViewer.
**Não usar matplotlib estático para isso** — escala mal para registros
longos.

**Dados:** LFP dos 3 .ns2 juntos (já disponíveis em
`concatena_sessao`). A navegação é sobre o **registro completo** da
sessão, não sobre 1 min central.

**Interface esperada:**

```
┌─────────────────────────────────────────────────────────────┐
│  EXPLORAÇÃO — MTESC04 S1 (3 arquivos, 150 s, 32 canais)  │
├─────────────┬──────────────────────────────┬───────────────┤
│ CANAIS      │ LFP BRUTO (navegável, zoom)   │ ESPECTROGRAMA │
│ (scroll/    │ Ver teta (4...8 Hz) como ondulação│ (freq vs t) │
│  clique)    │ Ver gama (30...80 Hz) como    │   Teta ↑ quando│
│             │   bursts modulados pelo teta│   gama bursts │
│             │                              │               │
├─────────────┴──────────────────────────────┴───────────────┤
│  PSD (linear + log) por canal: ver picos em teta/gama    │
├─────────────────────────────────────────────────────────────┤
│  CONTROLES                                                  │
│  [Marcar instante] → t_start = 47s, t_end = 57s, ch=5    │
│  [Enviar para triagem] → executa triagem_pac.py nessa     │
│     janela + canal, devolve z e FDR no terminal            │
│  [Salvar carimbo] → adiciona linha ao candidatos.csv     │
│                                                             │
│  [Parar] → encerra, gera candidatos.csv com todos os     │
│     carimbos marcados                                      │
└─────────────────────────────────────────────────────────────┘
```

**Fluxo correto:**

```
1. open exploracao_interativo.ipynb (Jupyter Lab)
2. carrega 3 .ns2 → concatena (concatena_sessao)
3. navegacao: scroll pelo registro, zoom, escolha de canal
4. ve teta (4...8 Hz) e gama (30...80 Hz) juntos?
5. clica em "Marcar instante" → salva (t_start, t_end, canal)
6. (opcional) clica em "Enviar para triagem" → recebe z + FDR
   imediatamente na célula do Jupyter
7. repete 4...6 para todos os instantes interessantes
8. clica "Salvar carimbos" → gera candidatos.csv
9. fecha notebook → inicia pipeline canônico com candidatos.csv
```

**Saída:**

```
candidatos.csv
rotulo,t_start,t_end,canal,observacao
ep1_rearing,47,57,5,"teta forte em 8Hz, gama em 70Hz"
ep2_grooming,82,92,5,"teta moderada, bursts gama curtos"
ep3_walking,110,120,16,"teta 6Hz, gama 50Hz"
```

Este arquivo alimenta `triagem_pac.py` (etapa 1) como alternativa à
varredura cega — ele processa **só as janelas marcadas**, não toda a
sessão.

**CLI mínima (Jupyter):**

```bash
cd C:\acoplamento_theta...gamma\SCRIPT
jupyter lab
# abrir notebooks/exploracao_interativo.ipynb
```

**Depende de:** `ns2_utils.py` (leitura), VisPy ou Bqplot (instalar:
`pip install bqplot` ou `pip install vispy`). MNE...Python já está no
requirements.

**Status:** a ser implementado (script atual `comodulogram_interativo.py`
é o protótipo estático; reescrever como Jupyter + VisPy/Bqplot é o
próximo passo).

.........

## 7. `audita_transientes.py` — Auditoria de transientes

**O que faz:** testa a hipótese "espigão broadband → MI falso" nos
vencedores. O mecanismo de transientes broadband *não* é coberto
pela nula de deslocamento circular — o deslocamento quebra a
coincidência dos transientes, então transientes coincidentes
produzem z alto. Esta auditoria discrimina com 5 blocos:

1. **Contexto de amplitude:** RMS da janela vs. janelas de 10 s do
   próprio arquivo (percentil); rank do RMS entre os 32 canais na
   mesma janela; kurtose do bruto e da banda gamma (transiente →
   kurtose alta em ambas);
2. **Despike:** interpolar linearmente amostras > k·σ robusta
   (mediana/MAD da janela, guarda ±50 ms) e recomputar z com a
   *mesma* nula (mesmos deslocamentos, mesma semente). Se z colapsa,
   o acoplamento morava nos transientes; se sobrevive, não era
   dirigido por espigões;
3. **Sub...janelas de 2 s:** acoplamento genuíno aparece na maioria
   das sub...janelas; dirigido por 1–2 episódios concentra todo o z
   nelas;
4. **Controles:** chan20 @ 85–95 e chan10 @ 30–40 (bons, validados
   em comportamento + respiração) vs. chan32 @ 95–105 (artefato
   conhecido — janela rejeitada #3); o veredito sobre o caso
   disputado só vale se o teste separar esses dois grupos;
5. **Vizinhos:** z na célula do pico nos canais vizinhos imediatos
   — fonte cortical local tem gradiente; cabo/EMG difuso não.

**Saída:** `auditoria/auditoria_transientes.csv`, `auditoria/
contexto_amplitude.csv`, PNG por caso.

**Núcleo reutilizável:** `audita_transientes_de_sinal(lfp, fs, fp, fa, ...)`
encapsula os blocos 2-3 (baseline → despike k6/k5 → sub-janelas →
histograma de fase) para um sinal já em memória — é o que
`audita_janela.py` chama internamente, sem reler o arquivo.

**CLI (sessão 08/07 default):** `python audita_transientes.py`
**Nova sessão (via CLI, casos nunca no código):**
```bash
python pipeline/auditorias/audita_transientes.py ......pasta_ns2 "<sessao>/<BASAL>" \
    ......saida_dir "<sessao>/RESULTADOS/auditoria" \
    ......casos "rotulo1,arq1,chan1,ini1,fim1,fp1,fa1;rotulo2,..." \
    ......vizinhos "chan24:z,chan26:z,..."
```

.........

## 8. `audita_segmentos.py` — Localização temporal do acoplamento

**O que faz:** complementa a auditoria de transientes — localiza,
*no tempo*, onde mora o acoplamento de um caso. Para cada sub...
segmento, recomputa o mapa z completo e reporta o z na célula do
pico original e o pico ΘΓ do próprio segmento.

**Uso (default sessão 08/07):** `python audita_segmentos.py`
**Nova sessão:**
```bash
python pipeline/auditorias/audita_segmentos.py ......pasta "<sessao>/<BASAL>" \
    ......arquivo <arq>.ns2 ......canal chanXX \
    ......segmentos "ini1...fim1,ini2...fim2,..." ......fp F ......fa A
```

**Exemplo:** `......segmentos "20...26,26...30,28...30,20...30"` divide a janela s para ver onde o z aparece (ou não).

---

## 9. `audita_footprint.py` — Pegada espacial do acoplamento

**O que faz:** mede o z na célula do pico em **todos os 32 canais** da mesma janela. O discriminador: fonte cortical local produz *gradiente suave* (poucos canais vizinhos significantes); artefato difuso (respiração, movimento, volume conduzido, cabo) aparece simultaneamente em muitos canais distantes.

**Leitura:** o caso "rearing" (chan18/20/30/32, 4 canais) é referência do que é pegada *local*. O grooming (chan22 + cluster 22/24/26/28) é *focal*. Se 26 dos 32 canais estão em z ≥ 3, é incompatível com cabo/EMG difuso.

**Núcleo reutilizável:** `footprint_de_janela(dados, fs, nomes, ini, fim, fp, fa)` roda o loop pelos 32 canais sobre dados já carregados — `audita_janela.py` chama isso, não relê o arquivo. É a auditoria mais cara das 4 do orquestrador (32× o custo de `mi_z_par` para um único canal); `--pula footprint` descarta.

**Uso:** `python audita_footprint.py` (default: sessão 08/07, 3 janelas pré-configuradas). Para outra sessão, edite a lista `JANELAS` no topo do script ou passe alvos via linha de comando.

---

## 10. `audita_grooming_robustez.py` — Robustez do núcleo de grooming *(script não encontrado no repositório atual — histórico da investigação preservado abaixo, mas o arquivo em si parece ter sido removido/renomeado antes da refatoração 2026-09; não é parte do trabalho de consolidação)*

**O que faz:** verifica se o *núcleo* de grooming do vencedor 2 (003 @ 20–27 s, par 5×35 Hz, re-ancorado no vídeo) sobrevive ao mesmo batizado dos outros vencedores: sweep de n_bins + MVL no chan22 + consistência nos vizinhos chan24/26/28.

**Contexto:** o z=9,4 da janela original 20–30 s era inflação por escolha de janela (mistura de estados + nula estreita). O núcleo robusto é z ≈ 4 durante o grooming puro.

**Uso:** `python audita_grooming_robustez.py` (sessão 08/07 default).

---

## 11. `diagnostico_janela.py` — Diagnóstico visual e espectral (5 painéis)

**O que faz:** para qualquer janela de interesse, plota 5 painéis integrados:
1. LFP bruto com filtros de banda lenta e rápida
2. Densidade espectral de potência (PSD / Welch)
3. Espectrograma tempo-frequência (STFT)
4. Comodulograma de Fase-Amplitude com Z-score
5. Distribuição de fase polar (histograma circular de Tort)

---

## 12. `audita_harmonico.py` — Teste de razão harmônica Θ→Γ (auditoria)

**O que faz:** usa FOOOF (Kühn et al. 2026) para estimar `cf_teta` numa janela de contexto longa (45 s) e testa se `amp_pico ≈ n × cf_teta` em frequência (razão inteira, tolerância 10%) e em fase (PLV entre `n × phi_theta` e `phi_gamma`).

**Três dimensões independentes:**
1. **Frequência**: razão inteira `amp/fase ≈ n` (tolerância 10% de `cf_teta`)
2. **Forma de onda**: skewness do teta (reusado de `audita_skewness.py`)
3. **Fase**: PLV entre `n × phi_theta` e `phi_gamma` — discriminador mais forte

**Veredito:** `CLEAN` / `REVISAR_RAZAO_INTEIRA` / `REVISAR_FASE_TRAVADA` / `SUSPEITO_HARMONICO_FORTE` / `SEM_REFERENCIA_TETA` / `AMBIGUO_MULTIPLOS_N` / `REVISAR_TETA_ASSIMETRICO`

**Núcleo reutilizável:** `avalia_harmonico(sinal_ctx, sinal_cand, fs, amp_pico, skew=...)` encapsula a árvore de decisão inteira (FOOOF → razão harmônica → PLV → veredito). **Atenção ao dtype:** `sinal_ctx`/`sinal_cand` precisam manter o dtype bruto de `le_ns2` (`int16`, sem `.astype(float)` prévio) — convertê-los antes faz o ajuste do FOOOF divergir na 9ª casa decimal (achado da migração, não é diferença de precisão trivial). `audita_janela.py` respeita isso.

**CLI (congelada — `processa_sessao.py` chama este script em produção via subprocess, não mudar nome de flag nem coluna de saída):**
```bash
python pipeline/auditorias/audita_harmonico.py --csv "<sessao>/RESULTADOS/vencedores.csv"     --pasta_ns2 "<sessao>/<BASAL>"     --saida "<sessao>/RESULTADOS/harmonico.csv"     --janela_contexto_s 45 --modo_preprocesso hibrido --f_linha 60.0
```

---

## 13. `audita_harmonico_hfo.py` — Teste de razão harmônica G→HFO *(novo 2026-09-05)*

Testa se HFO (150-250 Hz) é harmônico de Gamma (30-80 Hz) usando FOOOF para estimar `cf_gamma` e PLV(n × phi_gamma, phi_hfo). Análogo ao `audita_harmonico.py`.
**Vereditos:** `CLEAN` / `REVISAR_RAZAO_INTEIRA` / `REVISAR_FASE_TRAVADA` / `SUSPEITO_HARMONICO_FORTE` / `SEM_REFERENCIA_GAMMA`.

---

## 14. Módulo de Análise e Anotação Comportamental *(novo 2026-09)*

Para correlacionar os episódios de acoplamento detectados com o comportamento real do animal registrado em vídeo (.MPG, .mp4), o pipeline conta com três ferramentas complementares:

### 14.1 `gerar_template_comportamento.py` — Extração de Janelas Exclusivas
- **O que faz:** Varre o CSV mestre de resultados e extrai uma lista única e cronológica de janelas de 10s onde houve detecção de PAC em pelo menos um canal da sessão.
- **Por que é essencial:** Quando múltiplos canais (ex.: hipocampo CA1, córtex ou estriado) detectam acoplamento na mesma janela temporal, o rato executou exatamente o mesmo comportamento. Agrupar por janela exclusiva evita ter que anotar redundantemente os mesmos 10 segundos várias vezes.
- **Saída:** `template_comportamento.csv` com colunas `sessao`, `condicao`, `arquivo`, `janela_ini_s`, `janela_fim_s`, `video_tempo_ini`, `video_tempo_fim`, `pares_detectados`, `n_canais_pac`, `comportamento`, `observacoes`.

### 14.2 `anotador_comportamento.py` — Aplicativo Gráfico de Anotação Sincronizada
- **O que faz:** Interface gráfica interativa construída com **Tkinter + OpenCV + Pillow** para visualização do vídeo em sincronia com o LFP.
- **Sincronização com Gravação Contínua e Offsets por Arquivo (`offsets_arquivos`):**
  - Frequentemente a câmera grava um único vídeo longo contínuo (ex: 20 minutos em `.MPG`), enquanto o sistema neural particiona a aquisição em blocos de ~5 minutos (`001.ns2`, `002.ns2`, `003.ns2`), ou sofre pequenas interrupções (como perda de pacotes).
  - O anotador suporta offsets específicos por arquivo armazenados em `config_anotador.json`. Ao trocar de linha na tabela entre arquivos, o offset correspondente é aplicado e o vídeo salta instantaneamente para o ponto temporal correto.
- **Recursos de Alta Produtividade:**
  - **Loop contínuo de 10s:** Repete continuamente a janela de análise para inspeção sem precisar arrastar a barra de progresso.
  - **Autoplay ao avançar:** Ao salvar ou avançar a linha, o vídeo salta para a nova janela e já começa a tocar.
  - **Atalhos Rápidos de 1 Toque:** Teclas numéricas `[1]` a `[8]` preenchem categorias padrão (`Imóvel / Descanso`, `Exploração / Locomoção`, `Grooming / Limpeza`, `Rearing / Em pé`, `Sniffing / Farejando`, etc.), salvam silenciosamente e avançam para a próxima janela.
  - **Digitação Livre e Autocompletar:** Campo de texto com autocompletar de termos já utilizados e confirmação por `Enter`.
- **Integridade de Codificação (`utf-8-sig`):**
  - O aplicativo grava o CSV utilizando UTF-8 com BOM (`utf-8-sig`), garantindo que acentos da língua portuguesa não sejam corrompidos ao abrir no Microsoft Excel no Windows.

### 14.3 `junta_comportamento.py` — Mesclagem com o Dataset Mestre
- **O que faz:** Combina as anotações feitas no `template_comportamento.csv` de volta ao dataset mestre de resultados (`dataset_mestre_final_v2.csv`), propagando o comportamento anotado para todos os canais correspondentes àquela janela temporal.

---

## 15. `audita_skewness.py` — Assimetria do Teta (dente-de-serra)

**O que faz:** filtra a banda teta (4-8 Hz, Butterworth ordem 4, `filtfilt`) e calcula o coeficiente padronizado de Fisher-Pearson (skewness amostral). `|skewness| > limiar` (default 0.5) = `SUSPECT (Asymmetric)` — teta não-senoidal que pode gerar harmônicos espúrios na banda de amplitude (mesmo mecanismo do `audita_harmonico.py`, mas julgando pela forma de onda em vez de FOOOF+PLV).

**Núcleo reutilizável:** `skewness_de_sinal(sinal, fs, fs_banda=(4,8), order=4) -> (skew, n)` opera em array já em memória; `theta_skewness_for_window(file_path, ...)` é o wrapper que lê o arquivo e chama isso (mantido porque `audita_harmonico.py` importa esse nome). `classifica_skew(skew, n, limiar=0.5)` centraliza o veredito do caso "janela completa" (o caso "ilha reancorada" tem um sufixo de texto ligeiramente diferente e foi deixado como estava, de propósito, para não mudar comportamento).

**CLI (congelada — `processa_sessao.py` chama este script em produção via subprocess):**
```bash
python pipeline/auditorias/audita_skewness.py --csv "<sessao>/RESULTADOS/vencedores.csv" \
    --pasta_ns2 "<sessao>/<BASAL>" --saida "<sessao>/RESULTADOS/skewness.csv" --limiar 0.5
```

**Saída:** `skewness.csv` com `rotulo, arquivo, canal, janela_ini_s, janela_fim_s, janela, janela_tipo, n, skewness, veredito_skew, motivo`.

---

## 16. `audita_janela.py` — Orquestrador forense por janela *(novo 2026-09)*

**O que faz:** para UM canal/janela, roda numa só passada as 4 auditorias que compartilham o mesmo formato de entrada — skewness, transientes/despike, pegada espacial (footprint) e razão harmônica — lendo o `.ns2` **uma única vez** e derivando em memória as 3 variantes de sinal que cada uma precisa (bruto sem notch p/ skewness, com notch 60Hz p/ transientes/footprint, janela de contexto de 45s p/ FOOOF). Hoje, rodar os 4 scripts individualmente relê o mesmo arquivo pelo menos 3-4 vezes.

**Papel diferente de `processa_sessao.py`:** aquele é o orquestrador **em lote** (por canal, sobre um CSV de candidatos, alimenta `agrega_resultados.py`). `audita_janela.py` é forense **por janela única**: "tenho um caso suspeito, me diga tudo sobre ele numa chamada só". Os dois coexistem sem sobreposição.

**Fica de fora (pré-requisitos estruturais incompatíveis com "canal + janela genérica"):**
- `audita_held_out.py` — exige janela reancorada (ilha vs. resto).
- `audita_segmentos.py` — exige sub-segmentos explícitos definidos à mão; recomputa o comodulograma completo (~275 células × 200 surrogates POR segmento), ordens de magnitude mais caro.
- `audita_harmonico_hfo.py` — CSV de triagem HFO e fonte de dados (`.mat`) distintos; já integrado à rota separada de `processa_sessao.py`.

**Custo:** skewness ~0, harmônico ~0 em MI (1 fit FOOOF), transientes 8× `mi_z_par`, footprint 32× (o mais caro — use `--pula footprint` para descartar).

**Veredito consolidado:** regra conservadora e transparente — começa em `LIMPO`, acumula um motivo por sinal de alerta, reusando os limiares que cada auditoria individual já usa (`--limiar` do skewness, `z≥3` do footprint, `--limiar_plv`/`--tol_rel` do harmônico). Dois heurísticos de triagem novos (não critério científico, documentados como tal): razão `z_despike_k6/z_baseline < 0.5` e `≤1` de 5 sub-janelas com `z≥3`.

**CLI:**
```bash
# janela única
python pipeline/auditorias/audita_janela.py --pasta_ns2 "<sessao>/<BASAL>" \
    --arquivo <arquivo>.ns2 --canal chan22 --inicio 20 --fim 30 --fp 5 --fa 35

# modo lote, reaproveitando vencedores.csv
python pipeline/auditorias/audita_janela.py --pasta_ns2 "<sessao>/<BASAL>" --csv "<sessao>/RESULTADOS/vencedores.csv"
```

**Saída:** `auditoria_janela.csv` (1 linha por caso, colunas prefixadas `skew_*`/`trans_*`/`ctx_*`/`foot_*`/`harm_*` + `veredito_consolidado`/`motivos`) e `auditoria_janela_footprint.csv` (long format, `caso,canal,z`).

**Validado:** rodando sobre o caso de referência `chan22@20-30s, par 5×35 Hz` (sessão 08/07), todos os valores batem exatamente com os 4 scripts individuais rodados separadamente sobre o mesmo caso.

---

## 17. `pac_core/` e a cadeia do Dataset Mestre *(refatoração 2026-09)*

### 17.1 Núcleo matemático compartilhado
Antes da refatoração, `filtra_sinal`, `aplica_notch`, `_mi_de_bin_idx` e a lógica de surrogates existiam copiados em 4-6 scripts cada. Agora moram em `SCRIPT/pac_core/`:
- **`io.py`**: `le_ns2`, `le_bin_legado`, `le_mat`, dispatcher `carrega_dados()` por extensão, `fatia_janela`, `concatena_sessao`, `salva_csv` (utf-8-sig opcional, não aplicado retroativamente aos 26 `to_csv()` existentes — só disponível para quem migrar deliberadamente).
- **`filtering.py`**: `filtra_sinal` (Butterworth `order=3` canônico) e `aplica_notch` (multi-harmônico; aceita `linha_hz=` como alias legado de `freqs_notch=`).
- **`pac_metrics.py`**: `_mi_de_bin_idx`, `fase_para_bin_idx`, `gera_deslocamentos` (primitiva própria — o comodulograma sorteia UM conjunto de deslocamentos compartilhado por ~400 células, não um por célula), `mi_surrogates_de_deslocamentos`, `z_score_mi`/`z_score_mi_mapa`, `calcula_mi_com_surrogates`. `rng` é sempre parâmetro explícito — o núcleo nunca semeia sozinho por padrão, cada script mantém sua própria política (reaproveitar o gerador entre janelas, recriar a cada iteração, seed fixa 42, etc.).

**Deliberadamente fora do núcleo** (documentado nos próprios módulos): as 6 variantes de filtro com `order`/clamp diferentes (`deteccao_ripple.py`, `audita_held_out.py`, `utils_harmonico.py::narrow_band`, etc.), `audita_held_out.py` (métrica de MI diferente — informação mútua via histograma 2D), `mvl_z_par`/`mvl_bruto_e_rayleigh` de `robustez_parametros.py` (MVL de Canolty, não KL-MI).

**Shims de compatibilidade** (mesma CLI/nome antigo, delegam para `pac_core`): `ns2_utils.py`, `atualiza_fooof_mestre.py`, `aplica_portao_banda_larga_mestre.py`. (`adapta_lfp_mat.py`, `testa_metodos_matlab.py` e `compara_preprocesso_linha.py` foram removidos em 2026-09 por não terem mais uso.)

**Teste de paridade:** `tests/test_pac_metrics_parity.py` — duas camadas (cópias congeladas das implementações legadas + valores numéricos literais capturados antes de qualquer migração). Toda mudança em `pac_core/pac_metrics.py` deve manter esse teste verde.

### 17.2 `agrega_resultados.py` → `enriquece_dataset_mestre.py`
Cadeia de consolidação do dataset mestre (antes disso não estava documentada em lugar nenhum):
```bash
# 1. Constrói o dataset mestre bruto: merge de refinados.csv + skewness/comodulograma/harmônico(_hfo)
#    por canal (renomeia veredito→veredito_refino, adiciona sessao/condicao/canal)
python pipeline/dataset_mestre/agrega_resultados.py --resultados "<sessao>/RESULTADOS" --saida dataset_mestre.csv

# 2. Enriquece: FOOOF v2 (aperiodic_mode='knee', ajuste particionado 2-45Hz teta / 35-250Hz gama)
#    + portão de banda larga (rebaixa "Candidato robusto" com suspeito_banda_larga=True)
python pipeline/dataset_mestre/enriquece_dataset_mestre.py --entrada dataset_mestre.csv --saida dataset_mestre_v2.csv
```
As duas etapas do passo 2 (`--etapas fooof portao`, ambas por padrão) são independentes — `portao` só usa colunas que já vêm do passo 1, nenhuma delas é criada pelo `fooof`. Não-destrutivo por padrão; `--in_place` sobrescreve a entrada (comportamento do antigo `aplica_portao_banda_larga_mestre.py`). Guarda contra reenriquecer: se as 14 colunas FOOOF v2 já existirem, aborta com erro claro a menos que `--forca` seja passado (evita colunas `_x`/`_y` duplicadas silenciosas do `pandas.merge`).

---

## Convenções comuns a todos os scripts

- **Notch de Alta Frequência (Desvio da Literatura)**: Kühn et al. (2026) reportaram na prosa de seu artigo a aplicação de Notch apenas em 50Hz. Contudo, o script oficial deles (`rem_noise.m`) demonstra que aplicavam múltiplos Notches em todos os harmônicos (50, 100, 150, 200 Hz). Nosso pipeline faz o equivalente para 60Hz (60, 120, 180, 240 Hz) para limpar a contaminação harmônica da rede elétrica nas bandas HG e HFO.
- **O Paradoxo HFO vs Ripple (Passo 0.5)**: A triagem de coocorrência resolve a extrema permissividade da banda HFO ampla. O HFO (150-250 Hz) está presente na maior parte das janelas (devido à cauda ruidosa/aperiódica), mas um autêntico "Ripple" é definido por características transitórias em resolução de amostra (≥3 DP acima da mediana, ≥10-25 ms de duração contínua). Janelas com HFO mas sem Ripple são rejeitadas.
- **Ambiguidade Harmônica no HFO**: Ao estender `n_max` para testar se HFO é harmônico de Theta (ex: 200Hz / 8Hz = n=25), múltiplos harmônicos poderiam se sobrepor à mesma banda em `tol=10%` se o n_max fosse global. Adotamos o distanciamento exato (`n` dinâmico por candidato) resolvendo os casos onde o plv pudesse retornar ambíguo erroneamente.
- **Nula de surrogates:** deslocamento circular >= 1 s, 200 repetições, semente 42.
- **Parâmetros canônicos:** janela 10 s / passo 5 s; theta 4-8 Hz x gamma 30-80 Hz; n_bins 18; 200 surrogates.
- **filtra_sinal**: Butterworth bandpass, filtfilt (zero fase), ordem 3.
- **_mi_de_bin_idx**: núcleo vetorizado do KL-MI via np.bincount.
- **bh_fdr**: Benjamini-Hochberg; m_total = família completa de testes.

## Princípio revisor

> Casos específicos de sessão (canais, janelas, offsets, vencedores)
> **NUNCA** entram hardcoded no código. Entram por CLI/CSV/JSON.

---

## Pipeline multi-acoplamento: Theta-Gamma / Theta-HG / Theta-HFO

| Par | Banda amplitude | Substrato | Estado |
|---|---|---|---|
| `theta_gamma` | 30–80 Hz | Fast gamma: CA3 $\rightarrow$ CA1 | Exploração |
| `theta_hg` | 80–150 Hz | High gamma: EC $\rightarrow$ CA1 | Misto |
| `theta_hfo` | 150–250 Hz | HFO/ripple: CA1 local | Repouso/SWR |

**Eficiência:** 3 pares × 59 janelas × 200 surrogates $\approx$ 8 s para 300 s de sinal, 1 canal.

```bash
# .mat direto (sem .ns2):
python pipeline/etapa1_triagem/triagem_pac_mat.py --mat DADOS_EXEMPLO_LFP_HG_HFO/LFP_HG_HFO.mat

# .ns2 multi-par:
python pipeline/etapa1_triagem/triagem_pac.py --pasta <sessao>/BASAL     --pares theta_gamma theta_hg theta_hfo     --saida <sessao>/RESULTADOS/resultados_triplo.csv
```

### Resultado validado no LFP_HG_HFO.mat (2026-09-05)

| Par | z mediana | z std | Candidatos (z >= 3) |
|---|---|---|---|
| `theta_gamma` | 6.93 | 2.96 | 58/59 |
| `theta_hg` | 12.39 | 4.75 | 57/59 |
| `theta_hfo` | 0.50 | 1.52 | 7/59 |

- `teta_ok`: 17/59 (29%), std=0.457 varia por janela (verificação dinâmica calculada dentro do loop).
- `ratio_hfo_gamma`: mediana=0.013 (HFO não é harmônico de Gamma).
