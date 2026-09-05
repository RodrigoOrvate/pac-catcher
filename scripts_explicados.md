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
comodulogram.py         → etapa 3: mapas z-scoredos ±notch + FDR do mapa
diagnostico_janela.py   → etapa 5: 5 painéis + MI respiração vs θ×γ
robustez_parametros.py  → etapa 7 (validação): sweep de parâmetros + MVL
figura_apresentacao.py  → etapa 7 (opcional): figuras dos vencedores
audita_*.py             → auditorias pós-hoc (transientes, segmentos,
                          pegada espacial, robustez de grooming)
ns2_utils.py            → utilitários compartilhados (leitura .ns2 + fatia)
```

---

## 1. `triagem_pac.py` — Varredura estatística (etapa 1)

**O que faz:** varre todas as janelas de 10 s (passo 5 s) de todos os
arquivos .ns2 de uma sessão, calcula o **KL-MI** (informação mútua)
entre a *fase* do theta (4–8 Hz) e o *envelope* do gamma (30–80 Hz),
e reporta um **z-score** e p-valor empírico — não o MI bruto.

**Por que z-score e não MI bruto:**
- MI bruto não tem escala (0,03 é muito ou pouco sem contexto);
- Theta hipocampal é assimétrico (sawtooth) — gera MI espúrio
  mesmo sem acoplamento (Kramer et al. 2008; Cole & Voytek 2017);
- Ruído motor/contaminação EMG infla o gamma exatamente nas
  frequências que queremos olhar.

**A nula:** 200 *surrogates* por **deslocamento circular** do envelope
de gamma (≥ 1 s de deslocamento, semente 42). Isso quebra a relação
temporal fase-amplitude mas preserva o espectro de cada sinal,
simulando "nenhum acoplamento". O z-score é (MI_obs − MI_surr_mean)
/ MI_surr_std.

**Proxy de artefato motor:** potência relativa em banda larga de alta
frequência (150–450 Hz), onde o EMG costuma vazar no LFP. Não
substitui EMG real — serve só para *marcar* candidatos suspeitos
de contaminação motora (desconfie se z alto + proxy alto no
mesmo candidato).

**CLI principal:**
```bash
python triagem_pac.py --pasta "<sessao>/<BASAL>" \
    --saida "<sessao>/RESULTADOS/resultados.csv"
# janela 10 s / passo 5 s / 200 surrogates / z_corte 3.0 (padrão)
```

**Saída:** `resultados.csv` — uma linha por (arquivo, canal, janela):
`janela_ini_s`, `janela_fim_s`, `mi_observado`, `mi_surrogate_media`,
`mi_surrogate_dp`, `z_score`, `p_empirico`, `proxy_artefato_motor`.

**Notas:**
- `--demo` roda um teste com dados sintéticos (theta com jitter de
  frequência + acoplamento genuíno na metade) para validar o corte
  por z antes de usar em dados reais.
- A variável `filtra_sinal` (filtro butterworth bandpass, `filtfilt`)
  é interna e idêntica à dos outros scripts.

---

## 2. `refina_candidatos.py` — Refinamento FDR (etapa 2)

**O que faz:** pega os candidatos da triagem (z ≥ `--z_pre_filtro`,
default 2,0) e aplica filtros estatísticos mais rigorosos:

1. **p-valor paramétrico (Gama):** ajusta uma distribuição Gama aos
   200 surrogates de *cada janela* e tira o p-valor analítico da
   cauda. Resolve o problema do p-valor empírico travado em
   1/200 = 0,005 (resolução fina sem precisar de milhões de
   permutações). KL-MI é não-negativo e assimétrico → Gama é a
   família adequada.
2. **Correção FDR-Benjamini-Hochberg** sobre **todas** as janelas da
   triagem original (m = 5472 na sessão 1, 2832 na #4). Conservador:
   trata as que não foram refinadas como automaticamente não-
   significantes.
3. **Co-ocorrência entre canais:** quantos canais do mesmo arquivo
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
- `"Não significativo após correção FDR"` → não sobreviveu ao BH;
- `"Revisar: ..."` → significativo mas com alerta (kurtose alta /
  saturação / muitos canais simultâneos);
- `"Candidato robusto"` → significativo + sem alertas.

**CLI:**
```bash
python refina_candidatos.py --csv resultados.csv \
    --pasta_ns2 "<sessao>/<BASAL>" \
    --saida "<sessao>/RESULTADOS/resultados_refinados.csv"
# --z_pre_filtro 2.0 --n_surr 1000 --fdr_q 0.05 (padrão)
```

**Saída:** `resultados_refinados.csv` com colunas: `arquivo`, `canal`,
`janela_ini_s`, `janela_fim_s`, `mi_observado`, `z_score_refinado`,
`p_analitico`, `n_canais_simultaneos`, `kurtose_gamma`,
`proxy_saturacao`, `proxy_artefato_motor_150_450hz`, `significativo_fdr`,
`veredito`.

**Nota:** `filtra_sinal`, `_mi_de_bin_idx` e `bh_fdr` são funções
internas reutilizadas em quase todos os scripts do pipeline (mesma
matemática do triagem e do comodulograma).

---

## 3. `comodulogram.py` — Mapas z-scoredos (etapa 3)

**O que faz:** gera um **mapa de calor** (heatmap) do z-score do MI
para cada par (frequência de fase × frequência de amplitude),
antes e depois do notch 60 Hz. Opcionalmente aplica FDR por célula
do mapa e classifica o padrão.

**Dois modos:**
- **Lote (`--csv`):** lê `resultados_refinados.csv`, gera um PNG por
  candidato. Usa `--top_n N` para pegar só os N melhores por
  `z_score_refinado` (exploratório). Filtra por `--veredito_prefixo`.
- **Janela única (`--arquivo`):** exploração rápida de um trecho
  específico — `--canal`, `--inicio`, `--fim`.

**O que muda vs. a versão antiga:** MI agora é **z-scoredo** contra
surrogates (a mesma nula do triagem). Colormap `RdBu_r` centrado em
0 — azul = abaixo do acaso, vermelho = acima, branco = nada. Isso
permite distinguir acoplamento genuíno (pico focal no quadrante
θ×γ) de transientes (coluna inteira de 8 Hz, sem banda focal).

**Notch 60 Hz (`--notch 60`):** rejeita a rede elétrica e harmônicos
(iirnotch) antes de qualquer filtragem. Sempre rodar **com e sem**
notch e comparar: pico que cai > 1,5z com o notch é rede elétrica,
não acoplamento.

**FDR do mapa (`--fdr_q 0.05`):** Benjamini-Hochberg sobre as 275
células do mapa (fases 4–14 Hz × amplitudes 30–150 Hz, 5 Hz de passo).
Classifica cada janela em:
- **"concentrado em ΘΓ"** — ≥ 2 células sig, ≥ 50% dentro de ΘΓ
  (≥2 células sig, ≥ 50% dentro do quadrante theta-gamma) =
  acoplamento genuíno e estreito;
- **"esparso/fora de ΘΓ"** — poucas células significativas, fora do
  quadrante = transientes ritmados;
- **"nada sobrevive ao FDR"** — zero células significantes.

**Discriminador validado:** a janela rejeitada tem *mais* células
significantes mas menor % em ΘΓ; acoplamento genuíno concentra.
Janelas de 60 Hz têm 0% em ΘΓ.

**Saída:** PNGs em `<saida_dir>/<arquivo>_<canal>_<ini>-<fim>s_zcomodo.png`
+ `resumo_comodulogramas.csv` com: `canal`, `janela_ini_s`,
`janela_fim_s`, `z_pico_theta_gamma`, `fase_pico_hz`, `amp_pico_hz`,
`mi_pico`, `png`, e (com FDR) `classe_fdr`, `n_sig_fdr`, etc.

**Funções-chave internas:** `calcula_comodulograma_z` (mapa z por
célula), `p_valores_por_celula` (p por célula via Gama),
`bh_fdr_mapa` (BH sobre o mapa), `resume_cluster_fdr` (classifica
o padrão), `z_pico_theta_gamma` (pico no quadrante θ×γ).

---

## 4. `diagnostico_janela.py` — Figura diagnóstica (etapa 5)

**O que faz:** para um (arquivo, canal, janela) específico, plota 5
painéis alinhados no tempo + 4 números que discriminam acoplamento
genuíno de artefato de transiente ou respiração.

**Os 5 painéis:**
1. LFP bruto (com notch, se `--notch`);
2. Banda respiratória (0,5–3 Hz) — deflexões lentas grandes =
   artefato/potencial respiratório;
3. Theta (4–12 Hz) — visual;
4. Gamma (30–80 Hz) + envelope;
5. Espectrograma (STFT) com bandas θ e γ marcadas.

**Os 4 números:**
- **MI θ×γ z** — o acoplamento original (deve ser alto e estável);
- **MI resp×γ z** — se alto, a "fase" que organiza o gamma é
  respiratória (~1–2 Hz), não theta;
- **MI resp×θ z** — respiração modula o "theta"? (sobreposição
  sniffing 4–8 Hz ≈ theta);
- **Fator de crista do theta filtrado** — senoide pura ≈ 1,4;
  transientes / ondas agudas >> 2 (discrimina "oscilação" de
  "bursts").

**Limitação declarada:** o proxy respiratório (0,5–3 Hz) é **CEGO
para sniffing a 4–8 Hz**, que cai dentro da banda "theta". Para
desempate fino, o vídeo é o árbitro (etapa 4).

**CLI:**
```bash
python diagnostico_janela.py --arquivo "<BASAL>/<arquivo>.ns2" \
    --canal chan20 --inicio 85 --fim 95 --notch 60
```

**Saída:** PNG diagonal + 4 linhas no console.

---

## 5. `robustez_parametros.py` — Validação de robustez (etapa 7)

**O que faz:** testa se o acoplamento de cada *vencedor* sobrevive a
mudanças de parâmetro e a uma métrica alternativa — sem bins.

**O que sweepa para cada vencedor (no par de pico fixo):**
- **A. n_bins** ∈ {10, 12, 15, 18, 24, 30} — o parâmetro do KL-MI.
  Acoplamento genuíno não pode desaparecer com 12 ou 24 bins;
  o esperado é um **plató de z alto** em torno de n_bins canônico
  (±5 Hz), não imunidade aos extremos;
- **B. Largura do filtro de fase** ∈ {±0,7, ±1,0, ±1,5, ±2,5} Hz;
- **C. Largura do filtro de amplitude** ∈ {±2,5, ±3,5, ±5,0, ±7,5,
  ±10,0} Hz;
- **D. Mapa completo recomputado com n_bins=12 e n_bins=24** — o
  pico ΘΓ deve continuar no mesmo lugar (±1 Hz fase, ±5 Hz amp);
- **E. MVL (mean vector length, Canolty et al. 2006)** — métrica
  *sem bins*: |média(envelope · e^{i·fase})|. Confirma sem o
  binning do KL-MI. MVL z < z do KL-MI é normal (captura só o 1º
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
python robustez_parametros.py --pasta "<sessao>/<BASAL>" \
    --vencedores "<sessao>/vencedores.csv" \
    --resumo_fdr "<sessao>/RESULTADOS/comodulogramas_fdr/resumo_comodulogramas.csv" \
    --saida_csv "<sessao>/RESULTADOS/robustez_parametros.csv"
```

**Saída:** console + `robustez_parametros.csv` com linhas:
`janela`, `canal`, `par_pico`, `teste`, `parametro`, `z`.

---

## 6. `figura_apresentacao.py` — Figura dos vencedores (etapa 7, opcional)

**O que faz:** gera uma figura de apresentação para cada vencedor do
`vencedores.csv`. Cada figura tem 4 séries temporais + comodulograma
+ polar fase×amplitude (estilo Tort et al. 2010).

**Layout (4×2):**
- Esquerda (séries alinhadas, janela completa):
  1. LFP bruto (notch 60 Hz);
  2. theta filtrado no par de pico (±1 Hz);
  3. gamma filtrado no par de pico (±5 Hz) + envelope;
  4. espectrograma (STFT) com bandas θ/γ marcadas;
- Direita:
  5. comodulograma z-scoredo (recalculado aqui, com notch) com o
     par de pico marcado;
  6. distribuição polar fase×amplitude (18 bins), amplitude de γ
     normalizada por bin de fase do theta — o MI polar é reportado
     abaixo do círculo.

**Nota:** o comodulograma é recalculado aqui com o par de pico do
vencedor — garante consistência visual com os números da robustez.

**CLI:**
```bash
python figura_apresentacao.py --pasta_ns2 "<sessao>/<BASAL>" \
    --vencedores "<sessao>/vencedores.csv" \
    --saida_dir "<sessao>/RESULTADOS/figuras"
```

**Saída:** `<saida_dir>/<rotulo>_<canal>.png` (nome do PNG = rótulo +
canal, conforme convenção).

**Depende de:** `vencedores.csv` (formato: `rotulo,arquivo,canal,
inicio_s,fim_s,fase_pico_hz,amp_pico_hz[,comportamento]`).

---

## 0 (NOVO). `exploracao_interativo.py` — Passo 0: exploração visual em Jupyter

**Arquitetura:** navegador assíncrono para milhões de pontos,
integrado ao Jupyter. O pipeline canônico (etapas 1–7) varre
**cego** — processa todas as janelas e devolve números. O passo 0
inverte isso: o pesquisador **vê o sinal primeiro**, marca o que
interessa, e só então o pipeline processa. O pipeline não acha
acoplamento? Provavelmente a janela была errada — você viu o
acoplamento no LFP e precisa informar o carimbo.

**Tecnologia:** VisPy (CanvasSci, GPU-accelerated) ou Bqplot (d3.js
no browser) — lidam com milhões de pontos sem travar, zoom livre,
scroll, atualização em tempo real. Alternativa: Neuroglancer (se
houver voxels 3D) ou painel MNE-Python com TimeSeriesViewer.
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
│ (scroll/    │ Ver teta (4-8 Hz) como ondulação│ (freq vs t) │
│  clique)    │ Ver gama (30-80 Hz) como    │   Teta ↑ quando│
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
4. ve teta (4-8 Hz) e gama (30-80 Hz) juntos?
5. clica em "Marcar instante" → salva (t_start, t_end, canal)
6. (opcional) clica em "Enviar para triagem" → recebe z + FDR
   imediatamente na célula do Jupyter
7. repete 4-6 para todos os instantes interessantes
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
cd C:\acoplamento_theta-gamma\SCRIPT
jupyter lab
# abrir notebooks/exploracao_interativo.ipynb
```

**Depende de:** `ns2_utils.py` (leitura), VisPy ou Bqplot (instalar:
`pip install bqplot` ou `pip install vispy`). MNE-Python já está no
requirements.

**Status:** a ser implementado (script atual `comodulogram_interativo.py`
é o protótipo estático; reescrever como Jupyter + VisPy/Bqplot é o
próximo passo).

---

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
3. **Sub-janelas de 2 s:** acoplamento genuíno aparece na maioria
   das sub-janelas; dirigido por 1–2 episódios concentra todo o z
   nelas;
4. **Controles:** chan20 @ 85–95 e chan10 @ 30–40 (bons, validados
   em comportamento + respiração) vs. chan32 @ 95–105 (artefato
   conhecido — janela rejeitada #3); o veredito sobre o caso
   disputado só vale se o teste separar esses dois grupos;
5. **Vizinhos:** z na célula do pico nos canais vizinhos imediatos
   — fonte cortical local tem gradiente; cabo/EMG difuso não.

**Saída:** `auditoria/auditoria_transientes.csv`, `auditoria/
contexto_amplitude.csv`, PNG por caso.

**CLI (sessão 08/07 default):** `python audita_transientes.py`
**Nova sessão (via CLI, casos nunca no código):**
```bash
python audita_transientes.py --pasta_ns2 "<sessao>/<BASAL>" \
    --saida_dir "<sessao>/RESULTADOS/auditoria" \
    --casos "rotulo1,arq1,chan1,ini1,fim1,fp1,fa1;rotulo2,..." \
    --vizinhos "chan24:z,chan26:z,..."
```

---

## 8. `audita_segmentos.py` — Localização temporal do acoplamento

**O que faz:** complementa a auditoria de transientes — localiza,
*no tempo*, onde mora o acoplamento de um caso. Para cada sub-
segmento, recomputa o mapa z completo e reporta o z na célula do
pico original e o pico ΘΓ do próprio segmento.

**Uso (default sessão 08/07):** `python audita_segmentos.py`
**Nova sessão:**
```bash
python audita_segmentos.py --pasta "<sessao>/<BASAL>" \
    --arquivo <arq>.ns2 --canal chanXX \
    --segmentos "ini1-fim1,ini2-fim2,..." --fp F --fa A
```

**Exemplo:** `--segmentos "20-26,26-30,28-30,20-30"` divide a janela
em sub-blocos de 2–4 s para ver onde o z aparece (ou não).

---

## 9. `audita_footprint.py` — Pegada espacial do acoplamento

**O que faz:** mede o z na célula do pico em **todos os 32 canais**
da mesma janela. O discriminador: fonte cortical local produz
*gradiente suave* (poucos canais vizinhos significantes); artefato
difuso (respiração, movimento, volume conduzido, cabo) aparece
simultaneamente em muitos canais distantes.

**Leitura:** o caso "rearing" (chan18/20/30/32, 4 canais) é referência
do que é pegada *local*. O grooming (chan22 + cluster 22/24/26/28)
é *focal*. Se 26 dos 32 canais estão em z≥3, é incompatível com
cabo/EMG difuso.

**Uso:** `python audita_footprint.py` (default: sessão 08/07, 3
janelas pré-configuradas). Para outra sessão, edite a lista `JANELAS`
no topo do script — casos específicos *nunca* entram no código via
CLI neste script (use `--casos` no `audita_transientes.py`).

---

## 10. `audita_grooming_robustez.py` — Robustez do núcleo de grooming

**O que faz:** verifica se o *núcleo* de grooming do vencedor 2
(003 @ 20–27 s, par 5×35 Hz, re-ancorado no vídeo) sobrevive ao
mesmo batizado dos outros vencedores: sweep de n_bins + MVL no
chan22 + consistência nos vizinhos chan24/26/28.

**Contexto:** o z=9,4 da janela original 20–30 s era inflação por
escolha de janela (mistura de estados + nula estreita). O núcleo
robusto é z≈4 durante o grooming puro.

**Uso:** `python audita_grooming_robustez.py` (sessão 08/07 default).

---


## 12. `audita_harmonico.py` — Teste de razão harmônica Θ→Γ (auditoria)

**O que faz:** usa FOOOF (Kühn et al. 2026) para estimar `cf_teta` numa janela
de contexto longa (45 s) e testa se `amp_pico ≈ n × cf_teta` em frequência
(razão inteira, tolerância 10%) e em fase (PLV entre `n×phi_theta` e
`phi_gamma`).

**Três dimensões independentes:**
1. **Frequência**: razão inteira `amp/fase ≈ n` (tolerância 10% de `cf_teta`)
2. **Forma de onda**: skewness do teta (reusado de `audita_skewness.py`)
3. **Fase**: PLV entre `n×phi_theta` e `phi_gamma` — discriminador mais forte

**Veredito:** CLEAN / REVISAR_RAZAO_INTEIRA / REVISAR_FASE_TRAVADA / SUSPEITO_HARMONICO_FORTE / SEM_REFERENCIA_TETA

**CLI:**
```bash
python audita_harmonico.py --csv "<sessao>/RESULTADOS/vencedores.csv" \
    --pasta_ns2 "<sessao>/<BASAL>" \
    --saida "<sessao>/RESULTADOS/harmonico.csv" \
    --janela_contexto_s 45 --modo_preprocesso hibrido --f_linha 60.0
```

---

## 13. udita_harmonico_hfo.py � Teste de raz�o harm�nica G?HFO *(novo 2026-09-05)*

Testa se HFO (150�250 Hz) � harm�nico de Gamma (30�80 Hz) usando FOOOF para
estimar cf_gamma e PLV(n�phi_gamma, phi_hfo). An�logo ao udita_harmonico.py.
Vereditos: CLEAN / REVISAR_RAZAO_INTEIRA / REVISAR_FASE_TRAVADA / SUSPEITO_HARMONICO_FORTE / SEM_REFERENCIA_GAMMA.

---

## Conven��es comuns a todos os scripts

- **Notch 60 Hz** (--notch 60): padr�o recomendado. Comparar sempre com/sem.
- **Nula de surrogates:** deslocamento circular >= 1 s, 200 repeti��es, semente 42.
- **Par�metros canonicos:** janela 10 s / passo 5 s; theta 4-8 Hz x gamma 30-80 Hz; n_bins 18; 200 surrogates.
- **iltra_sinal**: Butterworth bandpass, iltfilt (zero fase), ordem 3.
- **_mi_de_bin_idx**: nucleo vetorizado do KL-MI via 
p.bincount.
- **h_fdr**: Benjamini-Hochberg; m_total = familia completa de testes.

## Princ�pio revisor

> Casos espec�ficos de sess�o (canais, janelas, offsets, vencedores)
> **NUNCA** entram hardcoded no c�digo. Entram por CLI/CSV.

---

## Pipeline multi-acoplamento: Theta-Gamma / Theta-HG / Theta-HFO

*(atualiza��o 2026-09-05)*

| Par | Banda amplitude | Substrato | Estado |
|---|---|---|---|
| 	heta_gamma | 30-80 Hz | Fast gamma: CA3->CA1 | Explora��o |
| 	heta_hg | 80-150 Hz | High gamma: EC->CA1 | Misto |
| 	heta_hfo | 150-250 Hz | HFO/ripple: CA1 local | Repouso/SWR |

**Efici�ncia:** 3 pares x 59 janelas x 200 surrogates ~ 8 s para 300 s de sinal, 1 canal.

```bash
# .mat direto (sem .ns2):
python pipeline/triagem_pac_mat.py --mat DADOS_EXEMPLO_LFP_HG_HFO/LFP_HG_HFO.mat

# .ns2 multi-par:
python pipeline/triagem_pac.py --pasta <sessao>/BASAL \
    --pares theta_gamma theta_hg theta_hfo \
    --saida <sessao>/RESULTADOS/resultados_triplo.csv
```

### Resultado validado no LFP_HG_HFO.mat (2026-09-05)

| Par | z mediana | z std | Candidatos (z>=3) |
|---|---|---|---|
| theta_gamma | 6.93 | 2.96 | 58/59 |
| theta_hg | 12.39 | 4.75 | 57/59 |
| theta_hfo | 0.50 | 1.52 | 7/59 |

	eta_ok: 17/59 (29%), std=0.457 � varia (bug de cache corrigido).
atio_hfo_gamma mediana=0.013 � HFO nao e harmonico de Gamma.

### Bugs corrigidos em 2026-09-05

| Bug | Sintoma | Causa | Correccao |
|---|---|---|---|
| Cache 	eta_ok | 58/58 = 1, variancia zero | Calculo fora do loop | Movido para 	eta_ok_por_janela() dentro do loop |
| z-score single-surrogate | z = 3e14 alternando com 0 | max(mi_s*0.1, 0.01) como "dp" | Substituido por mi_com_surrogates() 100 surrogates |
| Escala .mat | Diagnostico invisivel | Sem verificacao de std | dapta_lfp_mat.py com info_sinal() + --normaliza |
