# Scripts do SCRIPT central — explicação rápida

> Pasta `CONTEXTO` reúne documentação que não muda entre sessões e
> serve de referência para qualquer sessão futura do estudo
> (MTESC04/MTESC05 × NOCI/LAC). Nada de específico de sessão aqui
> — caminhos, offsets e vencedores entram por CLI/CSV, nunca no código.

## Visão geral do pipeline (Passos 1–5, com sub-passos 3.5–3.8)

**A exploração interativa NÃO é mais o Passo 0.** Até 2026-09 o plano era
"passo 0 e' primeiro, nao opcional": o pesquisador veria o sinal antes de
qualquer coisa, numa ferramenta assíncrona (VisPy/Neuroglancer/Bqplot)
para navegar milhões de pontos. Isso nunca foi implementado assim — a
ferramenta que existe hoje (`comodulogram_interativo.py` + notebook) faz
o oposto: roda **depois** de tudo, para dissecar visualmente os
vencedores já filtrados (Passo 5, ver seção dedicada abaixo).

```
processa_sessao.py      → orquestrador por sessão (raiz de pipeline/): estágio 0 (canais
                          vivos) → estágio 1 (triagem_coocorrencia.py decide se o canal
                          merece o pipeline pesado) → Passos 1, 2, 3 + auditorias de
                          skewness/harmônico, canal a canal. É o que roda em lote.
triagem_pac.py          → Passo 1 (etapa1_triagem): varre TODAS as janelas (10 s / 5 s)
refina_candidatos.py    → Passo 2 (etapa2_refinamento): FDR de janela + filtros de artefato
comodulogram.py         → Passo 3 (etapa3_comodulograma): mapas z...scoredos ±notch + FDR do mapa
gerar_template_comportamento.py,
anotador_comportamento.py,
junta_comportamento.py  → Passo 3.5 (pipeline/comportamento/): anotação comportamental
                          sincronizada com vídeo, mesclada ao dataset mestre
agrega_resultados.py,
enriquece_dataset_mestre.py → Passo 3.6 (pipeline/dataset_mestre/): constrói o dataset
                          mestre (merge por canal) + FOOOF v2 + portão de banda larga
audita_*.py             → Passo 3.7 (pipeline/auditorias/): auditorias pós-hoc
                          (transientes, segmentos, pegada espacial, held-out, harmônico...)
audita_janela.py        → orquestrador forense: as 4 auditorias
                          compatíveis (skew+transientes+footprint+
                          harmônico) numa só passada por canal/janela
diagnostico_janela.py   → pipeline/auditorias/: 5 painéis + MI respiração vs θ×γ
consolida_vencedores.py → Passo 3.8 (pipeline/dataset_mestre/, novo 2026-09): filtro final
                          (estatístico + FOOOF + harmônico + comportamental) →
                          candidatos_vencedores_consolidados.csv
robustez_parametros.py,
figura_apresentacao.py  → Passo 4 (etapa4_validacao): sweep de parâmetros + MVL,
                          e figuras finais dos vencedores (opcional)
comodulogram_interativo.py,
notebook exploracao_interativo.ipynb → Passo 5 (etapa5_exploracao): dissecação
                          interativa (CLI ou notebook) dos vencedores consolidados,
                          com zoom FOOOF + comodulograma
pac_core/               → núcleo matemático compartilhado (refatoração
                          2026-09): io.py, filtering.py, pac_metrics.py,
                          workspace.py. ns2_utils.py e demais "shims"
                          listados abaixo reexportam daqui — ver seção 17.
pac_studio/             → API de 1 linha sobre pac_core (seção 19)
dashboard_desktop/      → programa desktop ThetaGamma-Studio: roda as etapas
                          acima com log ao vivo e inspeciona resultados (seção 20)
preditor/               → protótipo de previsão de PAC p/ closed-loop, não
                          validado p/ hardware (seção 21)
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

## 4. `diagnostico_janela.py` — Figura diagnóstica (pipeline/auditorias/, Passo 3.7)

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

## 5. `robustez_parametros.py` — Validação de robustez (etapa 4, Passo 4)

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
python pipeline/etapa4_validacao/robustez_parametros.py ......pasta "<sessao>/<BASAL>" \
    ......vencedores "<sessao>/vencedores.csv" \
    ......resumo_fdr "<sessao>/RESULTADOS/comodulogramas_fdr/resumo_comodulogramas.csv" \
    ......saida_csv "<sessao>/RESULTADOS/robustez_parametros.csv"
```

**Saída:** console + `robustez_parametros.csv` com linhas:
`janela`, `canal`, `par_pico`, `teste`, `parametro`, `z`.

.........

## 6. `figura_apresentacao.py` — Figura dos vencedores (etapa 4, Passo 4, opcional)

**O que faz:** gera uma figura de apresentação para cada vencedor de
`candidatos_vencedores_consolidados.csv` (Passo 3.8). Cada figura tem
4 séries temporais + comodulograma + polar fase×amplitude (estilo Tort
et al. 2010) + 2 painéis FOOOF (teta e gama/HG, novo 2026-09).

**Layout (4×3):**
... Coluna 1 (séries alinhadas, janela completa):
  1. LFP bruto (notch 60/120/180/240 Hz);
  2. theta filtrado no par de pico (±1 Hz);
  3. gamma filtrado no par de pico (±5 Hz) + envelope;
  4. espectrograma (STFT) com bandas θ/γ marcadas;
... Coluna 2:
  5. comodulograma z...scoredo (recalculado aqui, com notch) com o
     par de pico marcado;
  6. distribuição polar fase×amplitude (18 bins), amplitude de γ
     normalizada por bin de fase do theta — o MI polar é reportado
     abaixo do círculo;
... Coluna 3 (novo 2026-09, reusa `painel_fooof`/`ajusta_fooof_teta_gamma`
   de `etapa5_exploracao/comodulogram_interativo.py`, não duplicado aqui):
  7. FOOOF banda baixa (2-45 Hz): fundo aperiódico (knee) + pico de teta;
  8. FOOOF banda alta (35 Hz-~0,95×Nyquist, limpeza de linha Kuhn):
     fundo aperiódico + picos de gama/HG.

**Nota:** o comodulograma é recalculado aqui com o par de pico do
vencedor — garante consistência visual com os números da robustez.

**Canal:** convenção 1-based do dataset mestre (mesma de
`comodulogram_interativo.py`), não o nome nativo do `.ns2`.

**CLI:**
```bash
python pipeline/etapa4_validacao/figura_apresentacao.py --pasta_ns2 "<sessao>/Basal antes da infusao" \
    --vencedores resultados/candidatos_vencedores_consolidados.csv \
    --saida_dir resultados/figuras
```

**Saída:** `<saida_dir>/<rotulo>_ch<canal>.png` (`rotulo` é opcional,
gerado automaticamente a partir de arquivo+canal+janela se ausente).

**Depende de:** um CSV com colunas `arquivo,canal,janela_ini_s,
janela_fim_s,fase_pico_hz,amp_pico_hz[,rotulo,comportamento,par]`
(schema atual de `candidatos_vencedores_consolidados.csv`; aceita
também o formato antigo `inicio_s`/`fim_s`).

.........

## 6.5. `comodulogram_interativo.py` + notebook — Dissecação Interativa dos Vencedores (etapa5_exploracao, Passo 5)

**Isto substitui o antigo plano de "Passo 0" descrito acima** (nunca
implementado como navegador assíncrono de milhões de pontos). A
ferramenta real, renomeada de `etapa0_exploracao` para
`etapa5_exploracao` em 2026-09, roda **depois** da consolidação
(Passo 3.8), não antes da triagem: consome
`candidatos_vencedores_consolidados.csv` e gera, por evento, os mesmos
4 painéis do `figura_apresentacao.py` (LFP+teta+gama, FOOOF teta,
FOOOF gama/HG, comodulograma) — mas de forma interativa/pontual em vez
de em lote.

**Duas formas de uso (mesma lógica por baixo):**
- **Notebook** (`pipeline/etapa5_exploracao/notebooks/exploracao_interativo.ipynb`):
  dashboard com dropdowns em cascata Rato → Comportamento → Janela
  campeã (ordenada por z-score), widget "Visualizar".
- **CLI** (`comodulogram_interativo.py --zoom_t_center <s>`): gera uma
  janela específica sem abrir GUI — útil para automação/testes.

**Canal:** convenção 1-based do dataset mestre (`--canal` = mesmo
número da coluna `canal` do CSV); o script converte internamente para
o índice 0-based do array. Funções reusáveis por outros scripts:
`ajusta_fooof_teta_gamma` (fit FOOOF teta+gama com limpeza de linha
Kuhn) e `painel_fooof` (plotagem do painel), ambas usadas também por
`figura_apresentacao.py` (etapa4_validacao) — ver seção 6.

**Tecnologia real:** matplotlib estático (headless-compatível) +
ipywidgets no notebook — não VisPy/Neuroglancer/Bqplot como o plano
original de "Passo 0" previa; a escala de dados (janelas de 10s, não a
sessão inteira ponto-a-ponto) não exigiu isso.

**CLI mínima (notebook):**

```bash
cd D:\acoplamento_theta-gamma\SCRIPT
jupyter notebook pipeline/etapa5_exploracao/notebooks/exploracao_interativo.ipynb
```

**Depende de:** `pac_core.io`/`pac_core.filtering` (leitura e filtros),
`fooof` (ajuste FOOOF), `ipywidgets` (controles do notebook).

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

**CLI:** `--pasta_ns2` é obrigatório (antes tinha padrão relativo `../Basal antes da infusao`, que só funcionava rodando de dentro de uma pasta de sessão).
**Uso (casos nunca no código):**
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

**Uso** (`--pasta` obrigatório; o padrão relativo antigo foi removido):
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
- **Filtro:** só `veredito_refino == "Candidato robusto"` (não a cadeia completa de `consolida_vencedores.py`, que já exige comportamento anotado — dependência circular).
- **Preserva anotações ao regenerar** (2026-09-14): rodar de novo depois de processar mais sessões (ex.: infusão, depois de já ter anotado o basal) casa pela chave da janela e traz de volta `comportamento`/`observacoes`/`video_tempo_ini`/`video_tempo_fim` do arquivo antigo — nunca apaga trabalho feito. Faz backup timestamped em `backups_comportamento/` antes de sobrescrever. Janela que não sobrevive no dataset mestre atual (deixou de ser robusta, ou pertencia a uma sessão reprocessada do zero, como o LAC renomeado) é descartada silenciosamente da nova lista.

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
- **Sessão = (sessao, condicao), não só sessao (2026-09-14):** uma sessão de infusão tem 4 vídeos na mesma pasta (basal + 0h/1h/2h pós), cada um cobrindo um conjunto diferente de `.ns2` — ao contrário do basal, onde 1 `sessao` = 1 vídeo. O seletor de sessão agora lista pares `sessao :: condicao`; `config_anotador.json["sessoes"]` usa a chave composta só para condição ≠ basal (entradas antigas de basal continuam com a chave simples, sem precisar reconfigurar). O autodetect de vídeo (`_tentar_autodetectar_video`) desempata pela palavra-chave da condição no nome do arquivo (`"0h"`, `"1h"`, `"2h"`, `"basal"`), já que os 4 vídeos da rodada moram na mesma pasta e o achado por rato sozinho é ambíguo.

### 14.3 `junta_comportamento.py` — Mesclagem com o Dataset Mestre
- **O que faz:** Combina as anotações feitas no `template_comportamento.csv` de volta ao dataset mestre de resultados (`resultados/dataset_mestre_final.csv`), propagando o comportamento anotado para todos os canais correspondentes àquela janela temporal. Salva em `resultados/dataset_mestre_COM_COMPORTAMENTO.csv`.

*Nota (reorg 2026-09):* as três ferramentas acima e seus dados (`template_comportamento.csv`, `config_anotador.json`, `backups_comportamento/`) vivem juntos em `pipeline/comportamento/` — pasta autocontida, sem depender de arquivos soltos na raiz de `SCRIPT/`.

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
#    por canal (renomeia veredito→veredito_refino, adiciona sessao/grupo/condicao/canal)
python pipeline/dataset_mestre/agrega_resultados.py --resultados "<sessao>/RESULTADOS" --saida dataset_mestre.csv

# 2. Enriquece: FOOOF v2 (aperiodic_mode='knee', ajuste particionado 2-45Hz teta / 35-250Hz gama)
#    + portão de banda larga (rebaixa "Candidato robusto" com suspeito_banda_larga=True)
python pipeline/dataset_mestre/enriquece_dataset_mestre.py --entrada dataset_mestre.csv --saida dataset_mestre_v2.csv
```
As duas etapas do passo 2 (`--etapas fooof portao`, ambas por padrão) são independentes — `portao` só usa colunas que já vêm do passo 1, nenhuma delas é criada pelo `fooof`. Não-destrutivo por padrão; `--in_place` sobrescreve a entrada (comportamento do antigo `aplica_portao_banda_larga_mestre.py`). Guarda contra reenriquecer: se as 14 colunas FOOOF v2 já existirem, aborta com erro claro a menos que `--forca` seja passado (evita colunas `_x`/`_y` duplicadas silenciosas do `pandas.merge`).

**Colunas de identificação (corrigidas em 2026-09-13).** `agrega_resultados.py` varre `--resultados` recursivamente e tira do caminho de cada pasta `chanN/`:
- `sessao`: a primeira parte do caminho que contém `MTESC` (o nome com o prefixo do rato, ex. `MTESC05_LAC_Rodada-1-04-05-2024`); só cai para uma parte com `Rodada` se nenhuma tiver `MTESC`. Antes, o laço ficava com a **última** parte contendo "Rodada" — numa estrutura `MTESC05_LAC_Rodada-X/Rodada-X/chan5` a pasta interna vencia e o rato sumia do rótulo.
- `grupo`: `NOCI`, `LAC` ou `VEH` (`grupo_farmacologico()`). O que foi infundido na sessão (`...-lac_hemdir` / `...-veh_hemesq` no nome da pasta) vence a pasta do rato — `MTESCnn_LAC` guarda sessões de lactato **e** de veículo, e antes toda sessão de veículo saía rotulada LAC. Sem esse token, procura `NOCI`, `VEH` ou `LAC` em qualquer parte do caminho (`indeterminado` se nenhuma).
- `condicao`: `basal`, `0h_pos`, `1h_pos` ou `2h_pos`, pelo nome da pasta da condição de infusão. Antes, `condicao` guardava o grupo farmacológico (e na prática ficava sempre `basal`, porque o nome da sessão não contém "NOCI"/"LAC").

Pastas cujo nome começa com `_` são **ignoradas** (desde 2026-09-14): `_CONTAMINADO_nao_usar/` (sessão antiga que processou arquivos de MTESC03_LAC e MTESC05_LAC juntos) e `_LAC_basal_antigo_nao_usar/` (saídas do basal LAC com os nomes antigos das pastas). É a convenção para isolar saída descartada sem apagar. Organização atual: basal em `basal/NOCI/`, `basal/LAC/` e `basal/VEH/` (nome `<pasta da sessão>_Basal antes da infusao`), infusão em `<RATO>_NOCI_<N>_<DATA>/<condição>/`.

**Pastas LAC com o rato no nome (2026-09-14):** as sessões de `LAC_NOCI/MTESC03_LAC/` e `MTESC05_LAC/` tinham nomes idênticos entre os dois ratos (`Rodada-1-02-05-2024-lac_hemdir`...), e só a pasta-mãe as separava. Foram renomeadas para `MTESC03 -- Rodada-1-02-05-2024-lac_hemdir` (mesmo padrão `MTESCnn -- ...` do NOCI); o mapeamento antigo → novo está em `LAC_NOCI/_renomeacoes_pastas_LAC.csv`.

**Lote do basal:** `RESULTADOS_MESTRADO/roda_lote_mtesc.ps1` acha sozinho toda pasta `LAC_NOCI/<rato>/<sessão>/Basal antes da infusao` (antes era uma lista fixa, com nomes LAC desatualizados e 3 entradas do MTESC05_LAC apontando para a rodada inteira em vez do basal). `-ratos "MTESC0?_LAC"` filtra, `-listar` só mostra pasta/grupo/saída sem rodar, `-sem_agregar` pula a agregação final. Pula (com aviso) sessão cujo nome não traga `MTESCnn`.

**Template de comportamento:** `gerar_template_comportamento.py` lista toda janela do CSV que receber, sem filtrar veredito. Para anotar só o que importa, gere o template a partir de um CSV filtrado por `veredito_refino == "Candidato robusto"` — **não** pela cadeia completa do `consolida_vencedores.py`, cujo portão comportamental exige `comportamento` preenchido e descartaria exatamente as janelas ainda não anotadas. Foi esse filtro que gerou o template atual (reprodução janela a janela confirmada em 2026-09-13).

---

## 18. `processa_sessao.py` + `triagem_coocorrencia.py` — Orquestrador por sessão

`processa_sessao.py` (raiz de `pipeline/`) é o que roda em lote: processa uma pasta de sessão (ou de condição de infusão) canal a canal.

```bash
python pipeline/processa_sessao.py --pasta "<pasta com .ns2>" --saida "<pasta base>" [--min_janelas 3]
```
Saída em `<saida>/<nome da pasta>/chanN/`, mais `resumo_canais.csv`. Rode de dentro de `SCRIPT/`: ele chama as etapas com caminhos relativos (`pipeline/...`).

- **Estágio 0 — canais vivos:** lê os primeiros 30 s do primeiro `.ns2`, descarta canais com amplitude morta e sinaliza (sem descartar) outliers de RMS (`rms_outlier` no resumo).
- **Estágio 1 — co-ocorrência (`triagem_coocorrencia.py`, o "Passo 0.5"):** para cada `.ns2`, gera `chanN/coocorrencia_<arq>.csv` com, por janela de 10 s / passo 5 s:
  - `teta_ok` (de `teta_ok_por_janela`, do `triagem_pac.py`);
  - `gamma_pot_rel` (30–80 Hz) e `hg_pot_rel` (80–150 Hz) — potência da banda sobre a média do PSD (Welch); `gamma_ok`/`hg_ok` = 1 quando passa de 1,5 × a mediana das janelas **do próprio arquivo** (limiar relativo, muda de arquivo pra arquivo);
  - `hfo_cru` e `ripple` — detecção de eventos no sinal inteiro, em resolução de amostra (150–250 Hz, `exigir_sharp_wave=False`). O console mostra o "diagnóstico do paradoxo HFO" (taxa de HFO cru sem ripple; alerta acima de 95% ou abaixo de 5%).

  O orquestrador soma as janelas `teta_ok & gamma_ok` e `teta_ok & hg_ok` entre os arquivos e **pula o canal** se as duas ficarem abaixo de `--min_janelas`. Teta + HFO e teta + ripple são registrados no resumo, mas **não entram no corte** — um canal só com teta + HFO é pulado, mesmo que o Estágio 2 rode o par `theta_hfo`.

  CLI avulsa: `python pipeline/etapa1_triagem/triagem_coocorrencia.py --origem <arq.ns2|.mat> --canal <idx 0-based | var do .mat> [--saida ...] [--limiar_dp 4.0] [--duracao_ms 25.0]`.
- **Estágio 2 — pipeline pesado:** `triagem_pac.py` (3 pares) → `refina_candidatos.py` → `comodulogram.py` (um por par) → `audita_skewness.py`, `audita_harmonico.py` (por par) e `audita_harmonico_hfo.py`.

**Sessão com zero candidatos:** olhe o `resumo_canais.csv` antes de suspeitar de bug. Em 2026-09-13, MTESC03_LAC Rodada-1 e MTESC05_LAC Rodada-2-06 deram zero porque nenhum canal passou do Estágio 1 (teta + gama e teta + HG ≤ 2 janelas em todos os canais), apesar de teta + HFO alto — comportamento esperado do corte, não falha.

---

## 19. `pac_studio/` — API de 1 linha *(novo 2026-09-13)*

Consolida sequências que existiam copiadas à mão em vários scripts, sem reimplementar nenhum cálculo (só orquestra `pac_core` e a config de bandas `BAND_PAIRS` de `triagem_pac.py`).

- **`analisa_janela(caminho, canal, t_ini_s, t_fim_s, par="theta_gamma", ...)`** — carrega → fatia → notch → filtra fase/amplitude → Hilbert → MI com surrogates, numa chamada. Devolve `mi_observado`, `z_score`, `p_empirico` e as bandas usadas; não decide "significativo". `rng=42` por padrão (convenção do projeto). `canal` aceita nome nativo (`chan5`) ou índice 1-based do dataset mestre.
- **`transicao_estado(df_mestre, gap_max_s=None)`** — reconstrói os pares consecutivos comportamento(N-1) → janela N com a mesma chave de agrupamento (`sessao, arquivo, canal, par`) e o mesmo `GAP_MAX_S` (15 s) de `preditor/preditor_estado_comportamental.py::monta_dataset()`, importado de lá. A análise χ² da Hipótese H2 não existia como código em nenhum script (só os números, nos docs e em `gerar_figuras_dissertacao.py`, digitados à mão); com `pd.crosstab` + `chi2_contingency` sobre o retorno desta função, reproduz n = 2715, χ² = 41,26, p = 2,57×10⁻⁷ e os percentuais por categoria da figura 2 da dissertação.

**Teste:** `tests/test_pac_studio.py` — par inválido levanta erro; paridade bit a bit de `analisa_janela` contra a mesma sequência chamada à mão via `pac_core`; contraste acoplado × não acoplado em sinal sintético; reprodução exata de H2 contra `dataset_mestre_COM_COMPORTAMENTO.csv`.

---

## 20. `dashboard_desktop/` — Programa desktop ThetaGamma-Studio *(novo 2026-09-14)*

Interface gráfica (Tkinter + ttkbootstrap, tema `bootstrap-light`) para rodar o pipeline e inspecionar resultados. Substituiu um protótipo em Streamlit (removido). Rodar com `python dashboard_desktop/app.py`; guia de uso em `docs/manual_usuario.md`.

| Arquivo | Papel |
|---|---|
| `app.py` | Monta as 3 seções (Pipeline, Análise de Dados, Anotador de Vídeo); ao fechar com etapa rodando, pergunta e encerra os processos (evita órfãos). |
| `runner.py` | `PainelExecucao`: roda um script via `subprocess` com stdout ao vivo num painel de log. Força `PYTHONUNBUFFERED=1` e `PYTHONIOENCODING=utf-8` (herdados pelos processos que `processa_sessao.py` dispara), senão o log chegaria em blocos e `print` de θ/× poderia quebrar em cp1252. **Parar** usa `taskkill /T` no Windows para derrubar a árvore inteira — `terminate()` só mataria o `processa_sessao.py` e deixaria o neto rodando órfão. Também tem `_rodar_em_thread` (função Python in-process em thread, usado pelas abas de análise). |
| `aba_pipeline.py` | 6 sub-abas que montam a linha de comando de scripts existentes: Sessão Completa (`processa_sessao.py`), Triagem (`triagem_pac.py`, com botão Demo), Refinamento e Filtros (`refina_candidatos.py` + gráfico de vereditos ao terminar), Comodulograma (lote), Auditorias (skewness / harmônico / harmônico HFO) e Agregação (agrega → enriquece → junta comportamento → consolida). Campos = flags argparse reais. `audita_transientes.py` e `audita_footprint.py` ficam fora de propósito (CLI presa a uma sessão, casos passados como texto). Saídas padrão da Agregação em `*_novo.csv`; toda etapa confirma antes de sobrescrever. |
| `aba_analise.py` | Inspeção de resultados em 6 sub-abas: canais (LFP multicanal com notch/passa-faixa/detrend/referência), comodulograma de uma janela (`calcula_comodulograma_z`), portões de qualidade (lê vereditos já computados; skewness reprova com `SUSPECT (Asymmetric)`), galeria dos 190 vencedores, comportamento (boxplot por categoria, χ² N-1→N, painel farmacológico) e replay do preditor (reusa `prever_pac_tempo_real.py`; aviso fixo de protótipo não validado para hardware). |
| `aba_anotador.py` | Abre `pipeline/comportamento/anotador_comportamento.py` como processo separado e mostra o progresso de anotação por sessão lendo o template. |
| `estilo.py` | Tema e helpers de layout (`_secao`, `_callout`, `_entrada_arquivo`, `_embute_figura`). Cores do callout vêm do tema ativo, não de hex fixo. |

Nenhuma métrica é recalculada com lógica nova: tudo passa por `pac_core`, `pac_studio`, `preditor` ou pelos próprios scripts de `pipeline/`.

---

## 21. `preditor/` — Previsão de PAC para closed-loop *(protótipo)*

> ⚠️ `preditor/README_preditor.md` (13/09/2026): **protótipo não validado — não usar para disparar hardware.** Detalhes, histórico e resultados completos estão lá; aqui fica só o papel de cada script. Nenhum recebe argumentos, exceto `prever_pac_tempo_real.py`.

| Script | O que faz | Status |
|---|---|---|
| `preditor_estado_comportamental.py` | Prevê se a janela N vira "Candidato robusto" a partir do comportamento e da potência teta da janela N-1 (sem vazamento). Ablation de 4 melhorias, `StratifiedKFold` de 5, limiar calibrado por precisão-recall. Lê `dataset_mestre_COM_COMPORTAMENTO.csv` + `.ns2`; grava `resultados/_preditor_estado_teta.csv` e `modelo_estado_comportamental.pkl`. `pac_studio` importa `GAP_MAX_S` daqui. | Melhor linha até agora: AUC-ROC 0,618. Protótipo. |
| `prever_pac_tempo_real.py` | Inferência em replay (`--ns2 <arq> --canal chan16 --replay`) ou streaming sintético (`--simular`): janela de 10 s, passo de 1 s, chama `dispara_ttl()` quando P ≥ limiar (padrão: `limiar_f1` do `.pkl`, 0,431). O comportamento do replay é **um valor fixo de CLI** (`--comportamento`) para todas as janelas. `dispara_ttl()` só imprime — é o ponto de ligar hardware (ver `docs/hardware_closed_loop.md`). Se não achar o `.pkl` novo, cai no `modelo_pac.pkl` antigo. Usado pela aba de replay do programa desktop. | O docstring fala em "modelo validado"; o README do preditor diz o contrário. Vale o README. |
| `validar_preditor.py` | Validação com **controles reais** (janelas do mesmo `.ns2`, longe do evento), LOOCV, 8 variáveis (5 espectrais + MI_z, MVL_z no par oficial + `Footprint_n_z3`). Lê `candidatos_vencedores_OURO_PURIFICADO_v2.csv` (cada `.ns2` via `localiza_ns2`). Grava as probabilidades LOOCV em `resultados/_oof_lfp_10s.csv` (curva ROC real da figura 1). **Sobrescreve** `modelo_pac.pkl`. | Acurácia ≈ 0,43 (chance): nenhum precursor no LFP bruto. |
| `analisar_janelas_ultracurtas.py` | Procura precursores ultracurtos (1–3 s antes) e detecção pelo início do evento nos 190 vencedores vs. controles; BH-FDR + classificadores, semente 42. Grava `resultados/_oof_lfp_ultracurto.csv` (probabilidades fora-da-amostra por cenário) e `resultados/_trajetoria_teta_pre_evento.csv` (teta a -3/-2/-1 s, figuras 1 e 4). | Exploratório (o código diz: p bruto < 0,05 não é confirmatório). |
| `analisar_pre_evento.py` | "Versão 3" do extrator dos 10 s antes de cada vencedor; procura `**/RESULTADOS/vencedores.csv` e grava em `ANALISE_PRE_EVENTO/`. | Schema legado; só alimenta `treinar_preditor.py`. |
| `treinar_preditor.py` | RandomForest com 5 variáveis espectrais, negativos = surrogates com fase embaralhada. Grava `modelo_pac.pkl`. | Legado e otimista (surrogates como negativo). |
| `gerar_figuras_dissertacao.py` | Gera 4 PNGs em `docs/figuras/` (`--saida_dir`) **só a partir de dados reais** (reescrito em 2026-09-14; antes as figuras 1 e 3 eram curvas simuladas com `rng.normal` e a 2 e a 4 tinham valores digitados à mão). Fig. 1: ROC fora-da-amostra do modelo de estado + LFP 10 s (`_oof_lfp_10s.csv`) + ultracurto (`_oof_lfp_ultracurto.csv`, `--cenario_ultracurto`). Fig. 2: taxa de acoplamento em N por comportamento em N-1, via `pac_studio.transicao_estado` + χ². Fig. 3: precisão-recall do modelo de estado, limiar ótimo por F1 e prevalência reais. Fig. 4: medianas da teta a -3/-2/-1 s (evento × controle) + Wilcoxon nas inclinações (`_trajetoria_teta_pre_evento.csv`). | Figura que depende de um CSV ausente é **pulada** com aviso (rode `validar_preditor.py` / `analisar_janelas_ultracurtas.py` antes); `--permitir_parcial` desenha a fig. 1 só com as curvas disponíveis. |

---

## 22. Scripts auxiliares (fora do fluxo automático)

Nenhum destes é chamado pelo `processa_sessao.py`.

| Script | O que faz | Observação |
|---|---|---|
| `pac_core/workspace.py` | Fonte única dos caminhos-base (`BASE_WORKSPACE`, `BASE_LAC_NOCI`, `BASE_RESULTADOS_MESTRADO`, `BASE_SCRIPT`, `BASE_RESULTADOS` = `SCRIPT/resultados`, `BASE_FIGURAS` = `SCRIPT/figuras`); o workspace é a pasta-mãe de `SCRIPT/`, sobrescrevível pela variável `ACOPLAMENTO_BASE`. `figuras(*subpastas)` cria e devolve uma subpasta de `SCRIPT/figuras/` (destino padrão dos PNGs de diagnóstico). `localiza_ns2(arquivo, dica=None)` acha um `.ns2` pelo nome em qualquer pasta de `LAC_NOCI/` (basal ou pós-infusão), desempata pelo rato (`MTESCnn` da `dica`, normalmente a coluna `sessao`) e levanta `ValueError` se continuar ambíguo — substitui os 5 resolvedores copiados que existiam no preditor/consolida e perdiam em silêncio os 30 vencedores LAC. Sem CLI. | Nunca escrever `C:\`/`D:\` direto num script novo — importar daqui. |
| `auditorias/linha_noise_kuhn.py` | Limpa ruído de linha (60 Hz e harmônicos) do PSD antes do FOOOF, no estilo de Kuhn et al. 2026: `remove_pico_gaussiana` (subtrai só a gaussiana do pico, em log10) e `remove_faixa_1f` (cópia fiel do `rem_noise.m`: troca ±2 Hz pela curva 1/f). `aplica_modo`: `gaussiana`, `cirurgica` (só 60 Hz) ou `hibrido`. Importado por `utils_harmonico.py`, `audita_harmonico.py` (`--modo_preprocesso`, padrão `hibrido`) e `comodulogram_interativo.py`. | Só para ruído de linha — nunca para tirar harmônico de teta (tornaria o `audita_harmonico.py` circular). O dicionário `MODOS` aplica 4 harmônicos (60–240 Hz); o docstring fala em 3. Regressão: `tests/test_fooof_linha_preprocess.py`. |
| `etapa1_triagem/preprocessa_referencia_diferencial.py` | Referência diferencial automática sem mapa anatômico: escolhe N canais "silenciosos" (`seleciona_pool_referencia`), tira a média (`constroi_referencia`) e subtrai (`aplica_referencia_diferencial`). | `deteccao_ripple.py` chama a subtração sempre, mas com `sinal_referencia=None` — e nenhum código de produção passa referência. **Na prática não é aplicada.** Experimental; só `tests/test_referencia_diferencial.py` usa o pool. |
| `etapa1_triagem/inspeciona_evento.py` | Figura de 3 painéis (bruto, banda de ripple 150–250 Hz, sharp-wave < 30 Hz) do primeiro candidato a ripple. `--arquivo --canal --limiar_dp --duracao_min_ms --saida_dir`; sem argumentos reproduz o caso original (MTESC04 sessão 1, arquivo 002, `chan1`, 3 DP, 15 ms). Grava em `SCRIPT/figuras/inspecao_ripple/`. | Diagnóstico avulso, citado no README para validar candidatos de 15–20 ms. |
| `etapa5_exploracao/exploracao_minuto.py` | Linha do tempo minuto a minuto: PSD com bandas, MI z dos 3 pares (passo 5 s, **100** surrogates), potência teta e `ratio_hfo_gamma`. `--pasta_ns2 <pasta> [--canal chan20] [--saida_dir ...] [--top_n 5]` (`--top_n 0` usa o canal fixo). | Supõe 5 min por `.ns2`. Exploração pontual. |
| `utilitarios/extrair_picos.py` | Reverifica o pico de cada vencedor: refaz o comodulograma (notch 60 Hz, 200 surrogates, semente 42) e só aceita o pico se a célula passar no BH-FDR do mapa inteiro, buscando dentro da banda (4–8 × 30–80 Hz). Compara pico antigo × novo. Corrige o argmax cru do antigo `atualizar_picos.py`. Aceita o schema atual (`janela_ini_s`/`janela_fim_s`) e o legado (`inicio_s`/`fim_s`); `--pasta_ns2` é opcional (sem ela, cada `.ns2` é achado por `localiza_ns2`). | Roda direto no `OURO_PURIFICADO_v2`. Custo: 1 comodulograma completo com 200 surrogates por linha. |
| `utilitarios/gera_plot_fooof.py` | Figura de 3 painéis: espectro bruto × ajuste FOOOF v2 (knee, teta e gama separados). `--arquivo --canal --ini --fim --ctx_s --saida`; sem argumentos reproduz o caso original (MTESC04 s1, 003, `chan10`, 40–50 s). Grava em `SCRIPT/figuras/fooof/`. "Knee válido" é calculado (`_knee_valido`), não mais texto fixo; eixo Y sai dos dados. Limpa o ruído de rede no PSD antes do FOOOF com `linha_noise_kuhn.aplica_modo` (`--modo_linha`, padrão `hibrido`, o mesmo de `audita_harmonico.py`: só repõe ±2 Hz em 60/120/180/240 Hz pela 1/f local, e só se o pico for estreito, 0,5–2 Hz — gama largo perto de um harmônico sobrevive). Faixas tratadas aparecem hachuradas; PSD cru em cinza no painel A; pico a ≤ 3 Hz de um harmônico sai em cinza "rede? conferir", não como gama. | Figura de apresentação. Antes da limpeza, o caso padrão marcava o harmônico de 120 Hz como pico gama (119,7 Hz). `--modo_linha nenhum` reproduz a figura antiga. |
| `utilitarios/gerar_relatorio_pdf.py` | PDF (reportlab) com tabela por comportamento, tabela de todos os vencedores (ordem de Z) e uma página por vencedor com o comodulograma (coluna `png`, localizado sob `RESULTADOS_MESTRADO`). `--csv` (padrão `resultados/candidatos_vencedores_OURO_PURIFICADO_v2.csv`), `--saida`, `--pasta_figuras`. Título, data, ratos e condições saem do CSV. | Reescrito para o schema atual em 2026-09-14 (o antigo esperava `z`/`pico`/`classificacao`/`figura`). |
| `utilitarios/plot_basal_results.py` | 4 PNGs descritivos do dataset mestre (contagem por canal/par, boxplot de z por canal, z × MVL, curtose). `--csv` (padrão `resultados/dataset_mestre_final.csv`), `--condicao` (padrão `basal`; `''` = todas), `--saida_dir` (padrão `SCRIPT/figuras/dataset_mestre/<condicao>/`). | — |

---

## 23. `tests/`

Todos são scripts standalone: rode `python tests/<arquivo>.py`. **Rodar `pytest` não testa nada** — só `test_synthetic_harmonico.py` tem uma função `test_*`, e ela não tem `assert`.

| Arquivo | O que verifica | Dados |
|---|---|---|
| `test_pac_metrics_parity.py` | Paridade bit a bit de `pac_core.pac_metrics` contra cópias congeladas das implementações legadas + valores literais pré-migração. **Gate de qualquer mudança no núcleo.** Sai com código 1 se falhar. | Sintéticos |
| `test_pac_studio.py` | `pac_studio`: paridade bit a bit com `pac_core`, contraste em sinal sintético e reprodução exata de H2 (n = 2715, χ² = 41,26). Sai com código 1 se falhar. | Sintéticos + dataset mestre real |
| `test_synthetic_harmonico.py` | Lógica razão harmônica + PLV do `audita_harmonico` em cenários sintéticos (harmônico puro, acoplamento genuíno, coincidência próxima, fundo aperiódico). Pula de propósito o filtro de erro do FOOOF. | Sintéticos; só imprime |
| `test_fooof_pico_budget.py` | Quanto o erro do FOOOF e a detecção do pico de teta mudam com `max_n_peaks` = 1, 2, 4. | Sintéticos; benchmark impresso |
| `test_fooof_linha_preprocess.py` | Regressão do `linha_noise_kuhn.aplica_modo` (spike de 50/60 Hz + artefato motor); esperado: erro do cenário D de 0,265 para < 0,15. | Sintéticos; só imprime |
| `test_falso_positivo_hfo.py` | Monte Carlo (100 mil sorteios, semente 42) da taxa de falso positivo de razão harmônica teta 8 Hz × HFO 150–250 Hz — origem dos ~31% citados no README. | Sintéticos |
| `smoke_comod_interativo.py` | Teste rápido do `comodulogram_interativo` com 30 s sintéticos (teta 8 Hz modulando gama 60 Hz); gera os PNGs numa pasta temporária. O docstring promete checar z > 3, mas só imprime. | Sintéticos |
| `test_integracao_pipeline.py` | Roda `processa_sessao.py --min_janelas 1000` na sessão 1 real do MTESC04 — o limiar impossível executa só os Estágios 0 e 1. | **Dados reais** (sai com erro se a pasta não existir) |
| `test_referencia_diferencial.py` | Ripples com e sem referência diferencial (pool de 4/8/16 canais + leave-one-out) nos `.ns2` 002 e 003 do MTESC04 sessão 1 (`--arquivos` troca). Grava os plots em `SCRIPT/figuras/referencia_diferencial/` (`--saida_dir`). | **Dados reais**; experimental |
| `sweep_ripple.py` | Varre duração (12–20 ms) × limiar (3,0–4,0 DP) chamando `triagem_coocorrencia.py` nos 3 `.ns2` do MTESC04 sessão 1. Roda no nível do módulo (sem `__main__`); escreve `tmp_triagem_fina.csv` na pasta atual. | **Dados reais** |

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
