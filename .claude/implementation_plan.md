# Plano Metodológico & Arquitetural: Pipeline PAC Catcher

Este documento consolida a **jornada completa da sua pesquisa**, os **fundamentos científicos sobre aperiodicidade e o filtro notch** (baseados em Kuhn et al., 2026), o **status dos pares de frequência analisados** (teta-gamma, teta-HG e teta-HFO), e o **plano de transição para o pipeline refatorado** que entrará em ação assim que a sua etapa atual de anotação comportamental for concluída.

---

## 🧭 1. A Sua Jornada até o Momento Presente

O seu projeto não nasceu pronto em um bloco rígido; ele evoluiu cientificamente de forma adaptativa, passando de uma estratégia inicial de **alta sensibilidade** para um refinamento rigoroso de **alta especificidade e contextualização ecológica**:

```
[Etapa 1: Triagem Ampla] 
   └── Varredura contínua de todo o LFP (janelas de 10s, passo 5s) para não perder nenhum candidato a PAC.
         │
         ▼
[Etapa 2: Especificidade & Blindagem Estatística]
   └── Filtragem rigorosa contra falsos positivos: FDR Benjamini-Hochberg, ajuste Gama, portão de banda larga e pegada espacial.
         │
         ▼
[Etapa 3: Decomposição Periódico vs. Aperiódico (FOOOF / SpecParam)]
   └── Introdução do modelo 1/f com joelho (knee) para garantir que aumentos de potência aperiódica (spiking/ruído) não fossem tomados como oscilações reais.
         │
         ▼
[Etapa 4: Auditoria de Harmônicos de Teta vs. Gamma/HFO]
   └── Discriminação de acoplamentos espúrios causados pela forma de onda não-senoidal (dente-de-serra) do teta através de razão harmônica + PLV (Phase-Locking Value).
         │
         ▼
[Etapa 5: Módulo Comportamental & Sincronização Fina (ESTÁGIO ATUAL)]
   └── Alinhamento neural-vídeo: descoberta do corte de 272s em 001.ns2 (glitch Blackrock), medição do gap de 30.4s pelo .nev, calibração dos offsets por bloco (001: 325s, 002: 627.5s, 003: 928.5s) e anotação manual das 39 janelas iniciais (100% preservadas).
```

### Onde você está exatamente agora:
Você está utilizando a interface gráfica do `anotador_comportamento.py` para classificar o comportamento do animal (Locomoção, Imobilidade Atenta, Sono/Descanso, Grooming, etc.) nas janelas dos blocos `002.ns2` e `003.ns2`, com os offsets temporais devidamente sincronizados ao vídeo contínuo.

---

## 🔬 2. O Artigo de Kuhn et al. (2026), o Filtro Notch e a Discussão Periódico vs. Aperiódico

### 2.1 O que Kuhn et al. abordam sobre o Filtro Notch? Eles aprovam ou não?
No artigo *"Aperiodicity in Mouse CA1 and DG Power Spectra"* (Kuhn et al., 2026) e no código-fonte original correspondente (`rem_noise.m` presente na pasta `SCRIPT/FOOOF/` do nosso projeto), os autores discutem um problema crítico:

1. **A Desaprovação do Notch IIR Tradicional na Análise Aperiódica:**
   - Filtros notch clássicos (como Butterworth ou Chebyshev aplicados no domínio do tempo) introduzem um corte extremamente abrupto e profundo ("notch") na densidade espectral de potência (PSD) em 50 Hz ou 60 Hz.
   - Para algoritmos de decomposição espectral como o **FOOOF (SpecParam)**, esse "buraco artificial" é destrutivo: o algoritmo tenta ajustar a reta/curva contínua de fundo $L(f) = b - \log(k + f^\chi)$, e a depressão artificial do notch puxa a estimativa do expoente aperiódico $\chi$ e do joelho (*knee*), distorcendo a curva de base e criando resíduos artificiais ao redor da frequência de corte.

2. **A Solução Proposta por Kuhn et al. (`rem_noise.m`):**
   - Os autores **não aplicam notch destrutivo no sinal contínuo temporal** quando o objetivo é parametrizar o espectro.
   - Em vez disso, eles realizam uma **remoção espectral cirúrgica**:
     - No espectro de potência (PSD), delimitam uma janela estreita ao redor da frequência da rede elétrica (ex: 46 a 55 Hz para 50 Hz, ou 54 a 66 Hz para 60 Hz).
     - Ajustam um modelo aperiódico local e detectam se há pico de rede.
     - **Substituem a faixa de ruído diretamente pela curva aperiódica $1/f$ pura interpolada**, ou subtraem exclusivamente o pico gaussiano em escala logarítmica ($\log_{10}$).
     - **Resultado:** O ruído da linha desaparece sem deixar a "vala" artificial que o filtro notch causaria.

### 2.2 O nosso pipeline tem essa preocupação?
**Sim, totalmente.** O repositório já conta com essa formulação nos módulos de auditoria:
- No processamento inicial de PAC contínuo (triagem rápida no domínio do tempo), usamos `aplica_notch_multiplo` (60, 120, 180, 240 Hz) para evitar que a rede module espuriamente os filtros passa-faixa de amplitude.
- No entanto, na etapa avançada de validação espectral, temos implementados em `SCRIPT/pipeline/auditorias/linha_noise_kuhn.py` e `compara_preprocesso_linha.py` os três modos inspirados diretamente no artigo:
  1. `gaussiana`: Subtrai apenas o pico gaussiano da rede em $\log_{10}$, preservando o chão aperiódico local.
  2. `cirurgica`: Réplica exata em Python do `rem_noise.m` do MATLAB de Kuhn et al. (reconstitui a banda de $\pm 2\text{ Hz}$ com a equação $1/f$ estimada).
  3. `hibrido`: Aplica a substituição cirúrgica em 60 Hz e nos harmônicos 120 Hz e 180 Hz.

---

## ⚡ 3. O Pipeline ainda está pegando Teta-HFO? Ou só Teta-Gamma e Teta-HG?

O pipeline **continua suportando e calculando Teta-HFO**. Ele foi projetado para registrar os três pares de frequências fundamentais da literatura de hipocampo de roedores:

| Identificador | Banda de Fase ($\theta$) | Banda de Amplitude ($f_\text{amp}$) | Contexto Fisiológico |
| :--- | :---: | :---: | :--- |
| **`theta_gamma`** | $4 - 8\text{ Hz}$ | $30 - 80\text{ Hz}$ | Gama clássico (CA3 $\to$ CA1), associado à navegação e exploração locomotora. |
| **`theta_hg`** | $4 - 8\text{ Hz}$ | $80 - 150\text{ Hz}$ | High-Gamma / Fast Gamma (Córtex Entorrinal $\to$ CA1), associado a processamento de memória e alta demanda cognitiva. |
| **`theta_hfo`** | $4 - 8\text{ Hz}$ | $150 - 250\text{ Hz}$ | High-Frequency Oscillations / Ripple bands, associados a estados de quietude, repouso e Sharp-Wave Ripples (SWR). |

### Onde o Teta-HFO está no código:
1. **`triagem_pac.py`**: O dicionário central `BAND_PAIRS` possui a chave `"theta_hfo"` ativa por padrão. O script calcula `mi_theta_hfo`, `z_theta_hfo` e gera a razão `ratio_hfo_gamma`.
2. **`refina_candidatos.py` e `processa_sessao.py`**: Aceitam o argumento `--pares theta_gamma theta_hg theta_hfo`.
3. **`audita_harmonico_hfo.py`**: Script dedicado especificamente a verificar se um pico de HFO ($150-250\text{ Hz}$) é uma oscilação biológica genuína ou um harmônico superior de ordem elevada (ex: $20\times$ a frequência fundamental de teta) ou vazamento de disparo unitário (*spike bleed-through*).

---

## 🚀 4. Como Você Continuará com o Pipeline Refatorado (Pós-Anotação)

Enquanto você termina de preencher os comportamentos dos blocos `002` e `003` no `anotador_comportamento.py`, preparamos a estrutura para que a retomada seja imediata e limpa:

```
[Você conclui a anotação dos vídeos]
               │
               ▼
[Passo A: Consolidação dos Comportamentos]
Executa: python pipeline/junta_comportamento.py
-> Cruza os comportamentos anotados com o dataset_mestre, gerando o arquivo com estados puros.
               │
               ▼
[Passo B: Ativação da Camada Compartilhada `pac_core`]
-> Centralização de `filtra_sinal`, `aplica_notch_multiplo` e `_mi_de_bin_idx` (Zero duplicação).
-> Códigos antigos passam a herdar do núcleo limpo sem mudar flags de comando.
               │
               ▼
[Passo C: Geração da Tabela de Vencedores Consolidados]
-> Filtro final: Candidato com PAC Significativo (FDR + Gama + FOOOF Aprovados) + Comportamento Definido.
               │
               ▼
[Passo D: Exploração Interativa (A "cereja do bolo")]
-> Executa as ferramentas interativas de PSD periódico (zoom interativo de teta e gamma isolados)
-> Timeline interativa correlacionando a trajetória do animal no vídeo com os pulsos de acoplamento.
```

### O que faremos na Fase 1 (Base Técnica):
- Criar a pasta `SCRIPT/pac_core/` com módulos testados:
  - `filtering.py`: Funções de filtro Butterworth e Notch multi-harmônico.
  - `pac_metrics.py`: Modulação por Entropia de Kullback-Leibler (Tort) e rotina de surrogates com shuffle temporal circular.
  - `io.py`: Leitura robusta de arquivos `.ns2` e salvamento blindado com `utf-8-sig`.
- Adicionar testes de integridade para garantir que os cálculos geram valores idênticos aos anteriores até a última casa decimal.

> [!NOTE]
> Você pode continuar a sua anotação com total tranquilidade. Todo o seu trabalho manual nos 39 eventos de `001` e os novos eventos de `002` e `003` estão protegidos, e a refatoração será ativada de forma modular e transparente.
