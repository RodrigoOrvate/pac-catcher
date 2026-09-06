# PAC Catcher 🎣

Bem-vindo ao **PAC Catcher**! Este repositório centraliza o pipeline computacional para detectar, auditar e validar **Acoplamento Fase-Amplitude (Phase-Amplitude Coupling - PAC)** em registros eletrofisiológicos (LFP).

## 🎯 Objetivo

O foco principal do pipeline é extrair e validar com rigor estatístico os acoplamentos **Teta-Gama (30-80 Hz)**, **Teta-High Gamma (HG, 80-150 Hz)** e **Teta-HFO (150-250 Hz)**. 

Devido à vasta quantidade de falsos positivos reportados na literatura de PAC (gerados por harmônicos de fase, transientes e ruído), este pipeline não se baseia apenas no tradicional *Modulation Index*. Ele cruza estatística de deslocamento circular, paramétrica, comportamento de rede espacial e decomposição espectral avançada para assegurar que o acoplamento detectado é um fenômeno neurofisiológico genuíno.

### 📚 Fundamentos e Citações

O desenvolvimento deste fluxo é fortemente embasado em literatura de referência:
- **Separação Aperiódica/Oscilatória:** Utilizamos o algoritmo **FOOOF (SpecParam)**. A pasta `FOOOF` presente neste repositório contém os arquivos originais disponibilizados pelos autores do artigo *"Aperiodicity in Mouse CA1 and DG Power Spectra"* (Kuhn et al., 2026), enquanto o nosso pipeline faz as adaptações necessárias para LFP. Essa abordagem é essencial para High-Gamma e HFO, garantindo que a elevação do *spiking* (disparo neuronal) não seja confundida com oscilação real de banda estreita. Repositório original: [LFP_FOOOF](https://github.com/huRashidy/LFP_FOOOF).
- **Modulação Teta-HG e Teta-HFO:** As bases e os métodos de distribuição de graus de fase (Tort's Modulation Index, Polar Plots) têm fortes raízes em trabalhos como *"Theta phase modulates multiple layer-specific oscillations in the CA1 region"* de Adriano Tort.

---

## 🚀 Como Rodar o Pipeline

Abaixo estão as etapas principais do pipeline. Todo o código reside aqui. Você apenas apontará os comandos para as pastas com os seus dados `.ns2`.

> **Nota:** Para uma explicação mais aprofundada do embasamento teórico, lógica de auditorias e o que cada função faz internamente, consulte o arquivo `scripts_explicados.md`.

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
python pipeline/exploracao_minuto.py ...
```

---

### Passo 1: Triagem Estatística
Varre todas as janelas do registro e calcula o Z-score do Modulation Index (KL-MI) contra *surrogates* de deslocamento circular, buscando os 3 pares simultaneamente.
```bash
python pipeline/triagem_pac.py \
    --pasta "../SESSAO_EXEMPLO/Basal" \
    --pares theta_gamma theta_hg theta_hfo \
    --saida "../SESSAO_EXEMPLO/RESULTADOS/resultados.csv"
```
*Gera uma tabela bruta com Z-scores, p-valores empíricos e SNR espectral de todas as janelas (ex: janelas de 10s com sobreposição).*

---

### Passo 2: Refinamento e Correção FDR
Aplica o *False Discovery Rate* (FDR - Benjamini-Hochberg) em todas as janelas, descarta falsos positivos por co-ocorrência multicanal e kurtose na banda alta.
```bash
python pipeline/refina_candidatos.py \
    --csv "../SESSAO_EXEMPLO/RESULTADOS/resultados.csv" \
    --pasta_ns2 "../SESSAO_EXEMPLO/Basal" \
    --saida "../SESSAO_EXEMPLO/RESULTADOS/resultados_refinados.csv"
```
*Quem sobrevive a esse funil estatístico ganha o selo "Candidato robusto".*

---

### Passo 3: Comodulogramas e Controle de Linha
Gera os mapas bidimensionais (fase x amplitude) para os candidatos robustos. **Crucial:** aplica filtros *notch* em 60 Hz e 120 Hz para evitar que harmônicos da rede elétrica forjem acoplamento no High-Gamma.
```bash
python pipeline/comodulogram.py \
    --csv "../SESSAO_EXEMPLO/RESULTADOS/resultados_refinados.csv" \
    --pasta_ns2 "../SESSAO_EXEMPLO/Basal" \
    --saida_dir "../SESSAO_EXEMPLO/RESULTADOS/comodulogramas" \
    --notch 60 --fdr_q 0.05
```

---

### Passo 4: Auditorias Específicas
O pipeline conta com auditorias separadas para blindar os resultados contra falhas físicas e matemáticas do sinal. Para aplicar, você executa scripts da pasta `auditorias`:
- **`audita_skewness.py`**: Avalia a assimetria (formato dente de serra) do Teta para evitar geração de harmônicos de fase.
- **`audita_harmonico.py` / `audita_harmonico_hfo.py`**: Usa o algoritmo **FOOOF** para verificar se o Gama é harmônico do Teta ou se o HFO é harmônico do Gama.
- **`audita_footprint.py`**: Garante que o gerador da oscilação é focal e descarta propagação por condução de volume (co-detecção massiva nos 32 canais).
- **`diagnostico_janela.py`**: Confirma se o que estamos chamando de Teta não é na verdade respiração/sniffing do roedor (ritmo olfatório).

---

### Passo 5: Geração de Figuras Finais
Ao final, você seleciona os verdadeiros vencedores e passa para a geração das imagens (STFT, distribuição polar e LFP bruto) para apresentação, *lab meetings* ou artigo.
```bash
python pipeline/figura_apresentacao.py \
    --vencedores "../SESSAO_EXEMPLO/RESULTADOS/vencedores.csv" \
    --pasta_ns2 "../SESSAO_EXEMPLO/Basal" \
    --saida_dir "../SESSAO_EXEMPLO/RESULTADOS/figuras"
```

*(Imagens de exemplo serão adicionadas aqui em breve para facilitar a visualização).*
