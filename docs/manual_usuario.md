# ThetaGamma-Studio — Manual do Usuário

Guia passo a passo do programa desktop, sem precisar programar: rodar as
etapas do pipeline, inspecionar resultados e abrir o anotador de vídeo.

## 1. Instalação (uma vez só)

1. Instale o [Python 3.11+](https://python.org) (marque "Add to PATH" no
   instalador do Windows).
2. Abra um terminal (PowerShell) em `D:\acoplamento_theta-gamma\SCRIPT`.
3. Instale as dependências:
   ```
   pip install -r requirements.txt
   ```
4. Confirme que os dados estão em `D:\acoplamento_theta-gamma\LAC_NOCI\` e
   `D:\acoplamento_theta-gamma\RESULTADOS_MESTRADO\` (ao lado da pasta
   `SCRIPT\`, não dentro dela). Se o workspace estiver em outro lugar,
   defina a variável de ambiente `ACOPLAMENTO_BASE` antes de abrir o
   programa (ver `pac_core/workspace.py`).

## 2. Abrir o programa

No terminal, dentro de `SCRIPT\`:
```
python dashboard_desktop/app.py
```

A janela abre maximizada. Se fechar a janela com alguma etapa do pipeline
rodando, o programa pergunta antes e encerra os processos junto — não
deixa nada rodando escondido.

## 3. Seção Pipeline

Cada sub-aba roda um script de `pipeline/` (os mesmos usados pelo terminal).
Preencha o formulário, clique no botão azul e acompanhe o log ao vivo no
painel escuro. **Parar** encerra a etapa e tudo que ela disparou.

| Sub-aba | O que roda | Quando usar |
|---|---|---|
| Sessão Completa | `processa_sessao.py` | Processar uma pasta de sessão inteira (triagem → refinamento → comodulograma → auditorias), canal a canal. Leva de 10 a 90 min. |
| Triagem | `triagem_pac.py` | Varredura rápida de uma pasta/canal com parâmetros próprios. O botão **Demo** roda um teste sintético em segundos. |
| Refinamento e Filtros | `refina_candidatos.py` | Testar outros limiares (z pré-filtro, FDR q, surrogates) sobre uma triagem já feita. Ao terminar, mostra o gráfico de vereditos. |
| Comodulograma (lote) | `comodulogram.py` | Gerar os PNGs de comodulograma dos candidatos de um CSV refinado. |
| Auditorias | `audita_skewness.py`, `audita_harmonico.py`, `audita_harmonico_hfo.py` | Rodar uma auditoria de qualidade sobre um CSV. A dica abaixo do tipo diz qual CSV cada uma espera. |
| Agregação | `agrega_resultados.py` → `enriquece_dataset_mestre.py` → `junta_comportamento.py` → `consolida_vencedores.py` | Reconstruir o dataset mestre em 4 passos, na ordem. |

Cuidados:
- As saídas padrão da Agregação gravam em arquivos `*_novo.csv`, pra não
  sobrescrever os CSVs de produção sem querer. Toda etapa pergunta antes de
  sobrescrever um arquivo que já existe.
- O agregador varre a pasta recursivamente e ignora pastas que começam com
  `_` (`_CONTAMINADO_nao_usar/`, `_LAC_basal_antigo_nao_usar/`). Pra agregar
  só um grupo, aponte pra subpasta (ex.: `basal\NOCI`).
- Transientes e footprint não estão na aba Auditorias: a linha de comando
  delas é presa a uma sessão específica. Rode pelo terminal.

## 4. Seção Análise de Dados

Inspeção dos resultados já processados, em 6 sub-abas:

- **Inspeção de Canais** — traçado de LFP multicanal com notch, passa-faixa,
  remoção de deriva, canal de referência e ganho (mexem só na exibição).
- **Comodulograma** — mapa fase×amplitude de uma janela (200 surrogates;
  pode levar alguns segundos).
- **Portões de Qualidade** — se um candidato passou nos portões de refino
  (estatístico, skewness, harmônico, transiente) e o footprint (informativo).
- **Galeria de Vencedores** — os 190 candidatos "ouro" com imagem e
  metadados, paginada.
- **Comportamento** — acoplamento por categoria comportamental, matriz de
  transição N-1→N com o qui-quadrado da Hipótese H2, e o painel Basal vs
  Nociceptina vs Lactato (aparece quando o dataset mestre tiver as colunas
  de grupo/condição da infusão).
- **Closed-Loop / Replay** — simula offline a previsão de P(PAC) ao longo de
  um `.ns2`. **Não liga hardware nenhum** (ver aviso na própria aba).

## 5. Seção Anotador de Vídeo

Abre a ferramenta de anotação de comportamento
(`pipeline/comportamento/anotador_comportamento.py`) numa janela separada e
mostra o progresso de anotação por sessão. Depois de anotar e fechar o
anotador, clique em **Atualizar progresso**.

## 6. Problemas comuns

| Sintoma | Causa provável |
|---|---|
| "Pasta não encontrada" / "Arquivo não encontrado" | Caminho digitado errado — use o botão Procurar... |
| Etapa termina com "Falhou (exit=2)" | Algum parâmetro com valor inválido (o próprio script explica no log, na última linha) |
| Log parado por muito tempo | Algumas etapas imprimem só ao fim de cada canal/arquivo; o status "Rodando..." indica que continua |
| Comodulograma muito lento | Janela de tempo grande demais — tente 10-20 s primeiro |
| Painel farmacológico sempre "aguardando" | O dataset mestre ainda não tem as colunas grupo/condição (regenere pela Agregação) |

## 7. Validação em máquina limpa (checklist)

Antes de apresentar pra banca ou compartilhar com outra pessoa do
laboratório, confira numa máquina sem nada do projeto instalado:

- [ ] `pip install -r requirements.txt` termina sem erro.
- [ ] `python dashboard_desktop/app.py` abre a janela sem erro no terminal.
- [ ] Pipeline → Triagem → **Demo** termina com "Concluído (exit=0)".
- [ ] Consegue carregar um `.ns2` real em Análise de Dados → Inspeção de Canais.
- [ ] `python tests/test_pac_studio.py` passa (paridade bit a bit com o
      pipeline em lote).
- [ ] `python tests/test_pac_metrics_parity.py` continua passando.
