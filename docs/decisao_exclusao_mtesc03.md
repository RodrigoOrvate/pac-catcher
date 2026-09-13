# Justificativa para exclusão do MTESC03 da análise final

## Parágrafo (pronto para a dissertação)

O animal MTESC03 foi excluído das análises finais de acoplamento fase-amplitude por apresentar um rendimento de candidatos ordens de grandeza abaixo dos demais sujeitos, mesmo após verificação exaustiva de que a totalidade dos seus dados brutos havia sido efetivamente processada. Das quatro sessões basais gravadas para este animal, o pré-filtro de triagem (Estágio 1, aplicado uniformemente a todos os sujeitos) produziu apenas 757 janelas-candidato por canal e par de frequências, contra 3.614 no MTESC04 e 3.022 no MTESC05 — uma fração de aproximadamente um quinto do rendimento médio dos demais animais, apesar de duração de gravação comparável. Desse total, apenas 77 janelas sobreviveram ao portão estatístico de refinamento ("candidato robusto"), e, ao final da cadeia completa de validação (correção de múltiplos harmônicos de ruído de rede, auditoria de robustez paramétrica e purificação por comodulograma), restaram apenas 2 eventos na lista final de 190 candidatos consolidados (1,1%, frente a 147 do MTESC04 e 41 do MTESC05). Ambos os eventos remanescentes do MTESC03 apresentaram desalinhamento entre a frequência oficialmente reportada e o pico real de força estatística no mapa bidimensional fase×amplitude (100% de desalinhamento, contra 20,4% no MTESC04 e 24,4% no MTESC05), e apenas um deles teve sua célula oficial confirmada pelo controle de taxa de falsa descoberta (FDR). Antes de se concluir pela exclusão, investigou-se a hipótese de que esse baixo rendimento decorresse de falha no processamento: duas das quatro sessões do MTESC03 estavam de fato ausentes do conjunto de dados agregado. A auditoria revelou que uma dessas ausências decorria de uma colisão de nomenclatura no script de processamento em lote, que fez com que os resultados dessa sessão fossem silenciosamente sobrescritos pelos de uma sessão homônima de outro animal (MTESC05); corrigido o script e reprocessada a sessão a partir dos dados brutos originais, o resultado confirmado foi de zero candidatos válidos, idêntico ao da segunda sessão ausente (cuja ausência já havia sido verificada como um resultado legítimo, não uma falha). Descartada a hipótese de perda de dados, conclui-se que o LFP do MTESC03 genuinamente não sustenta, nos critérios adotados neste estudo, um acoplamento teta-gama detectável, o que é consistente com problemas de qualidade de sinal ou posicionamento do eletrodo específicos deste animal, e não com um artefato do pipeline de análise.

## Evidência quantitativa de suporte

| Métrica | MTESC03 | MTESC04 | MTESC05 |
|---|---|---|---|
| Janelas candidatas (Estágio 1, canal×par×janela) | 757 | 3.614 | 3.022 |
| "Candidato robusto" (pós-refino estatístico) | 77 | 1.285 | 343 |
| Eventos na lista final (`OURO_PURIFICADO_v2.csv`, N=190) | 2 (1,1%) | 147 (77,4%) | 41 (21,6%) |
| % dos eventos finais com desalinhamento comodulograma | 100% (2/2) | 20,4% (30/147) | 24,4% (10/41) |
| % dos eventos finais com célula oficial significativa (FDR) | 50% (1/2) | 33,3% (49/147) | 41,5% (17/41) |

## Nota metodológica: a lacuna de processamento investigada e descartada

Das 20 sessões brutas de `LAC_NOCI` (MTESC03_LAC, MTESC04_NOCI, MTESC05_LAC, MTESC05_NOCI), 4 estavam ausentes do dataset mestre agregado. A auditoria completa (arquivo a arquivo, e depois lógica do script de lote `roda_lote_mtesc.ps1`) determinou:

- **1 caso de perda de dados real**: MTESC03 `Rodada-1-02-05-2024` — o script de lote gerava o nome da pasta de saída a partir apenas do nome imediato da pasta de origem (`Rodada-1-02-05-2024_Basal antes da infusao`), sem incluir o prefixo do animal. Como `MTESC05_LAC` tem uma pasta de mesmo nome (`Rodada-1-02-05-2024\Basal antes da infusao`) e roda depois no script, seus resultados sobrescreveram os do MTESC03. Corrigido em `RESULTADOS_MESTRADO/roda_lote_mtesc.ps1` (prefixo do animal agora sempre incluído) e a sessão reprocessada — resultado: zero candidatos.
- **3 casos de resultado legítimo, não bug**: MTESC03 `Rodada-2-09-05-2024`, MTESC05 `Rodada-2-06-05-2024` e MTESC04 sessão 4 (15-07-2024) rodaram integralmente e, em nenhum dos 32 canais, o pré-filtro de co-ocorrência teta-gama/teta-HG atingiu o limiar mínimo (`min_janelas=3`) para disparar a análise pesada — resultado de dados real, não falha de execução.

Essa verificação foi necessária porque, sem ela, o baixo rendimento do MTESC03 poderia ser confundido com um artefato de processamento incompleto, o que enfraqueceria a justificativa de exclusão apresentada acima.
