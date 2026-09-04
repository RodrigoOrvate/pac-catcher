# Lições Aprendidas: Geração de Relatórios PDF

Este documento detalha erros críticos de integridade identificados durante a geração dos relatórios PDF e como evitá-los em futuras iterações do pipeline.

## 1. Sincronia Absoluta Tabela-Figura
- **Erro**: O Z-score na tabela do PDF divergia do Z-score impresso dentro da imagem da figura.
- **Lição**: O valor no `vencedores_consolidado.csv` deve ser a transcrição exata do valor renderizado na figura. Divergências, mesmo pequenas (ex: 0.1), são interpretadas como erro de pareamento ou falta de rigor.
- **Como aplicar**: O script que gera a figura e o script que gera o CSV devem usar a mesma variável de arredondamento e a mesma métrica final.

## 2. Prevenção de Pseudoreplicação (Duplicatas)
- **Erro**: Um mesmo evento (mesmo arquivo, canal e janela) aparecia múltiplas vezes com valores ligeiramente diferentes.
- **Lição**: O CSV consolidado deve garantir a unicidade do evento. Se houver múltiplas detecções para a mesma janela, apenas o representante de maior Z deve ser mantido.
- **Como aplicar**: Implementar um check de unicidade `(arquivo, canal, janela)` no `generate_consolidated_csv.py`.

## 3. Filtro Estrito de Veredito
- **Erro**: Eventos marcados como "SUSPEITO" ou com observações de "artefato de cabo" foram incluídos na contagem de vencedores.
- **Lição**: O PDF deve ser estritamente "Driven by Validated Data".
- **Como aplicar**: O `gerar_relatorio_pdf.py` deve filtrar apenas linhas onde `veredito == 'VALIDO'`. Qualquer status "Suspeito" deve ser movido para o arquivo de rejeitados.

## 4. Transparência do Funil
- **Erro**: O relatório afirmava que o evento passou por 5 camadas de defesa, mas não mostrava evidência alguma dessas etapas (FDR, Footprint).
- **Lição**: A validade de um vencedor é opaca se a evidência reside apenas em logs escondidos.
- **Como aplicar**: Futuras versões do PDF devem incluir, ou linkar, os artefatos de auditoria (ex: mapa de FDR) para cada vencedor.

---
**Objetivo**: Evitar que a integridade científica do estudo seja questionada por erros de formatação ou pareamento de dados.
**Ação**: Revisar os scripts de consolidação de CSV para incluir validações de unicidade e sincronia numérica antes de disparar o gerador de PDF.
=== RESULTADO DO TESTE SINTETICO ===

Cenario A (Harmonico): PLV=0.991, ordem=3x, cf_fooof=7.98Hz -> PASSOU
Cenario B (Genuino): PLV=N/A, ordem=N/A, cf_fooof=7.98Hz -> PASSOU
Nota: limiar erro_ajuste=0.15 rejeitou sinais sinteticos (erro=0.2589).
Nota: PLV com bw=2Hz funciona corretamente (antes 0.5Hz dava PLV=0.078)
Nota: Tolerancia relativa (10% de f_theta) melhor que fixa 1.5Hz
Nota: Nenhum filtro automatico aplicado - auditoria rotula, nao exclui
