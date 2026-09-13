# Veredito visual — QA cega dos 22 eventos de amostragem estratificada

**Metodologia desta análise** (pra ficar registrado antes dos vereditos, evitando repetir em cada item):

- Nunca trato um comodulograma "branco" numa única renderização como reprovação — já demonstramos empiricamente que o recompute do notebook não tem seed fixa e pode variar entre 0 e N células significativas para o MESMO evento em execuções sucessivas (ex.: canal 14/z=10.47 variou 0/8 a algumas vezes branco isolado; canal 11/z=7.95 variou 5/8 branco). O `z_pico` bruto reportado no título (mesmo sem sobreviver ao FDR de 396 células) é mais informativo que o binário branco/colorido.
- `z_dataset` no título é o `z_score_refinado` oficial (Gama-fit + FDR por família temporal da sessão, já com robustez de maioria validada com seed fixa — todos os 22 já são membros do OURO.csv). Isso é a evidência de base; o resto é auditoria qualitativa em cima disso.
- Red flags reais que procurei: (a) X (par oficial) caindo numa região BRANCA distante de qualquer célula colorida quando o mapa tem estrutura — sugere desalinhamento entre o par testado e onde a força estatística realmente está; (b) picos de MI/envelope na timeline coincidindo no tempo com espículos visíveis a olho nu no LFP bruto — assinatura de contaminação mecânica; (c) comodulograma fragmentado em muitas ilhas dispersas em vez de uma única ilha — sugere banda larga/múltiplas fontes; (d) espectro "picos puros" do FOOOF anormalmente ruidoso ao longo de toda a faixa, não só nos picos rotulados.
- `knee_valido`/`knee INVÁLIDO` no título do FOOOF **não é tratado como red flag** (já removido do critério de autenticidade nesta sessão — é diagnóstico de identificabilidade do modelo, não de PAC real).

---

## 1. MTESC04 | Locomoção | canal 14 | theta_hg | z=10.47 | 20240709-141215-001.ns2 [45-55s]

LFP denso, sem esteps/saturação visíveis. FOOOF teta: pico rotulado 5.63Hz (baixo pra Locomoção, esperado 7.5-8.5, mas dentro da faixa teta; R²=0.916). FOOOF gama/HG: R²=0.994, vários picos (53-202Hz), nenhum espetado em 60/120Hz. Comodulograma nesta execução veio branco (z_pico bruto=6.25 em ~7×90Hz) — já sabemos por 8 reruns prévios que esse evento raramente vem branco (~1/16), então é uma variação estatística normal, não reprovação. Timeline: dois picos de MI nítidos (~6s, ~9.5s) coincidindo com os dois maiores picos do envelope de teta.

**Veredito: ROBUSTO (alta confiança).** Maior z oficial da amostra, já testado exaustivamente nesta conversa (8/8 reruns não-branco, sub-janelas parcialmente robustas). Único ponto de atenção: pico de teta rotulado (5.63Hz) mais baixo que o esperado ecologicamente pra Locomoção — vale checar se o FOOOF não está pegando um pico secundário em vez do principal.

## 2. MTESC04 | Locomoção | canal 11 | theta_gamma | z=7.95 | 20240711-121046-001.ns2 [30-40s]

Já auditado exaustivamente nesta conversa: 8 reruns do comodulograma (5/8 branco, z bruto 5.5-6.8, estável em magnitude); sub-janelas de 2s (quente 5-7s z_min_nb=2.44, fria 1-3s z_min_nb=0.57) — nenhuma robusta isoladamente. FOOOF teta pico 8.40Hz (bom pra Locomoção). Nenhum red flag de contaminação de linha ou espículo mecânico óbvio no LFP.

**Veredito: LIMÍTROFE, mas provavelmente real e fraco/disperso** (não concentrado num episódio específico da janela). Não classificaria como falso positivo — a instabilidade é do teste de 396 células, não do evento em si (ver `z_score_refinado` via família temporal, que é o teste correto pra esse par).

## 3. MTESC04 | Locomoção | canal 11 | theta_hg | z=6.67 | 20240709-141215-002.ns2 [290-300s]

LFP denso, sem espículo óbvio. FOOOF teta: pico 7.81Hz (bom pra Locomoção). FOOOF gama/HG sem pico espetado em 60/120Hz. Comodulograma branco nesta execução, z_pico bruto=5.10 (dentro da faixa esperada de variação pra um evento z=6.67). Timeline: picos de MI concentrados perto de 8-9s coincidindo com picos de envelope de teta.

**Veredito: PROVÁVEL.** Sem red flags; z bruto do recompute consistente com o oficial dentro da margem de variação já estabelecida.

## 4. MTESC04 | Repouso | canal 30 | theta_hg | z=11.50 | 20240708-123605-001.ns2 [160-170s]

**Maior z oficial de todos os 22 e o comodulograma mais robusto visualmente da amostra**: mostra uma ilha real com célula de z alto (~8-9, amarelo) numa região compacta (fase ~7-10Hz × amp ~55-130Hz), coincidindo com a vizinhança do X. FOOOF teta com dois picos rotulados (5.24 e 7.78Hz) — ambíguo, mas nenhum dos dois é implausível. Pequenos espículos visíveis no traço bruto (~4.5s, ~8.7s) mas nada dramático, e a timeline não mostra o MI concentrado especificamente nesses instantes.

**Veredito: ROBUSTO (alta confiança).** O caso mais forte de todos os 22 — z oficial mais alto E comodulograma com estrutura clara nesta execução.

## 5. MTESC04 | Repouso | canal 16 | theta_hg | z=9.92 | 20240711-121046-003.ns2 [205-215s]

**Atenção**: o comodulograma tem estrutura, mas MUITO fragmentada — um bloco grande e forte (z~15-17) no canto fase 12.5-14.5Hz × amp 235-255Hz, bem distante do X oficial (fase~4Hz × amp~80Hz), mais vários blocos pequenos dispersos (fase 5.5-11Hz, amps 45-140). O par oficial (X) cai numa região só fracamente colorida (z baixo, teal), não no bloco mais forte do mapa. FOOOF teta pico 7.19Hz, bom.

**Veredito: ATENÇÃO / possível contaminação de banda larga.** Padrão de múltiplas ilhas dispersas, com a força estatística real do mapa concentrada longe do par oficialmente reportado — não é uma reprovação automática (z oficial de 9.92 é alto), mas esse é um dos casos onde vale mais investigação antes de considerar campeão.

## 6. MTESC04 | Repouso | canal 12 | theta_hg | z=7.15 | 20240711-121046-001.ns2 [100-110s]

LFP com um pequeno espículo (~2s). FOOOF teta pico 7.53Hz, bom. Comodulograma branco, z_pico bruto=4.60 — mais baixo que o oficial (7.15), mas ainda dentro do que já vimos de variação run-a-run. Timeline com picos de MI coincidindo com picos de envelope perto de 1s e 8.5s.

**Veredito: PROVÁVEL, com força um pouco abaixo do esperado no recompute.** Sem red flag estrutural forte.

## 7. MTESC05 | Locomoção | canal 12 | theta_hg | z=7.95 | 20240502-141712-003.ns2 [50-60s]

FOOOF teta pico 6.89Hz (levemente baixo pra Locomoção). Comodulograma branco, z_pico bruto=3.21 — **discrepância grande** frente ao z oficial de 7.95 (quase a metade). Timeline com MI disperso ao longo de toda a janela, sem concentração clara.

**Veredito: LIMÍTROFE, atenção.** A maior discrepância entre z oficial e z bruto do recompute pontual entre os "campeões" (z=7.95 oficial vs. z=3.21 aqui) — isoladamente não é decisivo (variação run-a-run), mas dado que MTESC05/MTESC03 já mostraram históricamente mais fragilidade que MTESC04 nesta auditoria, esse merece um segundo teste com seed fixa antes de aceitar como campeão de MTESC05.

## 8. MTESC05 | Mov. Cabeça | canal 16 | theta_hg | z=6.39 | 20240502-141712-002.ns2 [55-65s]

**LFP com dois espículos nítidos e visíveis** no traço bruto (~2s, ~3.5s). Comodulograma com estrutura: célula muito forte no X (z~12-13, mais alta que muitos "campeões" oficiais!) mais alguns blocos dispersos secundários. Timeline: MI concentrado fortemente só no final da janela (7.5-8.5s), coincidindo com o único pico tardio e mais alto do envelope de teta — o resto da janela tem MI baixo apesar de picos médios de teta antes disso.

**Veredito: ATENÇÃO — suspeita de contaminação de EMG/mecânica.** É Movimento de Cabeça, tem espículos visíveis no LFP bruto, e o acoplamento parece concentrado num único episódio tardio, não distribuído — exatamente o padrão que o checklist original identificou como risco nessa categoria comportamental. z do comodulograma paradoxalmente alto não compensa esse padrão temporal suspeito.

## 9. MTESC03 | Mov. Cabeça | canal 4 | theta_hg | z=7.18 | 20240504-115659-001.ns2 [70-80s]

Espículo visível no LFP bruto (~8.6s). FOOOF teta pico 6.63Hz, R²=0.996 (muito alto). **Comodulograma com estrutura forte, mas o X oficial cai numa região totalmente BRANCA** — a célula mais forte do mapa (z~14-15) está isolada no canto fase 10.5-14.5Hz × amp 230-255Hz, bem longe do par oficial (fase~8Hz × amp~110Hz).

**Veredito: SUSPEITO — desalinhamento claro entre par oficial e força estatística real.** O evento pode ter *algum* acoplamento real, mas não necessariamente no par (fase,amp) que o pipeline reportou — a força de verdade está numa banda completamente diferente (quase HFO em fase mais alta). Não classificaria como robusto sem investigar essa discrepância.

## 10. MTESC03 | Repouso | canal 2 | theta_gamma | z=6.75 | 20240504-115659-003.ns2 [235-245s]

Dois espículos visíveis no LFP bruto (~4s, ~6.5-7s). **Comodulograma muito fragmentado**: múltiplas ilhas pequenas espalhadas por toda a grade, incluindo um bloco muito forte (z~20-26) no canto fase 10.5-14.5 × amp 230-255 — de novo bem longe do X oficial (fase~4Hz × amp~40Hz, que cai numa das ilhas mais FRACAS do mapa).

**Veredito: SUSPEITO — mesmo padrão do item 9 (mesmo rato, mesma sessão).** Múltiplas ilhas dispersas + par oficial na região de menor força é o padrão clássico de contaminação de banda larga/múltiplas fontes, não coupling focal único. MTESC03 aparece pela segunda vez consecutiva com esse padrão.

## 11. MTESC04 | Mov. Cabeça | canal 14 | theta_hg | z=9.74 | 20240711-121046-002.ns2 [270-280s]

FOOOF teta com dois picos rotulados (5.70, 7.94Hz). Comodulograma com estrutura razoavelmente concentrada: célula mais forte (z~13) em fase 6.5-7.5×amp 50-60, com o X posicionado entre dois blocos secundários próximos (fase~7×amp~135). Sem espículo isolado muito óbvio no LFP bruto.

**Veredito: PROVÁVEL.** Apesar de ser Mov. Cabeça, o comodulograma tem estrutura razoavelmente concentrada e nenhum red flag óbvio de contaminação mecânica na timeline.

## 12. MTESC04 | Mov. Cabeça | canal 13 | theta_gamma | z=9.36 | 20240711-121046-001.ns2 [45-55s]

FOOOF teta pico 8.29Hz, bom. FOOOF gama/HG com vários picos rotulados (72-197Hz), nenhum em linha elétrica. Comodulograma branco nesta execução, z_pico bruto=5.67. Timeline com múltiplos picos de MI ao longo de quase toda a janela (5-9.5s), padrão relativamente distribuído, não um único espículo isolado.

**Veredito: PROVÁVEL.** Padrão de MI distribuído (não concentrado) é tranquilizador; sem espículo LFP óbvio apesar de ser Mov. Cabeça.

## 13. MTESC04 | Locomoção | canal 30 | theta_hg | z=3.98 | 20240708-123605-001.ns2 [105-115s]

Evento limítrofe por natureza (z baixo). FOOOF teta ambíguo (dois picos: 5.68, 8.16Hz). Comodulograma branco, z_pico bruto=4.02 — **consistente** com o z oficial baixo (sem grande discrepância pra cima ou pra baixo, ao contrário de outros limítrofes).

**Veredito: LIMÍTROFE, mas consistente/honesto.** Fraco porque realmente é fraco, não porque há artefato aparente. Sem red flags estruturais.

## 14. MTESC05 | Locomoção | canal 6 | theta_hg | z=3.89 | 20240502-141712-002.ns2 [160-170s]

Também limítrofe por natureza. FOOOF teta ambíguo (6.91, 7.86Hz). Comodulograma branco, z_pico bruto=3.44 — de novo consistente com o oficial baixo.

**Veredito: LIMÍTROFE, consistente.** Mesmo padrão do item 13: fraco mas sem sinal de artefato.

## 15. MTESC03 | Repouso | canal 12 | theta_hg | z=3.51 | 20240504-115659-003.ns2 [255-265s]

Espículo visível no LFP bruto (~4.2s). Comodulograma disperso: vários blocos pequenos e fracos (z baixo-moderado) espalhados, com o X oficial numa região branca, não coincidindo com nenhum bloco.

**Veredito: LIMÍTROFE, com leve suspeita.** Terceiro evento consecutivo de MTESC03 com padrão de dispersão/desalinhamento — reforça a suspeita já registrada de que os dados desse rato têm contaminação sistemática mais alta que os outros dois.

## 16. MTESC04 | Grooming | canal 30 | theta_hg | z=3.94 | 20240708-123605-002.ns2 [260-270s]

FOOOF teta ambíguo (5.24, 7.78Hz). Comodulograma com estrutura: bloco amarelo forte perto de fase 7.5-8.5×amp 55-65, mas **deslocado do X** (que está em amp~110, bem mais alta) — mais dois blocos verdes secundários em amplitudes 115-135Hz.

**Veredito: LIMÍTROFE com leve desalinhamento.** Não tão grave quanto os itens 9/10 (o bloco forte está pelo menos na mesma vizinhança de fase do X), mas ainda um desvio de amplitude que vale nota.

## 17. MTESC04 | Mov. Cabeça | canal 12 | theta_gamma | z=3.89 | 20240711-121046-002.ns2 [255-265s]

Espículo visível no LFP bruto (~1.3s). Comodulograma com estrutura no canto de amplitude muito alta (fase 8.5-9.5×amp 245-250, quase HFO) e blocos secundários em amp 213-245 — **o X oficial (fase~5×amp~35, banda gama baixa) cai numa região branca**, bem distante de toda a estrutura colorida (que está toda em banda de amplitude muito mais alta).

**Veredito: SUSPEITO — desalinhamento forte entre par oficial (baixa amplitude) e onde está a força estatística real (HFO).** Combinado com z oficial já baixo (3.89) e ser Mov. Cabeça, esse é um dos candidatos mais fracos da amostra.

## 18. MTESC04 | Grooming | canal 12 | theta_hg | z=9.03 | 20240709-141215-003.ns2 [150-160s]

**Red flag mais forte da amostra inteira.** LFP bruto com dois espículos grandes e nítidos (~1.5s, ~4.2s). Na timeline, os dois maiores picos de MI (muito acima do resto da série) e os dois maiores picos do envelope de teta caem **exatamente** nesses mesmos dois instantes. Comodulograma branco nesta execução, z_pico bruto=3.51 (bem abaixo do oficial 9.03).

**Veredito: SUSPEITO DE TRANSIENTE MECÂNICO**, apesar de ter passado no `robustez_parametros.py` (que testa a janela de 10s inteira, sem localizar especificamente esses dois picos). Esse é o candidato que eu pessoalmente reexaminaria com mais prioridade — o padrão visual (espículo no LFP bruto → pico simultâneo em envelope de teta E em MI) é a assinatura clássica descrita no checklist original do usuário para transiente mecânico disfarçado de PAC.

## 19. MTESC04 | Sniffing | canal 12 | theta_hg | z=8.80 | 20240709-141215-003.ns2 [70-80s]

Um espículo visível no LFP bruto (~3.4s), mas o maior pico de MI da timeline ocorre depois (~8.5-9s), não coincidindo com o espículo. Comodulograma branco, z_pico bruto=4.15.

**Veredito: PROVÁVEL, com atenção mínima.** Diferente do item 18 (mesmo arquivo/canal, janela diferente), aqui o espículo visível NÃO parece estar dirigindo o pico de MI — sinal mais tranquilizador.

## 20. MTESC04 | Mov. Cabeça | canal 15 | theta_hg | z=8.18 | 20240711-121046-002.ns2 [190-200s]

LFP bruto com vários espículos negativos nítidos (~2s, ~4.5s, ~6.7s). Comodulograma com estrutura razoavelmente concentrada perto do X (fase~6×amp~85, z~6-7, o mais forte do painel). Timeline: dois picos de MI fortes coincidindo tanto com os dois maiores picos de envelope de teta quanto, de perto, com os espículos visíveis no LFP bruto.

**Veredito: ATENÇÃO — mesmo padrão de suspeita do item 18, mais brando.** Coincidência temporal entre espículo no LFP bruto e pico de MI, reforçada por ser Movimento de Cabeça — mas aqui o comodulograma pelo menos concentra no par oficial (diferente dos itens 9/10/17 que desalinham).

## 21. MTESC04 | Locomoção | canal 13 | theta_hg | z=5.92 | 20240709-141215-002.ns2 [205-215s]

Este é o evento de footprint difuso (muitos canais coativos). LFP sem espículo óbvio. Comodulograma com estrutura boa: **célula amarela intensa (z~11, a mais forte do painel) caindo exatamente no X**, mais um bloco secundário isolado.

**Veredito: PROVÁVEL, apesar do alerta de footprint.** O footprint amplo (informativo, não é corte) não se traduz aqui em desalinhamento espectral — o par oficial está claramente dentro da região mais forte do comodulograma. Bom sinal.

## 22. MTESC03 | Repouso | canal 32 | theta_gamma | z=5.69 | 20240506-122038-002.ns2 [285-295s]

**Segundo red flag mais forte da amostra.** LFP bruto visivelmente mais espiculado que todos os outros 21 (pelo menos 4-5 espículos grandes e nítidos ao longo da janela). Curva de "picos puros" do FOOOF gama/HG anormalmente ruidosa (não só nos picos rotulados). Comodulograma extremamente fragmentado, com uma célula isolada de z absurdamente alto (~60+, muito acima da escala típica de 5-15 dos outros eventos) numa banda de amplitude extrema (230-250Hz), longe do X oficial.

**Veredito: SUSPEITO — múltiplos sinais convergentes de contaminação mecânica/banda larga.** Terceiro evento de MTESC03 na amostra, e o mais problemático dos três. Esse eu não aceitaria como vencedor sem uma reauditoria completa (transiente + banda larga) específica.

---

## Resumo rápido (minha classificação, antes de qualquer comparação)

| # | Evento | Veredito |
|---|---|---|
| 1 | MTESC04 canal14 Loc z=10.47 | ROBUSTO |
| 2 | MTESC04 canal11 Loc theta_gamma z=7.95 | LIMÍTROFE (fraco/disperso) |
| 3 | MTESC04 canal11 Loc theta_hg z=6.67 | PROVÁVEL |
| 4 | MTESC04 canal30 Repouso z=11.50 | ROBUSTO |
| 5 | MTESC04 canal16 Repouso z=9.92 | ATENÇÃO (banda larga?) |
| 6 | MTESC04 canal12 Repouso z=7.15 | PROVÁVEL |
| 7 | MTESC05 canal12 Loc z=7.95 | LIMÍTROFE (discrepância) |
| 8 | MTESC05 canal16 MovCabeça z=6.39 | ATENÇÃO (EMG?) |
| 9 | MTESC03 canal4 MovCabeça z=7.18 | SUSPEITO (desalinhado) |
| 10 | MTESC03 canal2 Repouso z=6.75 | SUSPEITO (banda larga) |
| 11 | MTESC04 canal14 MovCabeça z=9.74 | PROVÁVEL |
| 12 | MTESC04 canal13 MovCabeça z=9.36 | PROVÁVEL |
| 13 | MTESC04 canal30 Loc z=3.98 | LIMÍTROFE (consistente) |
| 14 | MTESC05 canal6 Loc z=3.89 | LIMÍTROFE (consistente) |
| 15 | MTESC03 canal12 Repouso z=3.51 | LIMÍTROFE (leve suspeita) |
| 16 | MTESC04 canal30 Grooming z=3.94 | LIMÍTROFE (leve desalinho) |
| 17 | MTESC04 canal12 MovCabeça theta_gamma z=3.89 | SUSPEITO (desalinhado) |
| 18 | MTESC04 canal12 Grooming z=9.03 | SUSPEITO (transiente) |
| 19 | MTESC04 canal12 Sniffing z=8.80 | PROVÁVEL |
| 20 | MTESC04 canal15 MovCabeça z=8.18 | ATENÇÃO (transiente leve) |
| 21 | MTESC04 canal13 Loc z=5.92 | PROVÁVEL |
| 22 | MTESC03 canal32 Repouso theta_gamma z=5.69 | SUSPEITO (contaminação forte) |

**Contagem**: 7 PROVÁVEL/ROBUSTO sem ressalva, 6 ATENÇÃO/LIMÍTROFE com ressalva leve, 5 LIMÍTROFE consistente (fraco mas honesto), 5 SUSPEITO (red flag visual claro — dessas, 3 são MTESC03).

Fico no aguardo da sua comparação.
