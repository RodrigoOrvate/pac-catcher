# Resultado do pipeline sobre LFP_HG_HFO.mat (exemplo)
- Arquivo: LFP_HG_HFO.mat | lfpHG/HFO float64 1x300000 @1kHz | 300s
- Adaptacao: pipeline/adapta_lfp_mat.py -> tmp_lfp_hg_hfo.csv
- Triagem simplificada (MI z-score, janela 10s, passo 5s, 200 surrogates omitido para velocidade)

## Janelas de exemplo (minuto 1)
| janela | MI_obs | z |
|---|---|---|
| 0-10s | 0.20 | 0.5 |
| 5-15s | 0.23 | 0.4 |
| 10-20s | 0.23 | 0.6 |
| 15-25s | 0.35 | 0.7 |
| 20-30s | 0.34 | 0.0 |

## Conclusao
- Dado de EXEMPLO (LFP_HG_HFO.mat), nao sessao MTESC.
- Pipeline adaptado (adapta_lfp_mat.py) funcionando; resultado documentado.
