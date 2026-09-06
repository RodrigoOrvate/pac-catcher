import sys
import os
import re
import subprocess

duracoes = [12, 14, 16, 18, 20]
dps = [3.0, 3.2, 3.5, 3.8, 4.0]
origens = {
    "exploracao_001": r"C:\acoplamento_theta-gamma\MTESC04_NOCI\MTESC04 -- 1 - infusao - 08-07-2024\Basal antes da infusao\20240708-123605-001.ns2",
    "imobilidade_002": r"C:\acoplamento_theta-gamma\MTESC04_NOCI\MTESC04 -- 1 - infusao - 08-07-2024\Basal antes da infusao\20240708-123605-002.ns2",
    "imobilidade_003": r"C:\acoplamento_theta-gamma\MTESC04_NOCI\MTESC04 -- 1 - infusao - 08-07-2024\Basal antes da infusao\20240708-123605-003.ns2"
}
script_path = r"C:\acoplamento_theta-gamma\SCRIPT\pipeline\triagem_coocorrencia.py"

resultados = []

for origem_nome, origem_path in origens.items():
    print(f"\nRodando sweep fino para {origem_nome}...")
    for dur in duracoes:
        for dp in dps:
            cmd = ["python", script_path, "--origem", origem_path, "--canal", "chan1", 
                   "--saida", "tmp_triagem_fina.csv", "--limiar_dp", str(dp), "--duracao_ms", str(dur)]
            
            try:
                res = subprocess.run(cmd, capture_output=True, text=True, check=True)
                out = res.stdout
                
                m = re.search(r"Ripples detectados:\s*(\d+)", out)
                ripples = int(m.group(1)) if m else -1
                
                m2 = re.search(r"Sinal carregado:\s*(\d+)", out)
                samples = int(m2.group(1)) if m2 else 1
                minutes = samples / 1000.0 / 60.0
                
                rate = ripples / minutes if minutes > 0 else 0
                
                resultados.append({
                    "origem": origem_nome,
                    "dp": dp,
                    "dur": dur,
                    "ripples": ripples,
                    "rate": rate
                })
                print(f"  {origem_nome} | DP={dp} | Dur={dur}ms | Ripples={ripples} ({rate:.2f}/min)")
            except subprocess.CalledProcessError as e:
                print(f"Erro ao rodar {origem_nome} DP={dp} Dur={dur}ms: {e}")

print("\n\n=== TABELA FINAL (SWEEP FINO) ===")
print("Origem          | DP  | Dur(ms) | Ripples | Taxa (/min)")
print("-" * 60)
for r in resultados:
    print(f"{r['origem']:<15} | {r['dp']:<3} | {r['dur']:<7} | {r['ripples']:<7} | {r['rate']:.2f}")
