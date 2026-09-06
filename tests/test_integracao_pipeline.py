import os
import subprocess
import time

def run_cmd(cmd, desc):
    print(f"\n[{time.strftime('%H:%M:%S')}] Iniciando: {desc}")
    print(f"Comando: {' '.join(cmd)}")
    t0 = time.time()
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERRO! Comando falhou com codigo {e.returncode}")
        sys.exit(1)
    t1 = time.time()
    print(f"[{time.strftime('%H:%M:%S')}] Concluido em {t1-t0:.1f}s")

if __name__ == "__main__":
    import sys
    base_dir = r"C:\acoplamento_theta-gamma\SCRIPT"
    pasta = r"C:\acoplamento_theta-gamma\MTESC04_NOCI\MTESC04 -- 1 - infusao - 08-07-2024\Basal antes da infusao"
    canal = "0"
    
    # 1. Triagem (teta_ok, janelas)
    run_cmd(["python", "pipeline/triagem_pac.py", "--pasta", pasta, "--canal", canal, "--saida", "teste_triagem.csv"], "Triagem PAC")
    
    # 2. Refina (FOOOF)
    run_cmd(["python", "pipeline/refina_candidatos.py", "--csv", "teste_triagem.csv", "--saida", "teste_refinados.csv", "--pasta_ns2", pasta], "Refinamento FOOOF")
    
    # 3. Comodulogram (Tort MI)
    run_cmd(["python", "pipeline/comodulogram.py", "--csv", "teste_refinados.csv", "--pasta_ns2", pasta, "--canal", canal], "Comodulograma")
    
    # 4. Auditoria Harmônico
    run_cmd(["python", "pipeline/auditorias/audita_harmonico.py"], "Auditoria Harmonico (HG)")
    
    # 5. Auditoria Harmônico HFO
    run_cmd(["python", "pipeline/auditorias/audita_harmonico_hfo.py"], "Auditoria Harmonico (HFO)")
    
    print("\nSUCESSO! Teste de integração ponta a ponta concluído sem erros.")
