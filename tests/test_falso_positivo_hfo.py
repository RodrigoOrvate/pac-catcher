import numpy as np
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "pipeline", "auditorias"))

from utils_harmonico import testa_razao_harmonica, calcula_n_max

def run_monte_carlo_false_positives(n_iter=100000, f_teta=8.0, hfo_min=150.0, hfo_max=250.0, tol_rel=0.15):
    """
    Testa a probabilidade de um oscilador aleatorio na banda HFO cair na tolerancia
    harmonica de um f_teta, usando a regra dinamica de n_max.
    """
    rng = np.random.default_rng(42)
    
    # Sorteia frequencias uniformemente na banda HFO
    freqs_hfo = rng.uniform(hfo_min, hfo_max, n_iter)
    
    suspeitos = 0
    tol_abs_teta = tol_rel * f_teta
    
    for f_cand in freqs_hfo:
        # A nova arquitetura calcula n_max dinamicamente para O CANDIDATO
        n_max_dinamico = calcula_n_max(f_teta, f_cand, tol_abs_teta)
        
        # Testa se cai na razao harmonica (a auditoria usa tol_rel=0.15 por default)
        susp, _, _, _, _ = testa_razao_harmonica(f_teta, f_cand, tol_rel, n_max=n_max_dinamico)
        if susp:
            suspeitos += 1
            
    taxa = (suspeitos / n_iter) * 100
    print(f"Taxa de Falso Positivo Puro (n={n_iter}): {taxa:.2f}%")
    print(f"Banda testada: {hfo_min}-{hfo_max}Hz | Teta: {f_teta}Hz | tol_rel: {tol_rel}")

if __name__ == "__main__":
    print("Simulando tol_rel = 0.15 (default do audita_harmonico_hfo.py)")
    run_monte_carlo_false_positives(tol_rel=0.15)
    print("\nSimulando tol_rel = 0.10 (teste sintético clássico)")
    run_monte_carlo_false_positives(tol_rel=0.10)
