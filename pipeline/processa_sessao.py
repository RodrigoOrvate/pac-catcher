import os
import argparse
import subprocess
import numpy as np
import pandas as pd
from ns2_utils import carrega_dados

def run(cmd, desc):
    print(f"\n[ORQUESTRADOR] >>> {desc} <<<")
    cmd_str = " ".join(cmd)
    print(f"Executando: {cmd_str}")
    res = subprocess.run(cmd)
    if res.returncode != 0:
        print(f"ERRO ao executar: {cmd_str}")
        import sys
        sys.exit(1)

def filtra_canais_vivos(dados, fs, limiar_rms_minimo=1e-6):
    """
    Retorna lista de índices de canais com amplitude RMS acima do limiar.
    Reusa o mesmo critério de detecção de canal morto de
    preprocessa_referencia_diferencial.seleciona_pool_referencia.
    """
    vivos = []
    for c in range(dados.shape[1]):
        rms = np.sqrt(np.mean(dados[:, c]**2))
        if rms >= limiar_rms_minimo:
            vivos.append(c)
    return vivos

def processa_sessao(pasta_sessao, saida_base, notch=(60, 120, 180, 240),
                    min_janelas_coocorrencia=3):
    nome_sessao = os.path.basename(pasta_sessao.rstrip("/\\"))
    saida_sessao = os.path.join(saida_base, nome_sessao)
    os.makedirs(saida_sessao, exist_ok=True)
    
    # Encontrar o primeiro arquivo ns2 na pasta para ler dados básicos
    arquivos_ns2 = sorted([f for f in os.listdir(pasta_sessao) if f.endswith(".ns2")])
    if not arquivos_ns2:
        print(f"[{nome_sessao}] Nenhum arquivo .ns2 encontrado.")
        return []
    primeiro_arq = os.path.join(pasta_sessao, arquivos_ns2[0])
    
    print(f"\n[{nome_sessao}] ESTÁGIO 0: Mapeamento de Canais Vivos")
    try:
        dados, fs, canal_ids = carrega_dados(primeiro_arq)
        # Limita leitura inicial apenas aos primeiros 30s se houver muitos dados para agilizar RMS.
        if dados.shape[0] > fs * 30:
            dados = dados[:int(fs*30), :]
    except Exception as e:
        print(f"[{nome_sessao}] Erro ao ler {primeiro_arq}: {e}")
        return []
        
    canais_vivos = filtra_canais_vivos(dados, fs)
    log_pulados = [c for c in range(dados.shape[1]) if c not in canais_vivos]
    if log_pulados:
        print(f"[{nome_sessao}] Canais descartados por amplitude morta: {log_pulados}")
    
    resumo_canais = []
    
    for c in canais_vivos:
        print(f"\n[{nome_sessao}] Iniciando Canal {c+1}...")
        saida_canal = os.path.join(saida_sessao, f"chan{c+1}")
        os.makedirs(saida_canal, exist_ok=True)
        
        # ESTÁGIO 1: Coocorrência
        n_teta_gama_total, n_teta_hg_total = 0, 0
        n_teta_hfo_total, n_teta_ripple_total = 0, 0
        for arq in arquivos_ns2:
            coocorrencia_csv = os.path.join(saida_canal, f"coocorrencia_{arq}.csv")
            run(["python", "pipeline/triagem_coocorrencia.py",
                 "--origem", os.path.join(pasta_sessao, arq), "--canal", str(c),
                 "--saida", coocorrencia_csv], f"Estágio 1 - {arq} - chan{c+1}")
                 
            if os.path.exists(coocorrencia_csv):
                df_coo = pd.read_csv(coocorrencia_csv)
                n_teta_gama_total += int((df_coo['teta_ok'] & df_coo['gamma_ok']).sum())
                n_teta_hg_total += int((df_coo['teta_ok'] & df_coo['hg_ok']).sum())
                n_teta_hfo_total += int((df_coo['teta_ok'] & df_coo['hfo_cru']).sum())
                n_teta_ripple_total += int((df_coo['teta_ok'] & df_coo['ripple']).sum())
            else:
                print(f"[{nome_sessao}] chan{c+1}: Falha ao gerar {coocorrencia_csv}")
                
        n_teta_gama = n_teta_gama_total
        n_teta_hg = n_teta_hg_total
        n_teta_hfo = n_teta_hfo_total
        n_teta_ripple = n_teta_ripple_total
        resumo_canais.append({"sessao": nome_sessao, "canal": c+1,
                              "n_teta_gama": n_teta_gama, "n_teta_hg": n_teta_hg,
                              "n_teta_hfo": n_teta_hfo, "n_teta_ripple": n_teta_ripple})
                              
        if n_teta_gama < min_janelas_coocorrencia and n_teta_hg < min_janelas_coocorrencia:
            print(f"[{nome_sessao}] chan{c+1}: coocorrencia insuficiente (Gama={n_teta_gama}, HG={n_teta_hg}), pulando pipeline pesado")
            continue
            
        # ESTÁGIO 2: Pipeline pesado
        triagem_csv = os.path.join(saida_canal, "triagem.csv")
        refinados_csv = os.path.join(saida_canal, "refinados.csv")
        
        run(["python", "pipeline/triagem_pac.py",
             "--pasta", pasta_sessao, "--canal", str(c), 
             "--saida", triagem_csv], f"Estágio 2.1 - Triagem PAC - {nome_sessao} chan{c+1}")
             
        run(["python", "pipeline/refina_candidatos.py",
             "--csv", triagem_csv, "--pasta_ns2", pasta_sessao, 
             "--saida", refinados_csv], f"Estágio 2.2 - Refina Candidatos - {nome_sessao} chan{c+1}")
             
        # Se refinados.csv está vazio ou tem poucas linhas úteis (só cabeçalho), pular.
        # Mas vamos rodar os filtros independentemente.
        
        notch_str = [str(n) for n in notch]
        run(["python", "pipeline/comodulogram.py",
             "--csv", refinados_csv, "--pasta_ns2", pasta_sessao,
             "--canal", str(c), "--notch", *notch_str,
             "--fdr_q", "0.05"], f"Estágio 2.3 - Comodulogram (Notch) - {nome_sessao} chan{c+1}")
             
        auditorias = [
            ("audita_harmonico.py", "harmonico.csv"),
            ("audita_harmonico_hfo.py", "harmonico_hfo.csv"),
            ("audita_skewness.py", "skewness.csv"),
            ("audita_footprint.py", "footprint.csv")
        ]
        
        for script, saida_nome in auditorias:
            run(["python", f"pipeline/auditorias/{script}",
                 "--csv", refinados_csv, "--pasta_ns2", pasta_sessao,
                 "--saida", os.path.join(saida_canal, saida_nome)],
                 f"Estágio 2.4 - {script} - {nome_sessao} chan{c+1}")
                 
    pd.DataFrame(resumo_canais).to_csv(
        os.path.join(saida_sessao, "resumo_canais.csv"), index=False)
    return resumo_canais

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Orquestrador em lote de sessões PAC.")
    parser.add_argument("--pasta", required=True, help="Pasta da sessão contendo arquivos .ns2")
    parser.add_argument("--saida", required=True, help="Pasta base de saída de resultados")
    parser.add_argument("--min_janelas", type=int, default=3, help="Mínimo de janelas promissoras no Estágio 1")
    args = parser.parse_args()
    
    processa_sessao(args.pasta, args.saida, min_janelas_coocorrencia=args.min_janelas)
