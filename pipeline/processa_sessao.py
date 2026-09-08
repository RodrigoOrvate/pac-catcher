import os
import argparse
import subprocess
import numpy as np
import pandas as pd
from ns2_utils import carrega_dados

# Os 3 pares reconhecidos por refina_candidatos.py/comodulogram.py.
# processa_sessao.py precisa iterar sobre TODOS eles: comodulogram.py filtra
# internamente por --par (default "theta_gamma"), entao chama-lo uma unica vez
# sem essa flag descarta silenciosamente 100% dos candidatos de theta_hg e
# theta_hfo do Estagio 2.3 (nenhum comodulograma, nenhuma linha em
# resumo_comodulogramas.csv para esses dois pares).
PARES = ["theta_gamma", "theta_hg", "theta_hfo"]

def run(cmd, desc):
    print(f"\n[ORQUESTRADOR] >>> {desc} <<<")
    cmd_str = " ".join(cmd)
    print(f"Executando: {cmd_str}")
    res = subprocess.run(cmd)
    if res.returncode != 0:
        print(f"ERRO ao executar: {cmd_str}")
        import sys
        sys.exit(1)

def filtra_canais_vivos(dados, fs, limiar_rms_minimo=1e-6, limiar_rms_outlier=2.0):
    """
    Retorna lista de índices de canais com amplitude RMS acima do limiar mínimo,
    e uma lista de índices de canais suspeitos de serem outliers de ruído/referência.
    """
    rmss = []
    for c in range(dados.shape[1]):
        rms = np.sqrt(np.mean(dados[:, c]**2))
        rmss.append(rms)
        
    mediana_rms = np.median(rmss)
    vivos = []
    outliers = []
    
    for c, rms in enumerate(rmss):
        if rms < limiar_rms_minimo:
            continue
        if rms > limiar_rms_outlier * mediana_rms:
            outliers.append(c)
        vivos.append(c)
        
    return vivos, outliers

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
        
    canais_vivos, canais_outliers = filtra_canais_vivos(dados, fs)
    log_pulados = [c for c in range(dados.shape[1]) if c not in canais_vivos]
    if log_pulados:
        print(f"[{nome_sessao}] Canais descartados por amplitude morta: {log_pulados}")
    if canais_outliers:
        print(f"[{nome_sessao}] [ALERTA] Canais marcados como OUTLIER DE RMS (ruído/ref suspeita): {canais_outliers}")
    
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
        rms_outlier = (c in canais_outliers)
        resumo_canais.append({"sessao": nome_sessao, "canal": c+1,
                              "n_teta_gama": n_teta_gama, "n_teta_hg": n_teta_hg,
                              "n_teta_hfo": n_teta_hfo, "n_teta_ripple": n_teta_ripple,
                              "rms_outlier": rms_outlier})
                              
        if n_teta_gama < min_janelas_coocorrencia and n_teta_hg < min_janelas_coocorrencia:
            print(f"[{nome_sessao}] chan{c+1}: coocorrencia insuficiente (Gama={n_teta_gama}, HG={n_teta_hg}), pulando pipeline pesado")
            continue
            
        # ESTÁGIO 2: Pipeline pesado
        triagem_csv = os.path.join(saida_canal, "triagem.csv")
        refinados_csv = os.path.join(saida_canal, "refinados.csv")
        
        run(["python", "pipeline/triagem_pac.py",
             "--pasta", pasta_sessao, "--canal", str(c), 
             "--pares", *PARES,
             "--saida", triagem_csv], f"Estágio 2.1 - Triagem PAC - {nome_sessao} chan{c+1}")
             
        run(["python", "pipeline/refina_candidatos.py",
             "--csv", triagem_csv, "--pasta_ns2", pasta_sessao, 
             "--saida", refinados_csv], f"Estágio 2.2 - Refina Candidatos - {nome_sessao} chan{c+1}")
             
        # Se refinados.csv está vazio ou tem poucas linhas úteis (só cabeçalho), pular.
        # Mas vamos rodar os filtros independentemente.
        
        notch_str = [str(n) for n in notch]

        # ESTÁGIO 2.3: Comodulogram — um run POR PAR, senão theta_hg e theta_hfo
        # nunca chegam a ser processados (ver comentário em PARES acima).
        resumos_comod = []
        for par in PARES:
            saida_dir_par = os.path.join(saida_canal, "comodulogramas", par)
            run(["python", "pipeline/comodulogram.py",
                 "--csv", refinados_csv, "--pasta_ns2", pasta_sessao,
                 "--canal", str(c), "--notch", *notch_str,
                 "--par", par, "--saida_dir", saida_dir_par,
                 "--fdr_q", "0.05"],
                f"Estágio 2.3 - Comodulogram ({par}) - {nome_sessao} chan{c+1}")

            resumo_par_csv = os.path.join(saida_dir_par, "resumo_comodulogramas.csv")
            if os.path.exists(resumo_par_csv):
                df_resumo_par = pd.read_csv(resumo_par_csv)
                if len(df_resumo_par) > 0:
                    resumos_comod.append(df_resumo_par)

        resumo_comod_canal_csv = os.path.join(saida_canal, "resumo_comodulogramas.csv")
        if resumos_comod:
            pd.concat(resumos_comod, ignore_index=True).to_csv(resumo_comod_canal_csv, index=False)
        else:
            print(f"[{nome_sessao}] chan{c+1}: nenhum comodulograma gerado em nenhum par — "
                  f"pulando auditorias harmônicas (dependem de fase_pico_hz/amp_pico_hz).")

        # ESTÁGIO 2.4: Auditorias
        # audita_skewness.py depende só da janela (nao do par) — roda 1x sobre
        # refinados.csv deduplicado por janela para nao triplicar linhas
        # (refinados.csv tem 1 linha por par; a mesma janela aparece ate 3x).
        df_refinados = pd.read_csv(refinados_csv) if os.path.exists(refinados_csv) else pd.DataFrame()
        skew_input_csv = refinados_csv
        if len(df_refinados) > 0:
            chave_janela = ["arquivo", "canal", "janela_ini_s", "janela_fim_s"]
            df_dedup = df_refinados.drop_duplicates(subset=chave_janela)
            if len(df_dedup) < len(df_refinados):
                skew_input_csv = os.path.join(saida_canal, "_refinados_dedup_janela.csv")
                df_dedup.to_csv(skew_input_csv, index=False)

        run(["python", "pipeline/auditorias/audita_skewness.py",
             "--csv", skew_input_csv, "--pasta_ns2", pasta_sessao,
             "--saida", os.path.join(saida_canal, "skewness.csv")],
            f"Estágio 2.4 - audita_skewness.py - {nome_sessao} chan{c+1}")

        # audita_harmonico.py (teta -> par ativo) precisa de fase_pico_hz/amp_pico_hz,
        # que só existem em resumo_comodulogramas.csv (não em refinados.csv).
        # Roda 1x por par que de fato gerou comodulogramas.
        harmonicos = []
        for df_resumo_par in resumos_comod:
            par_atual = df_resumo_par["par"].iloc[0]
            resumo_par_path = os.path.join(saida_canal, f"_resumo_{par_atual}.csv")
            df_resumo_par.to_csv(resumo_par_path, index=False)
            harmonico_par_csv = os.path.join(saida_canal, f"harmonico_{par_atual}.csv")
            run(["python", "pipeline/auditorias/audita_harmonico.py",
                 "--csv", resumo_par_path, "--pasta_ns2", pasta_sessao,
                 "--saida", harmonico_par_csv],
                f"Estágio 2.4 - audita_harmonico.py ({par_atual}) - {nome_sessao} chan{c+1}")
            if os.path.exists(harmonico_par_csv):
                df_h = pd.read_csv(harmonico_par_csv)
                if len(df_h) > 0:
                    harmonicos.append(df_h)

            # audita_harmonico_hfo.py (gama/teta -> HFO) só faz sentido para o par teta-HFO
            if par_atual == "theta_hfo":
                run(["python", "pipeline/auditorias/audita_harmonico_hfo.py",
                     "--csv", resumo_par_path, "--pasta_ns2", pasta_sessao,
                     "--saida", os.path.join(saida_canal, "harmonico_hfo.csv")],
                    f"Estágio 2.4 - audita_harmonico_hfo.py - {nome_sessao} chan{c+1}")

        if harmonicos:
            pd.concat(harmonicos, ignore_index=True).to_csv(
                os.path.join(saida_canal, "harmonico.csv"), index=False)
                 
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
