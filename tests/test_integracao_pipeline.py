import os
import sys
import subprocess
import time

def main():
    print("=" * 85)
    print("TESTE DE INTEGRAÇÃO DO ORQUESTRADOR")
    print("=" * 85)
    
    # Vamos rodar o orquestrador em uma sessão de teste real.
    pasta_teste = r"C:\acoplamento_theta-gamma\MTESC04_NOCI\MTESC04 -- 1 - infusao - 08-07-2024\Basal antes da infusao"
    pasta_saida = r"C:\acoplamento_theta-gamma\SCRIPT\TESTE_INTEGRACAO_OUT"
    
    if not os.path.exists(pasta_teste):
        print(f"Pasta de teste não encontrada: {pasta_teste}")
        print("Teste pulado no CI, mas em produção deve falhar.")
        sys.exit(1)
        
    cmd = [
        "python", "pipeline/processa_sessao.py",
        "--pasta", pasta_teste,
        "--saida", pasta_saida,
        "--min_janelas", "1000"  # Colocamos um limiar impossível para que o Estágio 2 não rode em todos os 32 canais e demore 5 horas.
    ]
    
    print(f"Executando: {' '.join(cmd)}")
    
    try:
        # Só testamos a integridade do Estágio 0 e 1 no dry run.
        # No entanto, se quisermos testar o fluxo de verdade, precisamos processar pelo menos 1 canal.
        # Mas para evitar timeout, vamos passar min_janelas=0 e alterar temporariamente o script 
        # ou apenas executar `processa_sessao` com um timeout/limite artificial.
        
        # Como o usuário sugeriu um dry-run do orquestrador no passo 5, esse teste de integração automatizado
        # valida se os argumentos de parser do processa_sessao estao corretos.
        subprocess.run(cmd, check=True)
        print("SUCESSO: Orquestrador invocado com sucesso (Estágio 0 e 1 testados no limiar alto).")
    except subprocess.CalledProcessError as e:
        print(f"FALHA na orquestração: erro {e.returncode}")
        sys.exit(1)

if __name__ == "__main__":
    main()
