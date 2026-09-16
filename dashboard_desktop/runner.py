"""
dashboard_desktop/runner.py — execução de subprocessos com log ao vivo.

Não existia nenhum padrão desse tipo no projeto antes (confirmado por
exploração 2026-09-14: o único `subprocess.Popen` do repo, em
`aba_anotador.py`, é fire-and-forget, sem capturar stdout). `PainelExecucao`
é o componente novo que a seção Pipeline usa pra rodar
triagem/refinamento/comodulograma/auditorias/agregação de verdade, com log
linha a linha em vez de só um resultado final.

`_rodar_em_thread` é o helper mais simples (função Python in-process, sem
subprocess) já usado pelas abas de Análise de Dados (comodulograma pontual,
replay) -- mora aqui por ser a mesma família de "não travar a UI".
"""
import os
import queue
import subprocess
import threading
import tkinter as tk

import ttkbootstrap as ttk


def _rodar_em_thread(job, ao_terminar, ao_erro):
    """Roda `job` (sem args, retorna o resultado) numa thread separada pra
    não travar a janela do Tkinter, chamando `ao_terminar(resultado)` ou
    `ao_erro(excecao)` de volta na thread principal via widget.after()."""
    def alvo():
        try:
            resultado = job()
        except Exception as e:  # noqa: BLE001 -- repassa qualquer erro pra UI
            ao_erro(e)
            return
        ao_terminar(resultado)
    threading.Thread(target=alvo, daemon=True).start()


class PainelExecucao(ttk.Frame):
    """Log ao vivo + botão Parar pra um subprocess. Quem decide QUANDO e
    COM QUAIS argumentos rodar é a aba que usa o painel (cada uma tem seu
    próprio formulário) -- este widget só sabe rodar um `cmd` e mostrar a
    saída, não constrói comando nenhum sozinho.

    Uso:
        painel = PainelExecucao(parent)
        painel.pack(fill="both", expand=True)
        ...
        painel.executar([sys.executable, "script.py", "--flag", valor],
                         cwd=workspace.BASE_SCRIPT, ao_concluir=minha_funcao)
    """

    def __init__(self, parent, altura_log=8):
        # altura_log é a altura MÍNIMA pedida (linhas); o painel é empacotado
        # com expand=True e cresce pro espaço que sobrar na janela.
        super().__init__(parent)
        self._fila = queue.Queue()
        self._proc = None
        self._rodando = False
        self._ao_concluir = None

        barra = ttk.Frame(self)
        barra.pack(fill="x")
        self.var_status = tk.StringVar(value="Ocioso.")
        ttk.Label(barra, textvariable=self.var_status).pack(side="left")
        self.btn_parar = ttk.Button(barra, text="Parar", bootstyle="danger-outline",
                                    command=self._parar, state="disabled")
        self.btn_parar.pack(side="right")

        corpo = ttk.Frame(self)
        corpo.pack(fill="both", expand=True, pady=(6, 0))
        self.txt_log = tk.Text(corpo, height=altura_log, wrap="none", font=("Consolas", 9),
                               bg="#1E1E1E", fg="#D4D4D4", insertbackground="#D4D4D4")
        self.txt_log.pack(side="left", fill="both", expand=True)
        scroll = ttk.Scrollbar(corpo, orient="vertical", command=self.txt_log.yview)
        scroll.pack(side="left", fill="y")
        self.txt_log.configure(yscrollcommand=scroll.set)

    def executando(self):
        return self._rodando

    def executar(self, cmd, cwd=None, ao_concluir=None):
        """Roda `cmd` (lista de argumentos) via subprocess, stdout+stderr
        combinados streamados linha a linha pro log. `ao_concluir(exit_code)`
        (opcional) roda na thread principal quando o processo termina."""
        if self._rodando:
            return
        self._rodando = True
        self._ao_concluir = ao_concluir
        self.txt_log.delete("1.0", "end")
        self.txt_log.insert("end", f"$ {' '.join(str(c) for c in cmd)}\n\n")
        self.var_status.set("Rodando...")
        self.btn_parar.configure(state="normal")

        # Com stdout em pipe, o Python filho faz buffer em blocos (o log
        # chegaria aos trancos) e no Windows usa cp1252 (print de θ/×/→ pode
        # dar UnicodeEncodeError). As duas variáveis são herdadas também pelos
        # netos -- processa_sessao.py dispara cada etapa como subprocess.
        env = dict(os.environ, PYTHONUNBUFFERED="1", PYTHONIOENCODING="utf-8")

        def alvo():
            try:
                proc = subprocess.Popen(
                    cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    text=True, bufsize=1, encoding="utf-8", errors="replace", env=env,
                )
                self._proc = proc
                for linha in proc.stdout:
                    self._fila.put(("linha", linha))
                proc.wait()
                self._fila.put(("fim", proc.returncode))
            except Exception as e:
                self._fila.put(("erro", str(e)))

        threading.Thread(target=alvo, daemon=True).start()
        self.after(80, self._drena_fila)

    def _drena_fila(self):
        try:
            while True:
                tipo, valor = self._fila.get_nowait()
                if tipo == "linha":
                    self.txt_log.insert("end", valor)
                    self.txt_log.see("end")
                elif tipo == "fim":
                    self._finaliza(valor)
                    return
                elif tipo == "erro":
                    self.txt_log.insert("end", f"\n[ERRO ao iniciar processo] {valor}\n")
                    self._finaliza(-1)
                    return
        except queue.Empty:
            pass
        if self._rodando:
            self.after(80, self._drena_fila)

    def _finaliza(self, exit_code):
        self._rodando = False
        self._proc = None
        ok = exit_code == 0
        self.var_status.set(f"Concluído (exit={exit_code})." if ok else f"Falhou (exit={exit_code}).")
        self.btn_parar.configure(state="disabled")
        if self._ao_concluir:
            self._ao_concluir(exit_code)

    def parar(self):
        """Encerra o processo E seus filhos. No Windows, terminate() só mata o
        processo direto -- processa_sessao.py dispara cada etapa como
        subprocess, então o neto (ex. triagem_pac.py) ficaria órfão rodando
        (e segurando arquivos abertos). taskkill /T derruba a árvore toda."""
        proc = self._proc
        if proc is None:
            return
        if os.name == "nt":
            subprocess.run(["taskkill", "/F", "/T", "/PID", str(proc.pid)],
                           capture_output=True)
        else:
            proc.terminate()
        self.var_status.set("Interrompido pelo usuário.")

    _parar = parar
