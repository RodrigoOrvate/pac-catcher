"""
dashboard_desktop/app.py — ThetaGamma-Studio (versão desktop).

Ponto de entrada: monta as 3 seções de topo e aplica o tema ttkbootstrap.

    Pipeline           -- roda as etapas de pipeline/ com log ao vivo (aba_pipeline.py)
    Análise de Dados   -- visualização dos resultados já processados (aba_analise.py)
    Anotador de Vídeo  -- abre o anotador de comportamento existente (aba_anotador.py)

Rodar com: python dashboard_desktop/app.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
matplotlib.use("TkAgg")
import ttkbootstrap as ttk
from tkinter import messagebox

from pac_core import workspace
from dashboard_desktop.estilo import TEMA_PADRAO, EstiloMixin
from dashboard_desktop.aba_pipeline import PipelineMixin
from dashboard_desktop.aba_analise import AnaliseMixin
from dashboard_desktop.aba_anotador import AnotadorMixin


class ThetaGammaStudioApp(EstiloMixin, PipelineMixin, AnaliseMixin, AnotadorMixin, ttk.Window):
    def __init__(self):
        super().__init__(themename=TEMA_PADRAO)
        self.title("ThetaGamma-Studio")
        # Tamanho a partir da tela real (o notebook de uso é 1366x768 -- um
        # 1350x900 fixo deixava o log escondido atrás da barra de tarefas).
        larg = min(1350, self.winfo_screenwidth() - 40)
        alt = min(900, self.winfo_screenheight() - 80)
        self.geometry(f"{larg}x{alt}+10+10")
        self.minsize(960, 600)
        if os.name == "nt":
            self.state("zoomed")
        self._aplica_estilo()

        header = ttk.Frame(self, padding=(16, 8, 16, 6))
        header.pack(fill="x")
        ttk.Label(header, text="ThetaGamma-Studio", style="Titulo.TLabel").pack(anchor="w")
        ttk.Label(header, text=workspace.BASE_WORKSPACE, style="Subtitulo.TLabel").pack(anchor="w")
        ttk.Separator(self, orient="horizontal").pack(fill="x")

        nb = ttk.Notebook(self)
        nb.pack(fill="both", expand=True, padx=10, pady=10)
        self.tab_pipeline = ttk.Frame(nb, padding=6); nb.add(self.tab_pipeline, text="Pipeline")
        self.tab_analise = ttk.Frame(nb, padding=6); nb.add(self.tab_analise, text="Análise de Dados")
        self.tab_anotador = ttk.Frame(nb, padding=6); nb.add(self.tab_anotador, text="Anotador de Vídeo")

        self._build_secao_pipeline(self.tab_pipeline)
        self._build_secao_analise(self.tab_analise)
        self._build_tab_anotador()

        self.protocol("WM_DELETE_WINDOW", self._ao_fechar)

    def _ao_fechar(self):
        """Fechar a janela com uma etapa rodando deixaria o processo órfão
        (foi um processo órfão que travou a pasta em C: no início da
        migração) -- pergunta e derruba a árvore antes de sair."""
        rodando = [p for p in getattr(self, "_paineis", []) if p.executando()]
        if rodando:
            if not messagebox.askyesno(
                    "Processos em execução",
                    f"{len(rodando)} etapa(s) do pipeline ainda rodando. Encerrar e sair?"):
                return
            for painel in rodando:
                painel.parar()
        self.destroy()


if __name__ == "__main__":
    app = ThetaGammaStudioApp()
    app.mainloop()
