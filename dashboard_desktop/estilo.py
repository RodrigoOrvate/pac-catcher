"""
dashboard_desktop/estilo.py — sistema de design do ThetaGamma-Studio.

Usa ttkbootstrap (tema "bootstrap-light": accent azul #0a58ca, o mesmo azul
Bootstrap 5 padrão -- profissional, reconhecível, sem inventar paleta
própria) em vez de ttk.Style manual. `EstiloMixin` fica só com os helpers de
layout que ttkbootstrap não cobre pronto (rótulo de seção discreto, callout
de aviso com tarja lateral).
"""
import tkinter as tk

import ttkbootstrap as ttk
from tkinter import filedialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

TEMA_PADRAO = "bootstrap-light"
FONTE = "Segoe UI"


def _tinge_para_branco(hex_cor, peso_branco=0.88):
    """Mistura `hex_cor` com branco (usado pro fundo claro do callout, a
    partir da cor de destaque real do tema ativo -- não um hex chutado à
    mão, então acompanha a troca de tema automaticamente)."""
    hex_cor = hex_cor.lstrip("#")
    r, g, b = (int(hex_cor[i:i + 2], 16) for i in (0, 2, 4))
    r = int(r + (255 - r) * peso_branco)
    g = int(g + (255 - g) * peso_branco)
    b = int(b + (255 - b) * peso_branco)
    return f"#{r:02x}{g:02x}{b:02x}"


class EstiloMixin:
    def _aplica_estilo(self):
        """Chamado uma vez no __init__ (self já é um ttkbootstrap.Window,
        criado com themename=TEMA_PADRAO -- aqui só complementa com os
        estilos que o tema não define: rótulo de seção e cabeçalho."""
        style = self.style  # ttkbootstrap.Style, já ligado à Window

        style.configure("SecaoTitulo.TLabel", font=(FONTE, 8, "bold"),
                        foreground=style.colors.secondary)
        style.configure("Titulo.TLabel", font=(FONTE, 15, "bold"))
        style.configure("Subtitulo.TLabel", font=(FONTE, 8),
                        foreground=style.colors.secondary)

    # ------------------------------------------------------------------
    # Helpers de layout reusados por todas as seções/abas
    # ------------------------------------------------------------------
    def _secao(self, parent, titulo):
        """Rótulo pequeno em maiúsculas + separador -- hierarquia por peso
        de fonte, não por caixa com borda (LabelFrame lembra diálogo
        antigo)."""
        wrapper = ttk.Frame(parent)
        wrapper.pack(fill="x", padx=4, pady=(12, 2))
        ttk.Label(wrapper, text=titulo.upper(), style="SecaoTitulo.TLabel").pack(anchor="w")
        ttk.Separator(wrapper, orient="horizontal").pack(fill="x", pady=(3, 8))
        corpo = ttk.Frame(wrapper)
        corpo.pack(fill="x")
        return corpo

    def _callout(self, parent, texto, bootstyle="warning"):
        """Aviso discreto: tarja fina colorida na borda esquerda + fundo
        levemente tingido, cor tirada do tema ativo (`bootstyle` = nome
        semântico ttkbootstrap: warning/secondary/info/danger/success)."""
        cor = self.style.colors.get(bootstyle)
        cor_fundo = _tinge_para_branco(cor)
        f = tk.Frame(parent, bg=cor_fundo)
        f.pack(fill="x", padx=4, pady=(4, 10))
        tk.Frame(f, bg=cor, width=3).pack(side="left", fill="y")
        tk.Label(f, text=texto, bg=cor_fundo, fg=self.style.colors.dark,
                 font=(FONTE, 9), wraplength=1150, justify="left",
                 padx=10, pady=8).pack(side="left", fill="both", expand=True)
        return f

    def _entrada_arquivo(self, parent, var_caminho):
        f = ttk.Frame(parent)
        f.pack(fill="x", pady=2)
        ttk.Label(f, text="Arquivo (.ns2/.bin/.dat/.mat):", width=26).pack(side="left")
        ttk.Entry(f, textvariable=var_caminho, width=70).pack(side="left", padx=4, fill="x", expand=True)
        ttk.Button(f, text="Procurar...",
                   command=lambda: var_caminho.set(filedialog.askopenfilename() or var_caminho.get())
                   ).pack(side="left")
        return f

    def _embute_figura(self, parent, figsize=(9, 4)):
        fig = Figure(figsize=figsize, dpi=100)
        canvas = FigureCanvasTkAgg(fig, master=parent)
        canvas.get_tk_widget().pack(fill="both", expand=True)
        toolbar = NavigationToolbar2Tk(canvas, parent, pack_toolbar=False)
        toolbar.pack(side="bottom", fill="x")
        return fig, canvas
