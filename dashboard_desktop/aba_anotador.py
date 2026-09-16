"""
dashboard_desktop/aba_anotador.py — seção "Anotador de Vídeo".

Abre pipeline/comportamento/anotador_comportamento.py (ferramenta madura já
existente, com vídeo + atalhos de teclado) como processo separado -- não
reimplementa nada dela. Só mostra progresso de anotação lendo o
CSV/template direto.
"""
import os
import subprocess
import sys
import tkinter as tk
from tkinter import filedialog, messagebox

import pandas as pd
import ttkbootstrap as ttk

from pac_core import workspace

ANOTADOR_SCRIPT = os.path.join(workspace.BASE_SCRIPT, "pipeline", "comportamento",
                                "anotador_comportamento.py")
TEMPLATE_COMPORTAMENTO_PADRAO = os.path.join(workspace.BASE_SCRIPT, "pipeline", "comportamento",
                                              "template_comportamento.csv")


class AnotadorMixin:
    def _build_tab_anotador(self):
        aba = self.tab_anotador
        self._callout(
            aba,
            "Abre pipeline/comportamento/anotador_comportamento.py -- a ferramenta de "
            "anotação já existente (vídeo + atalhos de teclado), sem reimplementar nada "
            "aqui. Roda como janela separada; feche-a e clique em \"Atualizar progresso\" "
            "pra ver o que mudou.",
            bootstyle="secondary",
        )

        f_csv = self._secao(aba, "Template de anotação")
        self.var_anotador_csv = tk.StringVar(value=TEMPLATE_COMPORTAMENTO_PADRAO)
        row = ttk.Frame(f_csv); row.pack(fill="x")
        ttk.Entry(row, textvariable=self.var_anotador_csv, width=75).pack(
            side="left", fill="x", expand=True)
        ttk.Button(row, text="Procurar...",
                  command=lambda: self.var_anotador_csv.set(
                      filedialog.askopenfilename() or self.var_anotador_csv.get())
                  ).pack(side="left", padx=4)

        f_acoes = ttk.Frame(f_csv); f_acoes.pack(fill="x", pady=(8, 0))
        ttk.Button(f_acoes, text="Abrir Anotador de Vídeo", bootstyle="primary",
                  command=self._anotador_abrir).pack(side="left")
        ttk.Button(f_acoes, text="Atualizar progresso",
                  command=self._anotador_atualiza_progresso).pack(side="left", padx=8)
        self.var_anotador_resumo = tk.StringVar(value="")
        ttk.Label(f_acoes, textvariable=self.var_anotador_resumo,
                 bootstyle="secondary").pack(side="left", padx=8)

        f_tabela = self._secao(aba, "Progresso por sessão")
        cols = ("sessao", "total", "anotado", "pct")
        self.tree_anotador = ttk.Treeview(f_tabela, columns=cols, show="headings", height=14)
        for c, titulo, largura in [("sessao", "Sessão", 480), ("total", "Total", 70),
                                   ("anotado", "Anotado", 80), ("pct", "%", 60)]:
            self.tree_anotador.heading(c, text=titulo)
            self.tree_anotador.column(c, width=largura, anchor="w" if c == "sessao" else "center")
        scroll = ttk.Scrollbar(f_tabela, orient="vertical", command=self.tree_anotador.yview)
        self.tree_anotador.configure(yscrollcommand=scroll.set)
        self.tree_anotador.pack(side="left", fill="both", expand=True)
        scroll.pack(side="left", fill="y")

        self._anotador_atualiza_progresso()

    def _anotador_abrir(self):
        if not os.path.exists(ANOTADOR_SCRIPT):
            messagebox.showerror("Erro", f"Não encontrado: {ANOTADOR_SCRIPT}")
            return
        csv_path = self.var_anotador_csv.get()
        args = [sys.executable, ANOTADOR_SCRIPT]
        if csv_path and os.path.exists(csv_path):
            args.append(csv_path)
        try:
            subprocess.Popen(args, cwd=workspace.BASE_SCRIPT)
        except Exception as e:
            messagebox.showerror("Erro ao abrir", str(e))

    def _anotador_atualiza_progresso(self):
        caminho = self.var_anotador_csv.get()
        for item in self.tree_anotador.get_children():
            self.tree_anotador.delete(item)
        if not os.path.exists(caminho):
            self.var_anotador_resumo.set("Template não encontrado.")
            return
        try:
            df = pd.read_csv(caminho, encoding="utf-8-sig")
        except Exception as e:
            self.var_anotador_resumo.set(f"Erro ao ler: {e}")
            return

        total_geral = len(df)
        anotado_geral = df["comportamento"].notna().sum() if "comportamento" in df.columns else 0
        pct_geral = 100 * anotado_geral / total_geral if total_geral else 0
        self.var_anotador_resumo.set(
            f"{anotado_geral}/{total_geral} janelas anotadas ({pct_geral:.0f}%)")

        if "sessao" not in df.columns:
            return
        # Agrupa por (sessao, condicao) quando a coluna existe -- mesmo conceito de
        # "sessão" que o anotador usa desde 2026-09-14 (uma sessão de infusão tem
        # basal+0h+1h+2h com vídeos e progresso distintos, não um total só).
        chaves = ["sessao", "condicao"] if "condicao" in df.columns else ["sessao"]
        resumo = df.groupby(chaves).agg(
            total=("sessao", "size"),
            anotado=("comportamento", lambda s: s.notna().sum()),
        ).reset_index().sort_values(chaves)
        for _, r in resumo.iterrows():
            pct = 100 * r["anotado"] / r["total"] if r["total"] else 0
            rotulo = f"{r['sessao']} :: {r['condicao']}" if "condicao" in chaves else r["sessao"]
            self.tree_anotador.insert("", "end", values=(rotulo, int(r["total"]),
                                                          int(r["anotado"]), f"{pct:.0f}%"))


