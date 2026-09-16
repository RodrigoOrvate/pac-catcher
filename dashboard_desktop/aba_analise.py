"""
dashboard_desktop/aba_analise.py — seção "Análise de Dados" do ThetaGamma-Studio.

Consolida as 6 visualizações que antes eram abas de topo (Inspeção de
Canais, Comodulograma, Portões de Qualidade, Galeria de Vencedores,
Comportamento, Closed-Loop/Replay) sob uma única seção com sub-abas. Lógica
idêntica à versão anterior -- só a montagem de abas mudou (ver plano de
reestruturação 2026-09-14). Nenhuma métrica científica é recalculada com
lógica nova aqui.
"""
import glob
import os
import tkinter as tk
from tkinter import filedialog, messagebox

import numpy as np
import pandas as pd
import ttkbootstrap as ttk
from PIL import Image, ImageTk
from scipy.signal import detrend
from scipy.stats import chi2_contingency

from pac_core import workspace
from pac_core.filtering import aplica_notch, filtra_sinal
from pac_core.io import carrega_dados, fatia_janela, resolve_canal_idx
from pac_studio import BAND_PAIRS, analisa_janela, transicao_estado
from preditor.prever_pac_tempo_real import (
    CATEGORIAS_COMPORTAMENTO, JANELA_S, PASSO_S,
    carrega_modelo, extrair_features_estado, resolve_canal,
)
from pipeline.etapa3_comodulograma.comodulogram import (
    AMPS_DEFAULT, FASES_DEFAULT, calcula_comodulograma_z,
)

from dashboard_desktop.runner import _rodar_em_thread

VENCEDORES_PADRAO = os.path.join(workspace.BASE_SCRIPT, "resultados",
                                  "candidatos_vencedores_OURO_PURIFICADO_v2.csv")
DATASET_MESTRE_PADRAO = os.path.join(workspace.BASE_SCRIPT, "resultados",
                                      "dataset_mestre_COM_COMPORTAMENTO.csv")


class AnaliseMixin:
    def _build_secao_analise(self, parent):
        """Monta a seção 'Análise de Dados': um sub-notebook com as 6
        visualizações que antes eram abas de topo."""
        nb = ttk.Notebook(parent)
        nb.pack(fill="both", expand=True)

        self.tab_canais = ttk.Frame(nb, padding=4); nb.add(self.tab_canais, text="Inspeção de Canais")
        self.tab_comod = ttk.Frame(nb, padding=4); nb.add(self.tab_comod, text="Comodulograma")
        self.tab_portoes = ttk.Frame(nb, padding=4); nb.add(self.tab_portoes, text="Portões de Qualidade")
        self.tab_galeria = ttk.Frame(nb, padding=4); nb.add(self.tab_galeria, text="Galeria de Vencedores")
        self.tab_comport = ttk.Frame(nb, padding=4); nb.add(self.tab_comport, text="Comportamento")
        self.tab_closedloop = ttk.Frame(nb, padding=4); nb.add(self.tab_closedloop, text="Closed-Loop / Replay")

        self._imagens_vivas = []  # evita coleta de lixo dos PhotoImage da galeria
        self._indice_png = None

        self._build_tab_canais()
        self._build_tab_comod()
        self._build_tab_portoes()
        self._build_tab_galeria()
        self._build_tab_comport()
        self._build_tab_closedloop()

    def _build_tab_canais(self):
        aba = self.tab_canais
        self.var_canais_caminho = tk.StringVar()
        self._entrada_arquivo(self._secao(aba, "Arquivo"), self.var_canais_caminho)

        f_ctrl = self._secao(aba, "Controles")
        row1 = ttk.Frame(f_ctrl); row1.pack(fill="x", pady=2)
        ttk.Button(row1, text="Carregar canais", command=self._canais_carregar, bootstyle="primary").pack(side="left")
        self.var_canais_info = tk.StringVar(value="(nenhum arquivo carregado)")
        ttk.Label(row1, textvariable=self.var_canais_info).pack(side="left", padx=10)

        row2 = ttk.Frame(f_ctrl); row2.pack(fill="x", pady=2)
        ttk.Label(row2, text="Início (s):").pack(side="left")
        self.var_canais_ini = tk.DoubleVar(value=0.0)
        ttk.Entry(row2, textvariable=self.var_canais_ini, width=8).pack(side="left", padx=(2, 10))
        ttk.Label(row2, text="Fim (s):").pack(side="left")
        self.var_canais_fim = tk.DoubleVar(value=10.0)
        ttk.Entry(row2, textvariable=self.var_canais_fim, width=8).pack(side="left", padx=(2, 10))
        ttk.Label(row2, text="Ganho:").pack(side="left")
        self.var_canais_ganho = tk.DoubleVar(value=1.0)
        ttk.Scale(row2, from_=0.1, to=10.0, variable=self.var_canais_ganho,
                  orient="horizontal", length=120).pack(side="left", padx=4)

        row3 = ttk.Frame(f_ctrl); row3.pack(fill="x", pady=2)
        ttk.Label(row3, text="Canal de referência:").pack(side="left")
        self.var_canais_ref = tk.StringVar(value="(nenhum)")
        self.cb_canais_ref = ttk.Combobox(row3, textvariable=self.var_canais_ref,
                                          state="readonly", width=12)
        self.cb_canais_ref.pack(side="left", padx=4)
        self.var_canais_notch = tk.BooleanVar(value=True)
        ttk.Checkbutton(row3, text="Notch 60/120/180/240Hz",
                        variable=self.var_canais_notch).pack(side="left", padx=10)
        self.var_canais_detrend = tk.BooleanVar(value=False)
        ttk.Checkbutton(row3, text="Remover deriva",
                        variable=self.var_canais_detrend).pack(side="left", padx=10)
        self.var_canais_passafaixa = tk.BooleanVar(value=False)
        ttk.Checkbutton(row3, text="Passa-faixa", variable=self.var_canais_passafaixa).pack(side="left")
        ttk.Label(row3, text="lo:").pack(side="left", padx=(6, 0))
        self.var_canais_lo = tk.DoubleVar(value=1.0)
        ttk.Entry(row3, textvariable=self.var_canais_lo, width=6).pack(side="left")
        ttk.Label(row3, text="hi:").pack(side="left")
        self.var_canais_hi = tk.DoubleVar(value=100.0)
        ttk.Entry(row3, textvariable=self.var_canais_hi, width=6).pack(side="left")

        f_lista = ttk.Frame(f_ctrl); f_lista.pack(fill="x", pady=4)
        ttk.Label(f_lista, text="Canais (Ctrl/Shift p/ múltiplos):").pack(anchor="w")
        self.lst_canais = tk.Listbox(f_lista, selectmode="extended", height=5, exportselection=False)
        self.lst_canais.pack(fill="x")

        ttk.Button(f_ctrl, text="Plotar", command=self._canais_plotar, bootstyle="primary").pack(pady=4)

        f_plot = ttk.Frame(aba); f_plot.pack(fill="both", expand=True, padx=8, pady=4)
        self.fig_canais, self.canvas_canais = self._embute_figura(f_plot)
        self._canal_ids_canais = None
        self._dados_canais = None
        self._fs_canais = None

    def _canais_carregar(self):
        caminho = self.var_canais_caminho.get()
        if not caminho or not os.path.exists(caminho):
            messagebox.showerror("Erro", "Arquivo não encontrado.")
            return
        try:
            dados, fs, canal_ids = carrega_dados(caminho, n_canais_bin=16, fs_bin=1000.0)
        except Exception as e:
            messagebox.showerror("Erro ao carregar", str(e))
            return
        self._dados_canais, self._fs_canais, self._canal_ids_canais = dados, fs, canal_ids
        dur = dados.shape[0] / fs
        self.var_canais_info.set(f"{len(canal_ids)} canais · {dur:.1f}s · fs={fs:.0f}Hz")
        self.lst_canais.delete(0, "end")
        for c in canal_ids:
            self.lst_canais.insert("end", c)
        self.cb_canais_ref["values"] = ["(nenhum)"] + list(canal_ids)
        self.var_canais_fim.set(min(10.0, dur))

    def _canais_plotar(self):
        if self._dados_canais is None:
            messagebox.showerror("Erro", "Carregue um arquivo primeiro.")
            return
        sel = [self.lst_canais.get(i) for i in self.lst_canais.curselection()]
        if not sel:
            messagebox.showerror("Erro", "Selecione ao menos um canal.")
            return
        t_ini, t_fim = self.var_canais_ini.get(), self.var_canais_fim.get()
        if t_fim <= t_ini:
            messagebox.showerror("Erro", "Fim deve ser maior que início.")
            return

        dados, fs, canal_ids = self._dados_canais, self._fs_canais, self._canal_ids_canais
        trecho = fatia_janela(dados, fs, t_ini, t_fim)
        t = np.linspace(t_ini, t_fim, trecho.shape[0])
        ref = self.var_canais_ref.get()
        ref_idx = resolve_canal_idx(canal_ids, ref) if ref != "(nenhum)" else None
        ganho = self.var_canais_ganho.get()

        self.fig_canais.clear()
        ax = self.fig_canais.add_subplot(111)
        for i, nome in enumerate(sel):
            idx = resolve_canal_idx(canal_ids, nome)
            sinal = trecho[:, idx].astype(float)
            if ref_idx is not None and idx != ref_idx:
                sinal = sinal - trecho[:, ref_idx].astype(float)
            if self.var_canais_detrend.get():
                sinal = detrend(sinal)
            if self.var_canais_notch.get():
                sinal = aplica_notch(sinal, fs, freqs_notch=[60, 120, 180, 240])
            if self.var_canais_passafaixa.get():
                sinal = filtra_sinal(sinal, self.var_canais_lo.get(),
                                     min(self.var_canais_hi.get(), fs * 0.49), fs)
            offset = i * (np.std(sinal) * 6 + 1e-9)
            ax.plot(t, sinal * ganho + offset, linewidth=0.7, label=nome)
        ax.set_xlabel("Tempo (s)")
        ax.set_ylabel("Canais (empilhados)")
        ax.legend(loc="upper right", fontsize=7, ncol=min(4, len(sel)))
        self.canvas_canais.draw()


    def _build_tab_comod(self):
        aba = self.tab_comod
        self.var_comod_caminho = tk.StringVar()
        self._entrada_arquivo(self._secao(aba, "Arquivo"), self.var_comod_caminho)

        f_ctrl = self._secao(aba, "Parâmetros")
        row = ttk.Frame(f_ctrl); row.pack(fill="x", pady=2)
        ttk.Button(row, text="Carregar canais", command=self._comod_carregar, bootstyle="primary").pack(side="left")
        ttk.Label(row, text="Canal:").pack(side="left", padx=(10, 2))
        self.var_comod_canal = tk.StringVar()
        self.cb_comod_canal = ttk.Combobox(row, textvariable=self.var_comod_canal,
                                           state="readonly", width=10)
        self.cb_comod_canal.pack(side="left")
        ttk.Label(row, text="Início (s):").pack(side="left", padx=(10, 2))
        self.var_comod_ini = tk.DoubleVar(value=0.0)
        ttk.Entry(row, textvariable=self.var_comod_ini, width=8).pack(side="left")
        ttk.Label(row, text="Fim (s):").pack(side="left", padx=(10, 2))
        self.var_comod_fim = tk.DoubleVar(value=10.0)
        ttk.Entry(row, textvariable=self.var_comod_fim, width=8).pack(side="left")
        ttk.Label(row, text="Par (resumo MI):").pack(side="left", padx=(10, 2))
        self.var_comod_par = tk.StringVar(value=list(BAND_PAIRS.keys())[0])
        ttk.Combobox(row, textvariable=self.var_comod_par, state="readonly",
                    values=list(BAND_PAIRS.keys()), width=12).pack(side="left")

        self.btn_comod_calcular = ttk.Button(f_ctrl, text="Calcular comodulograma", bootstyle="primary",
                                             command=self._comod_calcular)
        self.btn_comod_calcular.pack(pady=4)
        self.var_comod_status = tk.StringVar(value="")
        ttk.Label(f_ctrl, textvariable=self.var_comod_status, bootstyle="primary").pack()

        f_plot = ttk.Frame(aba); f_plot.pack(fill="both", expand=True, padx=8, pady=4)
        self.fig_comod, self.canvas_comod = self._embute_figura(f_plot)
        self.txt_comod_resumo = tk.Text(aba, height=5, wrap="word")
        self.txt_comod_resumo.pack(fill="x", padx=8, pady=(0, 8))

        self._dados_comod = None

    def _comod_carregar(self):
        caminho = self.var_comod_caminho.get()
        if not caminho or not os.path.exists(caminho):
            messagebox.showerror("Erro", "Arquivo não encontrado.")
            return
        try:
            dados, fs, canal_ids = carrega_dados(caminho, n_canais_bin=16, fs_bin=1000.0)
        except Exception as e:
            messagebox.showerror("Erro ao carregar", str(e))
            return
        self._dados_comod = (dados, fs, canal_ids)
        self.cb_comod_canal["values"] = list(canal_ids)
        if canal_ids:
            self.var_comod_canal.set(canal_ids[0])
        self.var_comod_fim.set(min(10.0, dados.shape[0] / fs))

    def _comod_calcular(self):
        if self._dados_comod is None:
            messagebox.showerror("Erro", "Carregue um arquivo primeiro.")
            return
        dados, fs, canal_ids = self._dados_comod
        canal = self.var_comod_canal.get()
        t_ini, t_fim = self.var_comod_ini.get(), self.var_comod_fim.get()
        caminho = self.var_comod_caminho.get()
        par = self.var_comod_par.get()

        idx = resolve_canal_idx(canal_ids, canal)
        trecho = fatia_janela(dados, fs, t_ini, t_fim)[:, idx].astype(float)

        self.btn_comod_calcular.state(["disabled"])
        self.var_comod_status.set("Calculando mapa fase × amplitude (200 surrogates)...")

        def job():
            z_mapa = calcula_comodulograma_z(trecho, fs, FASES_DEFAULT, AMPS_DEFAULT,
                                             n_surr=200, n_bins=18, rng=42,
                                             notch_hz=[60, 120, 180, 240])
            resumo = analisa_janela(caminho, canal, t_ini, t_fim, par=par)
            return z_mapa, resumo

        def ao_terminar(resultado):
            z_mapa, resumo = resultado
            self.fig_comod.clear()
            ax = self.fig_comod.add_subplot(111)
            im = ax.imshow(z_mapa, aspect="auto", origin="lower", cmap="turbo",
                           extent=[FASES_DEFAULT[0], FASES_DEFAULT[-1],
                                   AMPS_DEFAULT[0], AMPS_DEFAULT[-1]])
            self.fig_comod.colorbar(im, ax=ax, label="z-score")
            ax.set_xlabel("Frequência de fase (Hz)")
            ax.set_ylabel("Frequência de amplitude (Hz)")
            self.canvas_comod.draw()
            self.txt_comod_resumo.delete("1.0", "end")
            self.txt_comod_resumo.insert("1.0", f"z máximo no mapa: {np.max(z_mapa):.2f}\n"
                                                  f"Resumo MI pontual ({par}): {resumo}")
            self.var_comod_status.set("Concluído.")
            self.btn_comod_calcular.state(["!disabled"])

        def ao_erro(e):
            messagebox.showerror("Erro no cálculo", str(e))
            self.var_comod_status.set("Falhou.")
            self.btn_comod_calcular.state(["!disabled"])

        _rodar_em_thread(job, lambda r: self.after(0, ao_terminar, r),
                         lambda e: self.after(0, ao_erro, e))


    def _build_tab_portoes(self):
        aba = self.tab_portoes
        f_csv = self._secao(aba, "CSV")
        self.var_portoes_csv = tk.StringVar(value=VENCEDORES_PADRAO)
        row = ttk.Frame(f_csv); row.pack(fill="x")
        ttk.Entry(row, textvariable=self.var_portoes_csv, width=80).pack(side="left", fill="x", expand=True)
        ttk.Button(row, text="Procurar...",
                  command=lambda: self.var_portoes_csv.set(filedialog.askopenfilename() or "")
                  ).pack(side="left", padx=4)
        ttk.Button(f_csv, text="Carregar", command=self._portoes_carregar_csv, bootstyle="primary").pack(pady=4)

        f_sel = self._secao(aba, "Candidato")
        row2 = ttk.Frame(f_sel); row2.pack(fill="x", pady=2)
        ttk.Label(row2, text="Sessão:").pack(side="left")
        self.var_portoes_sessao = tk.StringVar()
        self.cb_portoes_sessao = ttk.Combobox(row2, textvariable=self.var_portoes_sessao,
                                              state="readonly", width=55)
        self.cb_portoes_sessao.pack(side="left", padx=4)
        self.cb_portoes_sessao.bind("<<ComboboxSelected>>", lambda e: self._portoes_atualiza_candidatos())

        row3 = ttk.Frame(f_sel); row3.pack(fill="x", pady=2)
        ttk.Label(row3, text="Candidato:").pack(side="left")
        self.var_portoes_candidato = tk.StringVar()
        self.cb_portoes_candidato = ttk.Combobox(row3, textvariable=self.var_portoes_candidato,
                                                 state="readonly", width=55)
        self.cb_portoes_candidato.pack(side="left", padx=4)
        self.cb_portoes_candidato.bind("<<ComboboxSelected>>", lambda e: self._portoes_mostra_checklist())

        f_check = self._secao(aba, "Checklist")
        self.var_check_estatistico = tk.StringVar()
        self.var_check_skew = tk.StringVar()
        self.var_check_harmonico = tk.StringVar()
        self.var_check_transiente = tk.StringVar()
        self.var_check_footprint = tk.StringVar()
        for var in [self.var_check_estatistico, self.var_check_skew, self.var_check_harmonico,
                   self.var_check_transiente, self.var_check_footprint]:
            ttk.Label(f_check, textvariable=var, font=("Segoe UI", 10)).pack(anchor="w", pady=1)

        self.txt_portoes_linha = tk.Text(aba, height=10, wrap="word")
        self.txt_portoes_linha.pack(fill="both", expand=True, padx=8, pady=8)

        self._df_portoes = None

    def _portoes_carregar_csv(self):
        caminho = self.var_portoes_csv.get()
        if not os.path.exists(caminho):
            messagebox.showerror("Erro", "CSV não encontrado.")
            return
        self._df_portoes = pd.read_csv(caminho, encoding="utf-8-sig")
        sessoes = sorted(self._df_portoes["sessao"].dropna().unique())
        self.cb_portoes_sessao["values"] = sessoes
        if sessoes:
            self.var_portoes_sessao.set(sessoes[0])
            self._portoes_atualiza_candidatos()

    def _portoes_atualiza_candidatos(self):
        df_sessao = self._df_portoes[self._df_portoes["sessao"] == self.var_portoes_sessao.get()]
        self._df_portoes_sessao = df_sessao
        opcoes = [f"[{i}] {r['par']} | canal {r['canal']} | {r['janela_ini_s']}-{r['janela_fim_s']}s"
                 for i, r in df_sessao.iterrows()]
        self.cb_portoes_candidato["values"] = opcoes
        if opcoes:
            self.var_portoes_candidato.set(opcoes[0])
            self._portoes_mostra_checklist()

    def _portoes_mostra_checklist(self):
        escolha = self.var_portoes_candidato.get()
        if not escolha:
            return
        idx = int(escolha.split("]")[0].lstrip("["))
        linha = self._df_portoes_sessao.loc[idx]

        def icone(ok):
            return "✅" if ok else "❌"

        self.var_check_estatistico.set(
            f"{icone(linha.get('veredito_refino') == 'Candidato robusto')} Estatístico (FDR) "
            f"— {linha.get('veredito_refino', 'n/d')}")
        self.var_check_skew.set(
            f"{icone(not str(linha.get('veredito_skew', '')).startswith('SUSPECT'))} "
            f"Skewness — {linha.get('veredito_skew', 'n/d')}")
        harm_ok = linha.get("veredito_harmonico") in ("CLEAN",) or \
            str(linha.get("veredito_harmonico", "")).startswith("REVISAR_RAZAO_INTEIRA")
        self.var_check_harmonico.set(f"{icone(harm_ok)} Harmônico — {linha.get('veredito_harmonico', 'n/d')}")
        self.var_check_transiente.set(
            f"{icone(not bool(linha.get('transiente_detectado', False)))} Transiente/saturação "
            f"— frac={linha.get('frac_transiente', 'n/d')}")
        self.var_check_footprint.set(
            f"ℹ️ Footprint (informativo) — canais coativos: "
            f"{linha.get('footprint_n_z_ge3', linha.get('n_canais_pac', 'n/d'))} de "
            f"{linha.get('footprint_n_canais', 'n/d')}")

        self.txt_portoes_linha.delete("1.0", "end")
        self.txt_portoes_linha.insert("1.0", linha.to_string())


    def _build_tab_galeria(self):
        aba = self.tab_galeria
        f_csv = self._secao(aba, "CSV de vencedores")
        self.var_galeria_csv = tk.StringVar(value=VENCEDORES_PADRAO)
        row = ttk.Frame(f_csv); row.pack(fill="x")
        ttk.Entry(row, textvariable=self.var_galeria_csv, width=80).pack(side="left", fill="x", expand=True)
        ttk.Button(row, text="Carregar", command=self._galeria_carregar, bootstyle="primary").pack(side="left", padx=4)
        self.var_galeria_info = tk.StringVar()
        ttk.Label(f_csv, textvariable=self.var_galeria_info).pack(anchor="w", pady=2)

        f_nav = ttk.Frame(aba); f_nav.pack(fill="x", padx=8)
        ttk.Button(f_nav, text="◀ Anterior", command=lambda: self._galeria_muda_pagina(-1)).pack(side="left")
        self.var_galeria_pagina = tk.StringVar(value="Página 0 de 0")
        ttk.Label(f_nav, textvariable=self.var_galeria_pagina).pack(side="left", padx=10)
        ttk.Button(f_nav, text="Próxima ▶", command=lambda: self._galeria_muda_pagina(1)).pack(side="left")

        self.frame_galeria_grid = ttk.Frame(aba)
        self.frame_galeria_grid.pack(fill="both", expand=True, padx=8, pady=8)

        self._df_galeria = None
        self._pagina_galeria = 0
        self._por_pagina_galeria = 8

    def _indexa_pngs(self):
        if self._indice_png is None:
            self._indice_png = {}
            for caminho in glob.iglob(os.path.join(workspace.BASE_RESULTADOS_MESTRADO, "**", "*.png"),
                                      recursive=True):
                self._indice_png[os.path.basename(caminho)] = caminho
        return self._indice_png

    def _galeria_carregar(self):
        caminho = self.var_galeria_csv.get()
        if not os.path.exists(caminho):
            messagebox.showerror("Erro", "CSV não encontrado.")
            return
        self._df_galeria = pd.read_csv(caminho, encoding="utf-8-sig")
        rato_info = self._df_galeria["rato"].value_counts().to_dict() if "rato" in self._df_galeria.columns else {}
        self.var_galeria_info.set(f"{len(self._df_galeria)} vencedores — {rato_info}")
        self._pagina_galeria = 0
        self._galeria_desenha_pagina()

    def _galeria_muda_pagina(self, delta):
        if self._df_galeria is None:
            return
        n_paginas = max(1, (len(self._df_galeria) - 1) // self._por_pagina_galeria + 1)
        self._pagina_galeria = max(0, min(n_paginas - 1, self._pagina_galeria + delta))
        self._galeria_desenha_pagina()

    def _galeria_desenha_pagina(self):
        for w in self.frame_galeria_grid.winfo_children():
            w.destroy()
        self._imagens_vivas.clear()
        if self._df_galeria is None or len(self._df_galeria) == 0:
            return

        n_paginas = max(1, (len(self._df_galeria) - 1) // self._por_pagina_galeria + 1)
        self.var_galeria_pagina.set(f"Página {self._pagina_galeria + 1} de {n_paginas}")
        inicio = self._pagina_galeria * self._por_pagina_galeria
        pagina_df = self._df_galeria.iloc[inicio:inicio + self._por_pagina_galeria]
        indice = self._indexa_pngs()

        for i, (_, row) in enumerate(pagina_df.iterrows()):
            r, c = divmod(i, 4)
            cel = ttk.Frame(self.frame_galeria_grid, borderwidth=1, relief="solid", padding=4)
            cel.grid(row=r, column=c, padx=4, pady=4, sticky="n")

            caminho_img = indice.get(os.path.basename(str(row.get("png", ""))))
            if caminho_img:
                try:
                    img = Image.open(caminho_img)
                    img.thumbnail((220, 165))
                    foto = ImageTk.PhotoImage(img)
                    self._imagens_vivas.append(foto)
                    ttk.Label(cel, image=foto).pack()
                except Exception:
                    ttk.Label(cel, text="(erro ao abrir imagem)").pack()
            else:
                ttk.Label(cel, text="(imagem não localizada)", width=28).pack()

            z_val = row.get("z_score_refinado", row.get("z_pico_par", "n/d"))
            ttk.Label(cel, text=f"{row.get('sessao', '')}", font=("Segoe UI", 7),
                     wraplength=200).pack()
            ttk.Label(cel, text=f"canal {row.get('canal', '')} · {row.get('par', '')} · z={z_val}",
                     font=("Segoe UI", 7)).pack()


    def _build_tab_comport(self):
        aba = self.tab_comport
        f_csv = self._secao(aba, "Dataset mestre com comportamento")
        self.var_comport_csv = tk.StringVar(value=DATASET_MESTRE_PADRAO)
        row = ttk.Frame(f_csv); row.pack(fill="x")
        ttk.Entry(row, textvariable=self.var_comport_csv, width=80).pack(side="left", fill="x", expand=True)
        ttk.Button(row, text="Carregar", command=self._comport_carregar, bootstyle="primary").pack(side="left", padx=4)

        nb2 = ttk.Notebook(aba)
        nb2.pack(fill="both", expand=True, padx=8, pady=8)
        sub_cruz = ttk.Frame(nb2); nb2.add(sub_cruz, text="LFP × Comportamento")
        sub_transicao = ttk.Frame(nb2); nb2.add(sub_transicao, text="Transição N-1 → N (χ²)")
        sub_pharm = ttk.Frame(nb2); nb2.add(sub_pharm, text="Basal vs NOCI vs LAC")

        ttk.Label(sub_cruz, text="Métrica:").pack(anchor="w", padx=4, pady=(4, 0))
        self.var_comport_metrica = tk.StringVar(value="z_score_refinado")
        ttk.Combobox(sub_cruz, textvariable=self.var_comport_metrica, state="readonly",
                    values=["z_score_refinado", "mi_observado"], width=20).pack(anchor="w", padx=4)
        ttk.Button(sub_cruz, text="Plotar distribuição", bootstyle="primary",
                  command=self._comport_plota_cruzamento).pack(anchor="w", padx=4, pady=4)
        self.fig_comport_cruz, self.canvas_comport_cruz = self._embute_figura(sub_cruz, figsize=(9, 4))

        ttk.Button(sub_transicao, text="Calcular transição N-1→N", bootstyle="primary",
                  command=self._comport_calcula_transicao).pack(anchor="w", padx=4, pady=4)
        self.var_comport_chi2 = tk.StringVar(value="")
        ttk.Label(sub_transicao, textvariable=self.var_comport_chi2,
                 font=("Segoe UI", 10, "bold")).pack(anchor="w", padx=4)
        self.fig_comport_transicao, self.canvas_comport_transicao = self._embute_figura(
            sub_transicao, figsize=(9, 4))

        self.var_comport_pharm_msg = tk.StringVar(value="")
        ttk.Label(sub_pharm, textvariable=self.var_comport_pharm_msg, wraplength=900,
                 justify="left").pack(anchor="w", padx=4, pady=4)
        ttk.Button(sub_pharm, text="Atualizar painel farmacológico", bootstyle="primary",
                  command=self._comport_plota_pharm).pack(anchor="w", padx=4)
        self.fig_comport_pharm, self.canvas_comport_pharm = self._embute_figura(sub_pharm, figsize=(9, 4))

        self._df_comport = None

    def _comport_carregar(self):
        caminho = self.var_comport_csv.get()
        if not os.path.exists(caminho):
            messagebox.showerror("Erro", "CSV não encontrado.")
            return
        self._df_comport = pd.read_csv(caminho, encoding="utf-8-sig")
        messagebox.showinfo("OK", f"{len(self._df_comport)} linhas carregadas.")

    def _comport_plota_cruzamento(self):
        if self._df_comport is None:
            messagebox.showerror("Erro", "Carregue o CSV primeiro.")
            return
        metrica = self.var_comport_metrica.get()
        df = self._df_comport
        df_anotado = df[df["comportamento"].notna() & (df["comportamento"] != "Artefato / Cabo")]
        if metrica not in df_anotado.columns:
            messagebox.showerror("Erro", f"Coluna '{metrica}' não disponível.")
            return
        cats = sorted(df_anotado["comportamento"].dropna().unique())
        dados_box = [df_anotado.loc[df_anotado["comportamento"] == c, metrica].dropna() for c in cats]

        self.fig_comport_cruz.clear()
        ax = self.fig_comport_cruz.add_subplot(111)
        ax.boxplot(dados_box, tick_labels=cats, showfliers=False)
        ax.set_ylabel(metrica)
        ax.tick_params(axis="x", rotation=30, labelsize=7)
        self.canvas_comport_cruz.draw()

    def _comport_calcula_transicao(self):
        if self._df_comport is None:
            messagebox.showerror("Erro", "Carregue o CSV primeiro.")
            return
        try:
            consec = transicao_estado(self._df_comport)
            tabela = pd.crosstab(consec["comport_anterior"], consec["vencedor"])
            chi2, p_valor, dof, _ = chi2_contingency(tabela)
        except KeyError as e:
            messagebox.showerror("Erro", f"Coluna esperada ausente: {e}")
            return

        self.var_comport_chi2.set(f"n = {len(consec)} pares consecutivos · χ² = {chi2:.2f} "
                                  f"(dof={dof}) · p = {p_valor:.3e}")
        prop = (tabela[True] / tabela.sum(axis=1) * 100).sort_values(ascending=False) \
            if True in tabela.columns else pd.Series(dtype=float)
        self.fig_comport_transicao.clear()
        ax = self.fig_comport_transicao.add_subplot(111)
        ax.bar(prop.index, prop.values)
        ax.set_ylabel("% janelas com PAC significativo (N)")
        ax.set_xlabel("Estado comportamental (N-1)")
        ax.tick_params(axis="x", rotation=30, labelsize=7)
        self.canvas_comport_transicao.draw()

    def _comport_plota_pharm(self):
        if self._df_comport is None:
            messagebox.showerror("Erro", "Carregue o CSV primeiro.")
            return
        df = self._df_comport
        if "grupo" not in df.columns or "condicao" not in df.columns:
            self.var_comport_pharm_msg.set(
                "Este CSV ainda não tem as colunas 'grupo'/'condicao' (produzidas pelo "
                "agrega_resultados.py já corrigido). Regenere depois que a infusão e o LAC "
                "forem reprocessados para habilitar este painel.")
            self.fig_comport_pharm.clear()
            self.canvas_comport_pharm.draw()
            return

        vencedor_col = df["veredito_refino"] == "Candidato robusto"
        resumo = (df.assign(vencedor=vencedor_col).groupby(["grupo", "condicao"])["vencedor"]
                 .mean().mul(100).reset_index())
        if resumo["condicao"].nunique() < 2 and resumo["grupo"].nunique() < 2:
            self.var_comport_pharm_msg.set(
                "Ainda só há uma combinação grupo/condição neste CSV -- sem variação "
                "para comparar Basal vs Nociceptina vs Lactato.")
            return

        self.var_comport_pharm_msg.set("")
        self.fig_comport_pharm.clear()
        ax = self.fig_comport_pharm.add_subplot(111)
        grupos = sorted(resumo["grupo"].unique())
        condicoes = sorted(resumo["condicao"].unique())
        largura = 0.8 / max(1, len(grupos))
        x = np.arange(len(condicoes))
        for i, g in enumerate(grupos):
            sub = resumo[resumo["grupo"] == g].set_index("condicao").reindex(condicoes)
            ax.bar(x + i * largura, sub["vencedor"].values, width=largura, label=g)
        ax.set_xticks(x + largura * (len(grupos) - 1) / 2)
        ax.set_xticklabels(condicoes)
        ax.set_ylabel("% candidatos robustos")
        ax.legend()
        self.canvas_comport_pharm.draw()


    def _build_tab_closedloop(self):
        aba = self.tab_closedloop
        aviso = ("Protótipo não validado para hardware real (ver preditor/README_preditor.md, "
                 "AUC-ROC 0.618). Esta aba é um replay offline -- não aciona hardware nenhum, "
                 "dispara_ttl() só imprime no console.")
        self._callout(aba, aviso)

        self.var_cl_caminho = tk.StringVar()
        self._entrada_arquivo(self._secao(aba, "Arquivo"), self.var_cl_caminho)

        f_ctrl = self._secao(aba, "Parâmetros")
        row = ttk.Frame(f_ctrl); row.pack(fill="x", pady=2)
        ttk.Button(row, text="Carregar canais", command=self._cl_carregar, bootstyle="primary").pack(side="left")
        ttk.Label(row, text="Canal:").pack(side="left", padx=(10, 2))
        self.var_cl_canal = tk.StringVar()
        self.cb_cl_canal = ttk.Combobox(row, textvariable=self.var_cl_canal, state="readonly", width=10)
        self.cb_cl_canal.pack(side="left")
        ttk.Label(row, text="Comportamento (contexto):").pack(side="left", padx=(10, 2))
        self.var_cl_comport = tk.StringVar(value=CATEGORIAS_COMPORTAMENTO[0])
        ttk.Combobox(row, textvariable=self.var_cl_comport, state="readonly",
                    values=CATEGORIAS_COMPORTAMENTO, width=24).pack(side="left")

        row2 = ttk.Frame(f_ctrl); row2.pack(fill="x", pady=2)
        ttk.Label(row2, text="Início (s):").pack(side="left")
        self.var_cl_ini = tk.DoubleVar(value=0.0)
        ttk.Entry(row2, textvariable=self.var_cl_ini, width=8).pack(side="left", padx=(2, 10))
        ttk.Label(row2, text="Fim (s):").pack(side="left")
        self.var_cl_fim = tk.DoubleVar(value=120.0)
        ttk.Entry(row2, textvariable=self.var_cl_fim, width=8).pack(side="left", padx=(2, 10))
        ttk.Label(row2, text="Limiar TTL:").pack(side="left")
        self.var_cl_limiar = tk.DoubleVar(value=0.431)
        ttk.Scale(row2, from_=0.0, to=1.0, variable=self.var_cl_limiar,
                 orient="horizontal", length=150).pack(side="left", padx=4)
        ttk.Label(row2, textvariable=self.var_cl_limiar).pack(side="left")
        ttk.Label(f_ctrl, text="Limiar padrão (0.431) calibrado por F1: recall≈0.896, "
                              "precisão≈0.452 (ver preditor/README_preditor.md). Ajustar o "
                              "slider só recorta onde ESTE replay dispararia.",
                 font=("Segoe UI", 8), foreground="#666", wraplength=1000).pack(anchor="w")

        row3 = ttk.Frame(f_ctrl); row3.pack(fill="x", pady=4)
        self.btn_cl_rodar = ttk.Button(row3, text="Rodar replay", command=self._cl_rodar, bootstyle="primary")
        self.btn_cl_rodar.pack(side="left")
        self.pb_cl = ttk.Progressbar(row3, length=300, mode="determinate")
        self.pb_cl.pack(side="left", padx=10)
        self.var_cl_status = tk.StringVar()
        ttk.Label(row3, textvariable=self.var_cl_status).pack(side="left")
        ttk.Button(row3, text="Exportar relatório (HTML)", command=self._cl_exportar).pack(side="right")

        f_plot = ttk.Frame(aba); f_plot.pack(fill="both", expand=True, padx=8, pady=4)
        self.fig_cl, self.canvas_cl = self._embute_figura(f_plot)

        self._dados_cl = None
        self._ultimo_replay = None

    def _cl_carregar(self):
        caminho = self.var_cl_caminho.get()
        if not caminho or not os.path.exists(caminho):
            messagebox.showerror("Erro", "Arquivo não encontrado.")
            return
        dados, fs, canal_ids = carrega_dados(caminho, n_canais_bin=16, fs_bin=1000.0)
        self._dados_cl = (dados, fs, canal_ids)
        self.cb_cl_canal["values"] = list(canal_ids)
        if canal_ids:
            self.var_cl_canal.set(canal_ids[0])
        self.var_cl_fim.set(min(120.0, dados.shape[0] / fs))

    def _cl_rodar(self):
        if self._dados_cl is None:
            messagebox.showerror("Erro", "Carregue um arquivo primeiro.")
            return
        dados, fs, canal_ids = self._dados_cl
        canal = self.var_cl_canal.get()
        comport = self.var_cl_comport.get()
        t_ini, t_fim = self.var_cl_ini.get(), self.var_cl_fim.get()
        caminho = self.var_cl_caminho.get()

        self.btn_cl_rodar.state(["disabled"])
        self.pb_cl["value"] = 0
        self.var_cl_status.set("Rodando...")

        def job():
            modelo, features, _ = carrega_modelo()
            ch = resolve_canal(canal, canal_ids)
            n_amostras = JANELA_S * int(fs)
            passo = PASSO_S * int(fs)
            i_ini, i_fim = int(t_ini * fs), int(t_fim * fs)
            passos_totais = max(1, (i_fim - i_ini - n_amostras) // passo)

            ts, probs = [], []
            for k, inicio in enumerate(range(i_ini, i_fim - n_amostras + 1, passo)):
                janela = dados[inicio:inicio + n_amostras, ch].astype(float)
                X = extrair_features_estado(janela, fs, comport, features)
                p_pac = modelo.predict_proba(X)[0, 1]
                ts.append(inicio / fs)
                probs.append(float(p_pac))
                if k % 5 == 0:
                    pct = min(100, int(100 * k / passos_totais))
                    self.after(0, lambda p=pct: self.pb_cl.configure(value=p))
            return np.array(ts), np.array(probs)

        def ao_terminar(resultado):
            ts, probs = resultado
            limiar = self.var_cl_limiar.get()
            dispara = probs >= limiar
            self._ultimo_replay = dict(ts=ts, probs=probs, limiar=limiar, canal=canal,
                                       comport=comport, caminho=caminho, t_ini=t_ini, t_fim=t_fim,
                                       dispara=dispara)

            self.fig_cl.clear()
            ax = self.fig_cl.add_subplot(111)
            ax.plot(ts, probs, linewidth=1, label="P(PAC)")
            ax.scatter(ts[dispara], probs[dispara], color="red", s=20, zorder=3,
                      label="dispararia TTL")
            ax.axhline(limiar, linestyle="--", color="gray")
            ax.set_xlabel("Tempo (s)")
            ax.set_ylabel("P(PAC)")
            ax.legend(loc="upper right", fontsize=8)
            self.canvas_cl.draw()

            self.pb_cl["value"] = 100
            self.var_cl_status.set(f"{int(dispara.sum())}/{len(probs)} janelas "
                                   f"disparariam ({100 * dispara.mean():.1f}%).")
            self.btn_cl_rodar.state(["!disabled"])

        def ao_erro(e):
            messagebox.showerror("Erro no replay", str(e))
            self.var_cl_status.set("Falhou.")
            self.btn_cl_rodar.state(["!disabled"])

        _rodar_em_thread(job, lambda r: self.after(0, ao_terminar, r),
                         lambda e: self.after(0, ao_erro, e))

    def _cl_exportar(self):
        if self._ultimo_replay is None:
            messagebox.showerror("Erro", "Rode um replay primeiro.")
            return
        r = self._ultimo_replay
        destino = filedialog.asksaveasfilename(defaultextension=".html",
                                               filetypes=[("HTML", "*.html")])
        if not destino:
            return
        html = f"""<html><head><meta charset="utf-8"><title>Relatorio Replay</title></head>
<body>
<h1>ThetaGamma-Studio — Relatório de Replay</h1>
<p><b>Arquivo:</b> {os.path.basename(r['caminho'])} | <b>Canal:</b> {r['canal']} |
<b>Comportamento (contexto):</b> {r['comport']}</p>
<p><b>Janela:</b> {r['t_ini']:.0f}s–{r['t_fim']:.0f}s | <b>Limiar:</b> {r['limiar']:.3f} |
<b>Disparos:</b> {int(r['dispara'].sum())}/{len(r['probs'])} ({100 * r['dispara'].mean():.1f}%)</p>
<p style="color:#a15c00"><b>Aviso:</b> protótipo não validado para hardware real
(AUC-ROC 0.618, ver preditor/README_preditor.md). Replay offline apenas.</p>
</body></html>"""
        with open(destino, "w", encoding="utf-8") as f:
            f.write(html)
        messagebox.showinfo("OK", f"Relatório salvo em {destino}")

