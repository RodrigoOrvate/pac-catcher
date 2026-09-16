"""
dashboard_desktop/aba_pipeline.py — seção "Pipeline" do ThetaGamma-Studio.

Cada sub-aba é um formulário que monta a linha de comando de um script JÁ
existente em pipeline/ e o roda via PainelExecucao (log ao vivo). Nenhuma
lógica científica aqui: os campos espelham as flags argparse reais de cada
script (levantadas em 2026-09-14) e a validação de tipo continua sendo do
próprio argparse (erro aparece no log com exit != 0).

Fora daqui de propósito: audita_transientes.py e audita_footprint.py -- CLI
presa a uma sessão específica (--pasta_ns2 default "../Basal antes da
infusao"; --casos/--csv recebem string de casos inline, não caminho). São
ferramentas forenses de investigação pontual, rodar pelo terminal.
"""
import os
import sys
import tkinter as tk
from tkinter import filedialog, messagebox

import pandas as pd
import ttkbootstrap as ttk

from pac_core import workspace
from pac_studio import BAND_PAIRS
from dashboard_desktop.runner import PainelExecucao

_PIPE = os.path.join(workspace.BASE_SCRIPT, "pipeline")
SCRIPTS = {
    "sessao": os.path.join(_PIPE, "processa_sessao.py"),
    "triagem": os.path.join(_PIPE, "etapa1_triagem", "triagem_pac.py"),
    "refino": os.path.join(_PIPE, "etapa2_refinamento", "refina_candidatos.py"),
    "comod": os.path.join(_PIPE, "etapa3_comodulograma", "comodulogram.py"),
    "skewness": os.path.join(_PIPE, "auditorias", "audita_skewness.py"),
    "harmonico": os.path.join(_PIPE, "auditorias", "audita_harmonico.py"),
    "harmonico_hfo": os.path.join(_PIPE, "auditorias", "audita_harmonico_hfo.py"),
    "agrega": os.path.join(_PIPE, "dataset_mestre", "agrega_resultados.py"),
    "enriquece": os.path.join(_PIPE, "dataset_mestre", "enriquece_dataset_mestre.py"),
    "junta": os.path.join(_PIPE, "comportamento", "junta_comportamento.py"),
    "consolida": os.path.join(_PIPE, "dataset_mestre", "consolida_vencedores.py"),
}
_RES = os.path.join(workspace.BASE_SCRIPT, "resultados")
_TEMPLATE = os.path.join(_PIPE, "comportamento", "template_comportamento.csv")

# Parâmetros específicos de cada auditoria: (flag, rótulo, default).
# `combo` = (flag, rótulo, opções) pra flags com choices no argparse.
AUDITORIAS = {
    "Skewness (dente-de-serra)": dict(
        chave="skewness",
        dica="CSV com arquivo/canal/janela_ini_s/janela_fim_s -- processa_sessao.py "
             "usa o _refinados_dedup_janela.csv do canal.",
        params=[("--limiar", "Limiar |skewness|", "0.5")],
    ),
    "Harmônico (FOOOF + PLV)": dict(
        chave="harmonico",
        dica="O CSV precisa ter fase_pico_hz/amp_pico_hz -- use o resumo de "
             "comodulogramas do canal (_resumo_<par>.csv), não o refinados.csv.",
        params=[("--janela_contexto_s", "Contexto FOOOF (s)", "45.0"),
                ("--tol_rel", "Tolerância relativa", "0.1"),
                ("--limiar_plv", "Limiar PLV", "0.8"),
                ("--limiar_skew", "Limiar skewness", "0.5"),
                ("--f_linha", "Rede elétrica (Hz)", "60.0")],
        combo=("--modo_preprocesso", "Pré-processo", ["hibrido", "gaussiana", "cirurgica"]),
    ),
    "Harmônico HFO": dict(
        chave="harmonico_hfo",
        dica="CSV de triagem com a coluna z_theta_hfo (triagem rodada com o par theta_hfo).",
        params=[("--z_corte", "z mínimo theta_hfo", "2.0"),
                ("--janela_contexto_s", "Contexto FOOOF (s)", "30.0"),
                ("--tol_rel", "Tolerância relativa", "0.15"),
                ("--limiar_plv", "Limiar PLV", "0.8"),
                ("--limiar_ratio", "Razão HFO/Gama", "0.3"),
                ("--f_linha", "Rede elétrica (Hz)", "60.0")],
    ),
}


class PipelineMixin:
    def _build_secao_pipeline(self, parent):
        self._paineis = []
        self._callout(
            parent,
            "Cada sub-aba roda um script já existente de pipeline/ (os mesmos que o lote "
            "usa pelo terminal). O log aparece ao vivo; \"Parar\" encerra o processo e "
            "os que ele disparou.",
            bootstyle="info",
        )
        nb = ttk.Notebook(parent)
        nb.pack(fill="both", expand=True)
        for titulo, construtor in [
            ("Sessão Completa", self._pl_aba_sessao),
            ("Triagem", self._pl_aba_triagem),
            ("Refinamento e Filtros", self._pl_aba_refino),
            ("Comodulograma (lote)", self._pl_aba_comod),
            ("Auditorias", self._pl_aba_auditorias),
            ("Agregação", self._pl_aba_agregacao),
        ]:
            f = ttk.Frame(nb, padding=4)
            nb.add(f, text=titulo)
            construtor(f)

    # ------------------------------------------------------------------
    # Helpers de formulário
    # ------------------------------------------------------------------
    def _pl_linha(self, parent):
        f = ttk.Frame(parent)
        f.pack(fill="x", pady=3)
        return f

    def _pl_campo(self, linha, rotulo, var, largura=10, procurar=None):
        """procurar: None | 'arquivo' | 'pasta' | 'salvar'."""
        ttk.Label(linha, text=rotulo).pack(side="left", padx=(0, 4))
        largo = largura >= 40
        ttk.Entry(linha, textvariable=var, width=largura).pack(
            side="left", padx=(0, 4), **({"fill": "x", "expand": True} if largo else {}))
        if procurar:
            def escolher():
                if procurar == "pasta":
                    v = filedialog.askdirectory()
                elif procurar == "salvar":
                    v = filedialog.asksaveasfilename(defaultextension=".csv",
                                                     filetypes=[("CSV", "*.csv")])
                else:
                    v = filedialog.askopenfilename()
                if v:
                    var.set(v)
            ttk.Button(linha, text="Procurar...", bootstyle="secondary-outline",
                       command=escolher).pack(side="left", padx=(0, 12))
        else:
            ttk.Frame(linha, width=12).pack(side="left")

    def _pl_painel(self, parent, altura_log=8):
        painel = PainelExecucao(parent, altura_log=altura_log)
        self._paineis.append(painel)
        return painel

    def _pl_valida(self, campos):
        """campos: lista de (rótulo, valor, tipo) com tipo 'arquivo'|'pasta'|'texto'."""
        for rotulo, valor, tipo in campos:
            valor = (valor or "").strip()
            if not valor:
                messagebox.showerror("Campo obrigatório", f"Preencha: {rotulo}")
                return False
            if tipo == "arquivo" and not os.path.isfile(valor):
                messagebox.showerror("Arquivo não encontrado", f"{rotulo}:\n{valor}")
                return False
            if tipo == "pasta" and not os.path.isdir(valor):
                messagebox.showerror("Pasta não encontrada", f"{rotulo}:\n{valor}")
                return False
        return True

    def _pl_confirma_sobrescrita(self, caminho):
        if caminho and os.path.exists(caminho):
            return messagebox.askyesno(
                "Sobrescrever?", f"O arquivo de saída já existe:\n{caminho}\n\nSobrescrever?")
        return True

    def _pl_rodar(self, painel, chave, args, ao_concluir=None):
        if painel.executando():
            messagebox.showinfo("Aguarde", "Já tem um processo rodando neste painel.")
            return
        cmd = [sys.executable, SCRIPTS[chave]] + [str(a) for a in args]
        painel.executar(cmd, cwd=workspace.BASE_SCRIPT, ao_concluir=ao_concluir)

    # ------------------------------------------------------------------
    # 1. Sessão completa -- processa_sessao.py
    # ------------------------------------------------------------------
    def _pl_aba_sessao(self, aba):
        f = self._secao(aba, "Parâmetros")
        self.var_pl_ses_pasta = tk.StringVar()
        self.var_pl_ses_saida = tk.StringVar(value=workspace.BASE_RESULTADOS_MESTRADO)
        self.var_pl_ses_min = tk.StringVar(value="3")
        self._pl_campo(self._pl_linha(f), "Pasta da sessão (.ns2):", self.var_pl_ses_pasta, 70, "pasta")
        self._pl_campo(self._pl_linha(f), "Pasta base de saída:", self.var_pl_ses_saida, 70, "pasta")
        self._pl_campo(self._pl_linha(f), "Mínimo de janelas promissoras (estágio 1):",
                       self.var_pl_ses_min, 6)
        ttk.Label(f, bootstyle="secondary", text=(
            "Roda triagem → refinamento → comodulograma → auditorias (skewness e harmônico) "
            "canal a canal. Saída em <pasta base>/<nome da pasta da sessão>/chanN/. "
            "Uma sessão costuma levar de 10 a 90 minutos.")).pack(anchor="w", pady=(4, 0))

        painel = self._pl_painel(aba)
        ttk.Button(f, text="Rodar sessão completa", bootstyle="primary",
                   command=lambda: self._pl_sessao_rodar(painel)).pack(anchor="w", pady=(8, 0))
        painel.pack(fill="both", expand=True, padx=4, pady=(10, 4))

    def _pl_sessao_rodar(self, painel):
        pasta, saida = self.var_pl_ses_pasta.get(), self.var_pl_ses_saida.get()
        if not self._pl_valida([("Pasta da sessão", pasta, "pasta"),
                                ("Pasta base de saída", saida, "pasta")]):
            return
        self._pl_rodar(painel, "sessao", ["--pasta", pasta, "--saida", saida,
                                          "--min_janelas", self.var_pl_ses_min.get()])

    # ------------------------------------------------------------------
    # 2. Triagem -- triagem_pac.py
    # ------------------------------------------------------------------
    def _pl_aba_triagem(self, aba):
        f = self._secao(aba, "Parâmetros")
        self.var_pl_tri_pasta = tk.StringVar()
        self.var_pl_tri_canal = tk.StringVar()
        self.var_pl_tri_janela = tk.StringVar(value="10.0")
        self.var_pl_tri_passo = tk.StringVar(value="5.0")
        self.var_pl_tri_nsurr = tk.StringVar(value="200")
        self.var_pl_tri_zcorte = tk.StringVar(value="3.0")
        self.var_pl_tri_saida = tk.StringVar(value=os.path.join(_RES, "triagem_avulsa.csv"))
        self.var_pl_tri_pares = {p: tk.BooleanVar(value=(p == "theta_gamma")) for p in BAND_PAIRS}

        self._pl_campo(self._pl_linha(f), "Pasta com .ns2:", self.var_pl_tri_pasta, 70, "pasta")
        linha = self._pl_linha(f)
        self._pl_campo(linha, "Canal (índice 0-based, vazio = todos):", self.var_pl_tri_canal, 6)
        self._pl_campo(linha, "Janela (s):", self.var_pl_tri_janela, 6)
        self._pl_campo(linha, "Passo (s):", self.var_pl_tri_passo, 6)
        self._pl_campo(linha, "Surrogates:", self.var_pl_tri_nsurr, 6)
        self._pl_campo(linha, "z de corte:", self.var_pl_tri_zcorte, 6)
        linha = self._pl_linha(f)
        ttk.Label(linha, text="Pares:").pack(side="left", padx=(0, 6))
        for par, var in self.var_pl_tri_pares.items():
            ttk.Checkbutton(linha, text=par, variable=var).pack(side="left", padx=(0, 10))
        self._pl_campo(self._pl_linha(f), "CSV de saída:", self.var_pl_tri_saida, 70, "salvar")

        painel = self._pl_painel(aba)
        botoes = ttk.Frame(f)
        botoes.pack(anchor="w", pady=(8, 0))
        ttk.Button(botoes, text="Rodar triagem", bootstyle="primary",
                   command=lambda: self._pl_triagem_rodar(painel)).pack(side="left")
        ttk.Button(botoes, text="Demo (teste sintético)", bootstyle="secondary-outline",
                   command=lambda: self._pl_rodar(painel, "triagem", ["--demo"])
                   ).pack(side="left", padx=8)
        painel.pack(fill="both", expand=True, padx=4, pady=(10, 4))

    def _pl_triagem_rodar(self, painel):
        pasta, saida = self.var_pl_tri_pasta.get(), self.var_pl_tri_saida.get()
        pares = [p for p, v in self.var_pl_tri_pares.items() if v.get()]
        if not self._pl_valida([("Pasta com .ns2", pasta, "pasta"), ("CSV de saída", saida, "texto")]):
            return
        if not pares:
            messagebox.showerror("Campo obrigatório", "Marque ao menos um par.")
            return
        if not self._pl_confirma_sobrescrita(saida):
            return
        args = ["--pasta", pasta, "--janela", self.var_pl_tri_janela.get(),
                "--passo", self.var_pl_tri_passo.get(), "--n_surr", self.var_pl_tri_nsurr.get(),
                "--z_corte", self.var_pl_tri_zcorte.get(), "--saida", saida, "--pares", *pares]
        if self.var_pl_tri_canal.get().strip():
            args += ["--canal", self.var_pl_tri_canal.get().strip()]
        self._pl_rodar(painel, "triagem", args)

    # ------------------------------------------------------------------
    # 3. Refinamento e filtros -- refina_candidatos.py
    # ------------------------------------------------------------------
    def _pl_aba_refino(self, aba):
        f = self._secao(aba, "Parâmetros")
        self.var_pl_ref_csv = tk.StringVar()
        self.var_pl_ref_pasta = tk.StringVar()
        self.var_pl_ref_zpre = tk.StringVar(value="2.0")
        self.var_pl_ref_nsurr = tk.StringVar(value="1000")
        self.var_pl_ref_fdr = tk.StringVar(value="0.05")
        self.var_pl_ref_saida = tk.StringVar(value=os.path.join(_RES, "refinados_avulso.csv"))

        self._pl_campo(self._pl_linha(f), "CSV da triagem:", self.var_pl_ref_csv, 70, "arquivo")
        self._pl_campo(self._pl_linha(f), "Pasta com os .ns2:", self.var_pl_ref_pasta, 70, "pasta")
        linha = self._pl_linha(f)
        self._pl_campo(linha, "z pré-filtro:", self.var_pl_ref_zpre, 6)
        self._pl_campo(linha, "Surrogates:", self.var_pl_ref_nsurr, 7)
        self._pl_campo(linha, "FDR q (Benjamini-Hochberg):", self.var_pl_ref_fdr, 6)
        self._pl_campo(self._pl_linha(f), "CSV de saída:", self.var_pl_ref_saida, 70, "salvar")

        corpo = ttk.Frame(aba)
        esquerda = ttk.Frame(corpo)
        direita = ttk.Frame(corpo)
        painel = self._pl_painel(esquerda)
        ttk.Button(f, text="Rodar refinamento", bootstyle="primary",
                   command=lambda: self._pl_refino_rodar(painel)).pack(anchor="w", pady=(8, 0))
        corpo.pack(fill="both", expand=True, padx=4, pady=(10, 4))
        esquerda.pack(side="left", fill="both", expand=True)
        direita.pack(side="left", fill="both", padx=(8, 0))
        painel.pack(fill="both", expand=True)
        ttk.Label(direita, text="VEREDITO DOS CANDIDATOS", style="SecaoTitulo.TLabel").pack(anchor="w")
        self.fig_pl_ref, self.canvas_pl_ref = self._embute_figura(direita, figsize=(5, 3.6))

    def _pl_refino_rodar(self, painel):
        csv_in, pasta = self.var_pl_ref_csv.get(), self.var_pl_ref_pasta.get()
        saida = self.var_pl_ref_saida.get()
        if not self._pl_valida([("CSV da triagem", csv_in, "arquivo"),
                                ("Pasta com os .ns2", pasta, "pasta"),
                                ("CSV de saída", saida, "texto")]):
            return
        if not self._pl_confirma_sobrescrita(saida):
            return
        args = ["--csv", csv_in, "--pasta_ns2", pasta, "--z_pre_filtro", self.var_pl_ref_zpre.get(),
                "--n_surr", self.var_pl_ref_nsurr.get(), "--fdr_q", self.var_pl_ref_fdr.get(),
                "--saida", saida]
        self._pl_rodar(painel, "refino", args,
                       ao_concluir=lambda codigo: self._pl_refino_grafico(codigo, saida))

    def _pl_refino_grafico(self, codigo, saida):
        """Fecha o ciclo 'mudei o filtro, vi o efeito': distribuição de
        veredito do CSV que acabou de ser gerado."""
        self.fig_pl_ref.clear()
        if codigo != 0 or not os.path.exists(saida):
            self.canvas_pl_ref.draw()
            return
        df = pd.read_csv(saida)
        col = "veredito" if "veredito" in df.columns else "veredito_refino"
        if col not in df.columns or df.empty:
            self.canvas_pl_ref.draw()
            return
        contagem = df[col].fillna("(vazio)").value_counts().sort_values()
        ax = self.fig_pl_ref.add_subplot(111)
        cores = [self.style.colors.primary if v == "Candidato robusto" else self.style.colors.secondary
                 for v in contagem.index]
        ax.barh([str(v)[:38] for v in contagem.index], contagem.values, color=cores)
        ax.set_title(f"{len(df)} janelas testadas", fontsize=9)
        ax.tick_params(axis="y", labelsize=7)
        ax.tick_params(axis="x", labelsize=8)
        for lado in ("top", "right"):
            ax.spines[lado].set_visible(False)
        self.fig_pl_ref.tight_layout()
        self.canvas_pl_ref.draw()

    # ------------------------------------------------------------------
    # 4. Comodulograma em lote -- comodulogram.py --csv
    # ------------------------------------------------------------------
    def _pl_aba_comod(self, aba):
        f = self._secao(aba, "Parâmetros")
        self.var_pl_com_csv = tk.StringVar()
        self.var_pl_com_pasta = tk.StringVar()
        self.var_pl_com_par = tk.StringVar(value="theta_gamma")
        self.var_pl_com_vered = tk.StringVar(value="Candidato robusto")
        self.var_pl_com_topn = tk.StringVar()
        self.var_pl_com_nsurr = tk.StringVar(value="200")
        self.var_pl_com_fdr = tk.StringVar(value="0.05")
        self.var_pl_com_saida = tk.StringVar(value=os.path.join(_RES, "comodulogramas_avulso"))
        self.var_pl_com_notch = {hz: tk.BooleanVar(value=True) for hz in (60, 120, 180, 240)}

        self._pl_campo(self._pl_linha(f), "CSV refinado:", self.var_pl_com_csv, 70, "arquivo")
        self._pl_campo(self._pl_linha(f), "Pasta com os .ns2:", self.var_pl_com_pasta, 70, "pasta")
        linha = self._pl_linha(f)
        ttk.Label(linha, text="Par destacado:").pack(side="left", padx=(0, 4))
        ttk.Combobox(linha, textvariable=self.var_pl_com_par, values=list(BAND_PAIRS),
                     state="readonly", width=12).pack(side="left", padx=(0, 16))
        self._pl_campo(linha, "Prefixo de veredito (vazio = todos):", self.var_pl_com_vered, 18)
        self._pl_campo(linha, "Top N (vazio = todos):", self.var_pl_com_topn, 5)
        linha = self._pl_linha(f)
        self._pl_campo(linha, "Surrogates:", self.var_pl_com_nsurr, 6)
        self._pl_campo(linha, "FDR q do mapa (vazio = sem FDR):", self.var_pl_com_fdr, 6)
        ttk.Label(linha, text="Notch:").pack(side="left", padx=(8, 4))
        for hz, var in self.var_pl_com_notch.items():
            ttk.Checkbutton(linha, text=f"{hz} Hz", variable=var).pack(side="left", padx=(0, 6))
        self._pl_campo(self._pl_linha(f), "Pasta de saída dos PNGs:", self.var_pl_com_saida, 70, "pasta")

        painel = self._pl_painel(aba)
        ttk.Button(f, text="Gerar comodulogramas", bootstyle="primary",
                   command=lambda: self._pl_comod_rodar(painel)).pack(anchor="w", pady=(8, 0))
        painel.pack(fill="both", expand=True, padx=4, pady=(10, 4))

    def _pl_comod_rodar(self, painel):
        csv_in, pasta = self.var_pl_com_csv.get(), self.var_pl_com_pasta.get()
        saida = self.var_pl_com_saida.get()
        if not self._pl_valida([("CSV refinado", csv_in, "arquivo"),
                                ("Pasta com os .ns2", pasta, "pasta"),
                                ("Pasta de saída", saida, "texto")]):
            return
        notch = [str(hz) for hz, v in self.var_pl_com_notch.items() if v.get()]
        if not notch:
            # --notch é nargs="+" (não aceita lista vazia) e o projeto exige
            # notch da rede em tudo (CLAUDE.md, convenção 1).
            messagebox.showerror("Notch", "Marque ao menos uma frequência de notch "
                                          "(convenção do projeto: notch da rede em tudo).")
            return
        args = ["--csv", csv_in, "--pasta_ns2", pasta, "--par", self.var_pl_com_par.get(),
                "--veredito_prefixo", self.var_pl_com_vered.get(),
                "--n_surr", self.var_pl_com_nsurr.get(), "--saida_dir", saida, "--notch", *notch]
        if self.var_pl_com_topn.get().strip():
            args += ["--top_n", self.var_pl_com_topn.get().strip()]
        if self.var_pl_com_fdr.get().strip():
            args += ["--fdr_q", self.var_pl_com_fdr.get().strip()]
        self._pl_rodar(painel, "comod", args)

    # ------------------------------------------------------------------
    # 5. Auditorias
    # ------------------------------------------------------------------
    def _pl_aba_auditorias(self, aba):
        f = self._secao(aba, "Auditoria")
        ttk.Label(f, bootstyle="secondary", text=(
            "Transientes e footprint ficam fora de propósito: a linha de comando delas é presa a "
            "uma sessão específica (casos como texto, pasta padrão fixa) -- rode pelo terminal."
        )).pack(anchor="w", pady=(0, 4))
        self.var_pl_aud_tipo = tk.StringVar(value=list(AUDITORIAS)[0])
        self.var_pl_aud_csv = tk.StringVar()
        self.var_pl_aud_pasta = tk.StringVar()
        self.var_pl_aud_saida = tk.StringVar(value=os.path.join(_RES, "auditoria_avulsa.csv"))

        linha = self._pl_linha(f)
        ttk.Label(linha, text="Tipo:").pack(side="left", padx=(0, 4))
        cb = ttk.Combobox(linha, textvariable=self.var_pl_aud_tipo, values=list(AUDITORIAS),
                          state="readonly", width=28)
        cb.pack(side="left")
        cb.bind("<<ComboboxSelected>>", lambda e: self._pl_aud_monta_params())
        self.var_pl_aud_dica = tk.StringVar()
        ttk.Label(f, textvariable=self.var_pl_aud_dica, bootstyle="secondary",
                  wraplength=1100, justify="left").pack(anchor="w", pady=(2, 4))

        self._pl_campo(self._pl_linha(f), "CSV de entrada:", self.var_pl_aud_csv, 70, "arquivo")
        self._pl_campo(self._pl_linha(f), "Pasta com os .ns2:", self.var_pl_aud_pasta, 70, "pasta")
        self._pl_campo(self._pl_linha(f), "CSV de saída:", self.var_pl_aud_saida, 70, "salvar")
        self.frame_pl_aud_params = ttk.Frame(f)
        self.frame_pl_aud_params.pack(fill="x", pady=(4, 0))
        self._pl_aud_monta_params()

        painel = self._pl_painel(aba)
        ttk.Button(f, text="Rodar auditoria", bootstyle="primary",
                   command=lambda: self._pl_aud_rodar(painel)).pack(anchor="w", pady=(8, 0))
        painel.pack(fill="both", expand=True, padx=4, pady=(10, 4))

    def _pl_aud_monta_params(self):
        for w in self.frame_pl_aud_params.winfo_children():
            w.destroy()
        spec = AUDITORIAS[self.var_pl_aud_tipo.get()]
        self.var_pl_aud_dica.set(spec["dica"])
        self._pl_aud_vars = {}
        linha = self._pl_linha(self.frame_pl_aud_params)
        for flag, rotulo, default in spec["params"]:
            var = tk.StringVar(value=default)
            self._pl_aud_vars[flag] = var
            self._pl_campo(linha, f"{rotulo}:", var, 6)
        if "combo" in spec:
            flag, rotulo, opcoes = spec["combo"]
            var = tk.StringVar(value=opcoes[0])
            self._pl_aud_vars[flag] = var
            ttk.Label(linha, text=f"{rotulo}:").pack(side="left", padx=(0, 4))
            ttk.Combobox(linha, textvariable=var, values=opcoes, state="readonly",
                         width=10).pack(side="left")

    def _pl_aud_rodar(self, painel):
        csv_in, pasta = self.var_pl_aud_csv.get(), self.var_pl_aud_pasta.get()
        saida = self.var_pl_aud_saida.get()
        if not self._pl_valida([("CSV de entrada", csv_in, "arquivo"),
                                ("Pasta com os .ns2", pasta, "pasta"),
                                ("CSV de saída", saida, "texto")]):
            return
        if not self._pl_confirma_sobrescrita(saida):
            return
        spec = AUDITORIAS[self.var_pl_aud_tipo.get()]
        args = ["--csv", csv_in, "--pasta_ns2", pasta, "--saida", saida]
        for flag, var in self._pl_aud_vars.items():
            if var.get().strip():
                args += [flag, var.get().strip()]
        self._pl_rodar(painel, spec["chave"], args)

    # ------------------------------------------------------------------
    # 6. Agregação -- agrega → enriquece → junta comportamento → consolida
    # ------------------------------------------------------------------
    def _pl_aba_agregacao(self, aba):
        self._callout(
            aba,
            "Saídas padrão gravam em *_novo.csv (não sobrescrevem produção). O agregador varre "
            "a pasta recursivamente -- a raiz de RESULTADOS_MESTRADO inclui _CONTAMINADO_nao_usar/.",
            bootstyle="warning",
        )
        # 4 formulários empilhados não cabem numa tela 1366x768 com o log embaixo
        # (o botão Consolidar sumia): passos à esquerda com rolagem, log à
        # direita -- num notebook widescreen o que falta é altura, não largura.
        colunas = ttk.Frame(aba)
        colunas.pack(fill="both", expand=True)
        rolagem = ttk.ScrolledFrame(colunas, auto_hide=True)
        rolagem.pack(side="left", fill="both", expand=True)
        direita = ttk.Frame(colunas, width=470)
        direita.pack(side="left", fill="y", padx=(10, 0))
        direita.pack_propagate(False)
        mestre_novo = os.path.join(_RES, "dataset_mestre_final_novo.csv")
        enriq_novo = os.path.join(_RES, "dataset_mestre_enriquecido_novo.csv")
        comp_novo = os.path.join(_RES, "dataset_mestre_COM_COMPORTAMENTO_novo.csv")

        self.var_pl_ag_res = tk.StringVar(value=os.path.join(workspace.BASE_RESULTADOS_MESTRADO,
                                                              "basal", "NOCI"))
        self.var_pl_ag_saida = tk.StringVar(value=mestre_novo)
        self.var_pl_en_ent = tk.StringVar(value=mestre_novo)
        self.var_pl_en_saida = tk.StringVar(value=enriq_novo)
        self.var_pl_en_fooof = tk.BooleanVar(value=True)
        self.var_pl_en_portao = tk.BooleanVar(value=True)
        self.var_pl_en_forca = tk.BooleanVar(value=False)
        self.var_pl_en_workers = tk.StringVar(value="4")
        self.var_pl_ju_mestre = tk.StringVar(value=enriq_novo)
        self.var_pl_ju_comp = tk.StringVar(value=_TEMPLATE)
        self.var_pl_ju_saida = tk.StringVar(value=comp_novo)
        self.var_pl_co_ent = tk.StringVar(value=comp_novo)
        self.var_pl_co_saida = tk.StringVar(value=os.path.join(_RES, "candidatos_vencedores_consolidados_novo.csv"))
        self.var_pl_co_erro = tk.StringVar(value="0.15")
        self.var_pl_co_knee = tk.BooleanVar(value=False)
        self.var_pl_co_harm = tk.BooleanVar(value=False)

        painel = self._pl_painel(direita, altura_log=6)
        painel.txt_log.configure(width=1)  # largura vem da coluna, não de 80 colunas de texto

        f = self._secao(rolagem, "1. Agregar resultados por canal")
        self._pl_campo(self._pl_linha(f), "Pasta de resultados:", self.var_pl_ag_res, 60, "pasta")
        linha = self._pl_linha(f)
        self._pl_campo(linha, "Saída:", self.var_pl_ag_saida, 60, "salvar")
        ttk.Button(linha, text="Agregar", bootstyle="primary",
                   command=lambda: self._pl_ag_agregar(painel)).pack(side="left")

        f = self._secao(rolagem, "2. Enriquecer (FOOOF + portão de banda larga)")
        self._pl_campo(self._pl_linha(f), "Entrada:", self.var_pl_en_ent, 60, "arquivo")
        self._pl_campo(self._pl_linha(f), "Saída:", self.var_pl_en_saida, 60, "salvar")
        linha = self._pl_linha(f)
        ttk.Checkbutton(linha, text="FOOOF", variable=self.var_pl_en_fooof).pack(side="left", padx=(0, 10))
        ttk.Checkbutton(linha, text="Portão", variable=self.var_pl_en_portao).pack(side="left", padx=(0, 10))
        ttk.Checkbutton(linha, text="Forçar recálculo do FOOOF",
                        variable=self.var_pl_en_forca).pack(side="left", padx=(0, 16))
        self._pl_campo(linha, "Processos em paralelo:", self.var_pl_en_workers, 4)
        ttk.Button(linha, text="Enriquecer", bootstyle="primary",
                   command=lambda: self._pl_ag_enriquecer(painel)).pack(side="left")

        f = self._secao(rolagem, "3. Juntar anotações de comportamento")
        self._pl_campo(self._pl_linha(f), "Dataset mestre:", self.var_pl_ju_mestre, 60, "arquivo")
        self._pl_campo(self._pl_linha(f), "Template anotado:", self.var_pl_ju_comp, 60, "arquivo")
        linha = self._pl_linha(f)
        self._pl_campo(linha, "Saída:", self.var_pl_ju_saida, 60, "salvar")
        ttk.Button(linha, text="Juntar", bootstyle="primary",
                   command=lambda: self._pl_ag_juntar(painel)).pack(side="left")

        f = self._secao(rolagem, "4. Consolidar vencedores")
        self._pl_campo(self._pl_linha(f), "Entrada:", self.var_pl_co_ent, 60, "arquivo")
        self._pl_campo(self._pl_linha(f), "Saída:", self.var_pl_co_saida, 60, "salvar")
        linha = self._pl_linha(f)
        self._pl_campo(linha, "Erro máximo do FOOOF:", self.var_pl_co_erro, 6)
        ttk.Checkbutton(linha, text="Exigir knee válido (filtro antigo)",
                        variable=self.var_pl_co_knee).pack(side="left", padx=(0, 10))
        ttk.Checkbutton(linha, text="Exigir harmônico CLEAN (filtro antigo)",
                        variable=self.var_pl_co_harm).pack(side="left", padx=(0, 16))
        ttk.Button(linha, text="Consolidar", bootstyle="primary",
                   command=lambda: self._pl_ag_consolidar(painel)).pack(side="left")

        ttk.Label(direita, text="LOG", style="SecaoTitulo.TLabel").pack(anchor="w", pady=(12, 2))
        painel.pack(fill="both", expand=True)

    def _pl_ag_agregar(self, painel):
        res, saida = self.var_pl_ag_res.get(), self.var_pl_ag_saida.get()
        if not self._pl_valida([("Pasta de resultados", res, "pasta"), ("Saída", saida, "texto")]):
            return
        if not self._pl_confirma_sobrescrita(saida):
            return
        self._pl_rodar(painel, "agrega", ["--resultados", res, "--saida", saida])

    def _pl_ag_enriquecer(self, painel):
        ent, saida = self.var_pl_en_ent.get(), self.var_pl_en_saida.get()
        etapas = [nome for nome, v in (("fooof", self.var_pl_en_fooof), ("portao", self.var_pl_en_portao))
                  if v.get()]
        if not self._pl_valida([("Entrada", ent, "arquivo"), ("Saída", saida, "texto")]):
            return
        if not etapas:
            messagebox.showerror("Campo obrigatório", "Marque ao menos uma etapa (FOOOF ou portão).")
            return
        if not self._pl_confirma_sobrescrita(saida):
            return
        args = ["--entrada", ent, "--saida", saida, "--etapas", *etapas,
                "--n_workers", self.var_pl_en_workers.get()]
        if self.var_pl_en_forca.get():
            args.append("--forca")
        self._pl_rodar(painel, "enriquece", args)

    def _pl_ag_juntar(self, painel):
        mestre, comp, saida = (self.var_pl_ju_mestre.get(), self.var_pl_ju_comp.get(),
                               self.var_pl_ju_saida.get())
        if not self._pl_valida([("Dataset mestre", mestre, "arquivo"),
                                ("Template anotado", comp, "arquivo"), ("Saída", saida, "texto")]):
            return
        if not self._pl_confirma_sobrescrita(saida):
            return
        self._pl_rodar(painel, "junta", ["--mestre", mestre, "--comportamento", comp, "--saida", saida])

    def _pl_ag_consolidar(self, painel):
        ent, saida = self.var_pl_co_ent.get(), self.var_pl_co_saida.get()
        if not self._pl_valida([("Entrada", ent, "arquivo"), ("Saída", saida, "texto")]):
            return
        if not self._pl_confirma_sobrescrita(saida):
            return
        args = ["--entrada", ent, "--saida", saida, "--erro_fooof_max", self.var_pl_co_erro.get()]
        if self.var_pl_co_knee.get():
            args.append("--exigir_knee_valido")
        if self.var_pl_co_harm.get():
            args.append("--exigir_harmonico_clean")
        self._pl_rodar(painel, "consolida", args)
