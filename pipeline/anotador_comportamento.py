#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Anotador Comportamental Sincronizado com Vídeo
------------------------------------------------
Ferramenta para acelerar o processo de anotação comportamental a partir de
janelas temporais pré-definidas (ex: 10 segundos) no arquivo template_comportamento.csv.

Recursos:
- Seleção de Sessão e Arquivo Neural (.ns2)
- Associação com arquivo de vídeo (MP4, AVI, MKV, etc.) e memória automática de caminhos
- Sincronização temporal flexível via Offset em segundos (com ajuste fino)
- Salto automático para o momento exato do vídeo com reprodução contínua ou em loop
- Botões de atalho rápido para comportamentos comuns (ou digitação livre)
- Fluxo de trabalho ágil por teclado: Digite/escolha -> Enter -> Salva e pula para o próximo momento
- Backup de segurança automático ao iniciar e salvamento direto no CSV
"""

import os
import sys
import time
import json
import shutil
import datetime
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import pandas as pd
from PIL import Image, ImageTk
import cv2

def _encontrar_csv_padrao():
    candidatos = [
        os.path.join(os.getcwd(), "template_comportamento.csv"),
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "template_comportamento.csv"),
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "template_comportamento.csv"),
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "template_comportamento.csv"),
    ]
    for c in candidatos:
        if os.path.exists(c):
            return os.path.abspath(c)
    return os.path.abspath("template_comportamento.csv")

def _encontrar_config_padrao():
    candidatos = [
        os.path.join(os.getcwd(), "config_anotador.json"),
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "config_anotador.json"),
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "config_anotador.json"),
    ]
    for c in candidatos:
        if os.path.exists(c):
            return os.path.abspath(c)
    csv_dir = os.path.dirname(_encontrar_csv_padrao())
    return os.path.join(csv_dir, "config_anotador.json")

DEFAULT_CSV = _encontrar_csv_padrao()
CONFIG_FILE = _encontrar_config_padrao()

# Comportamentos sugeridos com atalhos rápidos
COMPORTAMENTOS_PADRAO = [
    ("1", "Imóvel / Descanso"),
    ("2", "Exploração / Locomoção"),
    ("3", "Grooming / Limpeza"),
    ("4", "Rearing / Em pé"),
    ("5", "Sniffing / Farejando"),
    ("6", "Movimento de Cabeça"),
    ("7", "Artefato / Cabo"),
    ("8", "Sono"),
]


def format_time(seconds: float) -> str:
    """Converte segundos para formato MM:SS ou MM:SS.s"""
    if pd.isna(seconds) or seconds is None:
        return "00:00"
    if seconds < 0:
        return f"-{format_time(-seconds)}"
    m = int(seconds // 60)
    s = seconds % 60
    if abs(s - round(s)) > 0.05:
        return f"{m:02d}:{s:04.1f}"
    return f"{m:02d}:{int(round(s)):02d}"


def reparar_mojibake(texto):
    """Corrige strings que sofreram dupla decodificação (ex: UTF-8 lido como Latin-1 no Excel)."""
    if not isinstance(texto, str) or not texto:
        return texto
    try:
        if any(bad in texto for bad in ["Ã§", "Ã£", "Ã³", "Ã©", "Ãª", "Ã¡", "Ã", "Â"]):
            return texto.encode("latin1").decode("utf-8")
    except Exception:
        pass
    return texto


class VideoPlayer:
    """Gerencia a captura e controle de vídeo com OpenCV."""
    def __init__(self):
        self.cap = None
        self.video_path = None
        self.total_frames = 0
        self.fps = 30.0
        self.duration_s = 0.0
        self.is_playing = False
        self.loop_window = True
        self.window_start_s = 0.0
        self.window_end_s = 10.0
        self.current_time_s = 0.0
        self.playback_speed = 1.0
        self.last_frame_pil = None

    def load_video(self, path: str) -> bool:
        if self.cap is not None:
            self.cap.release()
            self.cap = None

        if not os.path.exists(path):
            return False

        cap = cv2.VideoCapture(path)
        if not cap.isOpened():
            return False

        self.cap = cap
        self.video_path = path
        self.fps = cap.get(cv2.CAP_PROP_FPS) or 30.0
        if self.fps <= 0 or self.fps > 240:
            self.fps = 30.0
        self.total_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
        self.duration_s = self.total_frames / self.fps if self.fps > 0 else 0.0
        self.is_playing = False
        self.seek(0.0)
        return True

    def release(self):
        if self.cap is not None:
            self.cap.release()
            self.cap = None

    def seek(self, seconds: float):
        if self.cap is None:
            return
        clamped_s = max(0.0, min(seconds, self.duration_s))
        self.current_time_s = clamped_s
        self.cap.set(cv2.CAP_PROP_POS_MSEC, clamped_s * 1000.0)

    def get_current_frame(self):
        if self.cap is None:
            return None
        ret, frame = self.cap.read()
        if not ret:
            return None
        self.current_time_s = self.cap.get(cv2.CAP_PROP_POS_MSEC) / 1000.0
        frame_rgb = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
        self.last_frame_pil = Image.fromarray(frame_rgb)
        return self.last_frame_pil


class AnotadorApp(tk.Tk):
    def __init__(self, csv_path: str = DEFAULT_CSV):
        super().__init__()
        self.title("Anotador Comportamental - Vídeo & LFP")
        self.geometry("1300x820")
        self.minsize(1050, 700)

        # Configura tema e estilo
        self.style = ttk.Style(self)
        try:
            self.style.theme_use("clam")
        except Exception:
            pass

        self._configure_styles()

        self.csv_path = csv_path
        self.df = None
        self.config = self._load_config()

        self.player = VideoPlayer()
        self.current_session = None
        self.current_file = None
        self.session_files = []
        self.file_offsets = {}
        self.current_row_idx = None
        self.filtered_indices = []
        self._is_updating_slider = False

        # Variáveis de controle
        self.var_session = tk.StringVar()
        self.var_current_file = tk.StringVar(value="")
        self.var_video_path = tk.StringVar()
        self.var_offset = tk.DoubleVar(value=0.0)
        self.var_loop_window = tk.BooleanVar(value=True)
        self.var_autoplay = tk.BooleanVar(value=True)
        self.var_speed = tk.StringVar(value="1.0x")
        self.var_comportamento = tk.StringVar()
        self.var_observacoes = tk.StringVar()
        self.var_status_filter = tk.StringVar(value="Todos")
        self.var_current_time = tk.StringVar(value="00:00 / 00:00")
        self.var_progress_info = tk.StringVar(value="")
        self.var_status_save = tk.StringVar(value="Pronto para anotar")

        # Inicializa dados e interface
        self._carregar_csv()
        self._criar_backup_inicial()
        self._construir_interface()
        self._vincular_atalhos()

        # Inicia loop de atualização do vídeo
        self.after(30, self._video_loop)

        # Popula sessões
        self._atualizar_lista_sessoes()

    def _configure_styles(self):
        """Define cores e fontes com visual moderno e limpo."""
        self.configure(bg="#F4F6F9")
        self.style.configure(".", background="#F4F6F9", font=("Segoe UI", 9))
        self.style.configure("TLabel", background="#F4F6F9", foreground="#212529")
        self.style.configure("Header.TLabel", font=("Segoe UI", 11, "bold"), foreground="#1A202C")
        self.style.configure("TButton", font=("Segoe UI", 9), padding=5)
        self.style.configure("Primary.TButton", font=("Segoe UI", 9, "bold"), foreground="#FFFFFF", background="#2563EB")
        self.style.configure("Success.TButton", font=("Segoe UI", 10, "bold"), foreground="#FFFFFF", background="#059669")
        self.style.configure("Action.TButton", font=("Segoe UI", 9), padding=3)
        self.style.configure("Treeview", font=("Segoe UI", 9), rowheight=24)
        self.style.configure("Treeview.Heading", font=("Segoe UI", 9, "bold"))
        self.style.map("Treeview", background=[("selected", "#3B82F6")], foreground=[("selected", "#FFFFFF")])

    def _load_config(self) -> dict:
        if os.path.exists(CONFIG_FILE):
            try:
                with open(CONFIG_FILE, "r", encoding="utf-8") as f:
                    return json.load(f)
            except Exception:
                return {}
        return {}

    def _save_config(self):
        try:
            with open(CONFIG_FILE, "w", encoding="utf-8") as f:
                json.dump(self.config, f, indent=2, ensure_ascii=False)
        except Exception as e:
            print(f"Erro ao salvar config: {e}")

    def _carregar_csv(self):
        if not os.path.exists(self.csv_path):
            messagebox.showerror("Erro", f"Arquivo CSV não encontrado em:\n{self.csv_path}")
            sys.exit(1)

        # Detecta separador (, ou ;) e prioriza utf-8-sig (compatível com Excel)
        carregado = False
        for enc in ["utf-8-sig", "utf-8", "cp1252", "latin1"]:
            for sep in [",", ";", None]:
                try:
                    df_temp = pd.read_csv(self.csv_path, sep=sep, engine="python" if sep is None else "c", encoding=enc)
                    if "sessao" in df_temp.columns or "janela_ini_s" in df_temp.columns:
                        self.df = df_temp
                        carregado = True
                        break
                except Exception:
                    continue
            if carregado:
                break

        if not carregado or self.df is None:
            self.df = pd.read_csv(self.csv_path, sep=None, engine="python", encoding_errors="replace")

        # Garante colunas esperadas e repara eventuais caracteres bugados por Excel
        for col in ["comportamento", "observacoes"]:
            if col not in self.df.columns:
                self.df[col] = ""
        self.df["comportamento"] = self.df["comportamento"].fillna("").astype(str).apply(reparar_mojibake)
        self.df["observacoes"] = self.df["observacoes"].fillna("").astype(str).apply(reparar_mojibake)

    def _criar_backup_inicial(self):
        """Cria um backup timestamped para segurança dos dados originais."""
        try:
            ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            backup_dir = os.path.join(os.path.dirname(self.csv_path), "backups_comportamento")
            os.makedirs(backup_dir, exist_ok=True)
            backup_path = os.path.join(backup_dir, f"backup_template_{ts}.csv")
            shutil.copy2(self.csv_path, backup_path)
        except Exception as e:
            print(f"Aviso de backup: {e}")

    def _salvar_csv(self, mostrar_mensagem=True):
        """Salva as anotações no CSV com utf-8-sig (com BOM), para visualização 100% correta no Excel."""
        if self.df is None:
            return
        try:
            # Garante que textos não tenham resquícios de mojibake
            self.df["comportamento"] = self.df["comportamento"].fillna("").astype(str).apply(reparar_mojibake)
            self.df["observacoes"] = self.df["observacoes"].fillna("").astype(str).apply(reparar_mojibake)

            tmp_file = self.csv_path + ".tmp"
            # Salva sempre com utf-8-sig (UTF-8 com BOM) para que o Excel no Windows abra perfeitamente com acentos
            self.df.to_csv(tmp_file, index=False, sep=",", encoding="utf-8-sig")
            if os.path.exists(self.csv_path):
                os.replace(tmp_file, self.csv_path)
            else:
                shutil.move(tmp_file, self.csv_path)

            if mostrar_mensagem:
                messagebox.showinfo("Salvo", "Planilha salva com sucesso em formato UTF-8 (compatível com Excel)!")
        except Exception as e:
            messagebox.showerror("Erro ao Salvar", f"Não foi possível salvar o CSV:\n{e}")

    def _construir_interface(self):
        # 1. BARRA SUPERIOR: Sessão, Vídeo e Offset
        top_frame = ttk.LabelFrame(self, text=" 1. Sessão, Vídeo e Sincronização ", padding=10)
        top_frame.pack(fill="x", padx=10, pady=(8, 4))

        # Linha 1: Sessão
        f_sess = ttk.Frame(top_frame)
        f_sess.pack(fill="x", pady=2)
        ttk.Label(f_sess, text="Sessão:", font=("Segoe UI", 9, "bold"), width=9).pack(side="left")
        self.cb_sessao = ttk.Combobox(f_sess, textvariable=self.var_session, state="readonly", width=65)
        self.cb_sessao.pack(side="left", padx=5)
        self.cb_sessao.bind("<<ComboboxSelected>>", self._ao_selecionar_sessao)

        btn_prox_sess = ttk.Button(f_sess, text="Próxima Sessão ⏭", command=self._avancar_sessao)
        btn_prox_sess.pack(side="left", padx=5)

        self.lbl_progress = ttk.Label(f_sess, textvariable=self.var_progress_info, font=("Segoe UI", 9, "bold"), foreground="#2563EB")
        self.lbl_progress.pack(side="right", padx=10)

        # Linha 2: Arquivo de Vídeo
        f_vid = ttk.Frame(top_frame)
        f_vid.pack(fill="x", pady=2)
        ttk.Label(f_vid, text="Vídeo:", font=("Segoe UI", 9, "bold"), width=9).pack(side="left")
        self.entry_video = ttk.Entry(f_vid, textvariable=self.var_video_path, width=67)
        self.entry_video.pack(side="left", padx=5)
        btn_browse = ttk.Button(f_vid, text="Selecionar Vídeo...", command=self._selecionar_video)
        btn_browse.pack(side="left", padx=5)

        # Linha 3: Arquivo da Sessão & Offset
        f_off = ttk.Frame(top_frame)
        f_off.pack(fill="x", pady=2)

        ttk.Label(f_off, text="Arquivo:", font=("Segoe UI", 9, "bold"), width=9).pack(side="left")
        self.cb_arquivo = ttk.Combobox(f_off, textvariable=self.var_current_file, state="readonly", width=34)
        self.cb_arquivo.pack(side="left", padx=5)
        self.cb_arquivo.bind("<<ComboboxSelected>>", self._ao_selecionar_arquivo_dropdown)

        ttk.Label(f_off, text="Offset (s):", font=("Segoe UI", 9, "bold")).pack(side="left", padx=(8, 2))
        self.spin_offset = ttk.Spinbox(f_off, from_=-3600.0, to=7200.0, increment=0.5,
                                       textvariable=self.var_offset, width=9, command=self._ao_alterar_offset)
        self.spin_offset.pack(side="left", padx=5)
        self.spin_offset.bind("<Return>", lambda e: self._ao_alterar_offset())

        self.lbl_offset_time = ttk.Label(f_off, text="(00:00)", font=("Segoe UI", 9, "bold"), foreground="#2563EB", width=9)
        self.lbl_offset_time.pack(side="left", padx=2)

        for d, txt in [(-5.0, "-5s"), (-1.0, "-1s"), (+1.0, "+1s"), (+5.0, "+5s")]:
            btn = ttk.Button(f_off, text=txt, width=4, command=lambda delta=d: self._ajustar_offset(delta))
            btn.pack(side="left", padx=1)

        ttk.Button(f_off, text="Fixar Frame Atual no Início deste Arquivo",
                   command=self._definir_offset_pelo_frame).pack(side="left", padx=8)

        ttk.Label(f_off, text="(Vídeo = LFP + Offset do Arquivo)",
                  font=("Segoe UI", 8, "italic"), foreground="#64748B").pack(side="left", padx=3)

        # 2. ÁREA CENTRAL DIVIDIDA (Lista de Momentos na esquerda + Vídeo na direita)
        main_paned = ttk.PanedWindow(self, orient="horizontal")
        main_paned.pack(fill="both", expand=True, padx=10, pady=4)

        # 2.1 Painel Esquerdo: Lista de Momentos
        left_frame = ttk.LabelFrame(main_paned, text=" 2. Momentos da Sessão (Janelas PAC) ", padding=6)
        main_paned.add(left_frame, weight=1)

        # Filtro de status
        f_filt = ttk.Frame(left_frame)
        f_filt.pack(fill="x", pady=(0, 4))
        ttk.Label(f_filt, text="Exibir:").pack(side="left", padx=2)
        cb_filt = ttk.Combobox(f_filt, textvariable=self.var_status_filter,
                               values=["Todos", "Apenas Pendentes", "Apenas Concluídos"],
                               state="readonly", width=18)
        cb_filt.pack(side="left", padx=4)
        cb_filt.bind("<<ComboboxSelected>>", lambda e: self._atualizar_tabela_momentos())

        # Tabela Treeview
        cols = ("status", "arquivo", "janela_lfp", "tempo_vid", "pares", "comportamento")
        self.tree = ttk.Treeview(left_frame, columns=cols, show="headings", selectmode="browse")
        self.tree.heading("status", text="St")
        self.tree.heading("arquivo", text="Arquivo NS2")
        self.tree.heading("janela_lfp", text="LFP (s)")
        self.tree.heading("tempo_vid", text="Vídeo")
        self.tree.heading("pares", text="Pares PAC")
        self.tree.heading("comportamento", text="Comportamento")

        self.tree.column("status", width=35, anchor="center")
        self.tree.column("arquivo", width=140, anchor="w")
        self.tree.column("janela_lfp", width=85, anchor="center")
        self.tree.column("tempo_vid", width=95, anchor="center")
        self.tree.column("pares", width=110, anchor="w")
        self.tree.column("comportamento", width=120, anchor="w")

        tree_scroll = ttk.Scrollbar(left_frame, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=tree_scroll.set)
        self.tree.pack(side="left", fill="both", expand=True)
        tree_scroll.pack(side="right", fill="y")

        self.tree.bind("<<TreeviewSelect>>", self._ao_selecionar_momento_tabela)

        # 2.2 Painel Direito: Vídeo Player
        right_frame = ttk.LabelFrame(main_paned, text=" 3. Visualização do Vídeo ", padding=6)
        main_paned.add(right_frame, weight=2)

        # Canvas para exibição dos frames
        self.canvas_video = tk.Canvas(right_frame, bg="#111827", highlightthickness=0)
        self.canvas_video.pack(fill="both", expand=True)
        self.canvas_video.bind("<Configure>", lambda e: self._redesenhar_frame_atual())

        # Barra de Controles do Vídeo
        f_vid_ctrl = ttk.Frame(right_frame)
        f_vid_ctrl.pack(fill="x", pady=(6, 2))

        # Slider de tempo
        self.slider_time = ttk.Scale(f_vid_ctrl, from_=0.0, to=100.0, orient="horizontal", command=self._ao_arrastar_slider)
        self.slider_time.pack(fill="x", expand=True, side="left", padx=5)

        self.lbl_time = ttk.Label(f_vid_ctrl, textvariable=self.var_current_time, font=("Segoe UI", 9, "bold"), width=16)
        self.lbl_time.pack(side="right", padx=5)

        # Botões de Controle de Playback
        f_btn_ctrl = ttk.Frame(right_frame)
        f_btn_ctrl.pack(fill="x", pady=2)

        self.btn_play = ttk.Button(f_btn_ctrl, text="▶ Play (Espaço)", width=14, command=self._toggle_play)
        self.btn_play.pack(side="left", padx=3)

        ttk.Button(f_btn_ctrl, text="↺ Repetir Janela (R)", width=17, command=self._repetir_janela).pack(side="left", padx=3)
        ttk.Button(f_btn_ctrl, text="◀ -2s", width=6, command=lambda: self._avancar_tempo(-2.0)).pack(side="left", padx=2)
        ttk.Button(f_btn_ctrl, text="+2s ▶", width=6, command=lambda: self._avancar_tempo(+2.0)).pack(side="left", padx=2)

        # Opções de Loop e Velocidade
        chk_loop = ttk.Checkbutton(f_btn_ctrl, text="Loop na Janela (10s)", variable=self.var_loop_window,
                                   command=self._ao_alternar_loop)
        chk_loop.pack(side="left", padx=10)

        chk_auto = ttk.Checkbutton(f_btn_ctrl, text="Auto-Play ao Avançar", variable=self.var_autoplay)
        chk_auto.pack(side="left", padx=5)

        ttk.Label(f_btn_ctrl, text="Velocidade:").pack(side="left", padx=(10, 2))
        cb_speed = ttk.Combobox(f_btn_ctrl, textvariable=self.var_speed, values=["0.5x", "1.0x", "1.5x", "2.0x"],
                                state="readonly", width=6)
        cb_speed.pack(side="left", padx=2)
        cb_speed.bind("<<ComboboxSelected>>", self._ao_mudar_velocidade)

        # 3. PAINEL INFERIOR: Anotação Rápida e Produtiva
        bot_frame = ttk.LabelFrame(self, text=" 4. Anotação Rápida do Comportamento ", padding=8)
        bot_frame.pack(fill="x", padx=10, pady=(4, 8))

        # Linha de botões de atalho 1 toque
        f_quick = ttk.Frame(bot_frame)
        f_quick.pack(fill="x", pady=(0, 6))
        ttk.Label(f_quick, text="Atalhos Rápidos:", font=("Segoe UI", 9, "bold")).pack(side="left", padx=(0, 6))

        for tecla, nome in COMPORTAMENTOS_PADRAO:
            btn = ttk.Button(f_quick, text=f"[{tecla}] {nome}",
                             command=lambda comp=nome: self._aplicar_comportamento_rapido(comp))
            btn.pack(side="left", padx=2)

        # Linha de digitação e salvamento
        f_input = ttk.Frame(bot_frame)
        f_input.pack(fill="x", pady=2)

        ttk.Label(f_input, text="Comportamento:", font=("Segoe UI", 9, "bold")).pack(side="left", padx=2)
        self.cb_comportamento = ttk.Combobox(f_input, textvariable=self.var_comportamento, width=28)
        self.cb_comportamento.pack(side="left", padx=5)
        self._atualizar_historico_comportamentos()
        self.cb_comportamento.bind("<Return>", lambda e: self._salvar_e_proximo())

        ttk.Label(f_input, text="Observações:", font=("Segoe UI", 9)).pack(side="left", padx=(15, 2))
        self.entry_obs = ttk.Entry(f_input, textvariable=self.var_observacoes, width=32)
        self.entry_obs.pack(side="left", padx=5)
        self.entry_obs.bind("<Return>", lambda e: self._salvar_e_proximo())

        btn_salvar_prox = ttk.Button(f_input, text="✔ Salvar e Próximo (Enter)", style="Success.TButton",
                                     command=self._salvar_e_proximo)
        btn_salvar_prox.pack(side="left", padx=12)

        self.lbl_status_salvo = ttk.Label(f_input, textvariable=self.var_status_save,
                                          font=("Segoe UI", 9, "italic"), foreground="#059669")
        self.lbl_status_salvo.pack(side="left", padx=10)

        ttk.Button(f_input, text="💾 Salvar Planilha (Ctrl+S)", command=lambda: self._salvar_csv(True)).pack(side="right", padx=5)
        ttk.Button(f_input, text="➡ Pular Momento", command=self._avancar_momento).pack(side="right", padx=3)
        ttk.Button(f_input, text="⬅ Momento Anterior", command=self._voltar_momento).pack(side="right", padx=3)

    def _vincular_atalhos(self):
        """Atalhos globais de teclado para máxima agilidade."""
        self.bind("<Control-s>", lambda e: self._salvar_csv(True))
        self.bind("<Control-S>", lambda e: self._salvar_csv(True))
        self.bind("<space>", self._ao_pressionar_espaco)
        self.bind("<r>", lambda e: self._repetir_janela())
        self.bind("<R>", lambda e: self._repetir_janela())
        self.bind("<F5>", lambda e: self._repetir_janela())

        for tecla, nome in COMPORTAMENTOS_PADRAO:
            self.bind(f"<Alt-Key-{tecla}>", lambda e, comp=nome: self._aplicar_comportamento_rapido(comp))

    def _ao_pressionar_espaco(self, event):
        if isinstance(event.widget, (ttk.Entry, tk.Entry, ttk.Combobox)):
            return
        self._toggle_play()

    def _atualizar_historico_comportamentos(self):
        """Coleta comportamentos já existentes na planilha para autocompletar."""
        if self.df is not None:
            valores = sorted(list(set(
                [c.strip() for c in self.df["comportamento"].dropna().unique() if c and c.strip()]
                + [nome for _, nome in COMPORTAMENTOS_PADRAO]
            )))
            self.cb_comportamento["values"] = valores

    def _atualizar_lista_sessoes(self):
        """Carrega a lista de sessões únicas do CSV."""
        if self.df is None or len(self.df) == 0:
            return

        sessoes = self.df["sessao"].unique().tolist()
        sessoes_display = []
        for s in sessoes:
            sub = self.df[self.df["sessao"] == s]
            total = len(sub)
            anotados = (sub["comportamento"].str.strip() != "").sum()
            sessoes_display.append(f"{s} ({anotados}/{total})")

        self.cb_sessao["values"] = sessoes_display

        ultima_sessao = self.config.get("ultima_sessao")
        idx_escolhido = 0
        if ultima_sessao:
            for i, s in enumerate(sessoes):
                if s == ultima_sessao:
                    idx_escolhido = i
                    break

        if sessoes:
            self.cb_sessao.current(idx_escolhido)
            self._ao_selecionar_sessao()

    def _get_offset_for_file(self, arq_nome: str) -> float:
        if arq_nome in self.file_offsets:
            return float(self.file_offsets[arq_nome])
        try:
            return float(self.var_offset.get())
        except Exception:
            return 0.0

    def _ao_selecionar_sessao(self, event=None):
        raw_val = self.var_session.get()
        if not raw_val:
            return
        sessao_nome = raw_val.rsplit(" (", 1)[0]
        self.current_session = sessao_nome
        self.config["ultima_sessao"] = sessao_nome

        # Identifica todos os arquivos únicos pertencentes a esta sessão
        sub_df = self.df[self.df["sessao"] == sessao_nome]
        self.session_files = list(dict.fromkeys(sub_df["arquivo"].dropna().astype(str).tolist()))

        sess_cfg = self.config.get("sessoes", {}).get(sessao_nome, {})
        video_salvo = sess_cfg.get("video_path", "")
        base_offset = float(sess_cfg.get("offset", 0.0))
        saved_file_offsets = sess_cfg.get("offsets_arquivos", {})

        # Popula offsets para cada arquivo da sessão
        self.file_offsets = {}
        for i, f in enumerate(self.session_files):
            if f in saved_file_offsets:
                self.file_offsets[f] = float(saved_file_offsets[f])
            else:
                # Cada arquivo de registro neural subsequente tem offset padrão cumulativo de 5 min (300s)
                self.file_offsets[f] = round(base_offset + (i * 300.0), 2)

        # Atualiza combobox de arquivos
        self.cb_arquivo["values"] = self.session_files
        if self.session_files:
            self.current_file = self.session_files[0]
            self.var_current_file.set(self.current_file)
            off = self.file_offsets[self.current_file]
            self.var_offset.set(off)
            self.lbl_offset_time.configure(text=f"({format_time(off)})")
        else:
            self.current_file = None
            self.var_current_file.set("")
            self.var_offset.set(base_offset)
            self.lbl_offset_time.configure(text=f"({format_time(base_offset)})")

        if video_salvo and os.path.exists(video_salvo):
            self.var_video_path.set(video_salvo)
            self._carregar_video_arquivo(video_salvo)
        else:
            self.player.release()
            self.var_video_path.set("")
            self._tentar_autodetectar_video(sessao_nome)

        self._sincronizar_tempos_video_sessao()
        self._atualizar_tabela_momentos()

    def _ao_selecionar_arquivo_dropdown(self, event=None):
        arq = self.var_current_file.get()
        if not arq:
            return
        self.current_file = arq
        off = self._get_offset_for_file(arq)
        self.var_offset.set(off)
        self.lbl_offset_time.configure(text=f"({format_time(off)})")

        # Seleciona o primeiro momento deste arquivo na tabela
        for item in self.tree.get_children():
            idx = int(self.tree.item(item, "text"))
            if self.df is not None and idx in self.df.index:
                if str(self.df.at[idx, "arquivo"]) == arq:
                    self.tree.selection_set(item)
                    self.tree.see(item)
                    self._ao_selecionar_momento_tabela()
                    break

    def _sincronizar_tempos_video_sessao(self):
        """Garante que video_tempo_ini e video_tempo_fim reflitam os offsets de cada arquivo da sessão."""
        if self.df is None or not self.current_session:
            return
        mask_sess = self.df["sessao"] == self.current_session
        alterou = False
        for idx in self.df[mask_sess].index:
            arq = str(self.df.at[idx, "arquivo"])
            off = self._get_offset_for_file(arq)
            ini_lfp = float(self.df.at[idx, "janela_ini_s"])
            fim_lfp = float(self.df.at[idx, "janela_fim_s"])
            novo_ini = format_time(ini_lfp + off)
            novo_fim = format_time(fim_lfp + off)
            if self.df.at[idx, "video_tempo_ini"] != novo_ini or self.df.at[idx, "video_tempo_fim"] != novo_fim:
                self.df.at[idx, "video_tempo_ini"] = novo_ini
                self.df.at[idx, "video_tempo_fim"] = novo_fim
                alterou = True
        if alterou:
            self._salvar_csv(mostrar_mensagem=False)

    def _tentar_autodetectar_video(self, sessao_nome: str):
        """Procura o vídeo (.MPG, .mp4, etc.) da sessão nas pastas do workspace."""
        token = sessao_nome.split("_")[0].strip().lower()
        root_dir = os.path.dirname(os.path.abspath(self.csv_path))
        pastas_busca = [
            root_dir,
            os.path.join(root_dir, "LAC_NOCI"),
            os.path.join(root_dir, "EXPLORACAO_OBJETOS"),
            self.config.get("ultima_pasta_videos", "")
        ]
        exts = (".mpg", ".mpeg", ".mp4", ".avi", ".mkv", ".mov", ".wmv")
        for raiz in pastas_busca:
            if not raiz or not os.path.exists(raiz):
                continue
            for dp, _, fns in os.walk(raiz):
                if token in dp.lower() or os.path.basename(dp).lower() == token:
                    for f in fns:
                        if f.lower().endswith(exts):
                            full_p = os.path.join(dp, f)
                            self.var_video_path.set(full_p)
                            self._carregar_video_arquivo(full_p)
                            return

    def _selecionar_video(self):
        inicial = self.config.get("ultima_pasta_videos", os.path.dirname(self.csv_path))
        caminho = filedialog.askopenfilename(
            title="Selecione o arquivo de vídeo da sessão",
            initialdir=inicial,
            filetypes=[
                ("Todos os Vídeos (*.mpg, *.mp4, *.avi, ...)", "*.mpg *.MPG *.mpeg *.MPEG *.mp4 *.MP4 *.avi *.AVI *.mkv *.MKV *.mov *.MOV *.wmv *.WMV *.flv *.FLV"),
                ("Vídeos MPG / MPEG (*.mpg, *.mpeg)", "*.mpg *.MPG *.mpeg *.MPEG"),
                ("Vídeos MP4 (*.mp4)", "*.mp4 *.MP4"),
                ("Todos os arquivos (*.*)", "*.*")
            ]
        )
        if caminho:
            self.var_video_path.set(caminho)
            self.config["ultima_pasta_videos"] = os.path.dirname(caminho)
            self._carregar_video_arquivo(caminho)
            self._salvar_info_sessao()

    def _carregar_video_arquivo(self, caminho: str):
        if not self.player.load_video(caminho):
            messagebox.showwarning("Aviso", f"Não foi possível abrir o vídeo:\n{caminho}")
            return

        self.slider_time.configure(to=self.player.duration_s)
        self._redesenhar_frame_atual()
        self._salvar_info_sessao()

    def _salvar_info_sessao(self):
        if not self.current_session:
            return
        if "sessoes" not in self.config:
            self.config["sessoes"] = {}

        if self.current_file:
            try:
                self.file_offsets[self.current_file] = float(self.var_offset.get())
            except Exception:
                pass

        base_off = self.file_offsets.get(self.session_files[0], self.var_offset.get()) if self.session_files else self.var_offset.get()

        self.config["sessoes"][self.current_session] = {
            "video_path": self.var_video_path.get(),
            "offset": base_off,
            "offsets_arquivos": self.file_offsets
        }
        self._save_config()

    def _ajustar_offset(self, delta: float):
        try:
            novo = round(float(self.var_offset.get()) + delta, 2)
            self.var_offset.set(novo)
            self._ao_alterar_offset()
        except Exception:
            pass

    def _ao_alterar_offset(self, event=None):
        try:
            off = float(self.var_offset.get())
        except (ValueError, tk.TclError):
            return

        self.lbl_offset_time.configure(text=f"({format_time(off)})")

        if self.current_file:
            self.file_offsets[self.current_file] = off

        self._salvar_info_sessao()

        # Atualiza video_tempo_ini e video_tempo_fim no DataFrame e no CSV para o arquivo atual
        if self.df is not None and self.current_session:
            mask = (self.df["sessao"] == self.current_session)
            if self.current_file:
                mask = mask & (self.df["arquivo"] == self.current_file)

            for idx in self.df[mask].index:
                ini_lfp = float(self.df.at[idx, "janela_ini_s"])
                fim_lfp = float(self.df.at[idx, "janela_fim_s"])
                self.df.at[idx, "video_tempo_ini"] = format_time(ini_lfp + off)
                self.df.at[idx, "video_tempo_fim"] = format_time(fim_lfp + off)

            self._salvar_csv(mostrar_mensagem=False)
            self.var_status_save.set(f"✔ Offset salvo ({self.current_file or 'Sessão'}: {off:+.1f}s)")

        self._atualizar_tabela_momentos(preservar_selecao=True)
        if self.current_row_idx is not None:
            self._pular_para_momento(self.current_row_idx, dar_play=self.player.is_playing)

    def _definir_offset_pelo_frame(self):
        """Define o offset baseado no frame atual como marco zero do arquivo neural atual."""
        if self.player.cap is None:
            messagebox.showinfo("Vídeo", "Carregue um vídeo primeiro.")
            return

        arq_nome = self.current_file or "deste arquivo"
        resposta = messagebox.askyesno(
            "Fixar Offset",
            f"Deseja definir o tempo atual do vídeo ({format_time(self.player.current_time_s)}) "
            f"como o início do arquivo:\n{arq_nome}\n\n(Offset = {self.player.current_time_s:.2f}s)?"
        )
        if resposta:
            self.var_offset.set(round(self.player.current_time_s, 2))
            self._ao_alterar_offset()

    def _atualizar_tabela_momentos(self, preservar_selecao=False):
        """Atualiza a lista de momentos (janelas) para a sessão selecionada."""
        if self.df is None or not self.current_session:
            return

        selected_id = self.tree.selection()
        item_selecionado_idx = None
        if preservar_selecao and selected_id:
            try:
                item_selecionado_idx = int(self.tree.item(selected_id[0], "text"))
            except Exception:
                pass

        for item in self.tree.get_children():
            self.tree.delete(item)

        sub_df = self.df[self.df["sessao"] == self.current_session]
        filtro = self.var_status_filter.get()

        total = len(sub_df)
        anotados = 0
        self.filtered_indices = []

        for idx, row in sub_df.iterrows():
            comp = str(row.get("comportamento", "")).strip()
            ja_anotado = bool(comp and comp.lower() != "nan")
            if ja_anotado:
                anotados += 1

            if filtro == "Apenas Pendentes" and ja_anotado:
                continue
            if filtro == "Apenas Concluídos" and not ja_anotado:
                continue

            self.filtered_indices.append(idx)

            status_sym = "✔" if ja_anotado else "⏳"
            ini_lfp = float(row.get("janela_ini_s", 0.0))
            fim_lfp = float(row.get("janela_fim_s", 10.0))
            arq_nome = str(row.get("arquivo", ""))

            offset = self._get_offset_for_file(arq_nome)
            ini_vid = ini_lfp + offset
            fim_vid = fim_lfp + offset

            tempo_lfp_str = f"{ini_lfp:.1f} - {fim_lfp:.1f}"
            tempo_vid_str = f"{format_time(ini_vid)} - {format_time(fim_vid)}"
            pares = str(row.get("pares_detectados", ""))

            item_id = self.tree.insert("", "end", text=str(idx), values=(
                status_sym, arq_nome, tempo_lfp_str, tempo_vid_str, pares, comp
            ))

            if item_selecionado_idx == idx:
                self.tree.selection_set(item_id)
                self.tree.see(item_id)

        pct = (anotados / total * 100) if total > 0 else 0
        self.var_progress_info.set(f"Progresso: {anotados}/{total} ({pct:.0f}%)")

        if not self.tree.selection() and self.tree.get_children():
            primeiro = self.tree.get_children()[0]
            self.tree.selection_set(primeiro)
            self.tree.see(primeiro)
            self._ao_selecionar_momento_tabela()

    def _ao_selecionar_momento_tabela(self, event=None):
        sel = self.tree.selection()
        if not sel:
            return
        idx = int(self.tree.item(sel[0], "text"))
        self.current_row_idx = idx

        row = self.df.loc[idx]
        arq_row = str(row.get("arquivo", ""))

        # Se o momento pertencer a outro arquivo da sessão, sincroniza dropdown e offset
        if arq_row and arq_row != self.current_file:
            self.current_file = arq_row
            self.var_current_file.set(arq_row)
            off = self._get_offset_for_file(arq_row)
            self.var_offset.set(off)
            self.lbl_offset_time.configure(text=f"({format_time(off)})")

        self.var_comportamento.set(str(row.get("comportamento", "")).replace("nan", ""))
        self.var_observacoes.set(str(row.get("observacoes", "")).replace("nan", ""))

        self._pular_para_momento(idx, dar_play=self.var_autoplay.get())
        self.cb_comportamento.focus_set()

    def _pular_para_momento(self, row_idx: int, dar_play: bool = True):
        if self.df is None or row_idx not in self.df.index:
            return

        row = self.df.loc[row_idx]
        ini_lfp = float(row.get("janela_ini_s", 0.0))
        fim_lfp = float(row.get("janela_fim_s", ini_lfp + 10.0))
        arq_row = str(row.get("arquivo", ""))

        offset = self._get_offset_for_file(arq_row)
        t_ini_vid = max(0.0, ini_lfp + offset)
        t_fim_vid = max(0.0, fim_lfp + offset)

        self.player.window_start_s = t_ini_vid
        self.player.window_end_s = t_fim_vid
        self.player.seek(t_ini_vid)

        if dar_play and self.player.cap is not None:
            self.player.is_playing = True
            self.btn_play.configure(text="⏸ Pausar (Espaço)")
        self._redesenhar_frame_atual()

    def _repetir_janela(self):
        if self.current_row_idx is not None:
            self._pular_para_momento(self.current_row_idx, dar_play=True)

    def _toggle_play(self):
        if self.player.cap is None:
            return
        self.player.is_playing = not self.player.is_playing
        txt = "⏸ Pausar (Espaço)" if self.player.is_playing else "▶ Play (Espaço)"
        self.btn_play.configure(text=txt)

    def _avancar_tempo(self, delta_s: float):
        if self.player.cap is None:
            return
        self.player.seek(self.player.current_time_s + delta_s)
        self._redesenhar_frame_atual()

    def _atualizar_slider_pos(self, pos_s: float):
        """Atualiza a posição do slider programaticamente sem disparar evento de arrasto em loop."""
        self._is_updating_slider = True
        try:
            self.slider_time.set(pos_s)
        finally:
            self._is_updating_slider = False

    def _ao_arrastar_slider(self, val):
        if self._is_updating_slider or self.player.cap is None:
            return
        try:
            pos = float(val)
        except (ValueError, TypeError):
            return
        self.player.seek(pos)
        frame_img = self.player.get_current_frame()
        if frame_img is not None:
            self._desenhar_frame(frame_img)
        self.var_current_time.set(
            f"{format_time(self.player.current_time_s)} / {format_time(self.player.duration_s)}"
        )

    def _ao_mudar_velocidade(self, event=None):
        sp_str = self.var_speed.get().replace("x", "")
        try:
            self.player.playback_speed = float(sp_str)
        except Exception:
            self.player.playback_speed = 1.0

    def _ao_alternar_loop(self):
        self.player.loop_window = self.var_loop_window.get()

    def _aplicar_comportamento_rapido(self, comp_nome: str):
        """Aplica um comportamento de atalho e já salva/avança."""
        self.var_comportamento.set(comp_nome)
        self._salvar_e_proximo()

    def _salvar_e_proximo(self):
        """Salva a anotação atual no DataFrame e pula para a próxima janela."""
        if self.current_row_idx is None or self.df is None:
            return

        comp = self.var_comportamento.get().strip()
        obs = self.var_observacoes.get().strip()

        # Atualiza DataFrame
        self.df.at[self.current_row_idx, "comportamento"] = comp
        self.df.at[self.current_row_idx, "observacoes"] = obs

        # Atualiza também video_tempo_ini e fim com o offset real atual do arquivo
        ini_lfp = float(self.df.at[self.current_row_idx, "janela_ini_s"])
        fim_lfp = float(self.df.at[self.current_row_idx, "janela_fim_s"])
        arq = str(self.df.at[self.current_row_idx, "arquivo"])
        off = self._get_offset_for_file(arq)
        self.df.at[self.current_row_idx, "video_tempo_ini"] = format_time(ini_lfp + off)
        self.df.at[self.current_row_idx, "video_tempo_fim"] = format_time(fim_lfp + off)

        # Salva em disco silenciosamente
        self._salvar_csv(mostrar_mensagem=False)
        self.var_status_save.set(f"✔ Salvo no template_comportamento.csv ({datetime.datetime.now().strftime('%H:%M:%S')})")

        # Atualiza a linha na tabela visual
        sel = self.tree.selection()
        if sel:
            valores = list(self.tree.item(sel[0], "values"))
            valores[0] = "✔" if comp else "⏳"
            valores[5] = comp
            self.tree.item(sel[0], values=valores)

        # Atualiza histórico de comportamentos
        self._atualizar_historico_comportamentos()

        # Avança para o próximo momento
        self._avancar_momento()

    def _avancar_momento(self):
        """Move a seleção para a próxima linha da tabela."""
        filhos = self.tree.get_children()
        if not filhos:
            return

        sel = self.tree.selection()
        if not sel:
            self.tree.selection_set(filhos[0])
            self.tree.see(filhos[0])
            return

        cur_id = sel[0]
        idx_in_tree = filhos.index(cur_id)
        if idx_in_tree + 1 < len(filhos):
            prox_id = filhos[idx_in_tree + 1]
            self.tree.selection_set(prox_id)
            self.tree.see(prox_id)
            self._ao_selecionar_momento_tabela()
        else:
            resposta = messagebox.askyesno(
                "Sessão Concluída",
                f"Você chegou ao final dos momentos desta sessão!\n\n"
                f"Deseja avançar para a próxima sessão agora?"
            )
            if resposta:
                self._avancar_sessao()

    def _voltar_momento(self):
        """Move a seleção para a linha anterior da tabela."""
        filhos = self.tree.get_children()
        if not filhos:
            return

        sel = self.tree.selection()
        if not sel:
            return

        cur_id = sel[0]
        idx_in_tree = filhos.index(cur_id)
        if idx_in_tree > 0:
            ant_id = filhos[idx_in_tree - 1]
            self.tree.selection_set(ant_id)
            self.tree.see(ant_id)
            self._ao_selecionar_momento_tabela()

    def _avancar_sessao(self):
        """Avança para a próxima sessão na lista."""
        vals = self.cb_sessao["values"]
        if not vals:
            return
        cur_idx = self.cb_sessao.current()
        if cur_idx + 1 < len(vals):
            self.cb_sessao.current(cur_idx + 1)
            self._ao_selecionar_sessao()
        else:
            messagebox.showinfo("Parabéns!", "Você já está na última sessão!")

    def _video_loop(self):
        """Loop de renderização de quadros de vídeo a ~30 fps."""
        delay = 33
        if self.player.cap is not None and self.player.is_playing:
            speed = max(0.2, self.player.playback_speed)
            delay = int(max(10, 33 / speed))

            if self.player.current_time_s >= self.player.window_end_s:
                if self.player.loop_window:
                    self.player.seek(self.player.window_start_s)
                else:
                    self.player.is_playing = False
                    self.btn_play.configure(text="▶ Play (Espaço)")

            frame_img = self.player.get_current_frame()
            if frame_img is not None:
                self._desenhar_frame(frame_img)
                self._atualizar_slider_pos(self.player.current_time_s)
                self.var_current_time.set(
                    f"{format_time(self.player.current_time_s)} / {format_time(self.player.duration_s)}"
                )

        self.after(delay, self._video_loop)

    def _redesenhar_frame_atual(self):
        """Redesenha o frame no canvas quando o vídeo está pausado ou a janela é redimensionada."""
        if self.player.cap is not None:
            frame_img = self.player.get_current_frame()
            if frame_img is not None:
                self._desenhar_frame(frame_img)
            self._atualizar_slider_pos(self.player.current_time_s)
            self.var_current_time.set(
                f"{format_time(self.player.current_time_s)} / {format_time(self.player.duration_s)}"
            )
        else:
            self._desenhar_placeholder()

    def _desenhar_frame(self, pil_image: Image.Image):
        """Ajusta proporcionalmente a imagem ao tamanho do canvas e exibe."""
        cw = self.canvas_video.winfo_width()
        ch = self.canvas_video.winfo_height()
        if cw <= 10 or ch <= 10:
            return

        iw, ih = pil_image.size
        escala = min(cw / iw, ch / ih)
        nw, nh = max(1, int(iw * escala)), max(1, int(ih * escala))

        resized = pil_image.resize((nw, nh), Image.Resampling.BILINEAR)
        self.photo_tk = ImageTk.PhotoImage(resized)

        self.canvas_video.delete("all")
        x_pos = (cw - nw) // 2
        y_pos = (ch - nh) // 2
        self.canvas_video.create_image(x_pos, y_pos, anchor="nw", image=self.photo_tk)

    def _desenhar_placeholder(self):
        """Exibe mensagem amigável quando nenhum vídeo estiver carregado."""
        cw = self.canvas_video.winfo_width()
        ch = self.canvas_video.winfo_height()
        self.canvas_video.delete("all")
        if cw > 10 and ch > 10:
            self.canvas_video.create_text(
                cw // 2, ch // 2,
                text="Nenhum vídeo carregado.\nClique em 'Selecionar Vídeo...' acima para abrir o vídeo da sessão.",
                fill="#9CA3AF", font=("Segoe UI", 12), justify="center"
            )


if __name__ == "__main__":
    csv_alvo = DEFAULT_CSV
    if len(sys.argv) > 1:
        csv_alvo = sys.argv[1]

    app = AnotadorApp(csv_path=csv_alvo)
    app.mainloop()
