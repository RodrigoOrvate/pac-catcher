"""
ns2_utils.py
==========================================
Shim de compatibilidade: as funções de leitura de dados eletrofisiológicos
foram movidas para `pac_core/io.py`. Este módulo continua existindo para que
os importadores atuais (`from ns2_utils import le_ns2`, etc.) não precisem
mudar. Novo código deve importar diretamente de `pac_core.io`.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pac_core.io import (  # noqa: F401 - re-export para compatibilidade
    le_ns2,
    le_bin_legado,
    carrega_dados,
    fatia_janela,
    concatena_sessao,
)
