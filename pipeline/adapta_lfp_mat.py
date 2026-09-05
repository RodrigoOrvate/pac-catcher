"""Adapta pipeline para LFP_HG_HFO.mat (exemplo) — lê float64 1×300000, trata como 1 canal"""
import scipy.io, numpy as np, os, sys
sys.path.insert(0, "pipeline")

def carrega_lfp_mat(path):
    m = scipy.io.loadmat(path)
    # LFP_HG_HFO.mat tem lfpHG e lfpHFO (1 canal cada, 300000 pts, float64, ~300s @ 1kHz)
    for chave in ["lfpHG", "lfpHFO"]:
        if chave in m:
            return m[chave]  # retorna array (1, 300000)
    raise KeyError("Nem lfpHG nem lfpHFO encontrados")

def rodar_exemplo():
    dados = carrega_lfp_mat("LFP_HG_HFO.mat")
    fs = 1000.0
    n_pts = dados.shape[1] if dados.ndim > 1 else len(dados)
    print(f"LFP_HG_HFO: {dados.shape}, fs={fs}Hz, duracao={n_pts/fs:.0f}s")
    # Salva como CSV temporário para pipeline usar (simula .ns2)
    np.savetxt("tmp_lfp_hg_hfo.csv", dados.T, delimiter=",", fmt="%.6f")
    print("-> tmp_lfp_hg_hfo.csv salvo (1 canal, 300s @ 1kHz)")
    # Registra resultado como exemplo
    with open("README.md", "a", encoding="utf-8") as f:
        f.write("\n\n## Exemplo LFP_HG_HFO (2026-09-05)\n")
        f.write(f"- Arquivo: `LFP_HG_HFO.mat` | `lfpHG` / `lfpHFO` (float64, {dados.shape}, 1kHz)\n")
        f.write(f"- Resultado de exemplo: `{n_pts/fs:.0f}s` de LFP real, usado como dado de teste do pipeline\n")
        f.write(f"- Adaptacao: `carrega_lfp_mat()` le .mat e exporta CSV para pipeline (simula ns2)\n")
    print("-> Resultado documentado no README.md")

if __name__ == "__main__":
    rodar_exemplo()
