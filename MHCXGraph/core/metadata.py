from __future__ import annotations
from typing import Any, Dict, Optional
from statistics import mean

def _graph_mean(values):
    vals = [v for v in values if v is not None]
    return mean(vals) if vals else None

def rsa(G, **ctx):
    rsa_vals = [d.get("rsa") for _, d in G.nodes(data=True) if d.get("kind") != "water"]
    G.graph["rsa_mean"] = _graph_mean(rsa_vals)
    G.graph["rsa_count_nonnull"] = sum(1 for x in rsa_vals if x is not None)
    G.graph["rsa_prop_exposed_025"] = (
        sum(1 for x in rsa_vals if x is not None and x >= 0.25) / len(rsa_vals)
    ) if rsa_vals else None
    return G

# ------------------ PATCH: robustez para DSSP ------------------

def _dssp_exec_from(dssp_cfg: Any) -> str:
    """
    Extrai o executável do DSSP de configs em formatos diferentes:
    - Graphein-like: dssp_config.executable
    - Versão antiga que usávamos: dssp_config.dssp_path
    Fallback: "mkdssp"
    """
    if dssp_cfg is None:
        return "mkdssp"
    exe = getattr(dssp_cfg, "executable", None) or getattr(dssp_cfg, "dssp_path", None)
    return exe or "mkdssp"

def secondary_structure(G, **ctx):
    """
    Anota estrutura secundária via DSSP.

    Aceita dssp_config com ou sem atributo `.enabled`.
    Usa dssp_config.executable (Graphein) ou dssp_config.dssp_path (legado).
    """
    dssp_cfg = ctx.get("dssp_config")
    structure = ctx.get("structure")
    residue_map = ctx.get("residue_map")  # {node_id: Residue}
    pdb_path = ctx.get("pdb_path")

    # Se 'enabled' existir e for False, aborta. Se não existir, considera habilitado.
    if structure is None or not residue_map:
        return G
    if dssp_cfg is not None and hasattr(dssp_cfg, "enabled") and not bool(getattr(dssp_cfg, "enabled")):
        return G

    try:
        from Bio.PDB.DSSP import DSSP
    except Exception:
        return G

    exec_path = _dssp_exec_from(dssp_cfg)

    # DSSP por modelo 0
    try:
        model = structure[0]
    except Exception:
        return G

    try:
        dssp = DSSP(model, pdb_path, dssp=exec_path)
    except Exception:
        # Não quebra o pipeline se DSSP falhar
        return G

    # Mapear SS para cada nó (pula águas)
    for nid, data in G.nodes(data=True):
        if data.get("kind") == "water":
            continue
        res = residue_map.get(nid)
        if res is None:
            continue
        chain_id = res.get_parent().id
        hetflag, resseq, icode = res.id

        # dssp indexa por (chain_id, (' ', resseq, icode))
        ss_val = None
        for key in ((chain_id, res.id), (chain_id, (' ', resseq, icode))):
            try:
                ss_val = dssp[key][2]  # coluna de SS
                break
            except Exception:
                continue

        if ss_val is not None and ss_val == ' ':
            ss_val = 'C'  # coil como no costume

        G.nodes[nid]["ss"] = ss_val

    # Estatística simples
    ss_vals = [d.get("ss") for _, d in G.nodes(data=True) if d.get("kind") != "water"]
    counts = {}
    for s in ss_vals:
        counts[s] = counts.get(s, 0) + 1
    G.graph["ss_counts"] = counts
    return G
