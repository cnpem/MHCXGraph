from __future__ import annotations

import pydssp
import torch
from statistics import mean
from typing import Any

import numpy as np

# from MHCXGraph.core.pdb_graph_builder import ResidueList

BACKBONE_ATOMS = ("N", "CA", "C", "O")

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
    chains = ctx.get("chains")
    structure = ctx.get("structure")
    residue_map = ctx.get("residue_map")  # {node_id: Residue}
    pdb_path = ctx.get("pdb_path")
    
    def build_pydssp_input(chain: dict[str, str], *, include_nonstandard_aa=True):        
        residues = chain["canonical_aminoacid_residues"]
        if include_nonstandard_aa:
            residues += chain["noncanonical_aminoacid_residues"]

        used_residues = []
        coords = []

        for (res_label, res, _, _) in residues:
            for atom in BACKBONE_ATOMS:
                if atom not in res:
                    break
            else:
                bb = []
                for atom in BACKBONE_ATOMS:
                    bb.append(res[atom].get_coord().astype(np.float32))
                coords.append(bb)
                used_residues.append((res_label, res))

        if len(coords) == 0:
            raise ValueError("No residues with complete backbone atoms (N,CA,C,O) were found.")

        coord = torch.from_numpy(np.asarray(coords, dtype=np.float32))  # [L,4,3]
        return coord, used_residues, residues

    def assign_ss_to_chain(chain, *, include_nonstandard_aa=True, write_xtra=False):
        """
        Runs pydssp for a chain and maps SS back to residues.
        Returns:
          ss_used: string length L for used_residues
          used_residues: list of residues used (index aligns with ss_used)
          full_ss: string length len(all_residues) (omitted residues filled with '?')
          all_residues: list of all residues in chain
        """
        coord, used_residues, all_residues = build_pydssp_input(
            chain, include_nonstandard_aa=include_nonstandard_aa
        )
        
        ss = pydssp.assign(coord)  # per README: returns array like ['-','H',...]
        ss_used = "".join(ss.tolist() if hasattr(ss, "tolist") else list(ss))

        # Map back to residues
        if write_xtra:
            for res, s in zip(used_residues, ss_used):
                res.xtra["SS_PYDSSP"] = s

        # Optional: expand to a "full chain" string aligned to *all* residues in the chain
        # Fill residues not used by pydssp (ligands, missing backbone, etc.) with '?'
        used_set = set(used_residues)
        full_chars = []
        used_iter = iter(ss_used)
        # But we can’t just iterate ss_used linearly unless we check membership
        ss_map = {res[0]: s for res, s in zip(used_residues, ss_used)}
        for res in all_residues:
            full_chars.append(ss_map.get(res[0], "?"))
        full_ss = "".join(full_chars)

        return ss_used, used_residues, full_ss, all_residues, ss_map
    
    for chain in chains:
        ss_used, used_residues, full_ss, all_residues, ss_map = assign_ss_to_chain(chains[chain])
        print(chain)
        print(ss_map, ss_used, full_ss)
    # coord, sequence = pydssp.read_pdbtext(open(pdb_path, 'r').read(), return_sequence=True)
    # coord = torch.Tensor(coord).to("cpu")
    # donor_mask = sequence != 'PRO'
    # main calcuration
    # dsspline = ''.join(pydssp.assign(coord, donor_mask=donor_mask))


    # Se 'enabled' existir e for False, aborta. Se não existir, considera habilitado.
    if structure is None or not residue_map:
        return G
    if dssp_cfg is not None and hasattr(dssp_cfg, "enabled") and not bool(dssp_cfg.enabled):
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

def secondary_structure_backup(G, **ctx):
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
    if dssp_cfg is not None and hasattr(dssp_cfg, "enabled") and not bool(dssp_cfg.enabled):
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
