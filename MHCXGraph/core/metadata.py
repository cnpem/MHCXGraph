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
    chains = ctx.get("chains")
    include_noncanonical_residues = ctx.get("include_noncanonical_residues")
    if chains is None:
        return G

    def build_pydssp_input(chain: dict[str, str], include_noncanonical_residues):
        residues = chain["canonical_aminoacid_residues"]
        if include_noncanonical_residues:
            residues += chain["noncanonical_aminoacid_residues"]

        used_residues = []
        coords = []
        sequence = []
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
                sequence.append(res_label.split(":")[1])

        sequence = np.array(sequence)

        if len(coords) == 0:
            raise ValueError("No residues with complete backbone atoms (N,CA,C,O) were found.")

        coord = torch.from_numpy(np.asarray(coords, dtype=np.float32))  # [L,4,3]
        return coord, used_residues, sequence

    def assign_ss_to_chain(chain, include_noncanonical_residues):
        """
        Runs pydssp for a chain and maps SS back to residues.
        Returns:
          ss_used: string length L for used_residues
          used_residues: list of residues used (index aligns with ss_used)
          full_ss: string length len(all_residues) (omitted residues filled with '?')
          all_residues: list of all residues in chain
        """
        coord, used_residues, sequence = build_pydssp_input(
            chain, include_noncanonical_residues=include_noncanonical_residues
        )

        donor_mask = sequence != 'PRO'
        ss = pydssp.assign(coord, donor_mask=donor_mask)
        ss_used = "".join(ss)

        ss_map = {res[0]: s for res, s in zip(used_residues, ss_used, strict=True)}

        return ss_map

    ss_map = {}

    for chain in chains:
        ss_map = ss_map | assign_ss_to_chain(chains[chain], include_noncanonical_residues=include_noncanonical_residues)

    for nid, data in G.nodes(data=True):
        if data.get("kind") == "water":
            continue

        ss_val = ss_map[nid]
        if ss_val is not None and ss_val == ' ':
            ss_val = '-'
        G.nodes[nid]["ss"] = ss_val

    ss_vals = [d.get("ss") for _, d in G.nodes(data=True) if d.get("kind") != "water"]
    counts = {}

    for s in ss_vals:
        counts[s] = counts.get(s, 0) + 1
    G.graph["ss_counts"] = counts

    return G

