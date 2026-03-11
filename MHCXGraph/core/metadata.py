from __future__ import annotations

import pydssp
import re
import torch

import numpy as np

from MHCXGraph.utils.logging_utils import get_log

log = get_log()

BACKBONE_ATOMS = ("N", "CA", "C", "O")

def secondary_structure(G, **ctx):
    """
    Anota estrutura secundária via DSSP.
    """
    chains = ctx.get("chains")
    include_noncanonical_residues = ctx.get("include_noncanonical_residues")
    max_gap_helix = ctx.get("max_gap_helix", 0)
    if chains is None:
        return G

    def build_pydssp_input(chain: dict[str, list], chain_id, include_noncanonical_residues):
        residues = chain["canonical_aminoacid_residues"].copy()
        if include_noncanonical_residues:
            residues += chain["noncanonical_aminoacid_residues"]

        if not residues:
            log.info(f"The chain {chain_id} doesn't have any aminoacid residue.")
            return None

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
            raise ValueError(f"No residues with complete backbone atoms (N,CA,C,O) were found. Chain: {chain_id}")

        coord = torch.from_numpy(np.asarray(coords, dtype=np.float32))
        return coord, used_residues, sequence

    def assign_ss_to_chain(chain, chain_id, include_noncanonical_residues, max_gap_helix):
        """
        Runs pydssp for a chain and maps SS back to residues.
        Returns:
          ss_used: string length L for used_residues
          used_residues: list of residues used (index aligns with ss_used)
          full_ss: string length len(all_residues) (omitted residues filled with '?')
          all_residues: list of all residues in chain
        """
        pydsspInput = build_pydssp_input(
            chain, chain_id, include_noncanonical_residues=include_noncanonical_residues
        )

        if pydsspInput is None:
            return {}

        coord, used_residues, sequence = pydsspInput
        donor_mask = sequence != 'PRO'

        ss = pydssp.assign(coord, donor_mask=donor_mask)
        ss_used = "".join(ss)

        if max_gap_helix > 0:
            pattern = rf"(?<=[Hh])[-]{{1,{max_gap_helix}}}(?=[Hh])"
            ss_used = re.sub(pattern, lambda m: "H" * len(m.group()), ss_used)

        ss_map = {res[0]: s for res, s in zip(used_residues, ss_used, strict=True)}

        return ss_map

    ss_map = {}

    for chain in chains:
        ss_map = ss_map | assign_ss_to_chain(chains[chain], chain, include_noncanonical_residues=include_noncanonical_residues, max_gap_helix=max_gap_helix)

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
