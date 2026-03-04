# core/config.py
from __future__ import annotations

"""
Configuration primitives for structural graph construction.

This module defines a single configuration dataclass (`GraphConfig`)
and a helper factory (`make_default_config`). The configuration covers:
- chain selection and waters handling,
- geometric cutoffs and distance-matrix storage,
- ASA/RSA computation settings (DSSP or SR),
- centroid policy for node coordinates.

Tracking is external and not represented here.
"""

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Literal

Granularity = Literal["all_atoms", "backbone", "side_chain", "ca_only"]
# from enum import Enum

# class Granularity(str, Enum):
#     ALL_ATOMS = "all_atoms"
#     BACKBONE = "backbone"
#     SIDE_CHAIN = "side_chain"
#     CA_ONLY = "ca_only"

"""Centroid policy used to derive node coordinates.

- ``"all_atoms"``  — centroid of all heavy atoms per residue (fallback to any atom if needed).
- ``"backbone"``   — centroid over backbone atoms (N, CA, C, O, OXT when present).
- ``"side_chain"`` — centroid over side-chain heavy atoms; GLY falls back to backbone.
- ``"ca_only"``    — CA coordinate; if CA is missing, CB is used; if both are missing, fallback to first heavy atom.
"""


@dataclass
class DSSPConfig:
    """
    Compatibility shim for call paths expecting a DSSP configuration object.
    Only the executable path/name is retained.
    """
    executable: str = "mkdssp"


@dataclass
class GraphConfig:
    """
    Unified configuration for building residue-level graphs from PDB/mmCIF files.

    Parameters
    ----------
    chains
        Iterable of chain IDs to include. ``None`` includes all chains in the selected model.
    include_waters
        If ``True``, includes water molecules (HOH) as nodes and residue–water edges.
    residue_distance_cutoff
        Residue–residue and Water-residue centroid distance cutoff (Å) for adding edges.
    store_distance_matrix
        If ``True``, keeps the residue–residue distance matrix in the returned object.
    compute_rsa
        If ``True``, computes ASA and RSA per residue.
    model_index
        Model index to load from the structure file.
    allow_empty_chains
        If ``False``, raises when requested chains are not found.
    granularity
        Centroid policy for node coordinates. See ``Granularity`` for options.
    verbose
        Enables additional logging when supported by the caller.
    """
    # Selection and waters
    chains: Iterable[str] | None = None
    include_waters: bool = True
    include_noncanonical_residues: bool = True
    include_ligands: bool = True

    # Geometry / graph
    residue_distance_cutoff: float = 10.0
    store_distance_matrix: bool = True

    # ASA / RSA
    compute_rsa: bool = True

    # Structure selection
    model_index: int = 0
    allow_empty_chains: bool = False

    # Node coordinate policy
    granularity: Granularity = "all_atoms"

    # Misc
    verbose: bool = False

    make_virtual_cb_for_gly: bool = True


def make_default_config(
    *,
    edge_threshold: float = 10.0,
    granularity: Granularity = "all_atoms",
    chains: Iterable[str] | None = None,
    compute_rsa: bool = True,
    make_virtual_cb_for_gly: bool = True,
    include_noncanonical_residues: bool = True,
    include_ligands: bool = True,
    include_waters: bool = True,
    verbose: bool = False
) -> GraphConfig:
    """
    Helper for creating a ``GraphConfig`` with common defaults.

    Parameters
    ----------
    edge_threshold
        Residue–residue distance cutoff (Å).
    granularity
        Centroid policy for node coordinates.
    include_waters
        If ``True``, includes water molecules as nodes.
    chains
        Iterable of chain IDs to include; ``None`` includes all chains.
    compute_rsa
        If ``True``, computes ASA/RSA per residue.

    Returns
    -------
    GraphConfig
        Ready-to-use configuration instance.
    """
    return GraphConfig(
        chains=chains,
        include_waters=include_waters,
        residue_distance_cutoff=float(edge_threshold),
        store_distance_matrix=True,
        compute_rsa=compute_rsa,
        model_index=0,
        allow_empty_chains=False,
        granularity=granularity,
        verbose=verbose,
        make_virtual_cb_for_gly=make_virtual_cb_for_gly,
        include_noncanonical_residues=include_noncanonical_residues,
        include_ligands=include_ligands
    )
