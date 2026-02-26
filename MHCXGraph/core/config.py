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

from dataclasses import dataclass
from typing import Iterable, Optional, Literal, Union


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
        Residue–residue centroid distance cutoff (Å) for adding edges.
    water_distance_cutoff
        Water–residue centroid distance cutoff (Å) for adding edges when waters are included.
    store_distance_matrix
        If ``True``, keeps the residue–residue distance matrix in the returned object.
    compute_rsa
        If ``True``, computes ASA and RSA per residue.
    rsa_method
        RSA method: ``"dssp"`` uses an external DSSP executable; ``"sr"`` uses Shrake–Rupley with table normalization.
    dssp_exec
        Executable name or absolute path for DSSP (e.g., ``"mkdssp"`` or ``"dssp"``).
    dssp_acc_array
        Reference maximum-ASA table used by DSSP and SR normalization. One of ``"Sander"``, ``"Wilke"``, ``"Miller"``.
    probe_radius
        Probe radius (Å) for Shrake–Rupley ASA.
    n_points
        Number of surface points per atom for Shrake–Rupley ASA.
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
    chains: Optional[Iterable[str]] = None
    include_waters: bool = True

    # Geometry / graph
    residue_distance_cutoff: float = 10.0
    water_distance_cutoff: float = 6.0
    store_distance_matrix: bool = True

    # ASA / RSA
    compute_rsa: bool = True
    # rsa_method = "dssp"
    rsa_method: Literal["sr", "dssp"]  = "dssp"
    dssp_exec: str = "mkdssp"
    # dssp_acc_array = "Wilke"
    dssp_acc_array: Literal["Sander", "Wilke", "Miller"] = "Wilke",
    probe_radius: float = 1.4
    n_points: int = 960

    # Structure selection
    model_index: int = 0
    allow_empty_chains: bool = False

    # Node coordinate policy
    granularity: Granularity = "all_atoms"

    # Misc
    verbose: bool = False

    make_virtual_cb_for_gly: bool = True

    include_noncanonical_residues: bool = True
    include_ligands: bool = False


def make_default_config(
    *,
    edge_threshold: float = 10.0,
    granularity: Granularity = "all_atoms",
    chains: Optional[Iterable[str]] = None,
    compute_rsa: bool = True,
    # rsa_method = "dssp",
    rsa_method: Literal["sr", "dssp"] = "dssp",
    dssp_exec: str = "mkdssp",
    # dssp_acc_array = "Wilke",
    dssp_acc_array: Literal["Sander", "Wilke", "Miller"] = "Wilke",
    make_virtual_cb_for_gly: bool = True,
    include_noncanonical_residues: bool = True,
    include_ligands: bool = False,
    include_waters: bool = True,
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
    rsa_method
        RSA method: ``"dssp"`` or ``"sr"``.
    dssp_exec
        DSSP executable name or path.
    dssp_acc_array
        DSSP/SR max-ASA reference table.

    Returns
    -------
    GraphConfig
        Ready-to-use configuration instance.
    """
    return GraphConfig(
        chains=chains,
        include_waters=include_waters,
        residue_distance_cutoff=float(edge_threshold),
        water_distance_cutoff=6.0,
        store_distance_matrix=True,
        compute_rsa=compute_rsa,
        rsa_method=rsa_method,
        dssp_exec=dssp_exec,
        dssp_acc_array=dssp_acc_array,
        probe_radius=1.4,
        n_points=960,
        model_index=0,
        allow_empty_chains=False,
        granularity=granularity,
        verbose=False,
        make_virtual_cb_for_gly=make_virtual_cb_for_gly,
        include_noncanonical_residues=include_noncanonical_residues,
        include_ligands=include_ligands
    )
