from __future__ import annotations

import bisect
import json
import logging
import math
import os
import time
from collections.abc import Sequence  #, TypeVarTuple, Unpack
from itertools import chain, combinations, product
from os import path
from typing import Any
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from Bio.PDB.Atom import Atom
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Residue import Residue
from sklearn.metrics.pairwise import cosine_similarity

from MHCXGraph.core.residue_tracking import ResidueTracker, TrackCtx
from MHCXGraph.core.tracking import save
from MHCXGraph.utils.vis_tracer import TraversalTracer

log = logging.getLogger("CRSProtein")

#Ts = TypeVarTuple("Ts")

# Triad = Tuple[str, str, str, Unpack[Ts]]
# Combo = Tuple[Triad, ...]
# Filtered = Dict[Tuple, List[Combo]]

def save_pdb_with_spheres(atomic_data, selected_residues_data, pdb_filename):
    """
    Salva o arquivo PDB com esferas para os resíduos selecionados.

    Parameters
    ----------
    atomic_data : numpy.ndarray
        Dados atômicos da estrutura.
    selected_residues_data : list
        Lista com os dados dos resíduos selecionados, incluindo coordenadas e profundidade.
    pdb_filename : str
        Nome do arquivo PDB de saída.
    """


    with open(pdb_filename, 'w') as f:
        # Adiciona os átomos originais
        for row in atomic_data:
            atom_line = f"ATOM  {int(row[0]):5d}  {row[3]:<4} {row[2]:<3} {row[1]:<1} {int(row[4]):>4} {row[5]:>8.3f}{row[6]:>8.3f}{row[7]:>8.3f}{row[8]:>6.2f}{row[9]:>6.2f}\n"
            f.write(atom_line)

        # Adiciona as esferas para os resíduos selecionados
        for residue in selected_residues_data:
            residue_number = residue["ResidueNumber"]
            chain = residue["Chain"]
            coordinates = residue["Coordinates"]
            radius = 1.5  # Definindo um raio para a es fera

            sphere_line = "HETATM{:5d}  SPC {:<3} {:<1} {:>4}   {:>8.3f}{:>8.3f}{:>8.3f}  1.00  0.00           SPH\n".format(
                residue_number * 100, "SPC", chain, residue_number, coordinates[0], coordinates[1], coordinates[2]
            )
            f.write(sphere_line)

    print(f"PDB com esferas salvo em {pdb_filename}")

def convert_edges_to_residues(edges: set[frozenset], maps: dict) -> tuple[list, list, list]:
    """Convert the edges that contains tuple of indices to tuple of residues

    Args:
        edges (List[Tuple]): A list that contains tuple of edges that are made of tuples of indices
        maps (Dict): A map that relates the indice to residue

    Returns:
        convert_edge (List[Tuple]): Return edges converted to residues notation
    """
    original_edges = []
    edges_indices = []
    converted_edges = []
    residue_maps_unique = maps["residue_maps_unique"]
    possible_nodes_map = maps["possible_nodes"]
    for edge in edges:
        edge_list = list(edge)
        node1, node2 = possible_nodes_map[edge_list[0]], possible_nodes_map[edge_list[1]]
        converted_node1 = tuple(f"{residue_maps_unique[idx][0]}:{residue_maps_unique[idx][2]}:{residue_maps_unique[idx][1]}" for idx in node1)
        converted_node2 = tuple(f"{residue_maps_unique[idx][0]}:{residue_maps_unique[idx][2]}:{residue_maps_unique[idx][1]}" for idx in node2)

        converted_node1_indice = tuple(idx for idx in node1)
        converted_node2_indice = tuple(idx for idx in node2)

        original_edges.append(edge)
        edges_indices.append((converted_node1_indice, converted_node2_indice))
        converted_edges.append((converted_node1, converted_node2))

    return original_edges, edges_indices, converted_edges

def filter_maps_by_nodes(data: dict,
                        matrices_dict: dict,
                        distance_threshold: float = 10.0,
                    ) -> tuple[dict, dict]:

    logger = logging.getLogger("association.filter_maps_by_nodes")

    contact_maps = data["contact_maps"]
    rsa_maps = data["rsa_maps"]
    residue_maps = data["residue_maps"]
    nodes_graphs = data["nodes_graphs"]


    maps = {"full_residue_maps": [], "residue_maps_unique": {}}
    pruned_contact_maps = []
    adjacent_contact_maps = []
    thresholded_rsa_maps = []

    for contact_map, rsa_map, residue_map, nodes in zip(
            contact_maps, rsa_maps, residue_maps, nodes_graphs):

        indices = []
        full_res_map = {}
        for i, node in enumerate(nodes):
            parts = node.split(":")
            if len(parts) != 3:
                logger.warning(f"Node '{node}' does not have three parts separated by ':'")
                continue

            chain, res_name, res_num_str = parts
            key = (chain, res_num_str, res_name)

            # The residue_map seems to be equal to full_res_map, but it isn't. The residue_map considers all residues in the protein, even these ones that didn't pass in the filter. But the full_res_map considers only the residues that through all the filters, and now composes the filtered base graph.
            full_res_map[key] = i
            if key in residue_map:
                indices.append(residue_map[key])

        # The np.ix_ creates a submatrix from contact_map, where the order is the same of my ordered list of nodes, matching the order of full_res_map.
        pruned_map = contact_map[np.ix_(indices, indices)]
        np.fill_diagonal(pruned_map, np.nan)
        pruned_contact_maps.append(pruned_map)

        thresh_map = pruned_map.copy()

        # The thresh_map is the same as pruned_map, except for the np.nan for all values greater or equal than "edge_threshold"
        thresh_map[thresh_map >= distance_threshold] = np.nan
        adjacent_contact_maps.append(thresh_map)

        thresholded_rsa_maps.append(rsa_map)

        maps["full_residue_maps"].append(full_res_map)

    if matrices_dict is not None:
        matrices_dict["induced_contact_maps"] = pruned_contact_maps
        matrices_dict["adjacent_contact_maps"] = adjacent_contact_maps
        matrices_dict["thresholded_rsa_maps"] = thresholded_rsa_maps

    return matrices_dict, maps

def indices_graphs(nodes_lists: list[list]) -> list[tuple[int, int]]:
    """Make a list that contains indices that indicates the position of each protein in graph

    Args:
        nodes_list (List): A list of protein's resdiues. Each List has their own residues.

    Returns:
        ranges (List[Tuple]): A list of indicesthat indicates the position of each protein in matrix
    """

    current = 0
    ranges = []
    for nodes in nodes_lists:
        length = len(nodes)
        ranges.append((current, current + length))
        current += length
    return ranges

def value_to_class(
    value: float,
    bin_width: float,        # Fixed size of each division
    threshold: float,
    inverse: bool = False,
    upper_bound: float = 100.0,
    close_tolerance: float = 0.1,  # Absolute tolerance in 'value' units
    ) -> int | list[int] | None:
    """
    Non-inverse:
      - Domain: (0, threshold]
      - Bins: (L, R] spaced by bin_width. Last bin may be smaller than bin_width.
      - Uncertainty: half of the bin_width.
      - If |value - center_of_bin| <= close_tolerance, returns [class_index].
      - Otherwise, returns all classes whose bins intersect with [value - inc, value + inc].

    Inverse:
      - Domain: [threshold, upper_bound]
      - Same multi-class behavior as non-inverse, but bins live on [threshold, upper_bound].
    """

    if bin_width <= 0:
        return None

    if not inverse:
        lower = 0.0
        upper = threshold
    else:
        lower = threshold
        upper = upper_bound

    if value <= lower or value > upper:
        return None

    span = upper - lower
    local_value = value - lower

    n_divisions = math.ceil(span / bin_width)

    i = math.ceil(local_value / bin_width)
    i = max(1, min(i, n_divisions))

    left = (i - 1) * bin_width
    right = min(i * bin_width, span)
    width = right - left

    center = left + width / 2.0

    tol = max(0.0, close_tolerance)
    tol = min(tol, (width / 2.0) - 1e-12)

    if abs(local_value - center) <= abs(tol):
        return [i]


    inc = bin_width / 2.0
    low = max(0.0, local_value - inc)
    high = min(span, local_value + inc)

    classes: list[int] = []

    for j in range(1, n_divisions + 1):
        L = (j - 1) * bin_width
        R = min(j * bin_width, span)

        if low < R and high > L:
            classes.append(j)

    return classes or None



def backup_value_to_class(
    value: float,
    bin_width: float,        # Fixed size of each division
    threshold: float,
    inverse: bool = False,
    upper_bound: float = 100.0,
    close_tolerance: float = 0.1,  # Absolute tolerance in 'value' units
    ) -> int | list[int] | None:
    """
    Non-inverse:
      - Domain: (0, threshold]
      - Bins: (L, R] spaced by bin_width. Last bin may be smaller than bin_width.
      - Uncertainty: half of the bin_width.
      - If |value - center_of_bin| <= close_tolerance, returns [class_index].
      - Otherwise, returns all classes whose bins intersect with [value - inc, value + inc].

    Inverse:
      - Domain: [threshold, upper_bound]
      - Same multi-class behavior as non-inverse, but bins live on [threshold, upper_bound].
    """
    if bin_width <= 0:
        return None

    if not inverse:
        # Out of bounds check
        if value <= 0 or value > threshold:
            return None

        # Calculate total number of divisions
        n_divisions = math.ceil(threshold / bin_width)

        # Identify the primary bin index (1-based)
        i = math.ceil(value / bin_width)
        i = max(1, min(i, n_divisions))

        # Actual bin boundaries (last bin may be smaller)
        current_bin_left = (i - 1) * bin_width
        current_bin_right = min(i * bin_width, threshold)
        current_bin_actual_width = current_bin_right - current_bin_left

        bin_center = current_bin_left + (current_bin_actual_width / 2.0)

        # Clamp tolerance to at most half the actual bin width
        tol = max(0.0, float(close_tolerance))
        tol = min(tol, (current_bin_actual_width / 2.0) - 1e-12)

        if abs(value - bin_center) <= tol:
            return [i]

        inc = bin_width / 2.0
        low = value - inc
        high = value + inc

        classes: list[int] = []
        for j in range(1, n_divisions + 1):
            L = (j - 1) * bin_width
            R = min(j * bin_width, threshold)
            if (low < R) and (high > L):
                classes.append(j)

        return classes if classes else None

    # Inverse logic: mapping values from threshold up to upper_bound
    if value < threshold or value > upper_bound:
        return None

    span = upper_bound - threshold
    if span <= 0:
        return None

    n_divisions = math.ceil(span / bin_width)

    # Put the inverse domain on a local axis starting at 0
    rel_value = value - threshold

    # Primary bin index (1-based)
    i = math.ceil(rel_value / bin_width)
    i = max(1, min(i, n_divisions))

    # Actual bin boundaries in the inverse domain (last bin may be smaller)
    current_bin_left = (i - 1) * bin_width
    current_bin_right = min(i * bin_width, span)
    current_bin_actual_width = current_bin_right - current_bin_left

    bin_center = current_bin_left + (current_bin_actual_width / 2.0)

    tol = max(0.0, float(close_tolerance))
    tol = min(tol, (current_bin_actual_width / 2.0) - 1e-12)

    if abs(rel_value - bin_center) <= tol:
        return [i]

    inc = bin_width / 2.0
    low = rel_value - inc
    high = rel_value + inc

    # Clamp to inverse domain [0, span]
    low = max(0.0, low)
    high = min(span, high)

    classes: list[int] = []
    for j in range(1, n_divisions + 1):
        L = (j - 1) * bin_width
        R = min(j * bin_width, span)
        if (low < R) and (high > L):
            classes.append(j)

    return classes if classes else None

def create_classes_bin(total_length: float, bin_width: float):
    """Creates a dictionary of bins with their respective ranges."""
    n_bins = math.ceil(total_length / bin_width)
    return {
        str(n + 1): (n * bin_width, min((n + 1) * bin_width, total_length))
        for n in range(n_bins)
    }

def find_class(classes: dict[str, dict[str, float]], value: float):
    """
    Return all class names whose interval contains value.
    Assumes structure:
        classes = { "bin_1": {"low": ..., "high": ...}, ... }
    """
    hits = []
    for name, interval in classes.items():
        low = interval["low"]
        high = interval["high"]
        # choose your convention: [low, high] or [low, high)
        if low <= value <= high:
            hits.append(name)
    if not hits:
        return None
    # Return single class or list depending on your current semantics
    return hits[0] if len(hits) == 1 else hits


def residue_to_tuple(res):
    res_split = res.split(":")

    return (res_split[0], res_split[2], res_split[1])

def _as_list(x):
    """Converte None -> [], int -> [int], lista/tupla/conjunto -> lista."""
    if x is None:
        return []
    if isinstance(x, (list, tuple, set)):
        return list(x)
    return [x]

def _unit(v: np.ndarray) -> np.ndarray:
    """Safe unit vector.
    Takes a 3D vector and get its norm (length) - unit vector"""
    n = np.linalg.norm(v)
    return v / (n + 1e-12) #prevents division by 0

def _triangle_normal(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray: # maybe change this to consider the central node instead of the first node!!!!!!!!
    """
    Takes the node coordinates (Ca for instance) and return the normal of the plane - a vector perpendicular to the plane - and then normalize it - the unit vector
    """
    return _unit(np.cross(b - a, c - a))

def _sidechain_dirs(ca_list: Sequence[np.ndarray],
                    cb_list: Sequence[np.ndarray]) -> list[np.ndarray]:
    """
    Unit Cα→Cβ vectors per residue.

    NOTE: If you have GLY (no Cβ):
      - skip the triad fir now!!!!!
    """
    out = []
    for ca, cb in zip(ca_list, cb_list):
        out.append(_unit(cb - ca))
    return out

def triad_chirality_with_cb(
    ca_a: np.ndarray, ca_b: np.ndarray, ca_c: np.ndarray,
    cb_a: np.ndarray, cb_b: np.ndarray, cb_c: np.ndarray,
    *,
    weights: tuple[float, float, float] | None = None,
    outward_normal: np.ndarray | None = None,
    majority_only: bool = True,
) -> dict[str, Any]:
    """
    Compute a pose-invariant, mirror-variant chirality for a single triad using Cα + Cβ.

    Returns a dict with:
      - chi     : int in {-1,+1}, the handedness bit (flips under mirror, invariant to pose)
      - score   : float in [-1,1], margin = n · s̄  (near 0 ⇒ ambiguous)
      - kappa   : float in (0,1], side-chain consistency (higher = better)
      - n       : np.ndarray(3,), triangle unit normal (from Cα only)
      - sbar    : np.ndarray(3,), final averaged side-chain unit direction
      - details : dict with per-residue signs, which residues contributed, etc.

    Arguments:
      weights        : optional per-residue weights (e.g., RSA, 1/B-factor). If None, all 1.0.
      outward_normal : if provided (pose-invariant local "outward" direction), each side-chain
                       is first oriented outward: s_i <- sign(outward_normal·s_i) * s_i
      majority_only  : if True (default), use only the majority set of side-chains according to
                       sign(n·s_i). This is robust when you have "two outward, one inward".
                       If False, all three side-chains are averaged.
    """
    # Triangle normal from Cα only
    n = _triangle_normal(ca_a, ca_b, ca_c)

    # Side-chain unit directions
    sA, sB, sC = _sidechain_dirs([ca_a, ca_b, ca_c], [cb_a, cb_b, cb_c])
    S = [sA, sB, sC]

    # Optional: orient side-chains outward first - we are not using this
    if outward_normal is not None:
        u = _unit(outward_normal)
        S = [np.sign(np.dot(u, s)) * s for s in S]

    # Per-residue signs relative to triangle normal
    sigmas = [int(np.sign(np.dot(n, s))) or +1 for s in S]  # treat exact 0 as +1
    maj_sign = +1 if (sigmas.count(+1) >= 2) else -1

    # Which residues to keep in the average?
    if majority_only:
        idx_keep = [i for i, sg in enumerate(sigmas) if sg == maj_sign]  # typically 2 or 3
    else:
        idx_keep = [0, 1, 2]

    # Weights
    W = [1.0, 1.0, 1.0] if weights is None else list(weights)
    sum_w = sum(W[i] for i in idx_keep)

    # Weighted sum of selected side-chains
    vec = np.zeros(3)
    for i in idx_keep:
        vec += W[i] * S[i]

    # Final direction and consistency
    if np.linalg.norm(vec) > 0:
        sbar = _unit(vec)
        kappa = float(np.linalg.norm(vec) / (sum_w + 1e-12))
    else:
        # Degenerate: should not happen unless all S are zero/nan
        sbar = np.array([1.0, 0.0, 0.0])
        kappa = 0.0

    # Margin & bit
    score = float(np.dot(n, sbar))        # ∈ [-1, 1]
    chi   = int(np.sign(score)) or +1     # use +1 when exactly 0

    return {
        "chi": chi,
        "score": score,
        "kappa": kappa,
        "n": n,
        "sbar": sbar,
        "details": {
            "sigmas": sigmas,                 # per-residue sign(n·s_i)
            "kept_indices": idx_keep,         # which residues contributed to s̄
            "maj_sign": maj_sign,
            "weights_used": [W[i] for i in idx_keep],
        },
    }


def find_triads(graph_data, classes, config, checks, protein_index, tracker: ResidueTracker | None = None):
    # print(f"protein_index: {protein_index}")
    ctx = TrackCtx(run_id=config.get("run_id", "default"), stage="triads", protein_i=protein_index)

    G = graph_data["graph"]
    rsa = graph_data["rsa"]
    contact_map = graph_data["contact_map"]
    residue_map = graph_data["residue_map_all"]

    triads = {}

    if "residues" in classes.keys():
        residue_classes = {
            res: class_name
            for class_name, residues in classes["residues"].items()
            for res in residues
        }
    else:
        residue_classes = None

    rsa_classes = classes["rsa"] if "rsa" in classes.keys() else None
    distance_classes = classes["distance"] if "distance" in classes.keys() else None

    n_nodes = len(G.nodes())
    n_edges = len(G.edges())
    for center in G.nodes():
        neighbors = {n for n in G.neighbors(center) if n != center}
        for u, w in combinations(neighbors, 2):
            chains = config.get("filter_triads_by_chain", None)
            u_chain = u.split(":")[0]
            center_chain = center.split(":")[0]
            w_chain = w.split(":")[0]

            if isinstance(chains, list) and len(chains) > 0:
                if u_chain not in chains and center_chain not in chains and w_chain not in chains:
                    log.debug(f"({u}, {center}, {w}) was filtered.")
                    continue
            elif isinstance(chains, str) and chains.strip() != "":
                if u_chain != chains and center_chain != chains and w_chain != chains:
                    log.debug(f"({u}, {center}, {w}) was filtered.")
                    continue

            if u[:5] == w[:5]:
                _u_index, _center_index, _w_index = residue_map[residue_to_tuple(u)], residue_map[residue_to_tuple(center)], residue_map[residue_to_tuple(w)]
                _du = contact_map[_u_index, _center_index]
                _dw = contact_map[_center_index, _w_index]

                if _du <= _dw:
                    outer_sorted = tuple([u, w])
                else:
                    outer_sorted = tuple([w, u])
            else:
                outer_sorted = tuple(sorted([u, w]))

            u_split, center_split, w_split = outer_sorted[0].split(":"), center.split(":"), outer_sorted[1].split(":")
            u_res, center_res, w_res = u_split[1], center_split[1], w_split[1]
            u_tuple, center_tuple, w_tuple = residue_to_tuple(outer_sorted[0]), residue_to_tuple(center), residue_to_tuple(outer_sorted[1])
            u_index, center_index, w_index = residue_map[u_tuple], residue_map[center_tuple], residue_map[w_tuple]

            if residue_classes is not None:
                u_res_class, center_res_class, w_res_class = residue_classes[u_res], residue_classes[center_res], residue_classes[w_res]
            else:
                u_res_class, center_res_class, w_res_class = u_res, center_res, w_res


            u_resChain = str(u_split[2]) + u_split[0]
            center_resChain = str(center_split[2]) + center_split[0]
            w_resChain = str(w_split[2]) + w_split[0]

            d1 = contact_map[u_index, center_index]
            d2 = contact_map[u_index, w_index]
            d3 = contact_map[center_index, w_index]

            u_ca = np.array(G.nodes[outer_sorted[0]].get("ca_coord", (np.nan, np.nan, np.nan)))
            w_ca = np.array(G.nodes[outer_sorted[1]].get("ca_coord", (np.nan, np.nan, np.nan)))
            center_ca = np.array(G.nodes[center].get("ca_coord", (np.nan, np.nan, np.nan)))

            u_cb = np.array(G.nodes[outer_sorted[0]].get("cb_coord", (np.nan, np.nan, np.nan)))
            w_cb = np.array(G.nodes[outer_sorted[1]].get("cb_coord", (np.nan, np.nan, np.nan)))
            center_cb = np.array(G.nodes[center].get("cb_coord", (np.nan, np.nan, np.nan)))


            if not any(np.all(np.isnan(x)) for x in [u_ca, center_ca, w_ca, u_cb, center_cb, w_cb]):
                chi = triad_chirality_with_cb(u_ca, center_ca, w_ca, u_cb, center_cb, w_cb, majority_only=True)["chi"]
            else:
                chi = 0
                print(f"Nan: {u_cb, center_cb, w_cb}")
                print(u, center, w)

            rsa1 = rsa[outer_sorted[0]]*100
            rsa2 = rsa[center]*100
            rsa3 = rsa[outer_sorted[1]]*100

            def _rsa_opts(val):
                if val is None:
                    return [None]
                try:
                    if np.isnan(val):
                        return [None]
                except TypeError:
                    pass

                if rsa_classes is not None:
                    return _as_list(find_class(rsa_classes, val))
                return _as_list(
                    value_to_class(
                        val,
                        config["rsa_bin_width"],
                        config["rsa_filter"] * 100,
                        inverse=True,
                        close_tolerance=config["close_tolerance_rsa"],
                    )
                )

            if checks["rsa"]:
                rsa1_opts = _rsa_opts(rsa1)
                rsa2_opts = _rsa_opts(rsa2)
                rsa3_opts = _rsa_opts(rsa3)
            else:
                rsa1_opts, rsa2_opts, rsa3_opts = [None], [None], [None]

            if distance_classes is not None:
                d1_opts = _as_list(find_class(distance_classes, d1))
                d2_opts = _as_list(find_class(distance_classes, d2))
                d3_opts = _as_list(find_class(distance_classes, d3))
            else:
                d1_opts = _as_list(
                    value_to_class(
                        d1,
                        config["distance_bin_width"],
                        config["edge_threshold"],
                        close_tolerance=config["close_tolerance"],
                    )
                )
                d2_opts = _as_list(
                    value_to_class(
                        d2,
                        config["distance_bin_width"],
                        2 * config["edge_threshold"],
                        close_tolerance=config["close_tolerance"],
                    )
                )
                d3_opts = _as_list(
                    value_to_class(
                        d3,
                        config["distance_bin_width"],
                        config["edge_threshold"],
                        close_tolerance=config["close_tolerance"],
                    )
                )

            if d1_opts and d2_opts and d3_opts:
                for d1_c, d2_c, d3_c, rsa1_class, rsa2_class, rsa3_class in product(d1_opts, d2_opts, d3_opts, rsa1_opts, rsa2_opts, rsa3_opts):
                    full_describer = (chi, rsa1_class, rsa2_class, rsa3_class, d1_c, d2_c, d3_c)
                    full_describer_absolute = (chi, rsa1_class, rsa2_class, rsa3_class, d1, d2, d3)

                    triad_class = [u_res_class, center_res_class, w_res_class]
                    triad_abs = [outer_sorted[0], center, outer_sorted[1]]
                    triad_token = (*triad_class, *full_describer)
                    triad_full = (*triad_abs, *full_describer)
                    triad_absolute = (*triad_abs, *full_describer_absolute)

                    if triad_token not in triads:
                        triads[triad_token] = {
                            "count": 1,
                            "triads_full": [triad_full],
                            "triads_absolute": [triad_absolute]
                        }
                    else:
                        triads[triad_token]["count"] += 1
                        triads[triad_token]["triads_full"].append(triad_full)
                        triads[triad_token]["triads_absolute"].append(triad_absolute)
    n_triad = 0
    counters = {}
    for triad, data in triads.items():
        n_triad += triads[triad]["count"]

        triads_full = data["triads_absolute"]
        data["triads_absolute_d1"] = sorted(triads_full, key=lambda r: r[-3])

        if triads[triad]["count"] not in counters:
            counters[triads[triad]["count"]] = 1
        else:
            counters[triads[triad]["count"]] += 1

        if tracker is not None:
            tracker.triads_built(ctx=ctx, token=triad, triads_absolute=data["triads_absolute"])


    logging.info(f"N Nodes: {n_nodes} | N Edges: {n_edges} | N Triad: {n_triad} | Unique Triad: {len(triads.keys())}")
    logging.debug(f"Counters: {counters}")

    return triads

def create_residues_classes(path, residues_similarity_cutoff):

    atchley_factors = pd.read_csv(path, index_col = 0)

    sim = cosine_similarity(atchley_factors.values)
    aas = [convert_1aa3aa(aa) for aa in atchley_factors.index.tolist()]

    parent = {aa: aa for aa in aas}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        root_x, root_y = find(x), find(y)
        if root_x != root_y:
            parent[root_y] = root_x

    for i in range(len(aas)):
        for j in range(i+1, len(aas)):
            if sim[i, j] >= residues_similarity_cutoff:
                union(aas[i], aas[j])

    root_to_class = {}

    residue_to_class = {}
    class_id = 1
    for aa in aas:
        root = find(aa)
        if root not in root_to_class:
            root_to_class[root] = f"C{class_id}"
            class_id += 1
        residue_to_class[aa] = root_to_class[root]

    return residue_to_class

ROLE_TO_PAIR_IDX = {
    "UC": (0, 1),  # (u, c)
    "UW": (0, 2),  # (u, w)
    "CW": (1, 2),  # (c, w)
}

def _build_index(abs_list):
    """
    abs_list: lista de tuplas (r1, r2, r3, d1, d2, d3)
    Retorna estruturas de índice por d1/d2/d3 e vetor de registros.
    """
    records = list(abs_list)
    by_d1 = sorted((rec[-3], i) for i, rec in enumerate(records))  # (d1, id)
    by_d2 = sorted((rec[-2], i) for i, rec in enumerate(records))  # (d2, id)
    by_d3 = sorted((rec[-1], i) for i, rec in enumerate(records))  # (d3, id)

    keys_d1 = [k for k, _ in by_d1]
    keys_d2 = [k for k, _ in by_d2]
    keys_d3 = [k for k, _ in by_d3]

    return {
        "records": records,
        "by_d1": by_d1, "keys_d1": keys_d1,
        "by_d2": by_d2, "keys_d2": keys_d2,
        "by_d3": by_d3, "keys_d3": keys_d3,
    }


def _range_ids(keys, sorted_pairs, lo, hi):
    """
    sorted_pairs: lista ordenada de (chave, id)
    retorna conjunto de ids com chave em [lo, hi]
    """
    L = bisect.bisect_left(keys, lo)
    R = bisect.bisect_right(keys, hi)
    # coletar ids só do range
    return {sorted_pairs[j][1] for j in range(L, R)}

def cross_protein_triads(step_idx, chunk_idx, triads_per_protein, diff, check_distances=True):
    common_tokens = set.intersection(*(set(d.keys()) for d in triads_per_protein))
    cross = dict()

    for token in common_tokens:
        if check_distances:
            triads_len = [len(prot[token].get("triads_full", [])) for prot in triads_per_protein]
        else:
            triads_len = [len(prot[token].get("triads_full", [])) for prot in triads_per_protein]
        prod_triads = math.prod(triads_len) if triads_len else 0
        log.debug(f"Token: {token} | Qtd in each protein: {triads_len} | Total teoric combinations: {prod_triads}")


    # Função auxiliar para extrair/calcular limites da tríade original (passo 1)
    def _calc_initial_bounds(triad_abs):
        # triad_abs tem d1, d2, d3 nas posições -3, -2, -1
        d1, d2, d3 = triad_abs[-3], triad_abs[-2], triad_abs[-1]
        return (d1, d1, d2, d2, d3, d3) # min_d1, max_d1, min_d2, max_d2, min_d3, max_d3

    def _ensure_combo(item):
            if isinstance(item[0], str) and ":" in item[0]: return (item,)
            return item

    if check_distances:
        for token in common_tokens:
            prot_data = []
            valid_token = True

            for prot in triads_per_protein:
                abs_list = prot[token].get("triads_absolute", [])
                full_list = prot[token].get("triads_full", [])

                if not abs_list:
                    valid_token = False
                    break

                # Se os bounds já foram calculados num passo anterior, usamos. Senão, calculamos.
                if "bounds" in prot[token]:
                    bounds_arr = np.array(prot[token]["bounds"], dtype=np.float32)
                else:
                    bounds_arr = np.array([_calc_initial_bounds(t) for t in abs_list], dtype=np.float32)

                prot_data.append({
                    "abs": abs_list,
                    "full": full_list,
                    "bounds": bounds_arr
                })

            if not valid_token:
                continue

            global_mins = np.array([p["bounds"][:, [0, 2, 4]].min(axis=0) for p in prot_data])
            global_maxs = np.array([p["bounds"][:, [1, 3, 5]].max(axis=0) for p in prot_data])

            max_of_maxs = global_maxs.max(axis=0)
            min_of_mins = global_mins.min(axis=0)
            if np.any((max_of_maxs - min_of_mins) > diff):
                 pass

            ref_p = min(range(len(prot_data)), key=lambda p: len(prot_data[p]["abs"]))
            other_ps = [p for p in range(len(prot_data)) if p != ref_p]

            ref_data = prot_data[ref_p]

            combos_full = []
            combos_abs = []
            combos_bounds = []

            for i, ref_bound in enumerate(ref_data["bounds"]):
                ref_min1, ref_max1, ref_min2, ref_max2, ref_min3, ref_max3 = ref_bound

                candidates_by_p = []
                ok = True

                for p in other_ps:
                    p_bounds = prot_data[p]["bounds"]

                    # cand_max <= ref_min + diff  E  cand_min >= ref_max - diff
                    mask = (p_bounds[:, 1] <= ref_min1 + diff) & (p_bounds[:, 0] >= ref_max1 - diff) & \
                           (p_bounds[:, 3] <= ref_min2 + diff) & (p_bounds[:, 2] >= ref_max2 - diff) & \
                           (p_bounds[:, 5] <= ref_min3 + diff) & (p_bounds[:, 4] >= ref_max3 - diff)

                    valid_indices = np.where(mask)[0]

                    if len(valid_indices) == 0:
                        ok = False
                        break

                    # Armazena os índices válidos para esta proteína
                    candidates_by_p.append(valid_indices)

                if not ok:
                    continue

                # 5. Formação dos combos finais e merge dos Limites em O(1)
                ref_abs_val = ref_data["abs"][i]
                ref_full_val = ref_data["full"][i]

                for tail_indices in product(*candidates_by_p):
                    ordered_full = [None] * len(prot_data)
                    ordered_abs = [None] * len(prot_data)
                    ordered_bounds = [None] * len(prot_data)

                    ordered_full[ref_p] = ref_full_val
                    ordered_abs[ref_p] = ref_abs_val
                    ordered_bounds[ref_p] = ref_bound

                    for k, p_idx in enumerate(other_ps):
                        # tail_indices[k] é o índice escolhido para a proteína p_idx
                        idx_in_prot = tail_indices[k]

                        ordered_full[p_idx] = prot_data[p_idx]["full"][idx_in_prot]
                        ordered_abs[p_idx] = prot_data[p_idx]["abs"][idx_in_prot]
                        ordered_bounds[p_idx] = prot_data[p_idx]["bounds"][idx_in_prot]

                    # Validação Final de Distância (usando a lista já ordenada, embora p/ math tanto faz)
                    all_b_arr = np.array(ordered_bounds)
                    new_mins = all_b_arr[:, [0, 2, 4]].min(axis=0)
                    new_maxs = all_b_arr[:, [1, 3, 5]].max(axis=0)

                    if np.any((new_maxs - new_mins) > diff):
                        continue

                    new_bound = (new_mins[0], new_maxs[0], new_mins[1], new_maxs[1], new_mins[2], new_maxs[2])

                    # Flatten usando a lista ORDENADA
                    combo_abs = tuple(chain.from_iterable(_ensure_combo(x) for x in ordered_abs))
                    combo_full = tuple(chain.from_iterable(_ensure_combo(x) for x in ordered_full))

                    combos_abs.append(combo_abs)
                    combos_full.append(combo_full)
                    combos_bounds.append(new_bound)
            if not combos_full:
                continue

            cross[token] = {
                "count": len(combos_full),
                "triads_full": combos_full,
                "triads_absolute": combos_abs,
                "triads_absolute_d1": combos_abs,
                "bounds": combos_bounds # <-- CACHE SALVO PARA O PRÓXIMO PASSO
            }

    return cross

# def cross_protein_triads_backup(step_idx, chunk_idx, triads_per_protein, diff, check_distances=True):
#     """
#     triads_per_protein: list of dictionaries per protein.
#       For a token t:
#         triads_per_protein[p][t]["triads_absolute_d1"]  -> list (r1,r2,r3,d1,d2,d3)
#     diff: tolerance for the difference of distances
#     """
#     common_tokens = set.intersection(*(set(d.keys()) for d in triads_per_protein))
#     cross = dict()

#     # ---------- helpers ----------
#     def _index_by_abs_d1(triads_abs_d1_list):
#         """
#         triads_abs_d1_list: lista de triads_absolute já ordenada por d1
#         retorna índice -> posição na lista original (a própria ordem serve)
#         """
#         return list(range(len(triads_abs_d1_list)))

#     def _build_sorted_key_list(records, k):
#         """
#         records: lista de triads_absolute
#         k: -3, -2, -1 para d1,d2,d3
#         retorna keys e pares (key, idx)
#         """
#         pairs = sorted((rec[k], i) for i, rec in enumerate(records))
#         keys = [x for x, _ in pairs]
#         return keys, pairs

#     def _range_ids(keys, pairs, lo, hi):
#         L = bisect.bisect_left(keys, lo)
#         R = bisect.bisect_right(keys, hi)
#         return {pairs[j][1] for j in range(L, R)}

#     for token in common_tokens:
#         if check_distances:
#             triads_len = [len(prot[token].get("triads_full", [])) for prot in triads_per_protein]
#         else:
#             triads_len = [len(prot[token].get("triads_full", [])) for prot in triads_per_protein]
#         prod_triads = math.prod(triads_len) if triads_len else 0
#         log.debug(f"Token: {token} | Qtd in each protein: {triads_len} | Total teoric combinations: {prod_triads}")

#     for token in common_tokens:
#         if check_distances:
#             # obter listas absolutas por proteína para este token
#             abs_d1_lists = [prot[token].get("triads_absolute_d1", []) for prot in triads_per_protein]
#             if any(len(lst) == 0 for lst in abs_d1_lists):
#                 continue

#             full_lists = [prot[token].get("triads_full", []) for prot in triads_per_protein]
#             abs_lists  = [prot[token].get("triads_absolute", []) for prot in triads_per_protein]

#             look_full = []
#             look_abs  = []
#             for p in range(len(triads_per_protein)):
#                 m_full = {}
#                 for t_full in full_lists[p]:
#                     # t_full não tem d1,d2,d3 absolutos; não dá pra mapear por ele sozinho.
#                     # Então vamos mapear pelo prefixo (u,c,w, chi, rsa1,rsa2,rsa3) e deixar resolver com abs via triads_absolute
#                     key = t_full[:-3]
#                     m_full.setdefault(key, []).append(t_full)
#                 look_full.append(m_full)

#                 m_abs = {}
#                 for t_abs in abs_lists[p]:
#                     key = t_abs[:-3]
#                     m_abs.setdefault(key, []).append(t_abs)
#                 look_abs.append(m_abs)

#             idxs = []
#             for recs in abs_d1_lists:
#                 keys_d1, by_d1 = _build_sorted_key_list(recs, -3)
#                 keys_d2, by_d2 = _build_sorted_key_list(recs, -2)
#                 keys_d3, by_d3 = _build_sorted_key_list(recs, -1)
#                 idxs.append({
#                     "records": recs,
#                     "keys_d1": keys_d1, "by_d1": by_d1,
#                     "keys_d2": keys_d2, "by_d2": by_d2,
#                     "keys_d3": keys_d3, "by_d3": by_d3,
#                 })

#             ref_p = min(range(len(idxs)), key=lambda p: len(idxs[p]["records"]))
#             other_ps = [p for p in range(len(idxs)) if p != ref_p]
#             ref_records = idxs[ref_p]["records"]

#             combos_full = []
#             combos_abs  = []
#             combos_abs_d1 = []

#             w1 = w2 = w3 = diff

#             for ref_rec in ref_records:
#                 d1_ref, d2_ref, d3_ref = ref_rec[-3], ref_rec[-2], ref_rec[-1]

#                 candidates_by_p = []
#                 ok = True
#                 for p in other_ps:
#                     idx = idxs[p]
#                     cand1 = _range_ids(idx["keys_d1"], idx["by_d1"], d1_ref - w1, d1_ref + w1)
#                     if not cand1: ok = False; break
#                     cand2 = _range_ids(idx["keys_d2"], idx["by_d2"], d2_ref - w2, d2_ref + w2)
#                     if not cand2: ok = False; break
#                     cand3 = _range_ids(idx["keys_d3"], idx["by_d3"], d3_ref - w3, d3_ref + w3)
#                     if not cand3: ok = False; break

#                     cands = cand1 & cand2 & cand3
#                     if not cands: ok = False; break

#                     candidates_by_p.append([idx["records"][i] for i in cands])

#                 if not ok:
#                     continue

#                 for tail in product(*candidates_by_p):
#                     # montar combo_abs_d1 na ordem original das proteínas
#                     tmp_abs = [None] * len(idxs)
#                     tmp_abs[ref_p] = ref_rec
#                     k = 0
#                     for p in other_ps:
#                         tmp_abs[p] = tail[k]
#                         k += 1
#                     combo_abs_d1 = tuple(tmp_abs)

#                     combo_abs = combo_abs_d1  # mesmo tipo e shape

#                     tmp_full = []
#                     for p, t_abs in enumerate(combo_abs):
#                         key = t_abs[:-3]
#                         full_cands = look_full[p].get(key, [])
#                         if not full_cands:
#                             tmp_full = None
#                             break
#                         tmp_full.append(full_cands[0])
#                     if tmp_full is None:
#                         continue

#                     combos_abs_d1.append(combo_abs_d1)
#                     combos_abs.append(combo_abs)
#                     combos_full.append(tuple(tmp_full))

#             if not combos_full:
#                 continue

#             cross[token] = {
#                 "count": len(combos_full),
#                 "triads_full": combos_full,
#                 "triads_absolute": combos_abs,
#                 "triads_absolute_d1": combos_abs_d1,
#             }
#         else:
#             def _is_triad(x):
#                 # triad: (U, C, W, ..., d1,d2,d3) então x[0] é string tipo "A:GLU:19"
#                 return isinstance(x, tuple) and len(x) >= 3 and isinstance(x[0], str) and ":" in x[0]

#             def _is_combo(x):
#                 # combo: tuple of triads, então x[0] é uma triad
#                 return isinstance(x, tuple) and len(x) > 0 and _is_triad(x[0])

#             def _ensure_combo(item):
#                 # item pode ser triad ou combo
#                 if _is_combo(item):
#                     return item
#                 if _is_triad(item):
#                     return (item,)
#                 raise TypeError(f"Expected triad or combo, got {type(item)}: {item}")

#             lists_full = [[_ensure_combo(c) for c in prot[token].get("triads_full", [])] for prot in triads_per_protein]
#             lists_abs  = [[_ensure_combo(c) for c in prot[token].get("triads_absolute", [])] for prot in triads_per_protein]

#             if any(len(lst) == 0 for lst in lists_full) or any(len(lst) == 0 for lst in lists_abs):
#                 continue

#             combos_full = []
#             combos_abs  = []

#             for inner_full, inner_abs in zip(product(*lists_full), product(*lists_abs)):
#                 flat_full = tuple(chain.from_iterable(inner_full))
#                 flat_abs  = tuple(chain.from_iterable(inner_abs))

#                 combos_full.append(flat_full)
#                 combos_abs.append(flat_abs)


#             if not combos_full:
#                 continue

#             cross[token] = {
#                 "count": len(combos_full),
#                 "triads_full": combos_full,
#                 "triads_absolute": combos_abs,
#                 "triads_absolute_d1": combos_abs
#             }

#     return cross

def build_graph_from_cross_combos(cross_combos) -> set[tuple[tuple[str, ...], tuple[str, ...]]]:
    """
    Para cada combo aprovado, cria as arestas (C,U) e (C,W) nos nós alinhados.
    Retorna um set de arestas entre nós-tupla (um rótulo por proteína).
    """
    edges: set[tuple[tuple[str, ...], tuple[str, ...]]] = set()
    for _, data in cross_combos.items():
        combos = data.get("triads_full", [])
        for combo in combos:
            try:
                U = tuple(tri[0] for tri in combo)
            except Exception as e:
                raise Exception(f"{e}\nCombo: {combo}")
            C = tuple(tri[1] for tri in combo)
            W = tuple(tri[2] for tri in combo)
            edge = (C, U) if C <= U else (U, C)
            edges.add(edge)

            edge = (C, W) if C <= W else (W, C)
            edges.add(edge)

    return edges

def rebuild_cross_combos(cross_combos: dict[dict, list[tuple[tuple, ...]]], graph_nodes):
    new = {}
    graph_nodes = set(graph_nodes)

    for token, data in cross_combos.items():
        # Listas originais
        orig_full = data.get("triads_full", [])
        orig_abs = data.get("triads_absolute", [])
        orig_abs_d1 = data.get("triads_absolute_d1", [])
        orig_bounds = data.get("bounds", []) # Recupera os bounds calculados

        # Listas filtradas
        keep_full = []
        keep_abs = []
        keep_abs_d1 = []
        keep_bounds = []

        # Itera pelo índice para manter sincronia
        for i, combo in enumerate(orig_full):
            # Lógica de filtro: Verifica se os nós do combo ainda existem no grafo
            U = tuple(tri[0] for tri in combo)
            C = tuple(tri[1] for tri in combo)
            W = tuple(tri[2] for tri in combo)

            if {U, C, W} <= graph_nodes:
                keep_full.append(combo)

                # Sincroniza as outras listas se existirem
                if i < len(orig_abs):
                    keep_abs.append(orig_abs[i])
                if i < len(orig_abs_d1):
                    keep_abs_d1.append(orig_abs_d1[i])
                if i < len(orig_bounds):
                    keep_bounds.append(orig_bounds[i])

        if keep_full:
            new[token] = {
                "count": len(keep_full),
                "triads_full": keep_full,
                "triads_absolute": keep_abs,
                "triads_absolute_d1": keep_abs_d1,
                "bounds": keep_bounds  # <--- CRUCIAL: Passa os bounds para o próximo passo
            }

    return new

def parse_node(node: str) -> tuple[str, str, int]:
    chain, res, num = node.split(":")
    try:
        num_i = int(num)
    except ValueError:
        import re
        m = re.match(r"(-?\d+)", num)
        if not m:
            raise
        num_i = int(m.group(1))
    return chain, res, num_i

def _extract_chain_from_idx(idx, maps):
    res_tuple = maps["residue_maps_unique"].get(idx, None)
    if res_tuple is None:
        return None
    return res_tuple[0]  # cadeia

def _build_threshold_matrix(nodes, maps, threshold_cfg):
    K = len(nodes)
    dim = len(nodes[0])

    if not isinstance(threshold_cfg, dict):
        return np.full((K, K), float(threshold_cfg), dtype=float)

    default_thr = float(threshold_cfg.get("default", 2.0))
    mix_thr = float(threshold_cfg.get("mix", default_thr))

    # chains[p, i] = cadeia do nó i na proteína p
    chains = np.array(
        [[_extract_chain_from_idx(nodes[i][p], maps) for i in range(K)] for p in range(dim)],
        dtype=object
    )

    # all_eq[i, j] = True se para todo p chains[p,i] == chains[p,j]
    all_eq = np.all(chains[:, :, None] == chains[:, None, :], axis=0)

    # threshold por cadeia (usa cadeia da proteína 0 como referência de linha)
    row_chain = chains[0, :]
    row_thr = np.array([float(threshold_cfg.get(c, default_thr)) for c in row_chain], dtype=float)

    T = np.full((K, K), mix_thr, dtype=float)
    T[all_eq] = row_thr[:, None][all_eq]

    return T

def create_coherent_matrices(nodes, matrices: dict, maps: dict, threshold: float | dict = 3.0):
    # Verify the impact of std in distances. Suppose that 9 proteins have a simillar distance but one has a higher distance, how much does it impact?
    dim = len(nodes[0])
    K = len(nodes)

    maps_out = {}
    maps_out["possible_nodes"] = {}
    for i, node in enumerate(nodes):
        maps_out["possible_nodes"][i]      = node
        maps_out["possible_nodes"][str(node)] = i
    maps["possible_nodes"] = maps_out["possible_nodes"]

    stacked_induced  = np.empty((dim, K, K))
    stacked_adjacent  = np.empty_like(stacked_induced)

    for p in range(dim):
        idx = [node[p] for node in nodes]

        for i in range(K):
            stacked_induced[p,  i, :] = matrices["dm_induced"][idx[i],  idx]
            stacked_adjacent[p, i, :]  = matrices["dm_adjacent"][idx[i], idx]

    mask_invalid_induced = np.any((stacked_induced == 0) | np.isnan(stacked_induced), axis=0)
    mask_invalid_adjacent = np.any((stacked_adjacent == 0) | np.isnan(stacked_adjacent), axis=0)

    range_induced = np.max(stacked_induced, axis=0) - np.min(stacked_induced, axis=0)
    range_adjacent = np.max(stacked_adjacent, axis=0) - np.min(stacked_adjacent, axis=0)

    T = _build_threshold_matrix(nodes, maps, threshold)

    mask_valid_ind = (range_induced <= T) & (~mask_invalid_induced)
    final_induced = np.full((K, K), np.nan)
    final_induced[mask_valid_ind] = 1.0

    mask_valid_adj = (range_adjacent <= T) & (~mask_invalid_adjacent)
    final_adjacent = np.full((K, K), np.nan)
    final_adjacent[mask_valid_adj] = 1.0

    new_matrices = {
        "coherent_global_nodes": final_induced,
        "coherent_adjacent_nodes": final_adjacent
    }


    return new_matrices, maps

def get_memory_usage_mb():
    """
    Retorna uso de memória RSS em MB, se psutil estiver disponível.
    Caso contrário, retorna None.
    """
    try:
        import psutil
    except ImportError:
        return None

    process = psutil.Process(os.getpid())
    return process.memory_info().rss / (1024 ** 2)

def write_json_raw(path, data):
    """
    Salva dados em JSON bruto, indentado.
    """
    with open(path, "w") as f:
        json.dump(data, f, indent=2)

def execute_step(
    step_idx: int,
    graph_collection,
    max_chunks: int,
    current_filtered_cross_combos,
    graphs_data,
    global_state,
    residue_tracker
):

    profiling = global_state.get("profiling", None)
    debug = bool(profiling and profiling.get("debug", False))

    if step_idx == 1:
        source = graph_collection["triads"]
    else:
        if not current_filtered_cross_combos:
            raise RuntimeError("current_filtered_cross_combos is empty")
        source = current_filtered_cross_combos

    num_items = len(source)
    chunks = [source[i:i + max_chunks] for i in range(0, num_items, max_chunks)]

    new_filtered_cross_combos = []
    all_step_graphs = []

    step_mem_before = get_memory_usage_mb()
    step_time_start = time.perf_counter()

    step_profile = {
        "step": step_idx,
        "num_items": num_items,
        "num_chunks": len(chunks),
        "chunks": [],       # lista de dicts, um por chunk
        "time_sec": None,   # preenchido no final
        "mem_before_mb": None,
        "mem_after_mb": None,
    }

    log.debug(f"[DEBUG] Step {step_idx} started. Memory: {step_mem_before:.2f} MB")

    for chunk_idx, chunk_triads in enumerate(chunks):
        chunk_time_start = time.perf_counter()
        chunk_mem_before = get_memory_usage_mb()

        log.debug(f"[DEBUG]   Chunk {chunk_idx}: started. "
                  f"Memory before: {chunk_mem_before:.2f} MB. "
                  f"Chunk size: {len(chunk_triads)}")

        rebuilt_combos, final_graphs = process_chunk(
            step_idx=step_idx,
            chunk_idx=chunk_idx,
            chunk_triads=chunk_triads,
            graphs_data=graphs_data,
            global_state=global_state,
            residue_tracker=residue_tracker
        )

        if isinstance(rebuilt_combos, list):
            new_filtered_cross_combos.extend(rebuilt_combos)
        elif rebuilt_combos is not None:
            new_filtered_cross_combos.append(rebuilt_combos)

        all_step_graphs.extend(final_graphs)

        chunk_time_end = time.perf_counter()
        chunk_mem_after = get_memory_usage_mb()

        log.debug(f"[DEBUG]   Chunk {chunk_idx}: finished. "
              f"Time: {chunk_time_end - chunk_time_start:.3f}s. "
              f"Memory after: {chunk_mem_after:.2f} MB.")

        chunk_profile = {
            "chunk_idx": chunk_idx,
            "num_triads_in_chunk": len(chunk_triads),
            "time_sec": chunk_time_end - chunk_time_start,
            "mem_before_mb": chunk_mem_before,
            "mem_after_mb": chunk_mem_after,
        }
        step_profile["chunks"].append(chunk_profile)


    step_time_end = time.perf_counter()
    step_mem_after = get_memory_usage_mb()


    log.debug(f"[DEBUG] Step {step_idx} finished. "
          f"Total time: {step_time_end - step_time_start:.3f}s. "
          f"Final memory: {step_mem_after:.2f} MB\n")

    step_profile["time_sec"] = step_time_end - step_time_start
    step_profile["mem_before_mb"] = step_mem_before
    step_profile["mem_after_mb"] = step_mem_after

    profiling["steps"].append(step_profile)

    return new_filtered_cross_combos, all_step_graphs

def process_chunk(step_idx, chunk_idx, chunk_triads, graphs_data, global_state, residue_tracker):
    config = global_state["config"]
    dm_adjacent = global_state["dm_adjacent"]
    inv_maps = global_state["inv_maps"]
    metadata = global_state["metadata"]
    maps = global_state["maps"]
    matrices_dict = global_state["matrices_dict"]
    steps = global_state["steps"]

    max_chunks = config["max_chunks"]
    if len(chunk_triads) == 1:
        return chunk_triads, []
    if step_idx == 1:
        log.debug(f"[step {step_idx}] Creating combos | Chunk: {chunk_idx}")
        cross_combos = cross_protein_triads(step_idx, chunk_idx, chunk_triads, config["distance_diff_threshold"])
        log.info(f"[step {step_idx}] Created {len(cross_combos)}")
    else:
        log.debug(f"[step {step_idx}] Creating combos | Chunk: {chunk_idx}")
        cross_combos = cross_protein_triads(step_idx, chunk_idx, chunk_triads, config["distance_diff_threshold"])

    if residue_tracker is not None:
        ctx = TrackCtx(run_id=config.get("run_id","default"), stage="combos", step_id=step_idx, chunk_id=chunk_idx)
        for token, combos in cross_combos.items():
            residue_tracker.combos_built(ctx=ctx, token=token, combos=combos)

    triad_graph = build_graph_from_cross_combos(cross_combos)
    tuple_edges = [tuple(edge) for edge in triad_graph]
    G = nx.Graph()
    G.add_edges_from(tuple_edges)

    triads_isolated_count = 0
    triads_isolated = []
    for comp in list(nx.connected_components(G)):
        if len(comp) == 3:
            triads_isolated_count += 1
            triads_isolated.append(comp)
            G.remove_nodes_from(comp)
    log.debug(f"[step {step_idx} chunk_idx {chunk_idx}] I remove {triads_isolated} triads isolated from the graph")
    components = list(nx.connected_components(G))

    index_offset_base = chunk_idx * (max_chunks ** step_idx)
    maps["inv_maps"] = inv_maps

    Frames_comp = set()
    final_graphs = []

    comp_id = 1
    processed_status = 1
    num_components = len(components)

    for component in components:
        len_component = len(component)
        log.debug(
            f"[{step_idx}] {comp_id}: Processing component "
            f"{processed_status}/{num_components} with {len_component} nodes"
        )
        processed_status += 1

        if len_component <= 4:
            log.debug(
                f"[{step_idx}] {comp_id} Skipping component "
                f"{processed_status-1}, because it has just {len_component} nodes"
            )
            comp_id += 1
            continue


        subG = G.subgraph(component).copy()

        dm_adjacent_graph = np.zeros((metadata["total_length"], metadata["total_length"]))
        log.debug(f"[{step_idx}] Creating dm_adjacent_graph matrices...")
        for u, v in subG.edges():
            for p, (res_u, res_v) in enumerate(zip(u, v)):
                if res_u != res_v:
                    split_res_u, split_res_v = res_u.split(":"), res_v.split(":")
                    res_u_tuple = (split_res_u[0], split_res_u[2], split_res_u[1])
                    res_v_tuple = (split_res_v[0], split_res_v[2], split_res_v[1])
                    idx_u = inv_maps[p+index_offset_base][res_u_tuple]
                    try:
                        idx_v = inv_maps[p+index_offset_base][res_v_tuple]
                    except Exception as e:
                        raise Exception(
                            f"Problem at the protein {p}, offset: {index_offset_base}, "
                            f"v tuple: {res_v_tuple}, step: {step_idx}, "
                            f"chunk_idx: {chunk_idx}, comp_id: {comp_id}\nError: {e}"
                        )
                    dm_adjacent_graph[idx_u, idx_v] = dm_adjacent[idx_u, idx_v]
                    dm_adjacent_graph[idx_v, idx_u] = dm_adjacent[idx_v, idx_u]

        save(f"comp_id_{comp_id}", "dm_adjacent_graph", dm_adjacent_graph)
        matrices_dict["dm_adjacent"] = dm_adjacent_graph

        nodes = list(subG.nodes())
        nodes_indices = []

        for node in nodes:
            node_converted = []
            for k, res in enumerate(node):
                res_split = res.split(":")
                res_tuple = (res_split[0], res_split[2], res_split[1])
                res_indice = inv_maps[k+index_offset_base][res_tuple]
                node_converted.append(res_indice)
            nodes_indices.append(node_converted)

        save(f"comp_id_{comp_id}", "nodes_indices", nodes_indices)
        log.debug(f"[{step_idx}] {comp_id} Creating the std_matrix...")
        coherent_matrices, coherent_maps = create_coherent_matrices(
            nodes=nodes_indices,
            matrices=matrices_dict,
            maps=maps,
            threshold=config["distance_std_threshold"]
        )

        log.debug(f"[{step_idx}] {comp_id} Finished creating the std_matrix.")

        save(f"comp_id_{comp_id}", "coherent_matrices", coherent_matrices)
        save(f"comp_id_{comp_id}", "coherent_maps", coherent_maps)

        debugar_ = True if len(nodes) > 300 else False


        # Tracer for DFS debugging, not useful in real application
        tracer = TraversalTracer(
            out_dir=Path(config["output_path"]) / "viz_runs",
            fmt="mp4",
            fps=12,
            sample_every=50,
            max_frames=2000,
            dpi=110,
            enabled=False
        )

        if not debugar_:
            tracer.enabled = False

        steps_end = True if step_idx == steps else False
        frames, union_graph, error = generate_frames(
            component_graph=subG,
            matrices=coherent_matrices,
            maps=coherent_maps,
            len_component=len_component,
            chunk_id=chunk_idx,
            step=step_idx,
            config=config,
            debug=debugar_,
            tracer=tracer,
            nodes=nodes,
            steps_end=steps_end,
            residue_tracker=residue_tracker
        )

        if len(frames.keys()) > 1:
            if step_idx < steps:
                Frames_comp = Frames_comp.union(union_graph["edges_residues"])
            else:
                final_graphs.append(
                    (create_graph(frames, typeEdge="edges_residues", comp_id=comp_id), comp_id)
                )
        else:
            log.debug("Component Refused")

        comp_id += 1

    if step_idx < steps:
        tuple_edges_frames = [tuple(edge) for edge in Frames_comp]
        Graph_union_frames = nx.Graph()
        Graph_union_frames.add_edges_from(tuple_edges_frames)

        rebuilt_combos = rebuild_cross_combos(cross_combos, Graph_union_frames)
    else:
        rebuilt_combos = None
        final_graphs.insert(0, ([G], 0))

    return rebuilt_combos, final_graphs


def association_product(graphs_data: list,
                        config: dict,
                        debug: bool = True) -> dict[str, list] | None:
    logger = logging.getLogger("association.association_product")

    residue_tracker = config.get("watch_residues")

    if residue_tracker is not None:
        pdb_dfs = []
        for gd in graphs_data:
            G = gd["graph"]  # nx.Graph
            df = G.graph.get("pdb_df")
            if df is None:
                df = G.graph.get("raw_pdb_df")
            pdb_dfs.append(df)

        residue_tracker.resolve_from_pdb_dfs(pdb_dfs, stage="resolve_watch_residues")

    profiling = {
        "debug": debug,
        "steps": []
    }
    checks = config.get("checks", {"rsa": True})
    classes = config.get("classes", {})
    max_chunks = config.get("max_chunks", 5)


    graph_collection = {
        "graphs": [gd["graph"] for gd in graphs_data],
        "triads": [find_triads(gd, classes, config, checks, i, residue_tracker) for i, gd in enumerate(graphs_data)],
        "contact_maps": [gd["contact_map"] for gd in graphs_data],
        "residue_maps_all": [gd["residue_map_all"] for gd in graphs_data],
        "rsa_maps": [gd["rsa"] for gd in graphs_data],
        "nodes_graphs": [sorted(list(gd["graph"].nodes())) for gd in graphs_data]
    }

    save("association_product", "graph_collection", graph_collection)
    ranges_graph = indices_graphs(graph_collection["graphs"])
    total_length = sum(len(g.nodes()) for g in graph_collection["graphs"])
    metadata = {
        "total_length": total_length,
        "ranges_graph": ranges_graph
    }

    matrices_dict = {
        "type": 0,
        "neighbors": None,
        "rsa": None,
        "identity": None,
        "associated": None,
        "similarity": None,
        "dm_adjacent": None,
        "dm_induced": None,
        "metadata": metadata
    }

    filter_input = {
        "contact_maps": graph_collection["contact_maps"],
        "rsa_maps": graph_collection["rsa_maps"],
        "residue_maps": graph_collection["residue_maps_all"],
        "nodes_graphs": graph_collection["nodes_graphs"]
    }

    logger.info("Creating pruned and thresholded arrays...")
    matrices_dict, maps = filter_maps_by_nodes(filter_input,
                                            distance_threshold=config["edge_threshold"],
                                            matrices_dict=matrices_dict)
    logger.info("Arrays created successfully!")

    current_value = 0
    maps["residue_maps_unique_break"] = {}
    for i, res_map in enumerate(maps["full_residue_maps"]):
        maps["residue_maps_unique"].update({val + current_value: key for key, val in res_map.items()})
        maps["residue_maps_unique_break"][i] = {val + current_value: key for key, val in res_map.items()}
        current_value += len(res_map)

    inv_maps = {
        k: { res: idx for idx, res in br.items() }
        for k, br in maps["residue_maps_unique_break"].items()
    }

    current_index = 0

    dm_adjacent = np.zeros((metadata["total_length"], metadata["total_length"]))
    dm_induced = np.zeros((metadata["total_length"], metadata["total_length"]))

    save("association_product", "maps", maps)
    save("association_product", "inv_maps", inv_maps)

    for i, graph in enumerate(graph_collection["graphs"]):
        graph_length = len(graph.nodes())
        new_index = current_index + graph_length
        dm_adjacent[current_index:new_index, current_index:new_index] = matrices_dict["adjacent_contact_maps"][i]
        dm_induced[current_index:new_index, current_index:new_index] = matrices_dict["induced_contact_maps"][i]
        current_index = new_index

    save("association_product", "dm_adjacent", dm_adjacent)
    save("association_product", "dm_induced", dm_induced)
    matrices_dict["dm_adjacent"] = dm_adjacent
    matrices_dict["dm_induced"] = dm_induced

    log.debug("Criando cross combos")
    start_cross = time.perf_counter()

    Graphs = []

    qtd_graphs = len(graph_collection["graphs"])
    n = qtd_graphs

    steps = 0
    while n > 1:
        n = (n + max_chunks - 1) // max_chunks
        steps += 1

    filtered_cross_combos = []
    final_graphs = []

    for step_idx in range(1, steps + 1):
        log.debug(f"Executing step {step_idx}")
        filtered_cross_combos, step_graphs = execute_step(
            step_idx,
            graph_collection,
            max_chunks,
            filtered_cross_combos,
            graphs_data=graphs_data,
            global_state={
                "dm_adjacent": dm_adjacent,
                "dm_induced": dm_induced,
                "inv_maps": inv_maps,
                "metadata": metadata,
                "maps": maps,
                "matrices_dict": matrices_dict,
                "config": config,
                "steps": steps,
                "profiling": profiling
            },
            residue_tracker=residue_tracker
        )

        if step_idx == steps:
            final_graphs.extend(step_graphs)

    log.debug("Finished creating ALL Frames")
    save("association_product", "Graphs", final_graphs)

    if debug:
        write_json_raw("profiling_report.json", profiling)
        log.debug("[DEBUG] Profiling report saved to association_product/profiling_report.json")

    return {
        "AssociatedGraph": final_graphs
    }

def generate_frames(component_graph, matrices, maps, len_component, chunk_id, step, config, debug=False, debug_every=5000, tracer: TraversalTracer | None=None, nodes=None, steps_end=False, residue_tracker: ResidueTracker | None=None):
    """
    Build frames by branching on coherent groups of the frontier.

    dm  : coherence indicator (1 = coherent, else != 1 or NaN)
    adj : adjacency indicator for edges (1 = edge, NaN on diagonal)

    A branch is accepted when the induced subgraph (by adj) over the chosen
    nodes has >= 4 edges and all nodes have degree > 1 within that subgraph.
    """


    dm_raw  = matrices["coherent_global_nodes"].copy()
    adj_raw = matrices["coherent_adjacent_nodes"].copy()
    np.fill_diagonal(dm_raw, 1)
    np.fill_diagonal(adj_raw, np.nan)

    C = (dm_raw == 1)     # Matriz de coerência na forma booleana
    A = (adj_raw == 1)    # Matriz de adjacência na forma booleana

    np.fill_diagonal(C, True)
    np.fill_diagonal(A, False)
    K = A.shape[0]

    N_adj = [np.nonzero(A[u])[0] for u in range(K)]  # Vizinhos por adjacência

    if nodes is None:
        nodes = list(component_graph.nodes())

    if len(nodes) != K:
        raise ValueError(
            f"Mismatch between nodes ({len(nodes)}) and matrices size ({K}). "
            "Check node ordering passed to create_std_matrix vs generate_frames."
        )

    node_to_idx = {node: idx for idx, node in enumerate(nodes)}
    frames = {}

    edges_base = set()
    for u, v in component_graph.edges():
        iu = node_to_idx[u]
        iv = node_to_idx[v]
        if iu == iv:
            continue
        edges_base.add(frozenset((iu, iv)))

    _, edges_idx_0, edges_res_0 = convert_edges_to_residues(edges_base, maps)
    frames[0] = {"edges_indices": edges_idx_0, "edges_residues": edges_res_0}

    def frontier_for(chosen: np.ndarray, inter_mask: np.ndarray) -> tuple:
        """
        A função frontier_for constrói a nova fronteira a partir do conjunto de nós já escolhidos.
        A ideia é a seguinte: dado um conjunto de nós escolhidos (chosen), buscamos todos os vizinhos 
        imediatos desses nós no grafo de adjacência A. Porém, precisamos garantir duas coisas:
            1) Nenhum nó que já foi escolhido pode voltar para a fronteira.
            2) Apenas nós que permanecem coerentes com o inter_mask (máscara de coerência acumulada)
            são mantidos na fronteira.
        O resultado é retornado como uma tupla de índices inteiros representando os nós que formam a nova fronteira.
        """
        if chosen.size == 0:
            # Se não há nós escolhidos, a fronteira é vazia
            return ()
        mask = A[chosen].any(axis=0)   # pega todos os vizinhos dos nós escolhidos
        mask[chosen] = False           # remove os próprios escolhidos da fronteira
        mask &= inter_mask             # mantém apenas nós ainda coerentes
        return tuple(np.nonzero(mask)[0].tolist())  # retorna índices dos vizinhos coerentes


    def valid_subgraph(nodes):
        """
        A função valid_subgraph verifica se um conjunto de nós forma um subgrafo válido.
        As condições para validade são:
            1) O subgrafo induzido deve ter pelo menos 4 arestas.
            2) Todos os nós dentro do subgrafo devem ter grau > 1 (ou seja, sem nós pendurados).
        Caso seja válido, retorna (True, es) onde 'es' é o conjunto de arestas representado por pares de nós.
        Caso contrário, retorna (False, set()).
        """
        nodes = np.asarray(sorted(nodes), dtype=np.int32)  # ordena nós para consistência
        sub = A[np.ix_(nodes, nodes)]  # submatriz de adjacência restrita aos nós escolhidos
        deg = sub.sum(axis=1)          # grau de cada nó dentro do subgrafo

        # Verifica se o número de arestas é suficiente (deg.sum() // 2 = nº de arestas)
        # if int(deg.sum() // 2) < 4:
        #     return False, set()

        if np.any(deg < 1): # Colocar para reportar grafos com mais de 1 componente
            return False, set()

        # Extrai as arestas do subgrafo: usamos apenas a parte triangular superior para evitar duplicatas
        ii, jj = np.where(np.triu(sub, 1))
        es = {frozenset((int(nodes[i]), int(nodes[j]))) for i, j in zip(ii, jj)}

        return True, es

    # Bron–Kerbosch com pivô, restrito à fronteira usando a matriz de coerência C
    def maximal_cliques(frontier_set):
        """
        Dado um conjunto de nós da fronteira, queremos dividi-los em grupos
        de nós que são todos coerentes entre si (cliques no grafo de coerência).
        Para isso usamos o algoritmo clássico de Bron–Kerbosch com pivô,
        aplicado no subgrafo de coerência induzido pela fronteira.
        """

        # Transformamos a fronteira (set) em array ordenado, para termos consistência.
        F = np.array(sorted(frontier_set), dtype=np.int32)
        if F.size == 0:
            return []

        # fmask marca quais nós fazem parte da fronteira dentro do universo total K
        fmask = np.zeros(K, dtype=bool)
        fmask[F] = True

        # N é um dicionário de vizinhança no grafo de coerência restrito à fronteira:
        # para cada nó u, guardamos apenas os vizinhos coerentes que também estão na fronteira.
        N = {}
        for u in F:
            # nu = todos os nós coerentes com u (linha de C[u])
            nu = np.flatnonzero(C[u])
            # filtramos para manter só os que também estão na fronteira
            nu = nu[fmask[nu]]
            # salvamos como conjunto, excluindo o próprio u
            N[int(u)] = set(int(x) for x in nu if x != u)

        cliques = []

        def bk(R, P, X):
            """
            Bron–Kerbosch recursivo:
            - R = conjunto atual que estamos construindo (candidato a clique)
            - P = candidatos que ainda podem entrar em R
            - X = candidatos já explorados e que não podem mais entrar
            Quando P e X ficam vazios ao mesmo tempo, R é um clique maximal.
            """
            if not P and not X:
                # Encontramos um clique maximal
                cliques.append(tuple(sorted(R)))
                return

            # Escolhemos um pivô (nó com maior vizinhança em N) para reduzir ramificações
            U = P | X
            pivot = max(U, key=lambda v: len(N[v])) if U else None
            ext = P - (N[pivot] if pivot is not None else set())

            for v in list(ext):
                # Chamamos recursivamente adicionando v ao clique atual
                bk(R | {v}, P & N[v], X & N[v])
                # Após explorar v, movemos ele de P para X (já foi tratado)
                P.remove(v)
                X.add(v)

        # Iniciamos Bron–Kerbosch com:
        # - R vazio (nenhum nó escolhido ainda),
        # - P = todos os nós da fronteira,
        # - X vazio.
        bk(set(), set(int(x) for x in F), set())

        # Ordenamos os cliques encontrados: maiores primeiro
        cliques.sort(key=lambda c: (-len(c), c))
        return cliques


    def mask_signature(mask_bool: np.ndarray) -> bytes:
        # assinatura compacta da máscara boolean
        return np.packbits(mask_bool, bitorder="little").tobytes()

    # expansão unificada: não marca visited aqui
    def expand_groups(chosen, inter_mask, frontier):
        # Achar os grupos de nós coerentes em nossa fronteira é equivalente a converter a nossa fronteira
        # Em um grafo de coerência e então procrurarmos o clique máximo nele
        cliques = maximal_cliques(frontier) if frontier else []
        if not cliques:
            cliques = tuple((v,) for v in sorted(frontier))
        groups = []
        for clc in cliques:
            inter_g = inter_mask.copy()
            for node in clc:
                inter_g &= C[node]
            chosen_g = tuple(sorted(set(chosen) | set(clc)))
            frontier_g = frontier_for(np.asarray(chosen_g, dtype=np.int32), inter_g)
            sig = (chosen_g, mask_signature(inter_g))
            groups.append((sig, chosen_g, inter_g, frontier_g))
        return groups

    # ordem de nodes por grau
    deg_all = A.sum(axis=1)
    order_all = np.argsort(-deg_all)
    len_order = len(order_all)
    # order_all = list(nx.algorithms.core.degeneracy_ordering(nx.from_numpy_array(A)))

    seen_edge_sets    = set()   # chaves de arestas aceitas
    visited_states    = set()   # (chosen_tuple, mask_signature)
    checked_node_sets = set()   # frozenset(chosen)
    next_frame_id     = 1

    pushes = pops = accepts = 0
    pruned_dupe_state = 0

    log.debug("Traversing nodes in descending order of degree ")
    for node_id, node in enumerate(order_all): # Percorremos os nós na ordem do maior grau para o menor grau
        # if tracer and tracer.enabled:
        #     # passa a matriz A do componente que você está explorando aqui; se for o grafo inteiro A, serve
        #     tracer.start_anchor(anchor=int(node), A_sub=A, nodes_map=np.arange(A.shape[0], dtype=int))
        str_node = f"{node_id}/{len_order}"
        log.debug(f"Node: {str_node}")
        log.debug(f"[{str_node}] Making copy of inter list...")
        inter = C[int(node)].copy() # Primeiro começamos pegando a linha de coerência que representa o nó e fazemos uma cópia dela (para não modificar o original)
        log.debug(f"[{str_node}] Copy finished")
        accepted = (int(node),) # Iniciamos os nós que serão percorridos com o próprio nó inicial

        neigh = N_adj[int(node)]
        if neigh.size == 0:
            continue

        # filtra vizinhos pela coerência com o node
        # inter[neigh] é bool; pega apenas índices coerentes
        frontier_idx = neigh[inter[neigh]]
        if frontier_idx.size == 0:
            continue

        log.debug(f"[{str_node}] Transforming frontier to set...")
        frontier = set(int(v) for v in frontier_idx) # Lista da fronteira
        log.debug(f"[{str_node}] Finished transforming.")

        # Vamos criar um stack para realizar um busca em profundidade pelo LIFO (Last In First Out)
        stack = []
        stack_small = []

        """
        Sabemos que dado o nó raíz, nem todos os vizinhos são coerentes entre si, então precisamos separar a fronteira em grupos
        de nós que sejam coerentes entre si. Porém, coerência NÃO é transitiva: um nó u pode ser coerente com w e com v, mas w e v
        não serem coerentes entre si. Para cobrir TODAS as possibilidades válidas sem descartar combinações, construímos um “grafo
        de coerência” sobre a fronteira (arestas indicam dm==1 entre pares de vizinhos) e buscamos cliques máximos nele — cada
        clique representa um subconjunto de vizinhos mutualmente coerentes (par a par).
        
        A função expand_groups() encapsula exatamente isso: dado o conjunto 'accepted' (os nós já escolhidos), a máscara 'inter'
        que representa os nós que permanecem coerentes com TODOS os escolhidos, e a 'frontier' atual (vizinhos no grafo de adjacência
        que ainda sobreviveram na máscara), ela:
          1) monta o grafo de coerência restrito à fronteira atual;
          2) extrai cliques máximos (subgrupos de vizinhos mutuamente coerentes);
          3) para cada clique, produz um próximo estado:
             - 'chosen_g'  = accepted ∪ clique (isto é, avançamos escolhendo este subgrupo coerente);
             - 'inter_g'   = inter multiplicado pelo dm de cada nó do clique (filtrando apenas nós coerentes com TODOS os escolhidos + clique);
             - 'frontier_g'= união dos vizinhos de todos os nós em chosen_g, restrita a quem ainda sobrevive em 'inter_g' e que não está em chosen_g.
        
        O retorno é uma lista de assinaturas 'sig' para deduplicar estados (chosen, máscara, fronteira), e as tuplas (chosen_g, inter_g, frontier_g)
        que formam os filhos do estado atual no backtracking. Dessa forma exploramos ramos por “pacotes coerentes”, sem eliminar à força toda a fronteira.
        """
        log.debug(f"[{str_node}] Creating groups...")
        groups = expand_groups(accepted, inter, frontier)
        log.debug(f"[{str_node}] Groups created: {len(groups)}")
        # print(groups)
        # log.debug(f"[{str_node}] Adding to tracer...")
        # if tracer and tracer.enabled:
        #     tracer.tick(current=int(node), chosen=accepted, frontier=frontier)

        log.debug(f"[{str_node}] Adding to stack...")

        for i, (sig, chosen_g, inter_g, frontier_g) in enumerate(groups):
            log.debug(f"[{str_node}] {i}/{len(groups)} group")
            # Evitamos revisitar exatamente o mesmo estado (mesmo conjunto escolhido, mesma máscara e mesma fronteira),
            # o que previne ciclos e explosão desnecessária de busca.
            if sig in visited_states:
                # log.debug(f"[{node_id}/{len_order}] Sig already seen")
                pruned_dupe_state += 1
                continue
            visited_states.add(sig)
            added = tuple(sorted(set(chosen_g) - set(accepted)))

            stack.append((chosen_g, inter_g, frontier_g, added))
            stack_small.append((chosen_g, frontier_g))
            pushes += 1
        log.debug(f"[{str_node}] Stack created.")
        to_break = False
        # Log inicial para este nó-semente: quantos estados entraram na pilha e um snapshot compacto.
        # log.debug(f"[node {int(node)}] inicializado com {len(stack)} stack | {stack_small}")

        # Loop principal de busca em profundidade (DFS) sobre os estados empilhados.
        current_stack = 0
        while stack:
            current_stack+=1
            str_stack = f"{str_node}|{current_stack}"
            log.debug(f"[{str_node}] Current stack: {current_stack}")
            chosen_t, inter_m, frontier_t, added_t = stack.pop()
            # if tracer and tracer.enabled:
            #     tracer.tick(current=added_t if added_t else (),
            #                 chosen=chosen_t,
            #                 frontier=frontier_t)
            pops += 1
            chosen = list(chosen_t)
            frontier = set(frontier_t)

            log.debug(f"[{str_stack}] Expanding group")
            # primeiro tenta expandir
            children = []
            groups_s = expand_groups(chosen_t, inter_m, frontier)

            log.debug(f"[{str_stack}] Finished expanding groups: {groups_s} ")

            for sig, c_g, m_g, f_g in groups_s:
                if sig in visited_states:
                    pruned_dupe_state += 1
                    continue
                visited_states.add(sig)
                added = tuple(sorted(set(c_g) - set(chosen_t)))
                children.append((c_g, m_g, f_g, added))

            if children:
                log.debug(f"[{str_stack}] Found childrens, putting in the stack: {children}")
                # ainda há filhos: empilha e segue a busca
                stack.extend(children)
                pushes += len(children)
            else:
                # 2) nó-folha: avaliamos aceitação do subgrafo
                log.debug(f"[{str_stack}] Reached the leaf, avaliating...")
                if len(chosen) >= 4:
                    log.debug(f"[{str_stack}] Len of chose: {len(chosen)}")
                    cn = frozenset(chosen)
                    if cn not in checked_node_sets:
                        checked_node_sets.add(cn)
                        log.debug(f"[{str_stack}] Checking if it's a valid subraph...")
                        ok, es = valid_subgraph(chosen)
                        log.debug(f"[{str_stack}] Check finished")
                        if ok:
                            # chave canônica das arestas: pares (u<v), conjunto não-ordenado
                            edge_key = frozenset(tuple(sorted(e)) for e in es)
                            if edge_key not in seen_edge_sets:
                                seen_edge_sets.add(edge_key)
                                # if tracer and tracer.enabled:
                                #     tracer.tick(current=None,
                                #                 chosen=chosen,
                                #                 frontier=(),
                                #                 accepted_edges_latest=es)
                                log.debug(f"[{str_stack}] Converting edges...")
                                _, edges_idx, edges_res = convert_edges_to_residues(es, maps)
                                log.debug(f"[{str_stack}] Conversion finished.")
                                frames[next_frame_id] = {
                                    "edges_indices": edges_idx,
                                    "edges_residues": edges_res,
                                }
                                accepts += 1
                                next_frame_id += 1

                                if residue_tracker is not None:
                                    ctx = TrackCtx(run_id=config.get("run_id","default"), stage="frames", step_id=step, chunk_id=chunk_id, frame_id=next_frame_id)
                                    residue_tracker.frame_accepted(ctx=ctx, edges_residues=edges_res, edges_indices=edges_idx)

                else:
                    log.debug(f"Chose have len < 4: {len(chosen)}")

            if (pops % debug_every == 0):
                log.debug(
                    f"[progress] pops={pops:,} pushes={pushes:,} accepts={accepts:,} "
                    f"visited={len(visited_states):,} checked={len(checked_node_sets):,} "
                    f"pruned_dupe={pruned_dupe_state:,} stack={len(stack)}"
                )
                if pushes > 600_000:
                    print(f"The stack is above than 600.000 pushes | {pushes}. The component have {len_component} nodes, and it came from the chunk {chunk_id}, step {step}.")
                    to_break = False
                    # return _, _, True
            if to_break:
                break
        log.debug(f"[{str_node}] Finished while stack.")
        # if tracer and tracer.enabled:
        #     out_file = tracer.end_anchor(name_prefix=f"chunk{chunk_id}_step{step}_anchor{int(node)}")
        #     if debug and out_file:
        #         log.debug(f"[{str_node}][viz] salvo: {out_file}")


    log.debug("Sorting frames...")

    def canon_edge(e):
        """Retorna a aresta em forma canônica (u, v) sem alterar a ordem interna dos nós."""
        u, v = e
        u = tuple(u)
        v = tuple(v)
        return (u, v) if u <= v else (v, u)

    def canon_edges(edges_idx):
        """Retorna conjunto não-direcionado e sem duplicatas de arestas canônicas."""
        return frozenset(canon_edge(e) for e in edges_idx)


    def filter_subframes(frames, debug=False):
        """
        Remove frames que são subconjuntos de outros (mantém apenas os máximos por inclusão),
        mantendo o frame base (id=0) como o primeiro no resultado final.
        """
        log.debug("Starting subframe filtering...")

        # 1) separar frame base
        base0 = frames.get(0)
        items = [(fid, frames[fid]) for fid in frames if fid != 0]


        log.debug(f"Found {len(items)} non-base frames to process.")
        if base0 is not None:
            log.debug("Frame 0 detected and will be re-added at the end as base frame.")

        # 2) ordenar por tamanho canônico de arestas (decrescente)
        log.debug("Canonicalizing and sorting frames by edge count...")
        items_sorted = sorted(
            ((fid, fr, canon_edges(fr["edges_indices"])) for fid, fr in items),
            key=lambda t: -len(t[2])
        )

        if debug:
            log.debug("Frames sorted by descending edge count:")
            for idx, (fid, _, es) in enumerate(items_sorted[:10]):
                log.debug(f"  [{idx}] fid={fid} edges={len(es)}")

        # 3) filtro de subconjunto
        kept_sets = []
        kept_frames = []
        total = len(items_sorted)
        for idx, (fid, fr, es) in enumerate(items_sorted):
            if any(es.issubset(k) for k in kept_sets):
                log.debug(f"[{idx}/{total}] Frame {fid} skipped (subset of existing).")
                continue

            # remove qualquer conjunto estritamente menor já contido nesse novo
            to_drop = [k for k in kept_sets if k.issubset(es) and k != es]
            if to_drop:
                log.debug(f"[{idx}/{total}] Frame {fid} supersedes {len(to_drop)} smaller frames.")
                kept_sets = [k for k in kept_sets if k not in to_drop]

            kept_sets.append(es)
            kept_frames.append(fr)

            log.debug(f"[{idx}/{total}] Frame {fid} accepted. Total kept: {len(kept_frames)}")

        # 4) montar final_frames com frame 0 na frente
        final_frames = {}
        fid_out = 0
        if base0 is not None:
            final_frames[fid_out] = base0
            fid_out += 1

        for fr in kept_frames:
            final_frames[fid_out] = fr
            fid_out += 1

        log.debug(f"Subframe filtering complete: {len(final_frames)} frames kept total "
                      f"({len(kept_frames)} non-base + base0).")

        return final_frames

    if steps_end:
        final_frames = filter_subframes(frames, debug=debug)
    else:
        final_frames = frames

    # log.debug(f"[done] frames_out={len(final_frames)} "
                # f"(accepted_raw={len(others)}, kept_maximal={len(final_frames)-1})")

    edge_to_res = {}          # dict[tuple(edge)] -> tuple(residues) | None
    union_order = []          # ordem determinística: primeira ocorrência ao percorrer os frames

    for fid, fr in final_frames.items():
        if fid == 0:
            continue
        ei = fr.get("edges_indices", [])
        er = fr.get("edges_residues", [])
        for j, e in enumerate(ei):
            e_t = e if isinstance(e, tuple) else tuple(e)
            if e_t not in edge_to_res:
                union_order.append(e_t)
                edge_to_res[e_t] = tuple(er[j]) if (er and j < len(er)) else None

    union_edges_indices = union_order
    union_edges_residues = {frozenset(edge_to_res[e]) for e in union_order}

    union_graph = {
        "edges_indices": union_edges_indices,
        "edges_residues": union_edges_residues,
    }

    log.debug(
        f"[union] frames_used={len(final_frames)-1} edges={len(union_edges_indices)} "
        f"(excluded_fid=0)"
    )

    # log.debug(
        # f"[done] frames_out={len(final_frames)} "
        # f"(accepted_raw={len(others)}, kept_maximal={len(final_frames)-1})"
    # )

    return final_frames, union_graph, False

def create_graph(edges_dict: dict, typeEdge: str = "edges_indices", comp_id = 0):
    Graphs = []
    k = 0
    for frame in range(0, len(edges_dict.keys())):
        edges = edges_dict[frame][typeEdge]

        G_sub = nx.Graph()

        if len(edges) > 1:

            for sublist in edges:
                sublist = list(sublist)

                node_a = tuple(sublist[0]) if isinstance(sublist[0], np.ndarray) else sublist[0]
                node_b = tuple(sublist[1]) if isinstance(sublist[1], np.ndarray) else sublist[1]
                G_sub.add_edge(node_a, node_b)

            chain_color_map = {}
            color_palette = plt.cm.get_cmap('tab10', 20)
            color_counter = 1

            if typeEdge == "edges_residues":
                for nodes in G_sub.nodes:
                    chain_id = nodes[0][0]+nodes[1][0]

                    if chain_id not in chain_color_map and chain_id[::-1] not in chain_color_map:
                        chain_color_map[chain_id] = color_palette(color_counter)[:3]
                        chain_color_map[chain_id[::-1]] = chain_color_map[chain_id]  # RGB tuple
                        color_counter += 1

                    G_sub.nodes[nodes]['chain_id'] = chain_color_map[chain_id]

            G_sub.remove_nodes_from(list(nx.isolates(G_sub)))
            log.debug(f"{comp_id} Number of nodes graph {k}: {len(G_sub.nodes)}")
            k+= 1

            if k >= 100:
                break
            Graphs.append(G_sub)
    return Graphs

def create_sphere_residue(residue_name, residue_number, coord):
    atom = Atom("CA", coord, 1.0, 1.0, " ", "DUM", 1, element='CA')
    residue = Residue((' ', int(residue_number), ' '), residue_name, " ")
    residue.add(atom)
    return residue

def align_structures_by_chain(reference_pdb, target_pdb, chain_id):
    """
    Aligns a target structure to a reference structure based on a specified chain ID.

    Args:
        reference_pdb (str): Path to the reference PDB file.
        target_pdb (str): Path to the target PDB file.
        chain_id (str): The chain ID to use for alignment.

    Returns:
        Superimposer: The Superimposer object after alignment.
    """
    # Create a PDB parser
    parser = PDBParser(QUIET=True)

    # Load the structures
    reference_structure = parser.get_structure('reference', reference_pdb)
    target_structure = parser.get_structure('target', target_pdb)

    # Extract the chains
    reference_chain = reference_structure[0][chain_id]
    target_chain = target_structure[0][chain_id]

    # Extract the CA atoms from both chains for alignment
    reference_atoms = [atom for atom in reference_chain.get_atoms() if is_aa(atom.get_parent()) and atom.get_id() == 'CA']
    target_atoms = [atom for atom in target_chain.get_atoms() if is_aa(atom.get_parent()) and atom.get_id() == 'CA']

    if len(reference_atoms) != len(target_atoms):
        raise ValueError("The number of CA atoms in the chains do not match. Ensure both chains have the same number of residues.")

    # Create a Superimposer object
    super_imposer = Superimposer()

    # Set the reference and target atoms for alignment
    super_imposer.set_atoms(reference_atoms, target_atoms)

    # Apply the transformation to the target structure
    super_imposer.apply(target_structure.get_atoms())

    # Print the RMSD
    print(f"RMSD: {super_imposer.rms:.4f} Å")

    return super_imposer

def add_sphere_residues(graphs, list_node_names_mol, output_path, node_name):

    for graph, node_names_mol in zip(graphs, list_node_names_mol):
        mol_path = graph[1]

        # Read PDB file
        parser = PDBParser()
        structure = parser.get_structure('protein', mol_path)

        # Create a new structure to hold the spheres
        new_structure = Structure.Structure("spheres")

        # Create a new model and chain for each mol
        new_model = Model.Model(0)
        new_chain = Chain.Chain('X')
        new_model.add(new_chain)
        new_structure.add(new_model)

        # Keep track of added residues to avoid duplicates
        added_residues = set()

        # Add sphere residues to the new structure
        for residue_info in node_names_mol:
            chain_id, residue_name, residue_number = residue_info.split(':')
            residue_key = (chain_id, int(residue_number))
            if residue_key not in added_residues:
                residue = structure[0][chain_id][int(residue_number)]
                #ca_atom = residue['CA']
                atom_coords = [atom.coord for atom in residue]
                centroid_coords = np.mean(atom_coords, axis=0)
                sphere_residue = create_sphere_residue(residue_name, residue_number, centroid_coords)
                new_chain.add(sphere_residue)
                added_residues.add(residue_key)

        name = mol_path.replace("\\", "_").replace("/", "_")
        # Write the new PDB file
        io = PDBIO()
        io.set_structure(new_structure)
        io.save(path.join(output_path,f'spheres_{name}_{node_name}.pdb'))

def convert_1aa3aa(AA):
    amino_acid_codes = {
    'A': 'ALA',
    'R': 'ARG',
    'N': 'ASN',
    'D': 'ASP',
    'C': 'CYS',
    'E': 'GLU',
    'Q': 'GLN',
    'G': 'GLY',
    'H': 'HIS',
    'I': 'ILE',
    'L': 'LEU',
    'K': 'LYS',
    'M': 'MET',
    'F': 'PHE',
    'P': 'PRO',
    'S': 'SER',
    'T': 'THR',
    'W': 'TRP',
    'Y': 'TYR',
    'V': 'VAL'}

    return amino_acid_codes[AA]

def convert_3aa1aa(AA):
    amino_acid_codes = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLU': 'E',
    'GLN': 'Q',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V'}

    return amino_acid_codes[AA]

