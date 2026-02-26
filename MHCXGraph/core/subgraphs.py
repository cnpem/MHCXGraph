from __future__ import annotations
from typing import Dict, List, Optional, Tuple, Union, Set

import logging
import networkx as nx
import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

"""
Subgraph utilities for protein structure graphs.

Assumptions
-----------
- Nodes represent residues or atoms and may carry:
  'chain_id' (or 'chain'), 'residue_number' (or 'resseq'),
  'residue_name' (or 'resname'), and coordinates in 'coords' or 'centroid'.
- Graph-level metadata (G.graph) may include:
  'pdb_df', 'raw_pdb_df', 'rgroup_df', 'coords', 'distance_matrix',
  'dssp_df', 'residue_labels', 'water_labels', 'water_positions'.

The functions below select subsets by chain, residue type, spatial
radius, secondary structure, RSA, edge kind, k-hop, etc., and
propagate/update relevant graph metadata to the returned subgraph.
"""


# -------------------------- errors & small helpers --------------------------

class ProteinGraphConfigurationError(RuntimeError):
    """Raised when required graph/node annotations are missing."""


def _node_attr(d: dict, key: str, *fallbacks: str):
    """
    Get a node attribute with fallbacks.

    Parameters
    ----------
    d
        Node attributes mapping.
    key
        Primary key.
    *fallbacks
        Alternative keys to try.

    Returns
    -------
    Any or None
        First non-None value found.
    """
    for k in (key, *fallbacks):
        if k in d and d[k] is not None:
            return d[k]
    return None


def _node_coords(d: dict) -> Optional[np.ndarray]:
    """
    Return node coordinates.

    Tries 'coords', then 'centroid'.

    Parameters
    ----------
    d
        Node attributes mapping.

    Returns
    -------
    np.ndarray or None
        Array of shape (3,) or None.
    """
    c = _node_attr(d, "coords", "centroid")
    if c is None:
        return None
    return np.asarray(c, dtype=float)


def _ensure_set(x) -> set:
    """
    Normalize an edge-kind value to a set of strings.

    Parameters
    ----------
    x
        String, iterable, or None.

    Returns
    -------
    set
        Set of strings (possibly empty).
    """
    if x is None:
        return set()
    if isinstance(x, str):
        return {x}
    try:
        return set(str(v) for v in x)
    except Exception:
        return {str(x)}


def _update_coords_graph(g: nx.Graph) -> None:
    """
    Populate `g.graph['coords']` and `g.graph['residue_labels']` from nodes.

    Parameters
    ----------
    g
        Graph.
    """
    labels = list(g.nodes())
    coords = []
    missing = 0
    for n in labels:
        arr = _node_coords(g.nodes[n])
        if arr is None:
            missing += 1
            arr = np.array([np.nan, np.nan, np.nan], dtype=float)
        coords.append(arr)
    g.graph["residue_labels"] = labels
    g.graph["coords"] = np.vstack(coords) if coords else np.zeros((0, 3), dtype=float)
    if missing:
        log.debug(f"[subgraph] {missing} nodes without coords/centroid; filled with NaN.")


def compute_distmat(pdb_df: pd.DataFrame) -> np.ndarray:
    """
    Compute Euclidean distance matrix between nodes.

    Multiple rows per node_id are averaged first.

    Parameters
    ----------
    pdb_df : pandas.DataFrame
        Must contain: ['node_id', 'x_coord', 'y_coord', 'z_coord'].

    Returns
    -------
    np.ndarray
        Distance matrix (N, N) in the order of first occurrence of each node_id.
    """
    if pdb_df is None or len(pdb_df) == 0:
        return np.zeros((0, 0), dtype=float)

    seen, order = set(), []
    for nid in pdb_df["node_id"].tolist():
        if nid not in seen:
            seen.add(nid)
            order.append(nid)

    grouped = (
        pdb_df.groupby("node_id")[["x_coord", "y_coord", "z_coord"]]
        .mean()
        .reindex(order)
    )
    P = grouped.to_numpy(dtype=float)
    if P.size == 0:
        return np.zeros((0, 0), dtype=float)

    diff = P[:, None, :] - P[None, :, :]
    return np.sqrt((diff * diff).sum(axis=2))


# -------------------------- subgraph core --------------------------

def _filter_df_by_nodes(df, node_list):
    """
    Filter a DataFrame by node list.

    Parameters
    ----------
    df : pandas.DataFrame or None
        DataFrame to filter.
    node_list : list of str
        Node IDs to keep.

    Returns
    -------
    pandas.DataFrame or None
        Filtered copy (or None if input is None).
    """
    if df is None:
        return None
    try:
        if "node_id" in df.columns:
            return df[df["node_id"].isin(node_list)].copy()
    except Exception:
        pass
    try:
        if getattr(df.index, "name", None) == "node_id" or df.index.name is not None:
            return df[df.index.isin(node_list)].copy()
    except Exception:
        pass
    return df.copy()


def _carry_graph_level(
    parent: nx.Graph,
    child: nx.Graph,
    node_list: List[str],
    filter_dataframe: bool,
    update_coords: bool,
    recompute_distmat: bool,
):
    """
    Propagate graph-level metadata from parent to subgraph.

    Parameters
    ----------
    parent, child : nx.Graph
        Source and target graphs.
    node_list : list of str
        Nodes present in the subgraph.
    filter_dataframe : bool
        Whether to filter DataFrames to `node_list`.
    update_coords : bool
        Whether to recompute `child.graph['coords']` from node attributes.
    recompute_distmat : bool
        Whether to recompute a distance matrix for the subgraph.
    """
    for k in ("config", "path", "name"):
        if k in parent.graph:
            child.graph[k] = parent.graph[k]

    if "dssp_df" in parent.graph:
        child.graph["dssp_df"] = (
            _filter_df_by_nodes(parent.graph["dssp_df"], node_list)
            if filter_dataframe else parent.graph["dssp_df"]
        )

    for key in ("pdb_df", "raw_pdb_df", "rgroup_df"):
        if key in parent.graph:
            child.graph[key] = (
                _filter_df_by_nodes(parent.graph[key], node_list)
                if filter_dataframe else parent.graph[key]
            )
            if filter_dataframe and hasattr(child.graph[key], "reset_index"):
                child.graph[key] = child.graph[key].reset_index(drop=True)

    if "coords" in parent.graph and not update_coords:
        coords_map = {n: c for n, c in zip(parent.nodes(), parent.graph["coords"])}
        child.graph["coords"] = np.array([coords_map[n] for n in child.nodes()])

    if "distance_matrix" in parent.graph and not recompute_distmat:
        try:
            labels = parent.graph.get("residue_labels")
            if labels:
                pos = {nid: i for i, nid in enumerate(labels)}
                idx = [pos[n] for n in node_list if n in pos]
                dm = parent.graph["distance_matrix"]
                child.graph["distance_matrix"] = dm[np.ix_(idx, idx)]
        except Exception as e:
            log.debug(f"{e}")

    child.graph["residue_labels"] = list(child.nodes())
    if "water_labels" in parent.graph:
        child.graph["water_labels"] = [n for n in parent.graph["water_labels"] if n in child.nodes()]
    if "water_positions" in parent.graph:
        child.graph["water_positions"] = parent.graph["water_positions"]


# -------------------------- public API --------------------------

def extract_subgraph_from_node_list(
    g: nx.Graph,
    node_list: Optional[List[str]],
    filter_dataframe: bool = True,
    update_coords: bool = True,
    recompute_distmat: bool = False,
    inverse: bool = False,
    return_node_list: bool = False,
) -> Union[nx.Graph, List[str]]:
    """
    Build a subgraph from an explicit node list.

    Parameters
    ----------
    g : nx.Graph
        Input graph.
    node_list : list of str or None
        Nodes to keep. If None, returns `g`.
    filter_dataframe : bool, default=True
        Filter graph-level DataFrames to subgraph nodes.
    update_coords : bool, default=True
        Rebuild `graph['coords']` from node attributes.
    recompute_distmat : bool, default=False
        Recompute `graph['dist_mat']` from `pdb_df` if available.
    inverse : bool, default=False
        If True, keep the complement of `node_list`.
    return_node_list : bool, default=False
        If True, return the resolved node list instead of a subgraph.

    Returns
    -------
    nx.Graph or list of str
        Subgraph or node list.
    """
    if node_list is None:
        if return_node_list:
            return []
        return g

    if inverse:
        node_list = [n for n in g.nodes() if n not in node_list]

    if return_node_list:
        return node_list

    parent = g
    sub = parent.subgraph(node_list).copy()

    # Filter pdb_df if present
    if filter_dataframe and "pdb_df" in parent.graph:
        df_filtered = _filter_df_by_nodes(parent.graph["pdb_df"], node_list)
        if df_filtered is not None:
            sub.graph["pdb_df"] = df_filtered.reset_index(drop=True)

    # Update coordinates if requested
    if update_coords:
        coords = []
        for _, d in sub.nodes(data=True):
            arr = _node_coords(d)
            coords.append(arr if arr is not None else np.zeros(3, dtype=float))
        sub.graph["coords"] = np.vstack(coords) if coords else np.zeros((0, 3), dtype=float)

    # Optionally recompute distance matrix
    if recompute_distmat:
        if not filter_dataframe and "pdb_df" not in sub.graph and "pdb_df" in parent.graph:
            df_filtered = _filter_df_by_nodes(parent.graph["pdb_df"], node_list)
            if df_filtered is not None:
                sub.graph["pdb_df"] = df_filtered.reset_index(drop=True)
        try:
            if "pdb_df" in sub.graph and sub.graph["pdb_df"] is not None:
                sub.graph["dist_mat"] = compute_distmat(sub.graph["pdb_df"])
        except Exception as e:
            log.debug(f"Failed to recompute dist_mat: {e}")

    _carry_graph_level(parent, sub, node_list, filter_dataframe, update_coords, recompute_distmat)
    return sub

def extract_subgraph_from_point(
    g: nx.Graph,
    centre_point: Union[np.ndarray, Tuple[float, float, float]],
    radius: float,
    filter_dataframe: bool = True,
    update_coords: bool = True,
    recompute_distmat: bool = False,
    inverse: bool = False,
    return_node_list: bool = False,
) -> Union[nx.Graph, List[str], None]:
    """
    Select nodes within a sphere.

    Parameters
    ----------
    g : nx.Graph
        Input graph.
    centre_point : array-like of shape (3,)
        Sphere center.
    radius : float
        Sphere radius.
    filter_dataframe, update_coords, recompute_distmat, inverse, return_node_list
        See :func:`extract_subgraph_from_node_list`.

    Returns
    -------
    nx.Graph or list of str
        Subgraph or node list.
    """
    node_list: List[str] = []
    cp = np.asarray(centre_point, dtype=float)

    for n, d in g.nodes(data=True):
        coords = _node_coords(d)
        if coords is None:
            continue
        if np.linalg.norm(coords - cp) < float(radius):
            node_list.append(n)

    node_list = list(set(node_list))
    return extract_subgraph_from_node_list(
        g,
        node_list,
        filter_dataframe=filter_dataframe,
        inverse=inverse,
        return_node_list=return_node_list,
        recompute_distmat=recompute_distmat,
        update_coords=update_coords,
    )


def extract_subgraph_from_atom_types(
    g: nx.Graph,
    atom_types: List[str],
    filter_dataframe: bool = True,
    update_coords: bool = True,
    recompute_distmat: bool = False,
    inverse: bool = False,
    return_node_list: bool = False,
) -> Union[nx.Graph, List[str], None]:
    """
    Select nodes by atom type.

    Parameters
    ----------
    g : nx.Graph
        Input graph.
    atom_types : list of str
        Allowed atom types.
    filter_dataframe, update_coords, recompute_distmat, inverse, return_node_list
        See :func:`extract_subgraph_from_node_list`.

    Returns
    -------
    nx.Graph or list of str
        Subgraph or node list.
    """
    node_list: List[str] = [
        n for n, d in g.nodes(data=True) if _node_attr(d, "atom_type") in set(atom_types)
    ]
    node_list = list(dict.fromkeys(node_list))
    return extract_subgraph_from_node_list(
        g,
        node_list,
        filter_dataframe=filter_dataframe,
        inverse=inverse,
        return_node_list=return_node_list,
        update_coords=update_coords,
        recompute_distmat=recompute_distmat,
    )


def extract_subgraph_from_residue_types(
    g: nx.Graph,
    residue_types: Union[List[str], Set[str]],
    filter_dataframe: bool = True,
    update_coords: bool = True,
    recompute_distmat: bool = False,
    inverse: bool = False,
    return_node_list: bool = False,
) -> Union[nx.Graph, List[str], None]:
    """
    Select nodes by residue name.

    Parameters
    ----------
    g : nx.Graph
        Input graph.
    residue_types : list of str
        Allowed residue names (3-letter).
    filter_dataframe, update_coords, recompute_distmat, inverse, return_node_list
        See :func:`extract_subgraph_from_node_list`.

    Returns
    -------
    nx.Graph or list of str
        Subgraph or node list.
    """
    residue_types = set(residue_types)
    node_list: List[str] = []
    for n, d in g.nodes(data=True):
        rname = _node_attr(d, "residue_name", "resname")
        if rname in residue_types:
            node_list.append(n)

    node_list = list(dict.fromkeys(node_list))
    return extract_subgraph_from_node_list(
        g,
        node_list,
        filter_dataframe=filter_dataframe,
        inverse=inverse,
        return_node_list=return_node_list,
        update_coords=update_coords,
        recompute_distmat=recompute_distmat,
    )


def extract_subgraph_from_chains(
    g: nx.Graph,
    chains: Union[List[str], Set[str]],
    filter_dataframe: bool = True,
    update_coords: bool = True,
    recompute_distmat: bool = False,
    inverse: bool = False,
    return_node_list: bool = False,
) -> Union[nx.Graph, List[str], None]:
    """
    Select nodes by chain IDs.

    Parameters
    ----------
    g : nx.Graph
        Input graph.
    chains : list of str
        Chain IDs to include.
    filter_dataframe, update_coords, recompute_distmat, inverse, return_node_list
        See :func:`extract_subgraph_from_node_list`.

    Returns
    -------
    nx.Graph or list of str
        Subgraph or node list.
    """
    chains = set(chains)
    node_list: List[str] = []
    for n, d in g.nodes(data=True):
        cid = _node_attr(d, "chain_id", "chain")
        if cid in chains:
            node_list.append(n)

    node_list = list(dict.fromkeys(node_list))
    return extract_subgraph_from_node_list(
        g,
        node_list,
        filter_dataframe=filter_dataframe,
        update_coords=update_coords,
        recompute_distmat=recompute_distmat,
        inverse=inverse,
        return_node_list=return_node_list,
    )


def extract_subgraph_by_sequence_position(
    g: nx.Graph,
    sequence_positions: List[int],
    filter_dataframe: bool = True,
    update_coords: bool = True,
    recompute_distmat: bool = False,
    inverse: bool = False,
    return_node_list: bool = False,
) -> Union[nx.Graph, List[str], None]:
    """
    Select nodes by residue index.

    Parameters
    ----------
    g : nx.Graph
        Input graph.
    sequence_positions : list of int
        Residue numbers to include.
    filter_dataframe, update_coords, recompute_distmat, inverse, return_node_list
        See :func:`extract_subgraph_from_node_list`.

    Returns
    -------
    nx.Graph or list of str
        Subgraph or node list.
    """
    seqset = set(int(x) for x in sequence_positions)
    node_list: List[str] = []
    for n, d in g.nodes(data=True):
        rnum = _node_attr(d, "residue_number", "resseq")
        if rnum in seqset:
            node_list.append(n)

    node_list = list(dict.fromkeys(node_list))
    return extract_subgraph_from_node_list(
        g,
        node_list,
        filter_dataframe=filter_dataframe,
        inverse=inverse,
        return_node_list=return_node_list,
        recompute_distmat=recompute_distmat,
        update_coords=update_coords,
    )


def extract_subgraph_by_bond_type(
    g: nx.Graph,
    bond_types: Union[List[str], Set[str]],
    filter_dataframe: bool = True,
    update_coords: bool = True,
    recompute_distmat: bool = False,
    inverse: bool = False,
    return_node_list: bool = False,
) -> Union[nx.Graph, List[str], None]:
    """
    Select nodes incident to edges of specified kinds.

    Parameters
    ----------
    g : nx.Graph
        Input graph.
    bond_types : list of str
        Edge kinds to match.
    filter_dataframe, update_coords, recompute_distmat, inverse, return_node_list
        See :func:`extract_subgraph_from_node_list`.

    Returns
    -------
    nx.Graph or list of str
        Subgraph or node list.
    """
    bond_types = set(bond_types)
    node_list: List[str] = []

    for u, v, d in g.edges(data=True):
        kinds = _ensure_set(d.get("kind"))
        if kinds & bond_types:
            node_list.extend((u, v))

    node_list = list(dict.fromkeys(node_list))

    for _, _, d in g.edges(data=True):
        kinds = _ensure_set(d.get("kind"))
        if not inverse:
            kinds = {k for k in kinds if k in bond_types}
        else:
            kinds = {k for k in kinds if k not in bond_types}
        d["kind"] = kinds if len(kinds) != 1 else next(iter(kinds), None)

    return extract_subgraph_from_node_list(
        g,
        node_list,
        filter_dataframe=filter_dataframe,
        inverse=inverse,
        return_node_list=return_node_list,
        recompute_distmat=recompute_distmat,
        update_coords=update_coords,
    )


def extract_subgraph_from_secondary_structure(
    g: nx.Graph,
    ss_elements: List[str],
    inverse: bool = False,
    filter_dataframe: bool = True,
    recompute_distmat: bool = False,
    update_coords: bool = True,
    return_node_list: bool = False,
) -> Union[nx.Graph, List[str], None]:
    """
    Select nodes by secondary structure label.

    Parameters
    ----------
    g : nx.Graph
        Input graph. Nodes must carry 'ss'.
    ss_elements : list of str
        Allowed secondary structure labels.
    inverse : bool, default=False
        If True, exclude `ss_elements`.
    filter_dataframe, recompute_distmat, update_coords, return_node_list
        See :func:`extract_subgraph_from_node_list`.

    Returns
    -------
    nx.Graph or list of str
        Subgraph or node list.

    Raises
    ------
    ProteinGraphConfigurationError
        If any node lacks the 'ss' attribute.
    """


    name_to_codes = {
        "HELIX": {"H", "G", "I"},
        "STRAND": {"B", "E"},
        "LOOP": {"T", "S", "-"},

        "ALPHA HELIX": {"H"},
        "3-10 HELIX": {"G"},
        "PI HELIX": {"I"},
        "ISOLATED BETA-BRIDGE RESIDUE": {"B"},
        "BETA BRIDGE": {"B"},
        "BEND": {"S"},
        "TURN": {"T"},
        "NONE": {"-"},
        "COIL": {"-"},
        "LOOP STRUCTURE": {"T", "S", "-"},

        "H": {"H"},
        "G": {"G"},
        "I": {"I"},
        "B": {"B"},
        "E": {"E"},
        "T": {"T"},
        "S": {"S"},
        "-": {"-"},
    }

    normalized: set[str] = set()

    for ss in ss_elements:
        key = ss.strip().upper()

        if key not in name_to_codes:
            raise ValueError(
                f"Unknown secondary structure name {ss!r}. "
                f"Allowed: {sorted(name_to_codes.keys())}"
            )

        normalized |= name_to_codes[key]

    node_list: List[str] = []
    for n, d in g.nodes(data=True):
        if "ss" not in d:
            raise ProteinGraphConfigurationError(
                f"Secondary structure not set for node {n}. Annotate 'ss' first."
            )
        if d["ss"] in set(ss_elements):
            node_list.append(n)
    
    node_list = list(dict.fromkeys(node_list))

    return extract_subgraph_from_node_list(
        g,
        node_list,
        inverse=inverse,
        return_node_list=return_node_list,
        filter_dataframe=filter_dataframe,
        recompute_distmat=recompute_distmat,
        update_coords=update_coords,
    )


def extract_surface_subgraph_rsa(
    g: nx.Graph,
    rsa_threshold: float = 0.2,
    inverse: bool = False,
    filter_dataframe: bool = True,
    recompute_distmat: bool = False,
    update_coords: bool = True,
    return_node_list: bool = False,
    *,
    treat_water_as_surface: bool = True,
    unknown_policy: str = "skip",
    unknown_value: Optional[float] = None,
) -> Union[nx.Graph, List[str], None]:
    """
    Select nodes by relative solvent accessibility (RSA).

    Parameters
    ----------
    g : nx.Graph
        Input graph. Nodes may carry 'rsa' in [0, 1].
    rsa_threshold : float, default=0.2
        Minimum RSA to include.
    inverse : bool, default=False
        If True, include RSA < threshold.
    filter_dataframe : bool, default=True
        Filter graph-level DataFrames to subgraph nodes.
    recompute_distmat : bool, default=False
        Recompute `graph['dist_mat']` from `pdb_df` if available.
    update_coords : bool, default=True
        Rebuild `graph['coords']` from node attributes.
    return_node_list : bool, default=False
        If True, return the resolved node list instead of a subgraph.
    treat_water_as_surface : bool, default=True
        If True, nodes with residue name typical of water (e.g. HOH/WAT/DOD/TIP3)
        are treated as RSA=1.0 when 'rsa' is missing.
    unknown_policy : {'skip', 'value', 'error'}, default='skip'
        Behavior for nodes missing 'rsa' that are not water:
        - 'skip' : ignore node (do not include, do not raise);
        - 'value': use `unknown_value` as RSA;
        - 'error': raise ProteinGraphConfigurationError.
    unknown_value : float, optional
        RSA value to use when `unknown_policy='value'`.

    Returns
    -------
    nx.Graph or list of str
        Subgraph or node list.

    Raises
    ------
    ProteinGraphConfigurationError
        If `unknown_policy='error'` and a node lacks 'rsa'.
    """
    WATER_NAMES = {"HOH", "WAT", "H2O", "DOD", "TIP3", "TIP4", "SOL"}

    node_list: List[str] = []
    thr = float(rsa_threshold)

    for n, d in g.nodes(data=True):
        rsa = d.get("rsa", None)

        if rsa is None:
            rname = str(d.get("residue_name", d.get("resname", ""))).upper()
            is_water = rname in WATER_NAMES
            if is_water and treat_water_as_surface:
                rsa = 1.0
            else:
                if unknown_policy == "skip":
                    continue
                elif unknown_policy == "value":
                    if unknown_value is None:
                        continue
                    rsa = float(unknown_value)
                elif unknown_policy == "error":
                    raise ProteinGraphConfigurationError(
                        f"RSA not set for node {n}. Annotate 'rsa' first."
                    )
                else:
                    # Fallback: behave like 'skip'
                    continue

        try:
            rsa_f = float(rsa)
        except Exception:
            # Non-numeric RSA -> honor policy
            if unknown_policy == "error":
                raise ProteinGraphConfigurationError(
                    f"RSA not numeric for node {n}: {rsa!r}"
                )
            continue

        include = (rsa_f >= thr)
        if inverse:
            include = not include
        if include:
            node_list.append(n)

    node_list = list(dict.fromkeys(node_list))
    # Delegates subgraph construction; inverse already applied above.
    return extract_subgraph_from_node_list(
        g,
        node_list=node_list,
        filter_dataframe=filter_dataframe,
        update_coords=update_coords,
        recompute_distmat=recompute_distmat,
        inverse=False,
        return_node_list=return_node_list,
    )


def extract_surface_subgraph_asa(
    g: nx.Graph,
    asa_threshold: float,
    inverse: bool = False,
    filter_dataframe: bool = True,
    recompute_distmat: bool = False,
    update_coords: bool = True,
    return_node_list: bool = False,
) -> Union[nx.Graph, List[str], None]:
    """
    Select nodes by absolute solvent accessibility (ASA).

    Parameters
    ----------
    g : nx.Graph
        Input graph. Nodes are expected to carry 'asa' (float, in Å^2).
    asa_threshold : float
        Minimum ASA to include.
    inverse : bool, default=False
        If True, include ASA < threshold.
    filter_dataframe : bool, default=True
        Filter graph-level DataFrames to subgraph nodes.
    recompute_distmat : bool, default=False
        Recompute `graph['dist_mat']` from `pdb_df` if available.
    update_coords : bool, default=True
        Rebuild `graph['coords']` from node attributes.
    return_node_list : bool, default=False
        If True, return the resolved node list instead of a subgraph.

    Returns
    -------
    nx.Graph or list of str
        Subgraph or node list.
    """
    node_list: List[str] = []
    thr = float(asa_threshold)

    for n, d in g.nodes(data=True):
        asa = d.get("asa", None)
        if asa is None:
            continue

        try:
            asa_f = float(asa)
        except Exception:
            continue

        include = (asa_f >= thr)
        if inverse:
            include = not include

        if include:
            node_list.append(n)

    node_list = list(dict.fromkeys(node_list))

    return extract_subgraph_from_node_list(
        g,
        node_list=node_list,
        filter_dataframe=filter_dataframe,
        update_coords=update_coords,
        recompute_distmat=recompute_distmat,
        inverse=False,
        return_node_list=return_node_list,
    )

def extract_k_hop_subgraph(
    g: nx.Graph,
    central_node: str,
    k: int,
    k_only: bool = False,
    filter_dataframe: bool = True,
    update_coords: bool = True,
    recompute_distmat: bool = False,
    inverse: bool = False,
    return_node_list: bool = False,
) -> Union[nx.Graph, List[str], None]:
    """
    Select nodes by k-hop neighborhood.

    Parameters
    ----------
    g : nx.Graph
        Input graph.
    central_node : str
        Center node ID.
    k : int
        Number of hops.
    k_only : bool, default=False
        If True, include exactly k-hop nodes; otherwise include all <= k.
    filter_dataframe, update_coords, recompute_distmat, inverse, return_node_list
        See :func:`extract_subgraph_from_node_list`.

    Returns
    -------
    nx.Graph or list of str
        Subgraph or node list.
    """
    neighbours: Dict[int, List[str]] = {0: [central_node]}

    for i in range(1, k + 1):
        hop_set: set[str] = set()
        for node in neighbours[i - 1]:
            hop_set.update(g.neighbors(node))
        neighbours[i] = list(hop_set)

    if k_only:
        node_list: List[str] = neighbours[k]
    else:
        # flatten and remove duplicates
        all_nodes: set[str] = set()
        for nodes in neighbours.values():
            all_nodes.update(nodes)
        node_list = list(all_nodes)

    return extract_subgraph_from_node_list(
        g,
        node_list=node_list,
        filter_dataframe=filter_dataframe,
        inverse=inverse,
        return_node_list=return_node_list,
        recompute_distmat=recompute_distmat,
        update_coords=update_coords,
    )



def extract_interface_subgraph(
    g: nx.Graph,
    interface_list: Optional[List[str]] = None,
    chain_list: Optional[List[str]] = None,
    filter_dataframe: bool = True,
    update_coords: bool = True,
    recompute_distmat: bool = False,
    inverse: bool = False,
    return_node_list: bool = False,
) -> Union[nx.Graph, List[str], None]:
    """
    Select nodes at chain-chain interfaces.

    Parameters
    ----------
    g : nx.Graph
        Input graph.
    interface_list : list of str, optional
        Allowed chain pair labels (e.g., ["AB","BC"]). If None, any pairwise
        inter-chain contact qualifies.
    chain_list : list of str, optional
        Restrict to interactions among these chains.
    filter_dataframe, update_coords, recompute_distmat, inverse, return_node_list
        See :func:`extract_subgraph_from_node_list`.

    Returns
    -------
    nx.Graph or list of str
        Subgraph or node list.
    """
    node_list: List[str] = []

    for u, v in g.edges():
        u_chain = _node_attr(g.nodes[u], "chain_id", "chain")
        v_chain = _node_attr(g.nodes[v], "chain_id", "chain")
        if u_chain is None or v_chain is None:
            continue

        if chain_list is not None:
            if u_chain in chain_list and v_chain in chain_list and u_chain != v_chain:
                node_list.extend((u, v))
        if interface_list is not None:
            case_1 = f"{u_chain}{v_chain}"
            case_2 = f"{v_chain}{u_chain}"
            if case_1 in interface_list or case_2 in interface_list:
                node_list.extend((u, v))
        if chain_list is None and interface_list is None and u_chain != v_chain:
            node_list.extend((u, v))

    node_list = list(dict.fromkeys(node_list))
    return extract_subgraph_from_node_list(
        g,
        node_list,
        filter_dataframe=filter_dataframe,
        inverse=inverse,
        return_node_list=return_node_list,
        recompute_distmat=recompute_distmat,
        update_coords=update_coords,
    )


def extract_subgraph(
    g: nx.Graph,
    node_list: Optional[List[str]] = None,
    sequence_positions: Optional[List[int]] = None,
    chains: Optional[List[str]] = None,
    residue_types: Optional[List[str]] = None,
    atom_types: Optional[List[str]] = None,
    bond_types: Optional[List[str]] = None,
    centre_point: Optional[Union[np.ndarray, Tuple[float, float, float]]] = None,
    radius: Optional[float] = None,
    ss_elements: Optional[List[str]] = None,
    rsa_threshold: Optional[float] = None,
    asa_threshold: Optional[float] = None,
    k_hop_central_node: Optional[str] = None,
    k_hops: Optional[int] = None,
    k_only: Optional[bool] = None,
    filter_dataframe: bool = True,
    update_coords: bool = True,
    recompute_distmat: bool = False,
    inverse: bool = False,
    return_node_list: bool = False,
) -> Union[nx.Graph, List[str]]:
    """
    Aggregate subgraph selector with a unified API.

    Parameters
    ----------
    g : nx.Graph
        Input graph.
    node_list : list of str, optional
        Explicit nodes to include.
    sequence_positions : list of int, optional
        Residue numbers to include.
    chains : list of str, optional
        Chain IDs to include.
    residue_types : list of str, optional
        Residue names to include.
    atom_types : list of str, optional
        Atom types to include.
    bond_types : list of str, optional
        Edge kinds whose incident nodes to include.
    centre_point : array-like, optional
        Center for point-radius selection.
    radius : float, optional
        Radius for point-radius selection.
    ss_elements : list of str, optional
        Secondary structure labels to include.
    rsa_threshold : float, optional
        Minimum RSA to include.
    k_hop_central_node : str, optional
        Node ID for k-hop selection.
    k_hops : int, optional
        Number of hops for k-hop selection.
    k_only : bool, optional
        If True, include exactly k-hop nodes; else all <= k.
    filter_dataframe, update_coords, recompute_distmat, inverse, return_node_list
        See :func:`extract_subgraph_from_node_list`.

    Returns
    -------
    nx.Graph or list of str
        Subgraph or node list.
    """
    if node_list is None:
        node_list = []

    if sequence_positions is not None:
        node_list += _ensure_list_str(extract_subgraph_by_sequence_position(
            g, sequence_positions, return_node_list=True
        ))
    if chains is not None:
        node_list += _ensure_list_str(extract_subgraph_from_chains(
            g, chains, return_node_list=True
        ))
    if residue_types is not None:
        node_list += _ensure_list_str(extract_subgraph_from_residue_types(
            g, residue_types, return_node_list=True
        ))
    if atom_types is not None:
        node_list += _ensure_list_str(extract_subgraph_from_atom_types(
            g, atom_types, return_node_list=True
        ))
    if bond_types is not None:
        node_list += _ensure_list_str(extract_subgraph_by_bond_type(
            g, bond_types, return_node_list=True
        ))
    if centre_point is not None and radius is not None:
        node_list += _ensure_list_str(extract_subgraph_from_point(
            g, centre_point, radius, return_node_list=True
        ))
    if ss_elements is not None:
        node_list += _ensure_list_str(extract_subgraph_from_secondary_structure(
            g, ss_elements, return_node_list=True
        ))
    if rsa_threshold is not None:
        node_list += _ensure_list_str(extract_surface_subgraph_rsa(
            g, rsa_threshold, return_node_list=True
        ))

    if asa_threshold is not None:
        node_list += _ensure_list_str(extract_surface_subgraph_asa(
            g, asa_threshold, return_node_list=True
        ))
    if k_hop_central_node is not None and k_hops and k_only is not None:
        node_list += _ensure_list_str(extract_k_hop_subgraph(
            g, k_hop_central_node, k_hops, k_only, return_node_list=True
        ))

    seen, merged = set(), []
    for n in node_list:
        if n not in seen:
            seen.add(n)
            merged.append(n)

    return extract_subgraph_from_node_list(
        g,
        merged,
        filter_dataframe=filter_dataframe,
        inverse=inverse,
        return_node_list=return_node_list,
        recompute_distmat=recompute_distmat,
        update_coords=update_coords,
    )

def _ensure_list_str(value: Union[nx.Graph, List[str], None]) -> List[str]:
    if isinstance(value, list):
        return value
    return []
