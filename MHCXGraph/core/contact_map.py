from __future__ import annotations

import numpy as np
from typing import Iterable, List, Tuple, Dict, Optional


def _valid_coord(x) -> bool:
    try:
        a = np.asarray(x, dtype=float)
    except Exception:
        return False
    return a.shape == (3,) and np.isfinite(a).all()


def contact_map_from_graph(
    G,
    *,
    granularity: str,
    exclude_kinds: Iterable[str] = ("water",),
    fallback_to_centroid: bool = True,
) -> Tuple[
    np.ndarray,
    List[str],
    Dict[Tuple[str, str], int],
    Dict[Tuple[str, str, str], int]
]:
    """
    Return:
      contact_map (NxN)
      node_order (list of node_ids)
      residue_map_dict: (chain_id, residue_number) -> idx
      residue_map_dict_all: (chain_id, residue_number, residue_name) -> idx
    """

    gran = (granularity or "all_atoms").strip().lower()

    node_order: List[str] = []
    coords: List[np.ndarray] = []
    chain_ids: List[str] = []
    residue_numbers: List[str] = []
    residue_names: List[str] = []

    for nid, d in G.nodes(data=True):
        kind = d.get("kind") or ""
        if kind in exclude_kinds:
            continue

        cent = d.get("centroid", None)
        ca = d.get("ca_coord", None)
        cb = d.get("cb_coord", None)

        picked = None

        if gran == "ca_only":
            if _valid_coord(ca):
                picked = ca
            elif _valid_coord(cb):
                picked = cb
            elif fallback_to_centroid and _valid_coord(cent):
                picked = cent
        else:
            # other granularities: centroid already computed based on granularity
            if _valid_coord(cent):
                picked = cent
            elif fallback_to_centroid and _valid_coord(ca):
                picked = ca
            elif fallback_to_centroid and _valid_coord(cb):
                picked = cb

        if picked is None:
            continue

        node_order.append(nid)
        coords.append(np.asarray(picked, dtype=float))

        # Build residue maps
        chain = str(d.get("chain_id", d.get("chain")))
        resnum = str(d.get("residue_number", d.get("resseq")))
        resname = str(d.get("residue_name", d.get("resname")))

        chain_ids.append(chain)
        residue_numbers.append(resnum)
        residue_names.append(resname)

    if not coords:
        return (
            np.zeros((0, 0), dtype=float),
            [],
            {},
            {},
        )

    C = np.vstack(coords)
    diff = C[:, None, :] - C[None, :, :]
    contact_map = np.sqrt(np.sum(diff * diff, axis=2))

    # Build maps
    residue_map = list(zip(chain_ids, residue_numbers))
    residue_map_all = list(zip(chain_ids, residue_numbers, residue_names))

    residue_map_dict: Dict[Tuple[str, str], int] = {t: i for i, t in enumerate(residue_map)}
    residue_map_dict_all: Dict[Tuple[str, str, str], int] = {
        t: i for i, t in enumerate(residue_map_all)
    }

    return (
        contact_map,
        node_order,
        residue_map_dict,
        residue_map_dict_all,
    )
