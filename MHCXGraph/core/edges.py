# core/edges.py
from __future__ import annotations

import numpy as np
import pandas as pd
from typing import Optional, Iterable

def _pairs_under_threshold(coords: np.ndarray, thr: float):
    """Retorna índices (i,j) com i<j e 0<dist<=thr, mais as distâncias."""
    diff = coords[:, None, :] - coords[None, :, :]
    dist = np.sqrt(np.sum(diff * diff, axis=2))
    iu, ju = np.triu_indices(len(coords), k=1)
    mask = (dist[iu, ju] <= thr) & (dist[iu, ju] > 0.0)
    return iu[mask], ju[mask], dist[iu, ju][mask]

def _ensure_centroids_df(G, *, exclude_kinds: Iterable[str] = ("water",)) -> pd.DataFrame:
    """
    Constrói um DataFrame de nós a partir do grafo, contendo:
    node_id, chain_id, residue_number, residue_name, insertion, x_coord, y_coord, z_coord.
    Ignora nós cujo 'kind' esteja em exclude_kinds.
    """
    rows = []
    for n, d in G.nodes(data=True):
        if (d.get("kind") or "") in exclude_kinds:
            continue
        chain = d.get("chain_id", d.get("chain"))
        resnum = d.get("residue_number", d.get("resseq"))
        resname = d.get("residue_name", d.get("resname"))
        ins = d.get("insertion", d.get("icode", "")) or ""
        cent = d.get("centroid", None)
        if chain is None or resnum is None or resname is None or cent is None:
            continue
        x, y, z = map(float, cent)
        rows.append({
            "node_id": str(n),
            "chain_id": str(chain),
            "residue_number": int(resnum),
            "residue_name": str(resname),
            "insertion": str(ins),
            "x_coord": x, "y_coord": y, "z_coord": z,
        })
    return pd.DataFrame(rows)

def add_distance_threshold(
    G,
    *,
    threshold: float,
    long_interaction_threshold: int = 0,
    source_df: Optional[pd.DataFrame] = None,
    exclude_kinds: Iterable[str] = ("water",),
    **_
):
    """
    Adiciona arestas entre pares de nós cuja distância entre centróides <= threshold.
    Se forem da MESMA cadeia e |Δresidue_number| < long_interaction_threshold, NÃO adiciona.
    Por padrão, ignora nós com kind em exclude_kinds.
    """
    # Fonte de dados: prioriza DF em G.graph, senão monta a partir dos nós
    if source_df is None:
        pdb_df = G.graph.get("pdb_df", None)
        if (
            isinstance(pdb_df, pd.DataFrame)
            and {"node_id","x_coord","y_coord","z_coord","chain_id","residue_number"}.issubset(pdb_df.columns)
        ):
            df = pdb_df.loc[pdb_df["node_id"].isin(list(G.nodes()))].copy()
            # garantir que não entra água caso exista anotação
            if "residue_name" in df.columns:
                df = df[df["residue_name"] != "HOH"]
        else:
            df = _ensure_centroids_df(G, exclude_kinds=exclude_kinds)
    else:
        df = source_df.copy()

    if df.empty:
        return G

    coords = df[["x_coord","y_coord","z_coord"]].to_numpy(float)
    ixs, jxs, dists = _pairs_under_threshold(coords, float(threshold))
    node_ids = df["node_id"].to_list()
    chains = df["chain_id"].to_numpy()
    pos = df["residue_number"].to_numpy(int)

    added = 0
    lin = int(long_interaction_threshold)
    for i, j, d in zip(ixs, jxs, dists):
        same_chain = chains[i] == chains[j]
        near_in_seq = abs(pos[i] - pos[j]) < lin if lin > 0 else False
        if not (same_chain and near_in_seq):
            u, v = node_ids[i], node_ids[j]

            if not G.has_edge(u, v):
                G.add_edge(u, v, distance=float(d), kind="distance_threshold")
                added += 1

    return G