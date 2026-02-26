# core/pipeline.py
from __future__ import annotations
from dataclasses import asdict
from pathlib import Path
import networkx as nx

from MHCXGraph.core.pdb_graph_builder import PDBGraphBuilder
from MHCXGraph.core.config import GraphConfig

def build_graph_with_config(pdb_path: str, config: GraphConfig) -> nx.Graph:
    """
    Build the structural graph from a PDB/mmCIF using the unified GraphConfig.
    Returns a NetworkX graph with pickle-friendly artifacts in G.graph.
    """
    builder = PDBGraphBuilder(pdb_path, config)
    built = builder.build_graph()
    G = built.graph


    G.graph["path"]   = str(pdb_path)
    G.graph["name"]   = Path(pdb_path).stem
    G.graph["config"] = asdict(config)  # plain dict

    G.graph["residue_labels"]  = [nid for nid, _ in built.residue_index]
    G.graph["water_labels"]    = [nid for nid, _ in built.water_index]
    G.graph["distance_matrix"] = built.distance_matrix  # ndarray or None
    G.graph["coords"]          = built.residue_centroids  # ndarray
    G.graph["water_positions"] = built.water_centroids    # ndarray or None

    # DataFrames (great with pickle; not JSON)
    for key in ("raw_pdb_df", "pdb_df", "rgroup_df", "dssp_df"):
        val = getattr(built, key, None)
        if val is not None:
            G.graph[key] = val

    return G
