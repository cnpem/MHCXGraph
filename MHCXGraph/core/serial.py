# core/serial.py
from __future__ import annotations
from pathlib import Path
import json, numpy as np, pandas as pd, networkx as nx
from networkx.readwrite import json_graph
from types import SimpleNamespace
from datetime import datetime

BASIC = (str, int, float, bool, type(None))

def _np_save(arr: np.ndarray, outdir: Path, key: str):
    out = outdir / f"{key}.npy"
    np.save(out, arr)
    return {"__npy__": out.name}

def _df_save(df: pd.DataFrame, outdir: Path, key: str):
    out = outdir / f"{key}.parquet"
    df.to_parquet(out, index=True)
    return {"__parquet__": out.name}

def _sanitize_value(v, outdir: Path, key: str):
    # básicos
    if isinstance(v, BASIC):
        return v
    # numpy
    if isinstance(v, np.ndarray):
        return _np_save(v, outdir, key)
    # pandas
    if isinstance(v, pd.DataFrame):
        return _df_save(v, outdir, key)
    if isinstance(v, pd.Series):
        return _df_save(v.to_frame(), outdir, key)
    # paths
    if isinstance(v, Path):
        return {"__path__": str(v)}
    # simples containers
    if isinstance(v, (list, tuple)):
        return [_sanitize_value(x, outdir, f"{key}_{i}") for i, x in enumerate(v)]
    if isinstance(v, set):
        return {"__set__": [_sanitize_value(x, outdir, f"{key}_{i}") for i, x in enumerate(sorted(v, key=str))]}
    if isinstance(v, dict):
        return {str(k): _sanitize_value(val, outdir, f"{key}_{str(k)}") for k, val in v.items()}
    # SimpleNamespace → dict
    if isinstance(v, SimpleNamespace):
        return _sanitize_value(vars(v), outdir, key)
    # objetos arbitrários → JSON via repr (último recurso)
    return {"__repr__": repr(v), "__type__": type(v).__name__}

def dump_graph_bundle(G: nx.Graph, outdir: str | Path):
    """
    Salva:
      - graph.json (node-link + atributos sanetizados)
      - assets/* (npy/parquet dos atributos não-JSON)
    """
    outdir = Path(outdir)
    assets = outdir / "assets"
    outdir.mkdir(parents=True, exist_ok=True)
    assets.mkdir(exist_ok=True)

    # 1) node-link puro
    data = dict(json_graph.node_link_data(G))

    # 2) sanitizar atributos do graph (G.graph)
    meta = {}
    for k, v in G.graph.items():
        meta[k] = _sanitize_value(v, assets, f"graph_{k}")
    data["_graph_attrs"] = meta

    # 3) sanitizar atributos de nós
    node_attrs = []
    for i, n in enumerate(G.nodes()):
        d = {}
        for k, v in G.nodes[n].items():
            d[k] = _sanitize_value(v, assets, f"node{i}_{k}")
        node_attrs.append(d)
    data["_node_attrs"] = node_attrs  # paralelo à ordem de data["nodes"]

    # 4) sanitizar atributos de arestas
    edge_attrs = []
    for j, (u, v, d0) in enumerate(G.edges(data=True)):
        d = {}
        for k, v2 in d0.items():
            d[k] = _sanitize_value(v2, assets, f"edge{j}_{k}")
        edge_attrs.append(d)
    data["_edge_attrs"] = edge_attrs  # paralelo à ordem de data["links"]

    # 5) metadados do bundle
    data["_bundle"] = {
        "format": "nx-bundle-v1",
        "created": datetime.utcnow().isoformat() + "Z",
        "networkx": nx.__version__,
    }

    with open(outdir / "graph.json", "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)
    return outdir / "graph.json"

def load_graph_bundle(outdir: str) -> nx.Graph:
    """
    Recarrega o bundle. Para simplificar, os assets não são re-hidratados
    em objetos originais automaticamente: os marcadores (__npy__, __parquet__, …)
    ficam disponíveis para quem precisar carregar on-demand.
    """
    outdir = Path(outdir)
    with open(outdir / "graph.json", "r", encoding="utf-8") as f:
        data = json.load(f)

    G = json_graph.node_link_graph({k: v for k, v in data.items()
                                    if k not in {"_graph_attrs","_node_attrs","_edge_attrs","_bundle"}})

    for k, v in data.get("_graph_attrs", {}).items():
        G.graph[k] = v

    for (n, attrs), extra in zip(G.nodes(data=True), data.get("_node_attrs", [])):
        attrs.update(extra)

    for (u, v, attrs), extra in zip(G.edges(data=True), data.get("_edge_attrs", [])):
        attrs.update(extra)

    return G
