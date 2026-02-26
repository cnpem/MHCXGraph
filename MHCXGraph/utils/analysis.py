from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd



def _make_json_from_associated_graph(G: AssociatedGraph, out_json: Path) -> None:
    """
    Build a JSON payload from an AssociatedGraph instance.

    The payload contains:
    - "original_graphs": mapping from internal id to original graph data
      (name, nodes, peptide_nodes, mhc_nodes, edges, neighbors).
    - For each component:
      - frames with associated graph nodes grouped into peptide, mhc and mixed nodes,
        plus edges and neighbors.

    Parameters
    ----------
    G : AssociatedGraph
        AssociatedGraph instance already computed.
    out_json : Path
        Output path for the JSON file.
    """
    graphs_raw = G.graphs_data
    payload: Dict[str, Any] = {"original_graphs": {}}

    for graph_raw in graphs_raw:
        pdb_file = graph_raw["pdb_file"]
        _id = graph_raw["id"]
        name = Path(pdb_file).stem

        nodes = list(graph_raw["graph"].nodes)
        peptide_nodes = [
            node for node in nodes
            if isinstance(node, str) and node.split(":", 1)[0] == "C"
        ]
        mhc_nodes = [
            node for node in nodes
            if isinstance(node, str) and node.split(":", 1)[0] == "A"
        ]

        edges = list(graph_raw["graph"].edges)
        neighbors = {
            str(n): [str(nb) for nb in graph_raw["graph"].neighbors(n)]
            for n in nodes
        }

        payload["original_graphs"][_id] = {
            "name": name,
            "nodes": nodes,
            "peptide_nodes": peptide_nodes,
            "mhc_nodes": mhc_nodes,
            "edges": edges,
            "neighbors": neighbors,
        }

    for j, comps in enumerate(G.associated_graphs):
        payload[j] = {"comp": j, "frames": {}}
        for i in range(len(comps[0])):
            nodes = list(comps[0][i].nodes)

            peptide_nodes = [
                node for node in nodes
                if isinstance(node, tuple)
                and all(str(node_).startswith("C") for node_ in node)
            ]
            mhc_nodes = [
                node for node in nodes
                if isinstance(node, tuple)
                and all(str(node_).startswith("A") for node_ in node)
            ]
            mixed_nodes = [
                node for node in nodes
                if isinstance(node, tuple)
                and any(str(node_).startswith("C") for node_ in node)
                and not all(str(node_).startswith("C") for node_ in node)
                and not all(str(node_).startswith("A") for node_ in node)
            ]

            edges = list(comps[0][i].edges)
            neighbors = {
                str(n): [str(nb) for nb in comps[0][i].neighbors(n)]
                for n in nodes
            }

            payload[j]["frames"][i] = {
                "nodes": nodes,
                "peptide_nodes": peptide_nodes,
                "mhc_nodes": mhc_nodes,
                "mixed_nodes": mixed_nodes,
                "edges": edges,
                "neighbors": neighbors,
            }

    out_json.parent.mkdir(parents=True, exist_ok=True)
    with open(out_json, "w") as f:
        json.dump(payload, f, indent=4)


def get_protein_keys(original_graphs: Dict[Any, Any]) -> List[str]:
    """
    Return a list of protein keys sorted numerically if keys are numeric strings.

    Parameters
    ----------
    original_graphs : dict
        Mapping from internal id to original graph data.

    Returns
    -------
    list of str
        Ordered keys for proteins.
    """
    keys = list(original_graphs.keys())
    if all(isinstance(k, str) and k.isdigit() for k in keys):
        return [str(i) for i in sorted(map(int, keys))]
    return keys


def project_nodes_instances(frame_nodes: List[Any], p: int) -> List[str]:
    """
    Project associated nodes onto the p-th protein.

    Parameters
    ----------
    frame_nodes : list
        List of associated graph nodes as tuples of residue labels.
    p : int
        Index of the protein to project.

    Returns
    -------
    list of str
        Residue labels for protein p.
    """
    return [n[p] for n in frame_nodes]


def chain_signature(node_tuple: tuple) -> str:
    """
    Compute a chain signature for a tuple of node labels.

    Parameters
    ----------
    node_tuple : tuple
        Tuple of residue labels for an associated node.

    Returns
    -------
    str
        Concatenated chain identifiers.
    """
    chains = []
    for lab in node_tuple:
        s = lab if isinstance(lab, str) else str(lab)
        chains.append(s.split(":", 1)[0] if ":" in s else s)
    return "".join(chains)


def unique_chain_signatures(frame_nodes: List[tuple]) -> List[str]:
    """
    Compute sorted unique chain signatures for all nodes in a frame.

    Parameters
    ----------
    frame_nodes : list of tuple
        List of associated nodes as tuples of residue labels.

    Returns
    -------
    list of str
        Sorted unique chain signatures.
    """
    return sorted({chain_signature(n) for n in frame_nodes})


def node_similarity_for_protein(
    frame: Dict[str, Any],
    original_graphs: Dict[str, Any],
    protein_keys: List[str],
    p: int,
) -> Optional[Dict[str, Any]]:
    """
    Compute node coverage metrics for a single protein in one frame.

    Parameters
    ----------
    frame : dict
        Frame entry from the JSON payload built from AssociatedGraph.
    original_graphs : dict
        Mapping from protein key to original graph data.
    protein_keys : list of str
        Ordered protein keys.
    p : int
        Index of the protein to evaluate.

    Returns
    -------
    dict or None
        Coverage metrics for this protein and frame, or None if there are no nodes.
    """
    nodes_assoc = frame.get("nodes", [])
    if not nodes_assoc:
        return None

    prot_key = protein_keys[p]
    og = original_graphs[prot_key]
    prot_name = og.get("name", prot_key)

    Vp = set(og["nodes"])

    pep_orig_nodes = [
        node for node in Vp
        if isinstance(node, str) and node.startswith("C")
    ]
    mhc_orig_nodes = [
        node for node in Vp
        if isinstance(node, str) and node.startswith("A")
    ]

    inst = project_nodes_instances(nodes_assoc, p)
    Up = set(inst)

    total_orig = len(Vp) if Vp else 0
    node_coverage = (len(Up) / total_orig) if total_orig else 0.0

    total_inst = len(inst)
    unique_cnt = len(Up)
    duplication_ratio = (len(Up) / total_inst) if total_inst else 1.0
    duplication_rate = 1.0 - duplication_ratio
    avg_multiplicity = (total_inst / unique_cnt) if unique_cnt else float("inf")

    groups = defaultdict(set)
    for node_tuple in nodes_assoc:
        key = "".join(str(part).split(":", 1)[0] for part in node_tuple)
        groups[key].add(node_tuple[p])

    unique_nodes_per_chain = {k: len(v) for k, v in groups.items()}
    unique_nodes_per_chain_json = json.dumps(unique_nodes_per_chain, ensure_ascii=False)

    return dict(
        protein_index=p,
        protein_key=prot_key,
        protein_name=prot_name,
        total_nodes_associated=len(nodes_assoc),
        total_nodes_original=total_orig,
        pep_orig_nodes=len(pep_orig_nodes),
        mhc_orig_nodes=len(mhc_orig_nodes),
        frame_nodes_instances=total_inst,
        frame_nodes_unique=len(Up),
        node_coverage=node_coverage,
        duplication_ratio=duplication_ratio,
        duplication_rate=duplication_rate,
        avg_multiplicity=avg_multiplicity,
        unique_nodes_per_chain=unique_nodes_per_chain_json,
    )


def wmean(x, w):
    """
    Weighted mean.

    Parameters
    ----------
    x : array_like
        Data values.
    w : array_like
        Weights.

    Returns
    -------
    float
        Weighted mean.
    """
    x = np.asarray(x, float)
    w = np.asarray(w, float)
    s = w.sum()
    return float(np.sum(x * w) / s) if s > 0 else np.nan


def wstd(x, w):
    """
    Weighted standard deviation.

    Parameters
    ----------
    x : array_like
        Data values.
    w : array_like
        Weights.

    Returns
    -------
    float
        Weighted standard deviation.
    """
    x = np.asarray(x, float)
    w = np.asarray(w, float)
    m = wmean(x, w)
    s = w.sum()
    return float(np.sqrt(np.sum(w * (x - m) ** 2) / s)) if s > 0 else np.nan


def wmedian(x, w):
    """
    Weighted median.

    Parameters
    ----------
    x : array_like
        Data values.
    w : array_like
        Weights.

    Returns
    -------
    float
        Weighted median.
    """
    x = np.asarray(x, float)
    w = np.asarray(w, float)
    if w.sum() == 0:
        return np.nan
    order = np.argsort(x)
    x = x[order]
    w = w[order]
    cw = np.cumsum(w) / w.sum()
    return float(x[np.searchsorted(cw, 0.5)])


def wtrimmed_mean(x, w, trim=0.10):
    """
    Weighted trimmed mean removing tails.

    Parameters
    ----------
    x : array_like
        Data values.
    w : array_like
        Weights.
    trim : float, optional
        Fraction to trim at each tail, by default 0.10.

    Returns
    -------
    float
        Weighted trimmed mean.
    """
    x = np.asarray(x, float)
    w = np.asarray(w, float)
    if w.sum() == 0:
        return np.nan
    order = np.argsort(x)
    x = x[order]
    w = w[order]
    cw = np.cumsum(w) / w.sum()
    keep = (cw >= trim) & (cw <= 1.0 - trim)
    if not np.any(keep):
        keep = np.ones_like(cw, dtype=bool)
    return wmean(x[keep], w[keep])


def ivw_mean_proportions(cov, n):
    """
    Inverse variance weighted mean for proportions with shrinkage.

    Parameters
    ----------
    cov : array_like
        Coverage values between 0 and 1.
    n : array_like
        Sample sizes.

    Returns
    -------
    float
        Weighted mean proportion estimate.
    """
    cov = np.asarray(cov, float)
    n = np.asarray(n, float)
    p = ((cov * n) + 0.5) / (n + 1.0)
    var = p * (1.0 - p) / (n + 1.0) + 1e-12
    w = 1.0 / var
    return wmean(p, w)


def summarize_frame_nodes(df_fp_nodes_for_frame: pd.DataFrame) -> Dict[str, Any]:
    """
    Compute weighted summaries for node coverage across proteins in a frame.

    Parameters
    ----------
    df_fp_nodes_for_frame : pandas.DataFrame
        Per protein node coverage for one frame.

    Returns
    -------
    dict
        Summary statistics including weighted mean, median and dispersion.
    """
    if df_fp_nodes_for_frame.empty:
        return {}
    cov = df_fp_nodes_for_frame["node_coverage"].values
    n = df_fp_nodes_for_frame["total_nodes_original"].values
    w = n
    return {
        "node_cov_wmean": wmean(cov, w),
        "node_cov_wmedian": wmedian(cov, w),
        "node_cov_wtrimmed": wtrimmed_mean(cov, w, trim=0.10),
        "node_cov_ivw_meta": ivw_mean_proportions(cov, n),
        "node_cov_wstd": wstd(cov, w),
        "node_cov_p10": float(np.percentile(cov, 10)),
        "node_cov_p50": float(np.percentile(cov, 50)),
        "node_cov_p90": float(np.percentile(cov, 90)),
        "n_proteins": int(len(cov)),
        "mean_dup_rate": float(
            df_fp_nodes_for_frame.get(
                "duplication_rate", pd.Series([np.nan])
            ).mean()
        ),
        "mean_graph_size": float(np.mean(n)),
        "sum_graph_size": int(np.sum(n)),
    }


def evaluate_frame_nodes(
    component_id: Any,
    frame_id: Any,
    data: Dict[str, Any],
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """
    Evaluate node coverage metrics for one frame across all proteins.

    Parameters
    ----------
    component_id : hashable
        Component identifier.
    frame_id : hashable
        Frame identifier.
    data : dict
        JSON payload as built by _make_json_from_associated_graph.

    Returns
    -------
    df : pandas.DataFrame
        Per protein metrics for this frame.
    summary : dict
        Aggregated summary for the frame.
    """
    comp_key = str(component_id)
    frm_key = str(frame_id)
    frame = data[comp_key]["frames"][frm_key]

    original_graphs = data["original_graphs"]
    protein_keys = get_protein_keys(original_graphs)

    nodes_assoc = frame.get("nodes", [])
    chain_sigs = unique_chain_signatures(nodes_assoc)
    n_unique_chain_sigs = len(chain_sigs)

    rows = []
    if nodes_assoc:
        n_prot = len(nodes_assoc[0])
        for p in range(n_prot):
            r = node_similarity_for_protein(frame, original_graphs, protein_keys, p)
            if r is not None:
                rows.append(r)

    df = pd.DataFrame(rows)
    if df.empty:
        summary = {}
    else:
        summary = dict(
            component=component_id,
            frame=frame_id,
            n_proteins=len(df),
            node_cov_mean=float(df["node_coverage"].mean()),
            node_cov_median=float(df["node_coverage"].median()),
            node_cov_min=float(df["node_coverage"].min()),
            node_cov_std=float(
                df["node_coverage"].std(ddof=0) if len(df) > 1 else 0.0
            ),
            dup_ratio_mean=float(df["duplication_ratio"].mean()),
            unique_chain_signatures=chain_sigs,
            n_unique_chain_signatures=n_unique_chain_sigs,
        )
    return df, summary


def evaluate_all_frames_nodes(
    json_path: Path,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Evaluate node coverage for all components and frames from a JSON file.

    Parameters
    ----------
    json_path : Path
        Path to the JSON produced by _make_json_from_associated_graph.

    Returns
    -------
    df_fp_nodes : pandas.DataFrame
        Per frame and per protein metrics.
    df_frames_nodes : pandas.DataFrame
        Frame level aggregated metrics.
    """
    data = json.loads(json_path.read_text())
    component_ids = [k for k in data.keys() if k != "original_graphs"]

    try:
        component_ids = sorted(component_ids, key=lambda x: int(x))
    except Exception:
        pass

    all_fp: List[pd.DataFrame] = []
    summaries: List[Dict[str, Any]] = []
    for comp_id in component_ids:
        frames = data[comp_id]["frames"]
        frame_ids = list(frames.keys())
        try:
            frame_ids = sorted(frame_ids, key=lambda x: int(x))
        except Exception:
            pass

        for frm_id in frame_ids:
            df_fp, summ = evaluate_frame_nodes(comp_id, frm_id, data)
            if not df_fp.empty:
                df_fp.insert(0, "component_id", comp_id)
                df_fp.insert(1, "frame_id", frm_id)
                all_fp.append(df_fp)
            if summ:
                summaries.append(summ)

    df_fp_nodes = pd.concat(all_fp, ignore_index=True) if all_fp else pd.DataFrame()

    required = {"node_cov_mean", "node_cov_median", "node_cov_min", "node_cov_std"}

    summaries = [
        s for s in summaries
        if required.issubset(s.keys())
    ]

    df_frames_nodes = pd.DataFrame(summaries)

    if not df_frames_nodes.empty:
        df_frames_nodes = df_frames_nodes.sort_values(
            ["node_cov_mean", "node_cov_median", "node_cov_min", "node_cov_std"],
            ascending=[False, False, False, True],
        )

    return df_fp_nodes, df_frames_nodes


def evaluate_all_frames_nodes_weighted(
    json_path: Path,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Evaluate weighted node coverage summaries for all frames.

    Parameters
    ----------
    json_path : Path
        Path to the JSON produced by _make_json_from_associated_graph.

    Returns
    -------
    df_fp_nodes : pandas.DataFrame
        Per frame and per protein coverage metrics.
    df_frames_nodes_w : pandas.DataFrame
        Weighted summaries per frame.
    """
    df_fp_nodes, _ = evaluate_all_frames_nodes(json_path)

    summaries: List[Dict[str, Any]] = []
    if df_fp_nodes.empty:
        return df_fp_nodes, pd.DataFrame()

    for (comp_id, frame_id), g in df_fp_nodes.groupby(
        ["component_id", "frame_id"], dropna=False
    ):
        s = summarize_frame_nodes(g)
        s.update({"component_id": comp_id, "frame_id": frame_id})
        summaries.append(s)

    df_frames_nodes_w = pd.DataFrame(summaries)
    cols = ["component_id", "frame_id"] + [
        c for c in df_frames_nodes_w.columns if c not in ("component_id", "frame_id")
    ]
    df_frames_nodes_w = df_frames_nodes_w[cols]

    df_frames_nodes_w = df_frames_nodes_w.sort_values(
        ["node_cov_wmean", "node_cov_wmedian", "node_cov_p10", "node_cov_wstd"],
        ascending=[False, False, False, True],
    )
    return df_fp_nodes, df_frames_nodes_w


def _save_eval_tables(
    out_dir: Path,
    df_fp_nodes: pd.DataFrame,
    df_frames_nodes_w: pd.DataFrame,
) -> None:
    """
    Save evaluation tables as CSV files in a given directory.

    Parameters
    ----------
    out_dir : Path
        Output directory.
    df_fp_nodes : pandas.DataFrame
        Per frame and per protein metrics.
    df_frames_nodes_w : pandas.DataFrame
        Frame level weighted summaries.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    if not df_fp_nodes.empty:
        df_fp_nodes.to_csv(out_dir / "nodes_per_protein.csv", index=False)
    if not df_frames_nodes_w.empty:
        cols = list(df_frames_nodes_w.columns)
        for lead in ["component_id", "frame_id"]:
            if lead in cols:
                cols = [lead] + [c for c in cols if c != lead]
        df_frames_nodes_w[cols].to_csv(
            out_dir / "nodes_summary_weighted.csv", index=False
        )

