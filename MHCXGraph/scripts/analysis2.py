from __future__ import annotations

import json

from pathlib import Path
from collections import defaultdict
import numpy as np
import re
from graph.graph import AssociatedGraph
from utils.preprocessing import create_graphs
from itertools import combinations
import pandas as pd
from typing import Dict, Tuple, List, Literal, Optional
import logging
import sys, os

def get_protein_keys(original_graphs: dict):
    keys = list(original_graphs.keys())
    if all(isinstance(k, str) and k.isdigit() for k in keys):
        return [str(i) for i in sorted(map(int, keys))]
    return keys

# ---------- projection helpers ----------
def project_nodes_unique(frame_nodes, p):
    return set(n[p] for n in frame_nodes)

def project_nodes_instances(frame_nodes, p):
    return [n[p] for n in frame_nodes]

def chain_signature(node_tuple):
    """Map an associated node tuple -> chain signature string, e.g. ('A:ARG:23','A:ARG:34') -> 'AA'."""
    chains = []
    for lab in node_tuple:
        # tolerate tuples of tuples or plain strings
        s = lab if isinstance(lab, str) else str(lab)
        chains.append(s.split(":")[0] if ":" in s else s)
    return "".join(chains)

def unique_chain_signatures(frame_nodes):
    """Return a sorted, de-duplicated list of chain signatures present in this frame."""
    sigs = {chain_signature(n) for n in frame_nodes}
    return sorted(sigs)

def chain_combo_key(node_tuple) -> str:
    """e.g., ('A:ARG:23','A:ARG:34') -> 'AA'."""
    return ''.join(str(part).split(':', 1)[0] for part in node_tuple)

# ---------- per-protein node metrics ----------
def  node_similarity_for_protein(frame, original_graphs, protein_keys, p):
    nodes_assoc = frame.get("nodes", [])
    if not nodes_assoc:
        return None

    prot_key = protein_keys[p]
    og = original_graphs[prot_key]
    prot_name = og.get("name", prot_key)

    Vp = set(og["nodes"]) 

    inst = project_nodes_instances(nodes_assoc, p)
    Up   = set(inst)

    total_orig = len(Vp) if Vp else 0
    node_coverage = (len(Up) / total_orig) if total_orig else 0.0

    total_inst = len(inst)
    unique_cnt = len(Up)
    duplication_ratio = (len(Up) / total_inst) if total_inst else 1.0
    duplication_rate  = 1.0 - duplication_ratio
    avg_multiplicity  = (total_inst / unique_cnt) if unique_cnt else float('inf')

    groups = defaultdict(set)  # chain_key -> set of unique residues (for protein p)
    for node_tuple in nodes_assoc:
        key = ''.join(str(part).split(':', 1)[0] for part in node_tuple) # ('A:ARG:23','A:ARG:34') -> 'AA'.
        groups[key].add(node_tuple[p])  # only the residue from protein p

    unique_nodes_per_chain = {k: len(v) for k, v in groups.items()}
    # store as JSON so it round-trips through CSV cleanly
    unique_nodes_per_chain_json = json.dumps(unique_nodes_per_chain, ensure_ascii=False)

    return dict(
        protein_index=p,
        protein_key=prot_key,
        protein_name=prot_name,
        total_nodes_associated=len(nodes_assoc),
        total_nodes_original=total_orig,
        frame_nodes_instances=total_inst,
        frame_nodes_unique=len(Up),
        node_coverage=node_coverage,
        duplication_ratio=duplication_ratio,
        duplication_rate=duplication_rate,
        avg_multiplicity=avg_multiplicity,
        unique_nodes_per_chain=unique_nodes_per_chain_json
    )

def wmean(x, w):
    x = np.asarray(x, float); w = np.asarray(w, float)
    s = w.sum()
    return float(np.sum(x*w)/s) if s > 0 else np.nan

def wstd(x, w):
    x = np.asarray(x, float); w = np.asarray(w, float)
    m = wmean(x, w)
    s = w.sum()
    return float(np.sqrt(np.sum(w*(x-m)**2)/s)) if s > 0 else np.nan

def wmedian(x, w):
    x = np.asarray(x, float); w = np.asarray(w, float)
    if w.sum() == 0: return np.nan
    order = np.argsort(x); x = x[order]; w = w[order]
    cw = np.cumsum(w)/w.sum()
    return float(x[np.searchsorted(cw, 0.5)])

def wtrimmed_mean(x, w, trim=0.10):
    x = np.asarray(x, float); w = np.asarray(w, float)
    if w.sum() == 0: return np.nan
    order = np.argsort(x); x = x[order]; w = w[order]
    cw = np.cumsum(w)/w.sum()
    keep = (cw >= trim) & (cw <= 1.0-trim)
    if not np.any(keep): keep = np.ones_like(cw, dtype=bool)
    return wmean(x[keep], w[keep])

def ivw_mean_proportions(cov, n):
    cov = np.asarray(cov, float); n = np.asarray(n, float)
    p = ((cov*n) + 0.5) / (n + 1.0)
    var = p*(1.0-p) / (n + 1.0) + 1e-12
    w = 1.0/var
    return wmean(p, w)

# ---------- frame-level summary (nodes only) ----------
def summarize_frame_nodes(df_fp_nodes_for_frame):
    if df_fp_nodes_for_frame.empty:
        return {}
    cov = df_fp_nodes_for_frame["node_coverage"].values
    n   = df_fp_nodes_for_frame["total_nodes_original"].values
    w   = n

    return {
        "node_cov_wmean":       wmean(cov, w),
        "node_cov_wmedian":     wmedian(cov, w),
        "node_cov_wtrimmed":    wtrimmed_mean(cov, w, trim=0.10),
        "node_cov_ivw_meta":    ivw_mean_proportions(cov, n),
        "node_cov_wstd":        wstd(cov, w),
        "node_cov_p10":         float(np.percentile(cov, 10)),
        "node_cov_p50":         float(np.percentile(cov, 50)),
        "node_cov_p90":         float(np.percentile(cov, 90)),
        "n_proteins":           int(len(cov)),
        "mean_dup_rate":        float(df_fp_nodes_for_frame.get("duplication_rate", pd.Series([np.nan])).mean()),
        "mean_graph_size":      float(np.mean(n)),
        "sum_graph_size":       int(np.sum(n)),
    }



# ---------- evaluate a single frame ----------
def evaluate_frame_nodes(component_id, frame_id, data):
    comp_key = str(component_id)
    frm_key  = str(frame_id)
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
            node_cov_std=float(df["node_coverage"].std(ddof=0) if len(df)>1 else 0.0),
            dup_ratio_mean=float(df["duplication_ratio"].mean()),
            unique_chain_signatures=chain_sigs,                # NEW: list like ['AA','AC',...]
            n_unique_chain_signatures=n_unique_chain_sigs      # NEW: integer count
        )
    return df, summary

# ---------- evaluate ALL frames (nodes only) ----------
def evaluate_all_frames_nodes(json_path):
    data = json.loads(Path(json_path).read_text())
    component_ids = [k for k in data.keys() if k != "original_graphs"]
    try:
        component_ids = sorted(component_ids, key=lambda x: int(x))
    except Exception:
        pass

    all_fp, summaries = [], []
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
    df_frames_nodes = pd.DataFrame(summaries).sort_values(
        ["node_cov_mean","node_cov_median","node_cov_min","node_cov_std"],
        ascending=[False, False, False, True]
    )
    return df_fp_nodes, df_frames_nodes

def evaluate_all_frames_nodes_weighted(json_path):
    df_fp_nodes, _ = evaluate_all_frames_nodes(json_path)

    summaries = []
    if df_fp_nodes.empty:
        return df_fp_nodes, pd.DataFrame()

    for (comp_id, frame_id), g in df_fp_nodes.groupby(["component_id", "frame_id"], dropna=False):
        s = summarize_frame_nodes(g)
        s.update({"component_id": comp_id, "frame_id": frame_id})
        summaries.append(s)

    df_frames_nodes_w = pd.DataFrame(summaries)
    cols = ["component_id", "frame_id"] + [c for c in df_frames_nodes_w.columns if c not in ("component_id", "frame_id")]
    df_frames_nodes_w = df_frames_nodes_w[cols]

    df_frames_nodes_w = df_frames_nodes_w.sort_values(
        ["node_cov_wmean","node_cov_wmedian","node_cov_p10","node_cov_wstd"],
        ascending=[False, False, False, True]
    )
    return df_fp_nodes, df_frames_nodes_w

# --- Helpers ---------------------------------------------------------------


def _make_json_from_associated_graph(G, out_json: Path) -> None:
    """
    Serialize an AssociatedGraph into the expected JSON format.

    Parameters
    ----------
    G : AssociatedGraph
        Instance already built.
    out_json : pathlib.Path
        Output path (will be created/overwritten).
    """
    graphs_raw = G.graph_data
    payload: Dict = {"original_graphs": {}}

    for graph_raw in graphs_raw:
        pdb_file = graph_raw["pdb_file"]
        _id = graph_raw["id"]

        m = re.search(r'noTCR_([A-Za-z0-9]{4})\.trunc', pdb_file, re.IGNORECASE)
        name = m[1] if m else f"id{_id}"

        nodes = list(graph_raw["graph"].nodes)
        edges = list(graph_raw["graph"].edges)
        neighbors = {str(n): [str(nb) for nb in graph_raw["graph"].neighbors(n)] for n in nodes}

        payload["original_graphs"][_id] = {
            "name": name,
            "nodes": nodes,
            "edges": edges,
            "neighbors": neighbors,
        }

    for j, comps in enumerate(G.associated_graphs):
        payload[j] = {"comp": j, "frames": {}}
        for i in range(len(comps[0])):
            nodes = list(comps[0][i].nodes)
            edges = list(comps[0][i].edges)
            neighbors = {str(n): [str(nb) for nb in comps[0][i].neighbors(n)] for n in nodes}
            payload[j]["frames"][i] = {"nodes": nodes, "edges": edges, "neighbors": neighbors}

    out_json.parent.mkdir(parents=True, exist_ok=True)
    with open(out_json, "w") as f:
        json.dump(payload, f, indent=4)


def _save_eval_tables(out_dir: Path, df_fp_nodes: pd.DataFrame, df_frames_nodes_w: pd.DataFrame) -> None:
    """
    Save per-run evaluation tables.

    Parameters
    ----------
    out_dir : pathlib.Path
        Destination directory.
    df_fp_nodes : pandas.DataFrame
    df_frames_nodes_w : pandas.DataFrame
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    if not df_fp_nodes.empty:
        df_fp_nodes.to_csv(out_dir / "nodes_per_protein.csv", index=False)
    if not df_frames_nodes_w.empty:
        cols = list(df_frames_nodes_w.columns)
        for lead in ["component_id", "frame_id"]:
            if lead in cols:
                cols = [lead] + [c for c in cols if c != lead]
        df_frames_nodes_w[cols].to_csv(out_dir / "nodes_summary_weighted.csv", index=False)


def _build_associated_graph(files_name: str, run_name: str, out_dir: Path, manifest: Dict):
    """
    Build graphs + AssociatedGraph and persist JSON.

    Parameters
    ----------
    files_name : str
        Comma-separated file list passed to `args.files_name`.
    run_name : str
        Run identifier.
    out_dir : pathlib.Path
        Output directory.
    args : Any
        Argument namespace consumed by `create_graphs`.
    association_config : dict
        Passed to AssociatedGraph.

    Returns
    -------
    (G, json_path) : tuple
        AssociatedGraph instance and JSON path written.
    """
    manifest["inputs"] = [
        {
        "path": files_name,
        "enable_tui": False,
        "extensions": [".pdb"],
        "constrains": [
            { "name": "MHC1"}
        ]
        }
    ]
    manifest["settings"]["run_name"] = run_name
    manifest["settings"]["output_path"] = str(out_dir)

    graphs = create_graphs(manifest)
    S = manifest["settings"]
    checks = {
        "depth": S.get("check_depth"),
        "rsa":   S.get("check_rsa"),
    }
    association_config = {
        "centroid_threshold":          S.get("centroid_threshold"),
        "distance_diff_threshold":     S.get("distance_diff_threshold"),
        "rsa_filter":                  S.get("rsa_filter"),
        "depth_filter":                S.get("depth_filter"),
        "rsa_bins":                    S.get("rsa_bins"),
        "depth_bins":                  S.get("depth_bins"),
        "distance_bins":               S.get("distance_bins"),
        "checks":                      checks,
        "exclude_waters":              S.get("exclude_waters"),
        "classes": S.get("classes", {}),
    }
    G = AssociatedGraph(
        graphs=graphs,
        output_path=manifest["settings"]["output_path"],
        run_name=manifest["settings"]["run_name"],
        association_config=association_config
    )
    # Optional figures
    G.draw_graph_interactive(show=False, save=True)
    G.align_all_frames()
    G.create_pdb_per_protein()

    json_path = out_dir / f"graph_{run_name}.json"
    _make_json_from_associated_graph(G, json_path)
    return G


# --- Public API ------------------------------------------------------------

def run_allxall_per_group(
    cross_df: pd.DataFrame,
    manifest: dict,
    root: str = "Analysis/CrossGraphs"
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Run the "All × All" flow once per TCR_pair_id and aggregate outputs.

    Parameters
    ----------
    cross_df : pandas.DataFrame
        Must contain columns 'TCR_pair_id' and 'PDB_ID'.
    args : Any
        Namespace consumed by `create_graphs`.
    association_config : dict
        Passed to AssociatedGraph.
    root : str, default "Analysis/CrossGraphs"
        Root output directory.

    Returns
    -------
    (df_all_fp, df_all_frames) : tuple of DataFrame
        Global aggregates across groups (can be empty DataFrames).
    """
    root = Path(root)
    all_fp, all_frames = [], []

    for pair_id, group in cross_df.groupby("TCR_pair_id"):
        group_dir = root / str(pair_id) / "All"
        pdb_ids = [str(x).strip() for x in group["PDB_ID"]]
        files = [f"Analysis/selected_strs_renumber/without_TCR/noTCR_{pid.lower()}.trunc.fit_renum.pdb" for pid in pdb_ids]

        files_name = ",".join(files)
        run_name = str(pair_id)

        json_path = group_dir / f"graph_{run_name}.json"
        try:
            _ = _build_associated_graph(files_name, run_name, group_dir, manifest)
        except Exception as e:
            print(f"[SKIP] Could not build group {pair_id}: {e}")
            continue

        try:
            df_fp_nodes, df_frames_nodes_w = evaluate_all_frames_nodes_weighted(str(json_path))
            if "frame_nodes_unique" in df_fp_nodes.columns and "total_nodes_associated" not in df_fp_nodes.columns:
                df_fp_nodes["total_nodes_associated"] = df_fp_nodes["frame_nodes_unique"]

        except Exception as e:
            print(f"[SKIP] Evaluation failed for {pair_id}: {e}")
            continue


        out_dir = group_dir
        out_dir.mkdir(parents=True, exist_ok=True)
        if not df_fp_nodes.empty:
            df_fp_nodes.to_csv(out_dir / "nodes_per_protein.csv", index=False)
            df_fp_nodes.insert(0, "pair_id", pair_id)
            all_fp.append(df_fp_nodes)
        if not df_frames_nodes_w.empty:
            cols = list(df_frames_nodes_w.columns)
            for lead in ["component_id", "frame_id"]:
                if lead in cols:
                    cols = [lead] + [c for c in cols if c != lead]
            df_frames_nodes_w[cols].to_csv(out_dir / "nodes_summary_weighted.csv", index=False)
            df_frames_nodes_w.insert(0, "pair_id", pair_id)
            all_frames.append(df_frames_nodes_w)

        _save_eval_tables(group_dir, df_fp_nodes, df_frames_nodes_w)

    df_all_fp = pd.concat(all_fp, ignore_index=True) if all_fp else pd.DataFrame()
    df_all_frames = pd.concat(all_frames, ignore_index=True) if all_frames else pd.DataFrame()

    if not df_all_fp.empty:
        df_all_fp.to_csv(root / "ALL_nodes_per_protein.csv", index=False)
    if not df_all_frames.empty:
        lead = ["pair_id", "component_id", "frame_id"]
        cols = lead + [c for c in df_all_frames.columns if c not in lead]
        df_all_frames = df_all_frames[cols]
        df_all_frames.to_csv(root / "ALL_nodes_summary_weighted.csv", index=False)

    print("All×All completed.")
    return df_all_fp, df_all_frames


def run_pairwise_per_group(
    cross_df: pd.DataFrame,
    manifest: dict,
    root: str = "Analysis/CrossGraphs",
    score_column: str = "node_cov_wmean"
) -> Dict[str, Dict[str, pd.DataFrame]]:
    """
    Run the Pairwise flow (within each TCR_pair_id) and build similarity matrices.
    For each pair, also save the standard per-run evaluation tables.

    Parameters
    ----------
    cross_df : pandas.DataFrame
        Must contain columns 'TCR_pair_id' and 'PDB_ID'.
    args : Any
        Namespace consumed by `create_graphs`.
    association_config : dict
        Passed to AssociatedGraph.
    root : str, default "Analysis/CrossGraphs"
        Root output directory.
    score_column : str, default "node_cov_wmean"
        Frame-level column used to score a pair.

    Returns
    -------
    dict
        Mapping pair_id -> {"max": DataFrame, "mean": DataFrame}
    """
    out: Dict[str, Dict[str, pd.DataFrame]] = {}

    for pair_id, group in cross_df.groupby("TCR_pair_id"):
        refs = [str(x).strip() for x in group["PDB_ID"]]
        files_map = {r: f"Analysis/selected_strs_renumber/without_TCR/noTCR_{r.lower()}.trunc.fit_renum.pdb" for r in refs}

        idx = refs
        M_max = pd.DataFrame(0.0, index=idx, columns=idx)
        M_mean = pd.DataFrame(0.0, index=idx, columns=idx)

        group_dir = Path(root) / str(pair_id) / "Pairs"
        group_dir.mkdir(parents=True, exist_ok=True)

        for r1, r2 in combinations(refs, 2):
            files_name = ",".join([files_map[r1], files_map[r2]])
            run_name = f"{pair_id}_{r1}_{r2}"
            out_dir = group_dir / f"{r1}_{r2}"
            json_path = out_dir / f"graph_{run_name}.json"
            try:
                _ = _build_associated_graph(files_name, run_name, out_dir, manifest)
            except Exception as e:
                print(f"[SKIP] Could not build pair ({r1},{r2}): {e}")
                continue

            try:
                df_fp_nodes, df_frames_nodes_w = evaluate_all_frames_nodes_weighted(str(json_path))
                if "frame_nodes_unique" in df_fp_nodes.columns and "total_nodes_associated" not in df_fp_nodes.columns:
                    df_fp_nodes["total_nodes_associated"] = df_fp_nodes["frame_nodes_unique"]

            except Exception as e:
                print(f"[SKIP] Evaluation failed for pair ({r1},{r2}): {e}")
                continue

            # Save standard per-run outputs for the pair
            _save_eval_tables(out_dir, df_fp_nodes, df_frames_nodes_w)

            # Pair score from frames
            if not df_frames_nodes_w.empty and score_column in df_frames_nodes_w.columns:
                max_score = float(df_frames_nodes_w[score_column].max())
                mean_score = float(df_frames_nodes_w[score_column].mean())
            else:
                max_score = 0.0
                mean_score = 0.0

            M_max.loc[r1, r2] = M_max.loc[r2, r1] = max_score
            M_mean.loc[r1, r2] = M_mean.loc[r2, r1] = mean_score
            M_max.loc[r1, r1] = M_max.loc[r2, r2] = 1.0
            M_mean.loc[r1, r1] = M_mean.loc[r2, r2] = 1.0

        # Save matrices for the group
        M_max.to_csv(group_dir / f"matrix_{score_column}_MAX.csv")
        M_mean.to_csv(group_dir / f"matrix_{score_column}_MEAN.csv")
        print(f"Pairwise matrices saved in {group_dir}")

        out[str(pair_id)] = {"max": M_max, "mean": M_mean}

    print("Pairwise completed.")
    return out


def run_cross_analysis(
    mode: Literal["all", "pairwise", "both"],
    cross_df: pd.DataFrame,
    manifest: dict,
    root: str = "Analysis/CrossGraphs",
    score_column: str = "node_cov_wmean"
):
    """
    Unified entry point for running All×All, Pairwise, or both.

    Parameters
    ----------
    mode : {"all", "pairwise", "both"}
        Which analysis to run.
    cross_df : pandas.DataFrame
        Must contain columns 'TCR_pair_id' and 'PDB_ID'.
    args : Any
        Namespace consumed by `create_graphs`.
    association_config : dict
        Passed to AssociatedGraph.
    root : str, default "Analysis/CrossGraphs"
        Root output directory.
    score_column : str, default "node_cov_wmean"
        Frame-level column used to score a pair (Pairwise mode).

    Returns
    -------
    dict
        Results container. Keys present depend on the selected mode:
        - "all": {"df_all_fp": DataFrame, "df_all_frames": DataFrame}
        - "pairwise": {"matrices": dict(pair_id -> {"max": df, "mean": df})}
        - "both": union of the above.
    """
    results: Dict[str, object] = {}

    if mode in ("all", "both"):
        df_all_fp, df_all_frames = run_allxall_per_group(cross_df, manifest, root=root)
        results["all"] = {"df_all_fp": df_all_fp, "df_all_frames": df_all_frames}

    if mode in ("pairwise", "both"):
        matrices = run_pairwise_per_group(cross_df, manifest, root=root, score_column=score_column)
        results["pairwise"] = {"matrices": matrices}

    return results


manifest = {
  "settings": {
    "run_name": None,
    "output_path": None,
    "debug": True,
    "track_steps": False,
    "centroid_threshold": 10.0,
    "centroid_granularity": "ca_only",
    "exclude_waters": False,
    "check_rsa": True,
    "check_depth": False,
    "rsa_filter": 0.1,
    "depth_filter": 10.0,
    "distance_diff_threshold": 3.0,
    "rsa_table": "Wilke",

    "distance_bins": 3,
    "rsa_bins": 3,
    "depth_bins": 3,

    "serd_config": None,

    "classes_": {
      "residues": {
        "HID": ["ALA", "VAL", "LEU", "ILE", "MET"],
        "POL": ["SER", "THR", "ASN", "GLN", "CYS"],
        "POS": ["LYS", "ARG", "HIS"],
        "NEG": ["ASP", "GLU"],
        "ARO": ["PHE", "TYR", "TRP"],
        "ESP": ["GLY", "PRO"]
      }
    }
  },

  "inputs": [],
  "constrains": {
    "MHC1": {
      "chains": ["C"],
      "residues": {
        "A": [18, 19, 42, 43, 44, 54, 55, 56, 58, 59, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 79, 80, 81, 83, 84, 85, 89, 108, 109, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171]
      }
    },
    "MHC2": {
      "chains": ["C"],
      "residues": {
        "A": [39, 53, 54, 55, 57, 58, 60, 61, 62, 64, 65, 67, 68, 69, 71]
      }
    }
  }
}

diffMHC_diffPep = pd.read_csv("Analysis/crossreact_processed_diff_MHC_diff_pep_helder.csv")
diffMHC_SamePep = pd.read_csv("Analysis/crossreact_processed_diff_MHC_same_pep_helder.csv")
sameMHC_diffPep = pd.read_csv("Analysis/crossreact_processed_same_MHC_diff_pep_helder.csv")
rawCross = pd.read_csv("Analysis/crossreact_tcrs_v4.csv")

crossDf = rawCross[["TCR_ID", "TCR_pair_id", "MHC_allele_id", "peptide_id", "pMHC_id", "PDB_ID"]]

S = manifest["settings"]
logging.getLogger("matplotlib").setLevel(logging.ERROR)
logging.basicConfig(
    stream=sys.stdout,
    level=logging.DEBUG if S.get("debug", False) else logging.INFO
)
log = logging.getLogger("CRSProtein")
log.setLevel(logging.DEBUG if S.get("debug", False) else logging.INFO)


res_both = run_cross_analysis(
    mode="all",
    cross_df=crossDf,
    manifest=manifest,
    root="Analysis/CrossGraphs_10_Multiclass"
)
