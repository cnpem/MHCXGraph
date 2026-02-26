#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import logging
import sys
import re
from pathlib import Path
from collections import defaultdict
from itertools import combinations
from typing import Dict, Tuple, List, Literal, Optional, Any

import numpy as np
import pandas as pd

# Dependências do seu projeto
from graph.graph import AssociatedGraph
from utils.preprocessing import create_graphs


# =============================================================================
# Leitura do manifest no seu padrão
# =============================================================================
def load_manifest(manifest_path: str) -> Dict[str, Any]:
    if not manifest_path:
        return {}

    with open(manifest_path, "r") as f:
        data = json.load(f)

    data.setdefault("settings", {})
    data.setdefault("inputs", [])
    data.setdefault("constrains", {})

    S = data["settings"]
    S.setdefault("run_name", "test")
    S.setdefault("output_path", "./outputs")
    S.setdefault("debug", False)
    S.setdefault("track_steps", False)
    S.setdefault("rsa_table", "Wilke")

    S.setdefault("centroid_threshold", 8.5)
    S.setdefault("centroid_granularity", "all_atoms")
    S.setdefault("exclude_waters", True)

    S.setdefault("check_rsa", True)
    S.setdefault("check_depth", True)

    S.setdefault("rsa_filter", 0.1)
    S.setdefault("depth_filter", 10.0)

    S.setdefault("distance_diff_threshold", 2.0)

    S.setdefault("rsa_bins", 5)
    S.setdefault("depth_bins", 5)
    S.setdefault("distance_bins", 5)

    S.setdefault("serd_config", None)

    S.setdefault("classes", {})

    return data


# =============================================================================
# Scan de arquivos e construção de strings de input
# =============================================================================
def find_struct_files(dir_path: Path, patterns: List[str] = ["*.pdb", "*.cif"]) -> List[Path]:
    if not dir_path.exists() or not dir_path.is_dir():
        raise FileNotFoundError(f"Pasta não encontrada: {dir_path}")
    files: List[Path] = []
    for pat in patterns:
        files.extend(dir_path.rglob(pat))
    files = sorted(set(files))
    if not files:
        raise RuntimeError(f"Nenhum arquivo .pdb ou .cif encontrado em {dir_path}")
    return files

def build_label_map(files: List[Path]) -> Dict[str, str]:
    """
    Retorna label -> caminho, usando stem como base.
    Resolve colisões com sufixos _2, _3, ...
    """
    out: Dict[str, str] = {}
    counts: Dict[str, int] = defaultdict(int)
    for p in files:
        base = p.stem
        counts[base] += 1
        label = base if counts[base] == 1 else f"{base}_{counts[base]}"
        out[label] = str(p)
    return out

def files_to_comma_string(files: List[Path]) -> str:
    return ",".join(str(p) for p in files)


# =============================================================================
# Serialização do AssociatedGraph em JSON próprio do pipeline
# =============================================================================
def _make_json_from_associated_graph(G: AssociatedGraph, out_json: Path) -> None:
    """
    Sem noTCR, sem IDs mágicos. Usa stem do arquivo como 'name'.
    """
    graphs_raw = G.graph_data
    payload: Dict = {"original_graphs": {}}

    for graph_raw in graphs_raw:
        pdb_file = graph_raw["pdb_file"]
        _id = graph_raw["id"]  # chave interna mantida
        name = Path(pdb_file).stem

        nodes = list(graph_raw["graph"].nodes)
        peptide_nodes = [node for node in nodes if node.split(":")[0] == "C"]
        mhc_nodes = [node for node in nodes if node.split(":")[0] == "A"]

        edges = list(graph_raw["graph"].edges)
        neighbors = {str(n): [str(nb) for nb in graph_raw["graph"].neighbors(n)] for n in nodes}

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
            peptide_nodes = [node for node in nodes if all(node_.startswith("C") for node_ in node)]
            mhc_nodes = [node for node in nodes if all(node_.startswith("A") for node_ in node)]
            mixed_nodes = [node for node in nodes if any(node_.startswith("C") for node_ in node)] 
            edges = list(comps[0][i].edges)
            neighbors = {str(n): [str(nb) for nb in comps[0][i].neighbors(n)] for n in nodes}
            payload[j]["frames"][i] = {"nodes": nodes, "peptide_nodes": peptide_nodes, "mhc_nodes": mhc_nodes, "mixed_nodes": mixed_nodes, "edges": edges, "neighbors": neighbors}

    out_json.parent.mkdir(parents=True, exist_ok=True)
    with open(out_json, "w") as f:
        json.dump(payload, f, indent=4)


# =============================================================================
# Avaliação e métricas
# =============================================================================
def get_protein_keys(original_graphs: dict):
    keys = list(original_graphs.keys())
    if all(isinstance(k, str) and k.isdigit() for k in keys):
        return [str(i) for i in sorted(map(int, keys))]
    return keys

def project_nodes_instances(frame_nodes, p):
    return [n[p] for n in frame_nodes]

def chain_signature(node_tuple):
    chains = []
    for lab in node_tuple:
        s = lab if isinstance(lab, str) else str(lab)
        chains.append(s.split(":")[0] if ":" in s else s)
    return "".join(chains)

def unique_chain_signatures(frame_nodes):
    return sorted({chain_signature(n) for n in frame_nodes})

def node_similarity_for_protein(frame, original_graphs, protein_keys, p):
    nodes_assoc = frame.get("nodes", [])
    if not nodes_assoc:
        return None

    prot_key = protein_keys[p]
    og = original_graphs[prot_key]
    prot_name = og.get("name", prot_key)

    Vp = set(og["nodes"])
    
    pep_orig_nodes = [node for node in Vp if node.startswith("C")]
    mhc_orig_nodes = [node for node in Vp if node.startswith("A")]
    inst = project_nodes_instances(nodes_assoc, p)
    Up = set(inst)

    total_orig = len(Vp) if Vp else 0
    node_coverage = (len(Up) / total_orig) if total_orig else 0.0

    total_inst = len(inst)
    unique_cnt = len(Up)
    duplication_ratio = (len(Up) / total_inst) if total_inst else 1.0
    duplication_rate  = 1.0 - duplication_ratio
    avg_multiplicity  = (total_inst / unique_cnt) if unique_cnt else float('inf')

    groups = defaultdict(set)
    for node_tuple in nodes_assoc:
        key = ''.join(str(part).split(':', 1)[0] for part in node_tuple)
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
        unique_nodes_per_chain=unique_nodes_per_chain_json
    )

def wmean(x, w):
    x = np.asarray(x, float); w = np.asarray(w, float)
    s = w.sum()
    return float(np.sum(x*w)/s) if s > 0 else np.nan

def wstd(x, w):
    x = np.asarray(x, float); w = np.asarray(w, float)
    m = wmean(x, w); s = w.sum()
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
    keep = (cw >= trim) & (cw <= 1.0 - trim)
    if not np.any(keep): keep = np.ones_like(cw, dtype=bool)
    return wmean(x[keep], w[keep])

def ivw_mean_proportions(cov, n):
    cov = np.asarray(cov, float); n = np.asarray(n, float)
    p = ((cov*n) + 0.5) / (n + 1.0)
    var = p*(1.0-p) / (n + 1.0) + 1e-12
    w = 1.0/var
    return wmean(p, w)

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
            unique_chain_signatures=chain_sigs,
            n_unique_chain_signatures=n_unique_chain_sigs
        )
    return df, summary

def evaluate_all_frames_nodes(json_path: Path):
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

def evaluate_all_frames_nodes_weighted(json_path: Path):
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


# =============================================================================
# Construção do graph por string de arquivos
# =============================================================================
def _save_eval_tables(out_dir: Path, df_fp_nodes: pd.DataFrame, df_frames_nodes_w: pd.DataFrame) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    if not df_fp_nodes.empty:
        df_fp_nodes.to_csv(out_dir / "nodes_per_protein.csv", index=False)
    if not df_frames_nodes_w.empty:
        cols = list(df_frames_nodes_w.columns)
        for lead in ["component_id", "frame_id"]:
            if lead in cols:
                cols = [lead] + [c for c in cols if c != lead]
        df_frames_nodes_w[cols].to_csv(out_dir / "nodes_summary_weighted.csv", index=False)

def _build_associated_graph_from_string(file_string: str, run_name: str, out_dir: Path, manifest: Dict) -> Path:
    """
    Recebe a string única separada por vírgulas com caminhos de arquivos.
    Injeta no manifest e roda create_graphs + AssociatedGraph.
    Retorna o caminho do JSON gerado.
    """
    # deep copy simples
    mfest = json.loads(json.dumps(manifest))
    constrains = mfest["inputs"][0]["constrains"]
    mfest["inputs"] = [{
        "path": file_string,
        "enable_tui": False,
        "extensions": [".pdb", ".cif"],
        "constrains": constrains
    }]

    S = mfest["settings"]
    S["run_name"] = run_name
    S["output_path"] = str(out_dir)

    graphs = create_graphs(mfest)

    checks = {"depth": S.get("check_depth"), "rsa": S.get("check_rsa")}
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
        "classes":                     S.get("classes", {}),
        "rsa_table":                   S.get("rsa_table", "Wilke")
    }

    G = AssociatedGraph(
        graphs=graphs,
        output_path=S["output_path"],
        run_name=S["run_name"],
        association_config=association_config,
    )

    # Se quiser figuras e afins, pode habilitar aqui
    G.draw_graph_interactive(show=False, save=True)
    # G.draw_graph(show=False, save=True)
    G.align_all_frames()
    G.create_pdb_per_protein()

    out_json = out_dir / f"graph_{run_name}.json"
    _make_json_from_associated_graph(G, out_json)
    return out_json


# =============================================================================
# Modos de execução: ALL, PAIRWISE, BOTH
# =============================================================================
def run_all(dir_path: Path, out_root: Path, manifest: Dict) -> Tuple[pd.DataFrame, pd.DataFrame]:
    files = find_struct_files(dir_path)
    file_string = files_to_comma_string(files)

    out_dir = out_root / "ALL"
    out_dir.mkdir(parents=True, exist_ok=True)

    run_name = "ALL_IN_DIR"
    json_path = _build_associated_graph_from_string(file_string, run_name, out_dir, manifest)

    try:
        df_fp_nodes, df_frames_nodes_w = evaluate_all_frames_nodes_weighted(json_path)
        if "frame_nodes_unique" in df_fp_nodes.columns and "total_nodes_associated" not in df_fp_nodes.columns:
            df_fp_nodes["total_nodes_associated"] = df_fp_nodes["frame_nodes_unique"]
    except Exception as e:
        print(f"[SKIP] Falha na avaliação ALL: {e}")
        return pd.DataFrame(), pd.DataFrame()

    _save_eval_tables(out_dir, df_fp_nodes, df_frames_nodes_w)

    # agregados
    if not df_fp_nodes.empty:
        df_fp_nodes.to_csv(out_root / "ALL_nodes_per_protein.csv", index=False)
    if not df_frames_nodes_w.empty:
        lead = ["component_id", "frame_id"]
        cols = lead + [c for c in df_frames_nodes_w.columns if c not in lead]
        df_frames_nodes_w = df_frames_nodes_w[cols]
        df_frames_nodes_w.to_csv(out_root / "ALL_nodes_summary_weighted.csv", index=False)

    print("All×All concluído.")
    return df_fp_nodes, df_frames_nodes_w


def run_pairwise(dir_path: Path, out_root: Path, manifest: Dict, score_column: str = "node_cov_wmean") -> Dict[str, pd.DataFrame]:
    files = find_struct_files(dir_path)
    labels = build_label_map(files)  # label -> caminho
    refs = sorted(labels.keys())

    if len(refs) < 2:
        raise RuntimeError("Pairwise requer pelo menos 2 arquivos na pasta.")

    group_dir = out_root / "PAIRS"
    group_dir.mkdir(parents=True, exist_ok=True)

    M_max = pd.DataFrame(0.0, index=refs, columns=refs)
    M_mean = pd.DataFrame(0.0, index=refs, columns=refs)

    for r1, r2 in combinations(refs, 2):
        file_string = ",".join([labels[r1], labels[r2]])
        run_name = f"{r1}__{r2}"
        out_dir = group_dir / run_name
        json_path = _build_associated_graph_from_string(file_string, run_name, out_dir, manifest)

        try:
            df_fp_nodes, df_frames_nodes_w = evaluate_all_frames_nodes_weighted(json_path)
            if "frame_nodes_unique" in df_fp_nodes.columns and "total_nodes_associated" not in df_fp_nodes.columns:
                df_fp_nodes["total_nodes_associated"] = df_fp_nodes["frame_nodes_unique"]
        except Exception as e:
            print(f"[SKIP] Falha na avaliação do par ({r1},{r2}): {e}")
            continue

        _save_eval_tables(out_dir, df_fp_nodes, df_frames_nodes_w)

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

    M_max.to_csv(group_dir / f"matrix_{score_column}_MAX.csv")
    M_mean.to_csv(group_dir / f"matrix_{score_column}_MEAN.csv")
    print(f"Matrizes salvas em {group_dir}")
    print("Pairwise concluído.")
    return {"max": M_max, "mean": M_mean}


# =============================================================================
# CLI
# =============================================================================
def parse_args():
    ap = argparse.ArgumentParser(
        description="Roda análise a partir de um manifest e um diretório com .pdb/.cif. Modo all, pairwise ou both."
    )
    ap.add_argument("--manifest", required=True, type=str, help="Caminho do manifest JSON")
    ap.add_argument("--dir", required=True, type=str, help="Pasta com arquivos .pdb/.cif (scan recursivo)")
    ap.add_argument("--mode", choices=["all", "pairwise", "both"], default="both")
    ap.add_argument("--out", dest="out_root", default=None, help="Pasta de saída. Se vazio, usa settings.output_path do manifest")
    ap.add_argument("--score-column", default="node_cov_wmean", help="Coluna de score para pairwise")
    return ap.parse_args()

def main():
    args = parse_args()

    manifest = load_manifest(args.manifest)
    S = manifest["settings"]
    run_name = S["run_name"]
    output_path = args.out_root if args.out_root else S["output_path"]

    logging.getLogger("matplotlib").setLevel(logging.ERROR)
    logging.basicConfig(
        stream=sys.stdout,
        level=logging.DEBUG if S.get("debug", False) else logging.INFO
    )
    log = logging.getLogger("CRSProtein")
    log.setLevel(logging.DEBUG if S.get("debug", False) else logging.INFO)

    dir_path = Path(args.dir)
    out_root = Path(output_path)
    out_root.mkdir(parents=True, exist_ok=True)

    if args.mode in ("all", "both"):
        run_all(dir_path, out_root, manifest)

    if args.mode in ("pairwise", "both"):
        run_pairwise(dir_path, out_root, manifest, score_column=args.score_column)

    print("Concluído.")


if __name__ == "__main__":
    main()
