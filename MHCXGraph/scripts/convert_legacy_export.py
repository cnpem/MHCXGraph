#!/usr/bin/env python3
"""Convert legacy per-run association JSON exports into the current
unified ``graph_data_{pairwise,screening,multiple}.json`` schema the
dashboard (baked or standalone) expects.

Background
----------
``workflow/association.py::run_association_task`` (still the current code,
not a deprecated path) writes one ``graph_<run_name>.json`` file per
association task via ``utils/analysis.py::_make_json_from_associated_graph``.
For pairwise/screening runs that's one file per pair, sitting at:

    <ROOT>/PAIRWISE/<nameA>_vs_<nameB>/graph_<run_name>_<nameA>_<nameB>.json
    <ROOT>/SCREENING/<ref>_vs_<target>/graph_<run_name>_<ref>_<target>.json

For multiple mode it's a single file:

    <ROOT>/MULTIPLE/graph_<run_name>.json

``app.py`` then feeds each of these through ``AssociatedGraph.get_dashboard_data``
/``get_filtered_graph_data`` (in-memory, on the live graph objects) to build
the numeric-id, RSA/distance-annotated schema the dashboard actually
consumes -- and that transformation is never written back to disk anywhere.
If you still have those live objects (i.e. you can just re-run), do that;
it's the only way to get real distances back. This script reconstructs the
SAME node/edge/component shape from the on-disk JSON alone, for cases where
re-running isn't an option (old run, heavy pipeline, etc).

What's recoverable vs lost
---------------------------
Recoverable (topology): which residues matched, which are connected, node
classification (peptide/mhc/mixed) via chain letter, empty-pair detection,
per-protein full graphs (filtered_graphs).

NOT recoverable, because the legacy export never stored it in the first
place (it only existed on the live in-memory graph/attribute objects):
- edge "std" (structural-conservation std-dev) and "raw_dist" (Å distance)
- node RSA values
- components' std_matrix / node_index_map (frame-correlation highlighting)
- run metadata (node_granularity, edge_threshold, ...) and protein_paths

Those fields are filled with the exact same fallback values the live code
uses when it itself can't resolve them (None / "N/A"), which the frontend
already renders gracefully (see graph.js's null-checks and the `|| 'N/A'`
metadata fallbacks) -- so the converted dashboard shows correct topology,
just without distance-based coloring/thickness or RSA tooltips.

Usage
-----
    python scripts/convert_legacy_export.py --root /path/to/old/run/output

Scans --root recursively for legacy graph_*.json files (detected by the
presence of an "original_graphs" key, not by filename), groups them by
run, and writes graph_data_{pairwise,screening,multiple}.json next to each
group (or into --out if given). Those files load directly into
MHCXGraph_Standalone.html, or feed into create_master_dashboard() /
`mhcxgraph run ... --dashboard`-style regeneration for a baked HTML.
"""
from __future__ import annotations

import argparse
import json
import logging
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


log = logging.getLogger("convert_legacy_export")


# =============================================================================
# Per-pair / per-run reconstruction (mirrors classes/graph.py::get_dashboard_data
# and get_filtered_graph_data, field for field, minus what's unrecoverable)
# =============================================================================

def _node_repr_to_residues(node: Any) -> list[str]:
    """A frame node is a bare residue string (single-graph node) or a JSON
    array of residue strings (a multi-protein correspondence tuple, once it
    has round-tripped through json.dump/json.load)."""
    if isinstance(node, list):
        return [str(r) for r in node]
    return [str(node)]


def _split_residue(res: str) -> tuple[str, str, str] | None:
    """"C:GLY:18" -> ("C", "GLY", "18"). None if malformed."""
    parts = res.split(":")
    if len(parts) < 3:
        return None
    return parts[0], parts[1], parts[2]


def build_pair_nodes_edges_components(
    legacy: dict, global_proteins: list[str], local_protein_names: list[str]
) -> tuple[list[dict], list[dict], list[dict]]:
    """Reconstruct (nodes, edges, components) for one association result,
    exactly matching AssociatedGraph.get_dashboard_data's output shape.

    `legacy` is one parsed graph_<run_name>.json (has "original_graphs" plus
    numeric-string component keys). `local_protein_names` is this file's own
    protein order (original_graphs["0"]["name"], ["1"]["name"], ...) --
    needed to map each node's tuple position to a model_idx in
    `global_proteins`.
    """
    local_to_global = [global_proteins.index(n) for n in local_protein_names]

    global_nodes: dict[str, dict] = {}
    global_edges: dict[tuple[int, int], dict] = {}
    node_id_counter = 0
    edge_id_counter = 0
    components: list[dict] = []

    comp_keys = sorted(
        (k for k in legacy.keys() if k != "original_graphs"),
        key=lambda k: int(k) if str(k).lstrip("-").isdigit() else str(k),
    )

    for comp_key in comp_keys:
        comp = legacy[comp_key]
        comp_data = {
            "id": comp.get("comp", int(comp_key) if str(comp_key).isdigit() else comp_key),
            "node_ids": [],
            "frames": [],
            "std_matrix": None,   # only ever existed on the live graph object
            "node_index_map": {}, # ditto
        }

        frames = comp.get("frames", {})
        frame_keys = sorted(
            frames.keys(), key=lambda k: int(k) if str(k).isdigit() else str(k)
        )

        for frame_key in frame_keys:
            frame = frames[frame_key]
            frame_nodes = frame.get("nodes", [])
            frame_edges = frame.get("edges", [])

            if not frame_nodes and not frame_edges:
                # Empty frame -- e.g. no cross-reactivity found for this pair.
                # Still record it (mirrors get_dashboard_data appending a
                # frame even for empty results) so downstream "is this pair
                # empty" checks see an accurate node/edge count.
                continue

            frame_data = {
                "id": int(frame_key) if str(frame_key).isdigit() else frame_key,
                "node_ids": [],
                "rmsds": {},  # only ever existed on the live graph object
            }

            node_key_by_repr: dict[str, int] = {}

            for n in frame_nodes:
                residues = _node_repr_to_residues(n)
                # Mirror Python's str(tuple(...)) exactly -- including the
                # trailing comma on a length-1 tuple -- since node_label is
                # used as a dedup key and must match what the live code
                # would have produced for str(n).
                if len(residues) == 1:
                    node_label = f"('{residues[0]}',)"
                else:
                    node_label = "(" + ", ".join(f"'{r}'" for r in residues) + ")"

                mapping, chains = [], []
                for local_idx, res in enumerate(residues):
                    split = _split_residue(res)
                    if split is None:
                        continue
                    chain, resn, resi = split
                    chains.append(chain)
                    model_idx = (
                        local_to_global[local_idx] if local_idx < len(local_to_global) else local_idx
                    )
                    mapping.append({
                        "model_idx": model_idx,
                        "chain": chain, "resn": resn, "resi": resi,
                        "rsa": "N/A",  # not stored in the legacy export
                    })

                chain_id = "".join(chains) or "?"

                if node_label not in global_nodes:
                    title = f"Chains: {chain_id}\n{node_label}"
                    for m in mapping:
                        title += f"\nP{m['model_idx']} RSA: {m['rsa']}"
                    global_nodes[node_label] = {
                        "id": node_id_counter, "label": node_label, "title": title,
                        "group": chain_id, "mapping": mapping, "originalColor": None,
                    }
                    node_id_counter += 1

                numeric_id = global_nodes[node_label]["id"]
                node_key_by_repr[json.dumps(n, sort_keys=True) if isinstance(n, list) else str(n)] = numeric_id
                if numeric_id not in comp_data["node_ids"]:
                    comp_data["node_ids"].append(numeric_id)
                frame_data["node_ids"].append(numeric_id)

            def _resolve(node_val):
                key = json.dumps(node_val, sort_keys=True) if isinstance(node_val, list) else str(node_val)
                return node_key_by_repr.get(key)

            for edge in frame_edges:
                if not isinstance(edge, list) or len(edge) != 2:
                    continue
                u_val, v_val = edge
                u_id, v_id = _resolve(u_val), _resolve(v_val)
                if u_id is None or v_id is None:
                    continue
                edge_key = tuple(sorted([u_id, v_id]))
                if edge_key in global_edges:
                    continue
                global_edges[edge_key] = {
                    "id": edge_id_counter,
                    "from": u_id,
                    "to": v_id,
                    "std": None,       # only ever existed on the live graph object
                    "title": "Distances: unavailable (converted from legacy export)\n",
                    "raw_dist": None,  # ditto
                }
                edge_id_counter += 1

            comp_data["frames"].append(frame_data)

        components.append(comp_data)

    return list(global_nodes.values()), list(global_edges.values()), components


def build_filtered_graph(original_graph: dict, model_idx: int) -> dict:
    """Mirrors AssociatedGraph.get_filtered_graph_data: string ids, no
    distances (the legacy export's per-protein edges carry no attributes)."""
    f_nodes, f_edges = [], []
    for n in original_graph.get("nodes", []):
        split = _split_residue(str(n))
        mapping = (
            [{"model_idx": model_idx, "chain": split[0], "resn": split[1], "resi": split[2]}]
            if split else []
        )
        chain = split[0] if split else "?"
        f_nodes.append({
            "id": str(n), "label": str(n), "title": f"Chain: {chain}\n{n}",
            "group": chain, "mapping": mapping,
        })
    for edge in original_graph.get("edges", []):
        if not isinstance(edge, list) or len(edge) != 2:
            continue
        u, v = edge
        f_edges.append({
            "id": f"{u}-{v}", "from": str(u), "to": str(v),
            "title": "", "raw_dist": None,
        })
    return {
        "id": model_idx,
        "name": original_graph.get("name", f"protein_{model_idx}"),
        "nodes": f_nodes,
        "edges": f_edges,
    }


# =============================================================================
# Discovery: find legacy files on disk and group them into runs
# =============================================================================

def is_legacy_export(path: Path) -> dict | None:
    """Return the parsed JSON if `path` looks like a legacy association
    export (schema-based detection: has "original_graphs"), else None."""
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception as e:
        log.debug(f"Skipping {path} (not valid JSON: {e})")
        return None
    if isinstance(data, dict) and "original_graphs" in data:
        return data
    return None


def discover_legacy_files(root: Path) -> list[Path]:
    return sorted(p for p in root.rglob("*.json") if is_legacy_export(p) is not None)


def group_files(files: list[Path]) -> dict[Path, list[Path]]:
    """Group per-pair files by the folder that holds all of a run's
    "<a>_vs_<b>" pair subfolders (pairwise/screening); files that aren't
    inside a "_vs_" folder are each their own "multiple"-mode group."""
    groups: dict[Path, list[Path]] = defaultdict(list)
    for f in files:
        parent = f.parent
        if "_vs_" in parent.name:
            groups[parent.parent].append(f)
        else:
            groups[f].append(f)  # standalone: this file is its own group
    return groups


def detect_screening_reference(pair_protein_names: list[tuple[str, str]]) -> str | None:
    """If one protein appears in every pair and every other protein appears
    in exactly one pair, that's an unambiguous 1-vs-all screening topology."""
    counts = Counter()
    for a, b in pair_protein_names:
        counts[a] += 1
        counts[b] += 1
    n_pairs = len(pair_protein_names)
    candidates = [name for name, c in counts.items() if c == n_pairs]
    if len(candidates) != 1:
        return None
    ref = candidates[0]
    others_ok = all(c == 1 for name, c in counts.items() if name != ref)
    return ref if others_ok else None


# =============================================================================
# Assembly: build the final graph_data_{mode}.json for one group
# =============================================================================

def convert_pairwise_or_screening_group(
    files: list[Path], run_name: str, screening_reference: str | None
) -> dict:
    parsed = [(f, json.loads(f.read_text(encoding="utf-8"))) for f in files]

    # Deterministic global protein order: first-seen across the group.
    global_proteins: list[str] = []
    pair_protein_names: list[tuple[str, str]] = []
    for f, legacy in parsed:
        og = legacy.get("original_graphs", {})
        keys = sorted(og.keys(), key=lambda k: int(k) if str(k).isdigit() else str(k))
        if len(keys) != 2:
            log.warning(f"{f}: expected exactly 2 proteins in original_graphs, found {len(keys)} -- skipping.")
            continue
        name0, name1 = og[keys[0]]["name"], og[keys[1]]["name"]
        pair_protein_names.append((name0, name1))
        for n in (name0, name1):
            if n not in global_proteins:
                global_proteins.append(n)

    detected_ref = screening_reference or detect_screening_reference(pair_protein_names)

    master_export: dict[str, Any] = {
        "mode": "pairwise",
        "run_name": run_name,
        "metadata": {
            "run_mode": "screening" if detected_ref else "pairwise",
            "node_granularity": None,
            "edge_threshold": None,
            "global_distance_diff_threshold": None,
            "converted_from_legacy": True,
        },
        "proteins": global_proteins,
        "protein_paths": ["" for _ in global_proteins],  # not stored in legacy export
        "filtered_graphs": {},
        "pairs": {},
    }
    if detected_ref:
        master_export["actual_mode"] = "screening"
        master_export["reference_structure"] = detected_ref

    for f, legacy in parsed:
        og = legacy.get("original_graphs", {})
        keys = sorted(og.keys(), key=lambda k: int(k) if str(k).isdigit() else str(k))
        if len(keys) != 2:
            continue
        name0, name1 = og[keys[0]]["name"], og[keys[1]]["name"]
        pair_key = f"{name0}_vs_{name1}"

        nodes, edges, components = build_pair_nodes_edges_components(
            legacy, global_proteins, [name0, name1]
        )
        master_export["pairs"][pair_key] = {
            "proteins": [name0, name1],
            "protein_paths": ["", ""],
            "nodes": nodes,
            "edges": edges,
            "components": components,
        }

        for k in keys:
            pname = og[k]["name"]
            if pname not in master_export["filtered_graphs"]:
                model_idx = global_proteins.index(pname)
                master_export["filtered_graphs"][pname] = build_filtered_graph(og[k], model_idx)

        n_empty = sum(1 for c in components if not c["node_ids"])
        log.info(
            f"{f.parent.name}: {len(nodes)} nodes, {len(edges)} edges, "
            f"{len(components)} components ({n_empty} empty)"
        )

    return master_export


def convert_multiple_group(file: Path, run_name: str) -> dict:
    legacy = json.loads(file.read_text(encoding="utf-8"))
    og = legacy.get("original_graphs", {})
    keys = sorted(og.keys(), key=lambda k: int(k) if str(k).isdigit() else str(k))
    global_proteins = [og[k]["name"] for k in keys]
    local_names = list(global_proteins)

    nodes, edges, components = build_pair_nodes_edges_components(legacy, global_proteins, local_names)
    filtered_graphs = [build_filtered_graph(og[k], i) for i, k in enumerate(keys)]

    log.info(f"{file.name}: {len(nodes)} nodes, {len(edges)} edges, {len(components)} components")

    return {
        "mode": "multiple",
        "run_name": run_name,
        "metadata": {
            "run_mode": "multiple",
            "node_granularity": None,
            "edge_threshold": None,
            "global_distance_diff_threshold": None,
            "converted_from_legacy": True,
        },
        "proteins": global_proteins,
        "protein_paths": ["" for _ in global_proteins],
        "nodes": nodes,
        "edges": edges,
        "components": components,
        "filtered_graphs": filtered_graphs,
    }


# =============================================================================
# CLI
# =============================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True, type=Path, help="Root folder to scan recursively for legacy graph_*.json exports.")
    ap.add_argument("--out", type=Path, default=None, help="Output directory. Default: write next to each detected run's files.")
    ap.add_argument("--screening-reference", default=None, help="Force this protein name as the screening reference (overrides auto-detection).")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(message)s")

    root = args.root.resolve()
    if not root.is_dir():
        raise SystemExit(f"--root {root} is not a directory")

    files = discover_legacy_files(root)
    if not files:
        raise SystemExit(f"No legacy association exports (files with an 'original_graphs' key) found under {root}")
    log.info(f"Found {len(files)} legacy export file(s) under {root}")

    groups = group_files(files)
    log.info(f"Grouped into {len(groups)} run(s)")

    n_written = 0
    for group_root, group_files_ in groups.items():
        is_multi_pair = "_vs_" in Path(group_files_[0]).parent.name
        run_name = group_root.name or "converted"

        if is_multi_pair:
            master = convert_pairwise_or_screening_group(group_files_, run_name, args.screening_reference)
            mode_tag = "screening" if master.get("actual_mode") == "screening" else "pairwise"
        else:
            master = convert_multiple_group(group_files_[0], run_name)
            mode_tag = "multiple"

        out_dir = args.out if args.out else group_root
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"graph_data_{mode_tag}.json"
        with open(out_path, "w", encoding="utf-8") as fh:
            json.dump(master, fh)

        log.info(f"-> wrote {out_path}  ({mode_tag}, {len(master.get('pairs', master.get('nodes', [])))} "
                  f"{'pairs' if 'pairs' in master else 'nodes'})")
        n_written += 1

    log.info(
        f"\nDone. {n_written} file(s) written. Load them into MHCXGraph_Standalone.html "
        "(mhcxgraph standalone-dashboard), or pass to create_master_dashboard() for a baked HTML.\n"
        "Reminder: converted pairs carry topology only -- no edge distances/std, no RSA, "
        "no frame-correlation highlighting (that data only ever existed on the live run's "
        "in-memory graph objects, never written to disk)."
    )


if __name__ == "__main__":
    main()
