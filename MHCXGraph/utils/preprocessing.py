import logging
import os
import re
import time
from pathlib import Path
from typing import Any

import gemmi
import networkx as nx

from MHCXGraph.classes.graph import Graph
from MHCXGraph.core.config import make_default_config
from MHCXGraph.core.tracking import save
from MHCXGraph.io_utils.pdb_io import get_user_selection, list_pdb_files
from MHCXGraph.utils.logging_utils import get_log

log = get_log()

class LogicError(Exception):
    """Logic error while evaluating boolean set expression."""
    pass


def _eval_logic_expression(expr: str,
                           sets: dict[str, set],
                           universe: set) -> set:
    """Evaluate a boolean expression over named sets using &, |, ! and parentheses."""
    if not expr or not expr.strip():
        raise LogicError("Empty logic expression")

    tokens = re.findall(r"[A-Za-z_][A-Za-z0-9_:]*|[()!&|]", expr.replace(" ", ""))
    prec = {"!": 3, "&": 2, "|": 1}
    right_assoc = {"!"}
    output = []
    stack = []

    for t in tokens:
        if re.match(r"[A-Za-z_][A-Za-z0-9_]*", t):
            if t not in sets:
                raise LogicError(f"Unknown set name in logic: {t!r}")
            output.append(t)
        elif t in ("!", "&", "|"):
            while stack:
                top = stack[-1]
                if top in prec and (
                    prec[top] > prec[t] or
                    (prec[top] == prec[t] and t not in right_assoc)
                ):
                    output.append(stack.pop())
                else:
                    break
            stack.append(t)
        elif t == "(":
            stack.append(t)
        elif t == ")":
            while stack and stack[-1] != "(":
                output.append(stack.pop())
            if not stack:
                raise LogicError("Mismatched parentheses")
            stack.pop()
        else:
            raise LogicError(f"Invalid token in logic: {t!r}")

    while stack:
        top = stack.pop()
        if top in ("(", ")"):
            raise LogicError("Mismatched parentheses at end")
        output.append(top)

    st: list[set] = []
    for t in output:
        if t in ("!", "&", "|"):
            if t == "!":
                if not st:
                    raise LogicError("Missing operand for '!'")
                A = st.pop()
                st.append(universe - A)
            else:
                if len(st) < 2:
                    raise LogicError(f"Missing operands for {t}")
                B = st.pop()
                A = st.pop()
                if t == "&":
                    st.append(A & B)
                else:
                    st.append(A | B)
        else:
            st.append(set(sets[t]))

    if len(st) != 1:
        raise LogicError("Invalid logic expression evaluation")

    return st[0]

def remove_water_from_pdb(source_file, dest_file):
    """Remove water molecules from a PDB or mmCIF file and save the cleaned version safely."""

    if os.path.exists(dest_file) or "_nOH" in source_file:
        logger.debug(f"The file {dest_file} already exists.")
        return

    suffix = source_file.lower()
    is_cif = suffix.endswith((".cif", ".mmcif", ".mcif"))

    st = gemmi.read_structure(source_file)

    for model in st:
        model.remove_waters()

    if is_cif:
        doc = st.make_mmcif_document()
        with open(dest_file, "w") as f:
            f.write(doc.as_string())
    else:
        st.write_pdb(dest_file)

    logger.debug(f"Saved cleaned structure without waters: {dest_file}")

def get_exposed_residues(graph: Graph, rsa_filter = 0.1, asa_filter = 100.0, selection_params=None) -> nx.Graph:
    selection_params = selection_params or {}
    logic_expr = selection_params.get("logic")

    graph.create_subgraph(name="exposed_residues", rsa_threshold=rsa_filter, asa_threshold=asa_filter)
    exposed = graph.get_subgraph("exposed_residues")
    if exposed is None or exposed.number_of_nodes() == 0:
        raise Exception("No exposed residues found")

    sets: dict[str, set] = {}
    sets["exposed"] = set(exposed.nodes())

    universe = set(graph.graph.nodes())
    G = graph.graph

    if "chains" in selection_params:
        chains_cfg = selection_params["chains"]
        if not isinstance(chains_cfg, list):
            raise TypeError(f"`chains` must be list, got {type(chains_cfg)}")
        graph.create_subgraph(name="selected_chains", chains=chains_cfg)
        chains_sub = graph.get_subgraph("selected_chains")
        base_nodes = set(chains_sub.nodes()) if chains_sub is not None else set()
        sets["chains"] = set(chains_sub.nodes()) if chains_sub is not None else set()

        by_chain: dict[str, set] = {}
        for n in base_nodes:
            d = G.nodes[n]
            cid = d.get("chain_id") or d.get("chain")
            if cid is None:
                continue
            by_chain.setdefault(cid, set()).add(n)
        for cid, nodes in by_chain.items():
            sets[f"chains:{cid}"] = nodes

    if "residues" in selection_params:
        residues_cfg = selection_params["residues"]
        if not isinstance(residues_cfg, dict):
            raise TypeError(f"`residues` must be dict, got {type(residues_cfg)}")

        selected_nodes: list[str] = []
        per_chain: dict[str, set] = {}

        for n, d in G.nodes(data=True):
            cid = d.get("chain_id") or d.get("chain")
            rnum = d.get("residue_number") or d.get("resseq")
            if cid is None or rnum is None:
                continue
            if cid in residues_cfg and rnum in residues_cfg[cid]:
                selected_nodes.append(n)
                per_chain.setdefault(cid, set()).add(n)

        if selected_nodes:
            graph.create_subgraph(
                name="selected_residues_nodes",
                node_list=selected_nodes,
                return_node_list=False,
            )
        sets["residues"] = set(selected_nodes)

        for cid, nodes in per_chain.items():
            sets[f"residues:{cid}"] = nodes

    if "structures" in selection_params:
        structures_cfg = selection_params["structures"]

        if isinstance(structures_cfg, list):
            graph.create_subgraph(name="selected_structures", ss_elements=structures_cfg)
            s_sub = graph.get_subgraph("selected_structures")
            base_nodes = set(s_sub.nodes()) if s_sub is not None else set()
            sets["structures"] = base_nodes

            by_chain: dict[str, set] = {}
            for n in base_nodes:
                d = G.nodes[n]
                cid = d.get("chain_id") or d.get("chain")
                if cid is None:
                    continue
                by_chain.setdefault(cid, set()).add(n)
            for cid, nodes in by_chain.items():
                sets[f"structures:{cid}"] = nodes

        elif isinstance(structures_cfg, dict):
            allowed_by_chain = {
                ch: set(vals)
                for ch, vals in structures_cfg.items()
                if ch != "*"
            }
            default_structs = set(structures_cfg.get("*", []))

            selected_nodes: list[str] = []
            per_chain: dict[str, set] = {}

            for n, d in G.nodes(data=True):
                label = str(n)
                parts = label.split(":")
                if len(parts) < 3:
                    continue
                chain_id = parts[0]

                allowed_ss = allowed_by_chain.get(chain_id, default_structs)
                if not allowed_ss:
                    continue

                ss = d.get("ss", None)
                if ss in allowed_ss:
                    selected_nodes.append(n)
                    per_chain.setdefault(chain_id, set()).add(n)

            if selected_nodes:
                graph.create_subgraph(
                    name="selected_structures_nodes",
                    node_list=selected_nodes,
                    return_node_list=False,
                )
            sets["structures"] = set(selected_nodes)
            for cid, nodes in per_chain.items():
                sets[f"structures:{cid}"] = nodes
        elif structures_cfg is None:
            pass
        else:
            raise TypeError(
                f"`structures` must be list or dict, got {type(structures_cfg)}"
            )

    if not logic_expr:
        union_sets: list[set] = []
        for key in ("residues", "chains", "structures"):
            s = sets.get(key)
            if s:
                union_sets.append(s)

        if union_sets:
            combined = set().union(*union_sets)
            selected = sets["exposed"] & combined
        else:
            selected = sets["exposed"]
    else:
        selected = _eval_logic_expression(logic_expr, sets, universe)
        if "exposed" in sets and "exposed" not in logic_expr:
            selected &= sets["exposed"]

    if not selected:
        raise Exception("I did not find any nodes that pass your filter/logic")

    graph.create_subgraph(name="s_graph",
                          node_list=list(selected),
                          return_node_list=False)

    s_graph = graph.get_subgraph("s_graph")
    if s_graph is not None:
        return s_graph

    raise Exception("Unexpected error: s_graph is None")

def _is_within(child: Path, parent: Path) -> bool:
    """True se child está dentro de parent (ou igual)."""
    try:
        child.resolve().relative_to(parent.resolve())
        return True
    except Exception:
        return False

def _name_contains(fname: str, cond: Any) -> bool:
    """Suporta string ou lista de strings para 'file_name_contains'."""
    if cond is None:
        return True
    if isinstance(cond, str):
        return cond in fname
    try:
        return any(s in fname for s in cond)
    except Exception:
        return False

def _merge_constraints(base: dict[str, Any], add: dict[str, Any]) -> dict[str, Any]:
    """
    Merge chain, residue, and structure constraints from two selector blocks.

    Parameters
    ----------
    base : dict
        First constraint dictionary.
    add : dict
        Second constraint dictionary.

    Returns
    -------
    dict
        Unified constraint dictionary containing merged chains, residue
        positions per chain, and structure specifications. The ``structures``
        field is preserved either as a list or as a dict; mixing both
        representations across inputs is not allowed and raises TypeError.
    """
    out: dict[str, Any] = {
        "chains": [],
        "residues": {},
        "structures": None,
    }

    for src in (base or {}, add or {}):
        for c in src.get("chains") or []:
            if c not in out["chains"]:
                out["chains"].append(c)

        for ch, positions in (src.get("residues") or {}).items():
            lst = out["residues"].setdefault(ch, [])
            for p in positions:
                if p not in lst:
                    lst.append(p)

        structures = src.get("structures", None)
        if structures is None:
            continue

        if isinstance(structures, list):
            if out["structures"] is None:
                out["structures"] = []
            elif isinstance(out["structures"], dict):
                raise TypeError(
                    "Cannot merge list 'structures' with dict 'structures' in constraints."
                )
            for s in structures:
                if s not in out["structures"]:
                    out["structures"].append(s)

        elif isinstance(structures, dict):
            if out["structures"] is None:
                out["structures"] = {}
            elif isinstance(out["structures"], list):
                raise TypeError(
                    "Cannot merge dict 'structures' with list 'structures' in constraints."
                )
            for ch, vals in structures.items():
                lst = out["structures"].setdefault(ch, [])
                for v in vals:
                    if v not in lst:
                        lst.append(v)
        else:
            raise TypeError(
                f"'structures' must be list or dict, got {type(structures)}"
            )

    return out

def _infer_is_dir_or_file(path_str: str) -> str:
    p = Path(path_str).expanduser()
    if p.exists():
        return "dir" if p.is_dir() else "file"
    typical_exts = [".pdb", ".pdb.gz", ".cif", ".ent", ".mmcif"]
    return "file" if any(path_str.endswith(ext) for ext in typical_exts) else "dir"

def list_struct_files(folder: Path, extensions: list[str]) -> list[Path]:
    exts = set(extensions or [".pdb", ".pdb.gz", ".cif"])
    files: list[Path] = []
    for p in folder.rglob("*"):
        if p.is_file():
            for ext in exts:
                if str(p.name).endswith(ext):
                    files.append(p)
                    break
    files.sort()
    return files

def collect_selected_files_from_manifest(manifest):
    results = []

    for rule in manifest.get("inputs", []):
        path_val = rule.get("path")
        if not path_val:
            continue

        paths = path_val.split(",") if isinstance(path_val, str) else list(path_val)
        extensions = rule.get("extensions") or [".pdb", ".pdb.gz", ".cif"]
        enable_tui = bool(rule.get("enable_tui", False))

        for p_str in paths:
            kind = _infer_is_dir_or_file(p_str)
            p = Path(p_str).expanduser().resolve()

            if kind == "file":
                results.append({"input_path": str(p), "name": p.name})
            else:
                if enable_tui:
                    names = list_pdb_files(str(p), extensions=extensions)
                    selected = get_user_selection(names, str(p))
                    for full_path, fname in selected:
                        results.append({"input_path": str(Path(full_path).resolve()),
                                        "name": fname})
                else:
                    if not p.exists() or not p.is_dir():
                        continue
                    for fname in sorted(os.listdir(p)):
                        if any(fname.endswith(ext) for ext in extensions):
                            full = (p / fname).resolve()
                            results.append({"input_path": str(full), "name": fname})

    seen = set()
    out = []
    for it in results:
        if it["input_path"] not in seen:
            seen.add(it["input_path"])
            out.append(it)
    return out

def resolve_selection_params_for_file(file_path: Path, manifest: dict[str, Any]) -> dict[str, Any]:
    if not manifest:
        return {}
    fname = file_path.name
    fpath_str = str(file_path.resolve())
    merged: dict[str, Any] = {}
    merged_logic: str | None = None

    for rule in manifest.get("inputs", []):
        paths = rule.get("path")
        if not paths:
            continue
        if isinstance(paths, str):
            paths = paths.split(",")

        hit_path = False
        for p in paths:
            p_obj = Path(p).expanduser()
            if any(ch in p for ch in "*?[]"):
                if fnmatch.fnmatch(fpath_str, str(p)):
                    hit_path = True
                    break
            elif p_obj.exists() and p_obj.is_dir():
                if _is_within(file_path, p_obj):
                    hit_path = True
                    break
            else:
                try:
                    if file_path.resolve() == p_obj.resolve():
                        hit_path = True
                        break
                except Exception:
                    if str(p) in fpath_str:
                        hit_path = True
                        break
        if not hit_path:
            continue

        for c in (rule.get("selectors") or []):
            name = c.get("name")
            if not name:
                continue
            if not _name_contains(fname, c.get("file_name_contains")):
                continue
            spec = manifest.get("selectors", {}).get(name, {})
            merged = _merge_constraints(merged, spec)

            logic_expr = spec.get("logic")
            if logic_expr:
                if merged_logic and merged_logic != logic_expr:
                    logger.warning(
                        f"Conflicting logic for {file_path.name}: "
                        f"{merged_logic!r} vs {logic_expr!r}. Using last one."
                    )
                merged_logic = logic_expr

    if merged_logic:
        merged["logic"] = merged_logic

    return merged

def create_graphs(manifest: dict) -> list[tuple]:
    settings = manifest["settings"]

    output_path = Path(settings["output_path"]).expanduser().resolve()
    log.vinfo(f"Trying to create output directory in {output_path}", "Creating output diretory")
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Retrieve the list of files passed via manifest.
    selected_files = collect_selected_files_from_manifest(manifest)
    if not selected_files:
        msg = "None file was selected from manifest."
        log.warning(msg)
        raise Exception(msg)

    graph_config = make_default_config(
        edge_threshold=settings["edge_threshold"],
        granularity=settings["node_granularity"],
        include_waters=settings["include_waters"],
        dssp_acc_array=settings["rsa_table"],
        include_ligands=settings["include_ligands"],
        include_noncanonical_residues=settings["include_noncanonical_residues"]
    )

    graphs: list[tuple] = []
    start = time.perf_counter()
    for file_info in selected_files:
        orig_path = Path(file_info["input_path"]).resolve()
        dest_path = Path(settings["output_path"]).resolve()
        if not settings["include_waters"]:
            cleaned_name = file_info["name"]
            if cleaned_name.endswith(".pdb.gz"):
                cleaned_name = cleaned_name[:-7] + "_nOH.pdb" if "_nOH.pdb" not in cleaned_name else cleaned_name
            elif cleaned_name.endswith(".pdb"):
                cleaned_name = cleaned_name[:-4] + "_nOH.pdb" if "_nOH.pdb" not in cleaned_name else cleaned_name
            elif cleaned_name.endswith(".cif"):
                cleaned_name = cleaned_name[:-4] + "_nOH.cif" if "_nOH.cif" not in cleaned_name else cleaned_name
            else:
                cleaned_name = cleaned_name + "_nOH.pdb" if "_nOH.pdb" not in cleaned_name else cleaned_name
            cleaned_path = (dest_path.parent / cleaned_name).resolve()
            remove_water_from_pdb(str(orig_path), str(cleaned_path))
            graph_path = cleaned_path
        else:
            graph_path = orig_path

        graph_instance = Graph(config=graph_config, graph_path=str(graph_path))
        selection_params = resolve_selection_params_for_file(orig_path, manifest)

        subgraph = get_exposed_residues(
            graph=graph_instance,
            rsa_filter=settings.get("rsa_filter"),
            asa_filter=settings.get("asa_filter"),
            selection_params=selection_params or {},
        )

        sub_dir = output_path / "filtered_graphs"
        base_name = Path(orig_path).stem
        graph_instance.save_subgraph_view(
            g=subgraph,
            output_dir=sub_dir,
            name=f"{base_name}_filtered",
            with_html=True,
        )

        graph_instance.save_filtered_pdb(
            g=subgraph,
            output_path=sub_dir,
            name=f"{base_name}_filtered",
            use_cif=True
        )


        save("create_graphs", f"{graph_path.stem}_subgraph", subgraph)
        graphs.append((subgraph, str(orig_path)))

    end = time.perf_counter()

    log.vinfo(f"Took {end - start:.6f} seconds to create graphs")
    return graphs
