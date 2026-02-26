# immunograph/core/residue_tracking.py

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Set, Union
import json
import os
import re
import time

def _normalize_watch_item(item):
    # Accept:
    #   - (protein_i, chain, resseq, icode) tuple (already ok)
    #   - "0:A:ARG:44" or "0:A:44" (protein index explicit)
    #   - "A:ARG:44" or "A:44" (assume protein_i=0)
    if isinstance(item, tuple) and len(item) == 4:
        p_i, chain, resseq, icode = item
        return (int(p_i), str(chain), int(resseq), str(icode or ""))
    if isinstance(item, str):
        s = item.strip()

        # optional leading "p_i:"
        p_i = 0
        if re.match(r"^\d+:", s):
            p_str, s = s.split(":", 1)
            p_i = int(p_str)

        # allow "CHAIN:RES:NUM" or "CHAIN:NUM"
        parts = s.split(":")
        if len(parts) == 3:
            chain, _res, num = parts
            chain2, _res2, resseq, icode = parse_node_label(f"{chain}:{_res}:{num}")
            return (p_i, chain2, resseq, icode)
        if len(parts) == 2:
            chain, num = parts
            # make a fake resname field for parse_node_label
            chain2, _res2, resseq, icode = parse_node_label(f"{chain}:XXX:{num}")
            return (p_i, chain2, resseq, icode)

    raise ValueError(f"Invalid watch residue format: {item!r}")


# -----------------------------
# Context: where in the pipeline
# -----------------------------
@dataclass(frozen=True)
class TrackCtx:
    run_id: str
    stage: str                     # "triads" | "combos" | "component" | "frames" | etc.
    protein_i: Optional[int] = None
    step_id: Optional[int] = None
    chunk_id: Optional[int] = None
    comp_id: Optional[int] = None
    frame_id: Optional[int] = None


# -----------------------------
# Residue normalization helpers
# -----------------------------
# Your nodes look like: "A:TYR:42" or sometimes "A:TYR:42A" (insertion)
_NODE_RE = re.compile(r"^(?P<chain>[^:]+):(?P<res>[^:]+):(?P<num>.+)$")
_NUM_RE  = re.compile(r"^(?P<resseq>-?\d+)(?P<icode>[A-Za-z]?)$")


def parse_node_label(node: str) -> Tuple[str, str, int, str]:
    """
    Parse "CHAIN:RES:NUM" -> (chain, resname, resseq_int, icode_str)
    Accepts "42", "42A", "-1", "-1B".
    """
    m = _NODE_RE.match(node)
    if not m:
        raise ValueError(f"Invalid node label: {node!r} (expected 'CHAIN:RES:NUM')")
    chain = m.group("chain")
    res   = m.group("res")
    num   = m.group("num")

    m2 = _NUM_RE.match(num)
    if not m2:
        # fallback: strip leading digits
        mm = re.match(r"(-?\d+)", num)
        if not mm:
            raise ValueError(f"Invalid residue number field: {num!r} in {node!r}")
        resseq = int(mm.group(1))
        icode  = ""
    else:
        resseq = int(m2.group("resseq"))
        icode  = m2.group("icode") or ""

    return chain, res, resseq, icode


def residue_key(protein_i: int, chain: str, resseq: int, icode: str = "") -> Tuple[int, str, int, str]:
    """
    Canonical watched key. We do NOT include resname, because in PDBs
    you can get weird naming and you usually care about position.
    """
    return (protein_i, chain, resseq, icode)


def triad_residues_from_absolute(triad_absolute: Tuple[Any, ...]) -> List[str]:
    """
    triad_absolute = (*triad_abs, *full_describer_absolute)
    triad_abs = [u, center, w] (strings "A:TYR:42")
    """
    if len(triad_absolute) < 3:
        return []
    u = triad_absolute[0]
    c = triad_absolute[1]
    w = triad_absolute[2]
    out = []
    if isinstance(u, str):
        out.append(u)
    if isinstance(c, str):
        out.append(c)
    if isinstance(w, str):
        out.append(w)
    return out


def combo_residues(combo: Tuple[Tuple[Any, ...], ...]) -> List[List[str]]:
    """
    combo = tuple(per_protein_triad_absolute)
    Each triad absolute begins with 3 residue labels.
    Returns list per protein: [[u,c,w],[u,c,w],...]
    """
    out: List[List[str]] = []
    for tri in combo:
        if not isinstance(tri, tuple) or len(tri) < 3:
            out.append([])
            continue
        labels = []
        for j in range(3):
            if isinstance(tri[j], str):
                labels.append(tri[j])
        out.append(labels)
    return out


# -----------------------------
# Main tracker
# -----------------------------
class ResidueTracker:
    """
    Watches specific residues across the pipeline.

    Design goals:
      - provenance rich (ctx includes step/chunk/comp/frame)
      - bounded memory (stores only when watched residue appears)
      - easy export to JSON

    Typical use:
      tracker = ResidueTracker(watch_residues={ (0,"A",42,""), (1,"B",15,"") }, out_dir="tracking")
      ...
      tracker.triads_built(ctx, token, triads_absolute)
      tracker.combos_built(ctx, token, combos)
      tracker.frame_accepted(ctx, edges_residues, edges_indices)
      ...
      tracker.dump_json()
    """

    def __init__(
        self,
        watch_residues: Optional[Iterable[Tuple[int, str, int, str]]] = None,
        *,
        out_dir: str = "residue_tracking",
        max_examples_per_event: int = 25,
        keep_event_log: bool = True,
    ) -> None:
        self.out_dir = out_dir
        os.makedirs(self.out_dir, exist_ok=True)

        raw = list(watch_residues or [])
        self.watch: Set[Tuple[int, str, int, str]] = set(_normalize_watch_item(x) for x in raw)
        self.max_examples_per_event = int(max_examples_per_event)
        self.keep_event_log = bool(keep_event_log)

        # events keyed by residue_key -> list of event dicts
        self.events_by_residue: Dict[Tuple[int, str, int, str], List[Dict[str, Any]]] = {
            rk: [] for rk in self.watch
        }

        # global summary counters
        self.summary: Dict[str, Any] = {
            "created_at": time.time(),
            "num_watched": len(self.watch),
            "counts": {
                "triads_built": 0,
                "combos_built": 0,
                "component_selected": 0,
                "component_skipped": 0,
                "frame_accepted": 0,
                "triad_filtered": 0,
            },
        }

    # -------------------------
    # Optional resolution helper
    # -------------------------
    def resolve_from_pdb_dfs(self, pdb_dfs: Sequence[Any], stage: str = "resolve") -> None:
        """
        You call this already in association_product.
        Here we just record that resolution happened.
        If you want mapping from (chain,resseq,icode)->resname etc,
        you can extend this later.
        """
        self._note_global({"event": "resolve_from_pdb_dfs", "stage": stage, "num_proteins": len(pdb_dfs)})

    # -------------------------
    # Event helpers
    # -------------------------
    def _note_global(self, payload: Dict[str, Any]) -> None:
        if not self.keep_event_log:
            return
        path = os.path.join(self.out_dir, "global_events.jsonl")
        with open(path, "a") as f:
            f.write(json.dumps(payload, ensure_ascii=False) + "\n")

    def _add_event(self, rk: Tuple[int, str, int, str], event: Dict[str, Any]) -> None:
        if rk not in self.events_by_residue:
            # only store watched residues
            return
        self.events_by_residue[rk].append(event)

    def _ctx_dict(self, ctx: TrackCtx) -> Dict[str, Any]:
        return asdict(ctx)

    def _hit_keys_from_labels(self, protein_i: int, labels: Iterable[str]) -> Set[Tuple[int, str, int, str]]:
        hits: Set[Tuple[int, str, int, str]] = set()
        for lab in labels:
            try:
                chain, _res, resseq, icode = parse_node_label(lab)
            except Exception:
                continue
            rk = residue_key(protein_i, chain, resseq, icode)
            if rk in self.watch:
                hits.add(rk)
        return hits

    # -------------------------
    # Public API events
    # -------------------------
    def triads_built(self, ctx: TrackCtx, token: Any, triads_absolute: Sequence[Tuple[Any, ...]]) -> None:
        self.summary["counts"]["triads_built"] += 1

        # find watched hits in this token list
        if ctx.protein_i is None:
            return

        protein_i = int(ctx.protein_i)

        # scan triads_absolute but keep only examples up to limit
        examples = []
        hits: Set[Tuple[int, str, int, str]] = set()

        for tri in triads_absolute:
            labels = triad_residues_from_absolute(tri)
            tri_hits = self._hit_keys_from_labels(protein_i, labels)
            if tri_hits:
                hits |= tri_hits
                if len(examples) < self.max_examples_per_event:
                    examples.append(tri)

        if not hits:
            return

        ev = {
            "event": "triads_built",
            "ctx": self._ctx_dict(ctx),
            "token": repr(token),
            "num_triads_token": len(triads_absolute),
            "examples": examples,
        }
        for rk in hits:
            self._add_event(rk, ev)

    def triad_filtered(self, ctx: TrackCtx, reason: str, *, triad_abs: Optional[Tuple[Any, ...]] = None, token: Any = None) -> None:
        self.summary["counts"]["triad_filtered"] += 1
        # record only if triad has watched residue
        if ctx.protein_i is None or triad_abs is None:
            return
        labels = triad_residues_from_absolute(triad_abs)
        hits = self._hit_keys_from_labels(int(ctx.protein_i), labels)
        if not hits:
            return
        ev = {
            "event": "triad_filtered",
            "ctx": self._ctx_dict(ctx),
            "reason": reason,
            "token": repr(token) if token is not None else None,
            "triad": triad_abs,
        }
        for rk in hits:
            self._add_event(rk, ev)

    def combos_built(self, ctx: TrackCtx, token: Any, combos: Sequence[Tuple[Tuple[Any, ...], ...]]) -> None:
        self.summary["counts"]["combos_built"] += 1

        # combos contain per-protein triads. We check each protein positionally.
        examples = []
        hits: Set[Tuple[int, str, int, str]] = set()

        for combo in combos:
            per_prot = combo_residues(combo)
            for p_i, labels in enumerate(per_prot):
                hits |= self._hit_keys_from_labels(p_i, labels)
            if hits and len(examples) < self.max_examples_per_event:
                examples.append(combo)

        if not hits:
            return

        ev = {
            "event": "combos_built",
            "ctx": self._ctx_dict(ctx),
            "token": repr(token),
            "num_combos_token": len(combos),
            "examples": examples,
        }
        for rk in hits:
            self._add_event(rk, ev)

    def component_selected(self, ctx: TrackCtx, *, component_nodes: Sequence[Any], component_edges: Optional[Sequence[Any]] = None) -> None:
        self.summary["counts"]["component_selected"] += 1
        # nodes are tuples, one residue label per protein, e.g. ("A:TYR:42","B:GLY:9",...)
        hits: Set[Tuple[int, str, int, str]] = set()
        for node in component_nodes:
            if not isinstance(node, tuple):
                continue
            for p_i, lab in enumerate(node):
                if not isinstance(lab, str):
                    continue
                hits |= self._hit_keys_from_labels(p_i, [lab])

        if not hits:
            return

        ev = {
            "event": "component_selected",
            "ctx": self._ctx_dict(ctx),
            "num_nodes": len(component_nodes),
            "num_edges": len(component_edges) if component_edges is not None else None,
        }
        for rk in hits:
            self._add_event(rk, ev)

    def component_skipped(self, ctx: TrackCtx, *, reason: str, component_size: int) -> None:
        self.summary["counts"]["component_skipped"] += 1
        # no residue detail here unless you pass nodes; keep as global only
        self._note_global({"event": "component_skipped", "ctx": self._ctx_dict(ctx), "reason": reason, "component_size": component_size})

    def frame_accepted(self, ctx: TrackCtx, *, edges_residues: Sequence[Any], edges_indices: Optional[Sequence[Any]] = None) -> None:
        self.summary["counts"]["frame_accepted"] += 1

        # edges_residues look like: ((("A:TYR:42","B:..."),(...)), ...)
        # In your convert_edges_to_residues you produce converted_edges as tuples of tuples of residue strings.
        hits: Set[Tuple[int, str, int, str]] = set()
        examples = []

        for e in edges_residues:
            # e is (nodeA, nodeB) where each node is a tuple of residue labels per protein
            if not isinstance(e, tuple) or len(e) != 2:
                continue
            a, b = e
            if isinstance(a, tuple):
                for p_i, lab in enumerate(a):
                    if isinstance(lab, str):
                        hits |= self._hit_keys_from_labels(p_i, [lab])
            if isinstance(b, tuple):
                for p_i, lab in enumerate(b):
                    if isinstance(lab, str):
                        hits |= self._hit_keys_from_labels(p_i, [lab])

            if hits and len(examples) < self.max_examples_per_event:
                examples.append(e)

        if not hits:
            return

        ev = {
            "event": "frame_accepted",
            "ctx": self._ctx_dict(ctx),
            "num_edges": len(edges_residues),
            "examples": examples,
        }
        for rk in hits:
            self._add_event(rk, ev)

    # -------------------------
    # Export
    # -------------------------
    def dump_json(self, filename: str = "residue_tracking_report.json") -> str:
        """
        Writes a single JSON with:
          - summary
          - events_by_residue (only watched residues)
        Returns output path.
        """
        out_path = os.path.join(self.out_dir, filename)
        payload = {
            "summary": self.summary,
            "watched": [list(rk) for rk in sorted(self.watch)],
            "events_by_residue": {
                str(rk): self.events_by_residue[rk] for rk in sorted(self.events_by_residue.keys())
            },
        }
        with open(out_path, "w") as f:
            json.dump(payload, f, indent=2, ensure_ascii=False)
        return out_path

