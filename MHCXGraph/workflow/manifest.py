import json
import os
import shutil

from typing import Dict, Any

def load_manifest(manifest_path: str) -> Dict[str, Any]:
    if not manifest_path:
        return {}

    with open(manifest_path, "r") as f:
        data = json.load(f)

    data.setdefault("settings", {})
    data.setdefault("inputs", [])
    data.setdefault("selectors", {})

    S = data["settings"]
    S.setdefault("run_name", "test")
    S.setdefault("run_mode", "all")
    S.setdefault("output_path", "./outputs")

    os.makedirs(S["output_path"], exist_ok=True)
    
    shutil.copy2(manifest_path, S["output_path"]+"/manifest.json")
    S.setdefault("debug", False)
    S.setdefault("track_steps", False)
    S.setdefault("rsa_table", "Wilke")

    S.setdefault("edge_threshold", 8.5)
    S.setdefault("close_tolerance", 1.0)
    S.setdefault("node_granularity", "all_atoms")

    S.setdefault("include_ligands", True)
    S.setdefault("include_noncanonical_residues", True)
    S.setdefault("include_waters", True)

    S.setdefault("triad_rsa", False)
    S.setdefault("rsa_filter", 0.1)
    S.setdefault("asa_filter", 100)
    S.setdefault("close_tolerance_rsa", 0.1)
    S.setdefault("distance_std_threshold", 3.0)
    S.setdefault("distance_diff_threshold", 1.0)

    S.setdefault("rsa_bin_width", 20)
    S.setdefault("distance_bin_width", 2.0)

    S.setdefault("max_chunks", 5)

    S.setdefault("filter_triads_by_chain", None)
    S.setdefault("classes", {})
    
    S.setdefault("watch_residues", None)

    return data

def build_association_config(settings: Dict[str, Any], run_mode: str, tracker_residues) -> Dict[str, Any]:
    """
    Build the association_config dict passed to AssociatedGraph.

    This is just a thin adapter that takes the manifest 'settings'
    and produces the config structure expected by the core.
    """
    checks = {
        "rsa": settings.get("triad_rsa"),
    }

    return {
        "run_mode":                 run_mode,
        "edge_threshold":           settings.get("edge_threshold"),
        "distance_std_threshold":   settings.get("distance_std_threshold"),
        "distance_diff_threshold":  settings.get("distance_diff_threshold"),
        "rsa_filter":               settings.get("rsa_filter"),
        "rsa_bin_width":            settings.get("rsa_bin_width"),
        "distance_bin_width":       settings.get("distance_bin_width"),
        "close_tolerance":          settings.get("close_tolerance"),
        "close_tolerance_rsa":      settings.get("close_tolerance_rsa"),
        "checks":                   checks,
        "include_waters":           settings.get("include_waters"),
        "include_ligands":          settings.get("include_ligands"),
        "include_noncanonical_residues": settings.get("include_noncanonical_residues"),
        "classes":                  settings.get("classes", {}),
        "max_chunks":               settings.get("max_chunks"),
        "rsa_table":                settings.get("rsa_table", "Wilke"),
        "filter_triads_by_chain":   settings.get("filter_triads_by_chain", None),
        "watch_residues":           tracker_residues
    }

