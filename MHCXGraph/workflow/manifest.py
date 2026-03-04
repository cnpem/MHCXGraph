import json
import os
import shutil
from typing import Any


def load_manifest(manifest_path: str) -> dict[str, Any]:
    if not manifest_path:
        return {}

    with open(manifest_path) as f:
        data = json.load(f)

    data.setdefault("settings", {})
    data.setdefault("inputs", [])
    data.setdefault("selectors", {})

    settings = data["settings"]

    settings.setdefault("run_name", "test")
    settings.setdefault("run_mode", "all")
    settings.setdefault("output_path", "./outputs")

    os.makedirs(settings["output_path"], exist_ok=True)

    shutil.copy2(manifest_path, settings["output_path"] + "/manifest.json")
    settings.setdefault("debug_logs", False)
    settings.setdefault("debug_tracking", False)
    settings.setdefault("verbose", False)


    settings.setdefault("edge_threshold", 8.5)
    settings.setdefault("close_tolerance", 1.0)
    settings.setdefault("node_granularity", "all_atoms")

    settings.setdefault("include_ligands", True)
    settings.setdefault("include_noncanonical_residues", True)
    settings.setdefault("include_waters", True)

    settings.setdefault("triad_rsa", False)
    settings.setdefault("rsa_filter", 0.1)
    settings.setdefault("asa_filter", 5)
    settings.setdefault("close_tolerance_rsa", 0.1)
    settings.setdefault("local_distance_threshold", 1.0)
    settings.setdefault("global_distance_threshold", 2.0)

    settings.setdefault("rsa_bin_width", 0.2)
    settings.setdefault("distance_bin_width", 2.0)

    settings.setdefault("max_chunks", 5)

    settings.setdefault("filter_triads_by_chain", None)
    settings.setdefault("classes", {})

    settings.setdefault("watch_residues", None)

    return data


def build_association_config(settings: dict[str, Any], run_mode: str, tracker_residues) -> dict[str, Any]:
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
        "local_distance_diff_threshold":  settings.get("local_distance_diff_threshold"),
        "global_distance_diff_threshold":  settings.get("global_distance_diff_threshold"),
        "rsa_filter":               settings.get("rsa_filter"),
        "rsa_bin_width":            settings.get("rsa_bin_width"),
        "distance_bin_width":       settings.get("distance_bin_width"),
        "close_tolerance":          settings.get("close_tolerance"),
        "close_tolerance_rsa":      settings.get("close_tolerance_rsa"),
        "checks":                   checks,
        "include_waters":           settings.get("include_waters"),
        "include_ligands":          settings.get("include_ligands"),
        "include_noncanonical_residues": settings.get("include_noncanonical_residues"),
        "classes":                  settings.get("classes"),
        "max_chunks":               settings.get("max_chunks"),
        "filter_triads_by_chain":   settings.get("filter_triads_by_chain"),
        "watch_residues":           tracker_residues,
        "debug_logs":               settings.get("debug_logs"),
        "verbose":                  settings.get("verbose")
    }

