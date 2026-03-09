import webbrowser
from itertools import combinations
from pathlib import Path

from MHCXGraph.cli.cli_parser import parse_args
from MHCXGraph.core.residue_tracking import ResidueTracker
from MHCXGraph.core.tracking import init_tracker
from MHCXGraph.utils.logging_utils import setup_logging
from MHCXGraph.utils.preprocessing import create_graphs
from MHCXGraph.workflow.association import run_association_task
from MHCXGraph.workflow.manifest import build_association_config, load_manifest


def setup_trackers(output_dir, settings):
    tracker_residues = (
        ResidueTracker(settings.get("watch_residues"))
    ) if settings.get("watch_residues") else None

    init_tracker(
        root="CrossSteps",
        outdir=output_dir,
        enabled=settings.get("debug_tracking"),
        prefer_npy_for_ndarray=True,
        add_timestamp_prefix=False,
    )

    return tracker_residues


def run_all_mode(graphs, base_output, run_name, config, log):
    target_dir = base_output / "ALL"

    run_association_task(
        graphs=graphs,
        output_path=target_dir,
        run_name=run_name,
        association_config=config,
        log=log,
    )


def clean_graph_name(graph):
    """Extract cleaned stem name from graph tuple."""
    name = Path(graph[1]).stem
    return name.replace("_nOH", "")


def run_pair_mode(graphs, base_output, run_name, config, log):
    pair_base_dir = base_output / "PAIR"

    for g1, g2 in combinations(graphs, 2):
        name1 = clean_graph_name(g1)
        name2 = clean_graph_name(g2)

        pair_folder = f"{name1}_vs_{name2}"
        pair_run_name = f"{run_name}_{name1}_{name2}"

        run_association_task(
            graphs=[g1, g2],
            output_path=pair_base_dir / pair_folder,
            run_name=pair_run_name,
            association_config=config,
            log=log,
        )


def main():
    args = parse_args()
    manifest = load_manifest(args.manifest)
    settings = manifest["settings"]

    run_name = settings["run_name"]
    run_mode = settings.get("run_mode", "all")

    if run_mode not in {"all", "pair"}:
        raise ValueError("run_mode must be 'all' or 'pair'")

    base_output = Path(settings["output_path"])
    output_dir = base_output / run_name

    log = setup_logging(
        outdir=output_dir,
        debug=settings.get("debug_logs"),
        verbose=settings.get("verbose"),
    )

    tracker_residues = setup_trackers(output_dir=output_dir, settings=settings)
    association_config = build_association_config(settings, run_mode, tracker_residues)

    graphs = create_graphs(manifest)

    if run_mode == "all":
        run_all_mode(graphs, base_output, run_name, association_config, log)
    else:
        run_pair_mode(graphs, base_output, run_name, association_config, log)

    if tracker_residues:
        out_path = tracker_residues.dump_json()
        log.info(f"Residue tracking report saved to: {out_path}")

    if args.dashboard:
        log.info("Opening dashboard in the default web browser...")
        if run_mode == "all":
            dash_path = base_output / "ALL" / "Dashboard.html"
            if dash_path.exists():
                webbrowser.open(f"file://{dash_path.resolve()}")
        else:
            for dash_path in (base_output / "PAIR").rglob("Dashboard.html"):
                webbrowser.open(f"file://{dash_path.resolve()}")


if __name__ == "__main__":
    main()
