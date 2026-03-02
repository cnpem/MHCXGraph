import logging
import coloredlogs
import sys
from itertools import combinations
from pathlib import Path

from MHCXGraph.cli.cli_parser import parse_args
from MHCXGraph.core.residue_tracking import ResidueTracker
from MHCXGraph.core.tracking import init_tracker
from MHCXGraph.utils.preprocessing import create_graphs
from MHCXGraph.workflow.association import run_association_task
from MHCXGraph.workflow.manifest import build_association_config, load_manifest
from MHCXGraph.utils.logging_utils import setup_logging, get_log

def main():
    args = parse_args()
    manifest = load_manifest(args.manifest)
    settings = manifest["settings"]

    base_run_name = settings["run_name"]
    base_output_path = Path(settings["output_path"])
    run_mode = settings.get("run_mode", "all")

    tracker_residues = (
        ResidueTracker(settings.get("watch_residues")) if settings["watch_residues"] else None
    )

    init_tracker(
        root="CrossSteps",
        outdir=base_output_path / base_run_name,
        enabled=settings.get("debug_tracking", False),
        prefer_npy_for_ndarray=True,
        add_timestamp_prefix=False,
    )

    logging.getLogger("matplotlib").setLevel(logging.ERROR)
    logging.basicConfig(
        stream=sys.stdout,
        level=logging.DEBUG if settings.get("debug_logs", False) else logging.INFO,
    )

    setup_logging(
        outdir=Path(base_output_path) / base_run_name,
        debug=bool(settings.get("debug_logs", False)),
        verbose=bool(settings.get("verbose", False)),
    )

    log = get_log()
    log.info("Initializing MHCXGraph...")
    log.debug("Debug file logging initialized")

    graphs = create_graphs(manifest)

    base_association_config = build_association_config(settings, run_mode, tracker_residues)

    if run_mode == "all":
        target_dir = base_output_path / "ALL"
        run_association_task(
            graphs=graphs,
            output_path=target_dir,
            run_name=base_run_name,
            association_config=base_association_config,
            log=log,
        )

    elif run_mode == "pair":
        pair_base_dir = base_output_path / "PAIR"

        for g1, g2 in combinations(graphs, 2):
            name1 = Path(g1[1]).stem
            name2 = Path(g2[1]).stem

            name1 = name1.replace("_nOH", "")
            name2 = name2.replace("_nOH", "")

            pair_folder_name = f"{name1}_vs_{name2}"
            target_dir = pair_base_dir / pair_folder_name

            pair_run_name = f"{base_run_name}_{name1}_{name2}"

            run_association_task(
                graphs=[g1, g2],
                output_path=target_dir,
                run_name=pair_run_name,
                association_config=base_association_config,
                log=log,
            )

    else:
        log.error(f"Unknown run_mode: {run_mode}. Please use 'all' or 'pair'.")

    if tracker_residues is not None:
        out_path = tracker_residues.dump_json()
        log.info(f"Residue tracking report saved to: {out_path}")


if __name__ == "__main__":
    main()
