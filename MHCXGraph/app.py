import logging
import os
import webbrowser
from itertools import combinations
from pathlib import Path

from MHCXGraph.cli.cli_parser import parse_args
from MHCXGraph.core.residue_tracking import ResidueTracker
from MHCXGraph.core.tracking import init_tracker
from MHCXGraph.scripts.create_heatmaps import create_heatmap
from MHCXGraph.scripts.renumber_MHCI_imgt import load_mhci_templates, process_structure_file_mhci
from MHCXGraph.scripts.renumber_MHCII_imgt import load_mhcii_templates, process_structure_file_mhcii
from MHCXGraph.utils.logging_utils import setup_logging
from MHCXGraph.utils.preprocessing import create_graphs
from MHCXGraph.workflow.association import run_association_task
from MHCXGraph.workflow.manifest import build_association_config, load_manifest


def setup_trackers(output_dir, settings):
    """
    Initialize runtime tracking utilities.

    This function configures the global tracking system used to store
    intermediate artifacts produced during execution and optionally
    creates a :class:`ResidueTracker` to monitor selected residues.

    Parameters
    ----------
    output_dir : pathlib.Path
        Directory where tracking artifacts and debug files will be stored.

    settings : dict[str, Any]
        Runtime configuration dictionary loaded from the manifest.
        Relevant keys include ``watch_residues`` and ``debug_tracking``.

    Returns
    -------
    tracker_residues : ResidueTracker or None
        Residue tracker instance if residue monitoring is enabled,
        otherwise ``None``.
"""
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
    """
    Execute the association workflow in all-graphs mode.

    In this mode all graphs are processed together in a single
    association task.

    Parameters
    ----------
    graphs : list
        Collection of graph objects produced by the preprocessing stage.

    base_output : pathlib.Path
        Base directory where output results are written.

    run_name : str
        Identifier for the current execution run.

    config : dict[str, Any]
        Association configuration dictionary controlling the
        graph association algorithm.

    log : logging.Logger
        Logger instance used to record runtime messages.

    Returns
    -------
    None
    """
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
    """
    Execute the association workflow in pairwise mode.

    Each unique pair of graphs is processed independently and
    written to a dedicated output directory.

    Parameters
    ----------
    graphs : list
        Collection of graph objects produced by preprocessing.

    base_output : pathlib.Path
        Root directory where pairwise comparison results will be saved.

    run_name : str
        Base identifier for the run.

    config : dict[str, Any]
        Association configuration dictionary controlling the
        graph association algorithm.

    log : logging.Logger
        Logger instance used to record runtime messages.

    Returns
    -------
    None
    """
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


def run(args):
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


def renumber(args):
    if args.mhc_class.upper() == "MHCI":
        load_templates = load_mhci_templates
        process_structure_file = process_structure_file_mhci
    elif args.mhc_class.upper() == "MHCII":
        load_templates = load_mhcii_templates
        process_structure_file = process_structure_file_mhcii
    else:
        raise ValueError(f"{args.mhc_class} is an invalid class of MHC. Please, choose between MHCI and MHCII.")


    log = logging.getLogger("MHCXGraph")

    os.makedirs(args.output_dir, exist_ok=True)
    assets_dir = Path(__file__).resolve().parent / "assets"
    display_csv_path = assets_dir / "imgt_display_all.csv"
    numbering_csv_path = assets_dir / "imgt_numbering_mapping_all.csv"

    templates = load_templates(display_csv_path, numbering_csv_path)

    valid_ext = {".pdb", ".cif", ".mmcif"}
    files = sorted(
        f for f in os.listdir(args.input_dir)
        if os.path.isfile(os.path.join(args.input_dir, f))
        and os.path.splitext(f)[1].lower() in valid_ext
    )

    if not files:
        raise RuntimeError("No .pdb, .cif, or .mmcif files found in input directory.")

    n_ok = 0
    n_fail = 0

    for fname in files:
        input_path = os.path.join(args.input_dir, fname)
        stem, ext = os.path.splitext(fname)
        out_name = f"{stem}{args.suffix}{ext}" if args.suffix else fname
        output_path = os.path.join(args.output_dir, out_name)

        log.info(f"Processing: {fname}")
        try:
            process_structure_file(
                input_path=input_path,
                output_path=output_path,
                templates=templates,
                debug=args.debug,
                warn_score=args.warn_score
            )
            log.info(f"  OK -> {output_path}")
            n_ok += 1
        except Exception as e:
            log.info(f"  FAILED -> {fname}: {e}")
            n_fail += 1

    log.info("\nFinished.")
    log.info(f"  Success: {n_ok}")
    log.info(f"  Failed : {n_fail}")


def main():
    """
    Run the MHCXGraph command-line pipeline.

    This function orchestrates the full workflow:

    1. Parse command-line arguments.
    2. Load the execution manifest.
    3. Configure logging and runtime tracking.
    4. Generate graph representations from input structures.
    5. Execute the association workflow.

    The workflow can operate in two modes defined in the manifest:

    ``all``
        Process all graphs together in a single association task.

    ``pair``
        Perform pairwise comparisons between all graph combinations.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If the configured ``run_mode`` is not ``"all"`` or ``"pair"``.
    """
    args = parse_args()

    if args.command == "run":
        run(args)

    elif args.command == "renumber":
        renumber(args)
    
    elif args.command == "heatmap":
        create_heatmap(args) 

if __name__ == "__main__":
    main()
