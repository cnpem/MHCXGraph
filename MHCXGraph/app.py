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

def create_master_dashboard(export_data, output_dir, log):
    """Helper to inject the master aggregated JSON into the dashboard template."""
    import base64
    import json
    import re
    from pathlib import Path

    assets_dir = Path(__file__).resolve().parent / "assets"
    html_files = assets_dir / "dashboard"
    js_files = html_files / "js"

    def _load_asset(folder, filename, default=""):
        p = folder / filename
        return p.read_text(encoding="utf-8") if p.exists() else default

    def split_html(raw_html):
        css_match = re.search(r'<style>(.*?)</style>', raw_html, re.DOTALL)
        css = css_match.group(1) if css_match else ""
        clean_html = re.sub(r'<style>.*?</style>', '', raw_html, flags=re.DOTALL).strip()
        return css, clean_html

    def inject_js(html, name, code):
        placeholder = f"__{name.upper()}_JS_INJECTION__"
        return html.replace(placeholder, code)

    vis_local = js_files / "vis-network.min.js"
    if vis_local.exists():
        vis_injection = f'<script>\n{vis_local.read_text(encoding="utf-8")}\n</script>'
    else:
        vis_injection = '<script type="text/javascript" src="https://unpkg.com/vis-network/standalone/umd/vis-network.min.js"></script>'

    mol3d_local = js_files / "3Dmol-min.js"
    if mol3d_local.exists():
        mol3d_injection = f'<script>\n{mol3d_local.read_text(encoding="utf-8")}\n</script>'
    else:
        mol3d_injection = '<script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>'

    fabric_local = js_files / "fabric.min.js"
    if fabric_local.exists():
        fabric_injection = f'<script>\n{fabric_local.read_text(encoding="utf-8")}\n</script>'
    else:
        fabric_injection = '<script src="https://cdnjs.cloudflare.com/ajax/libs/fabric.js/5.3.1/fabric.min.js"></script>'

    plotly_local = js_files / "plotly.min.js"
    if plotly_local.exists():
        plotly_injection = f'<script>\n{plotly_local.read_text(encoding="utf-8")}\n</script>'
    else:
        plotly_injection = '<script src="https://cdnjs.cloudflare.com/ajax/libs/plotly.js/2.32.0/plotly.min.js"></script>'

    mhcx_logo_path = assets_dir / "images/MHCXGraph logo.png"
    mhcx_logo_injection = "<b>MHCXGraph</b>\n"
    favicon_injection = ""
    if mhcx_logo_path.exists():
        with open(mhcx_logo_path, "rb") as image_file:
            encoded = base64.b64encode(image_file.read()).decode("utf-8")
            mhcx_logo_injection = f'<img src="data:image/png;base64,{encoded}" alt="MHCXGraph Logo" style="width: 100%; height: 100%; margin-bottom: -10px;">'
            favicon_injection = f'<link rel="icon" type="image/png" href="data:image/png;base64,{encoded}">'

    logo_dark_path = assets_dir / "images/LNBio white.png"
    logo_light_path = assets_dir / "images/LNBio.png"
    logo_injection = ""
    if logo_light_path.exists():
        with open(logo_light_path, "rb") as image_file:
            encoded = base64.b64encode(image_file.read()).decode("utf-8")
            logo_injection += f'<img src="data:image/png;base64,{encoded}" alt="LNBio Logo" class="logo-light" style="height: 8rem; width: auto;">'
    if logo_dark_path.exists():
        with open(logo_dark_path, "rb") as image_file:
            encoded = base64.b64encode(image_file.read()).decode("utf-8")
            logo_injection += f'<img src="data:image/png;base64,{encoded}" alt="LNBio Logo" class="logo-dark" style="height: 8rem; width: auto;">'
    if not logo_injection:
        log.debug("LNBio logos not found in assets/images. Skipping logo injection.")

    html_template = _load_asset(html_files, "base.html")
    if not html_template:
        log.error(f"Template base.html not found at {html_files}.")
        return

    sidebar_html = _load_asset(html_files, "sidebar.html")
    modal_html = _load_asset(html_files, "export_modal.html")
    main_js = _load_asset(js_files, "main.js")
    modal_js = _load_asset(js_files, "export_modal.js")
    grid_js = _load_asset(js_files, "grid.js")
    data_js = _load_asset(js_files, "data.js")
    theme_js = _load_asset(js_files, "theme.js")
    viewer_js = _load_asset(js_files, "viewer.js")
    structures_js = _load_asset(js_files, "structures.js")
    init_js = _load_asset(js_files, "init_functions.js")
    analysis_js = _load_asset(js_files, "analysis.js")
    graph_js = _load_asset(js_files, "graph.js")

    sidebar_css, sidebar_dom = split_html(sidebar_html)
    modal_css, modal_dom = split_html(modal_html)

    # Inject everything
    final_html = html_template.replace("__FAVICON_INJECTION__", favicon_injection)
    final_html = inject_js(final_html, "vis", vis_injection)
    final_html = inject_js(final_html, "3Dmol", mol3d_injection)
    final_html = inject_js(final_html, "fabric", fabric_injection) # ADD THIS
    final_html = inject_js(final_html, "plotly", plotly_injection)
    final_html = final_html.replace("__SIDEBAR_CSS_INJECTION__", sidebar_css)
    final_html = final_html.replace("__MODAL_CSS_INJECTION__", modal_css)
    final_html = final_html.replace("__SIDEBAR_HTML_INJECTION__", sidebar_dom)
    final_html = final_html.replace("__MHCXGRAPH_LOGO_INJECTION__", mhcx_logo_injection)
    final_html = final_html.replace("__LNBIO_LOGO_INJECTION__", logo_injection)

    final_html = final_html.replace("__MODAL_HTML_INJECTION__", modal_dom)

    # Javascript injection
    final_html = inject_js(final_html, "main", main_js)
    final_html = inject_js(final_html, "data", data_js)
    final_html = inject_js(final_html, "graph_data", json.dumps(export_data))
    final_html = inject_js(final_html, "init", init_js)
    final_html = inject_js(final_html, "viewer", viewer_js)
    final_html = inject_js(final_html, "theme", theme_js)
    final_html = inject_js(final_html, "grid", grid_js)
    final_html = inject_js(final_html, "analysis", analysis_js)
    final_html = inject_js(final_html, "graph", graph_js)
    final_html = inject_js(final_html, "structures", structures_js)
    final_html = inject_js(final_html, "modal", modal_js)

    # Dynamically name the output file based on the mode
    mode = export_data.get("mode")
    file_name = "Dashboard_Pairwise.html" if mode == "pairwise" else "Dashboard_Multiple.html"
    
    full_path = output_dir / file_name
    with open(str(full_path), "w+", encoding="utf-8") as out:
        out.write(final_html)
    log.info(f"Interactive Dashboard saved to {full_path}")

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


def run_multiple_mode(graphs, base_output, run_name, config, log):
    """
    Execute the association workflow in multiple-graphs mode.

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
    target_dir = base_output / "MULTIPLE"

    G = run_association_task(
        graphs=graphs,
        output_path=target_dir,
        run_name=run_name,
        association_config=config,
        log=log,
    )

    if G and G.associated_graphs is not None:
        global_proteins = [clean_graph_name(g) for g in graphs]
        
        # G.get_dashboard_data correctly formats nodes, edges, components & filtered_graphs
        master_export = G.get_dashboard_data(global_proteins)
        
        # Append the top-level parameters required by the JS frontend
        master_export["mode"] = "multiple"
        master_export["run_name"] = run_name
        master_export["metadata"] = config
        
        create_master_dashboard(master_export, target_dir, log)


def clean_graph_name(graph):
    """Extract cleaned stem name from graph tuple."""
    name = Path(graph[1]).stem
    return name.replace("_nOH", "")


def run_pairwise_mode(graphs, base_output, run_name, config, log):
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
    pair_base_dir = base_output / "PAIRWISE"

    global_proteins = [clean_graph_name(g) for g in graphs]

    master_export = {
        "mode": "pairwise",
        "run_name": run_name,
        "metadata": config,
        "proteins": global_proteins,
        "protein_paths": [str(Path(g[1]).resolve()) for g in graphs],
        "pairs": {}
    }

    for g1, g2 in combinations(graphs, 2):
        name1 = clean_graph_name(g1)
        name2 = clean_graph_name(g2)

        pair_folder = f"{name1}_vs_{name2}"
        pair_key = f"{name1}_vs_{name2}"
        pair_run_name = f"{run_name}_{name1}_{name2}"

        G = run_association_task(
            graphs=[g1, g2],
            output_path=pair_base_dir / pair_folder,
            run_name=pair_run_name,
            association_config=config,
            log=log,
        )
        if G and G.associated_graphs is not None:
            master_export["pairs"][pair_key] = G.get_dashboard_data(global_proteins)

    create_master_dashboard(master_export, pair_base_dir, log)

def run(args):
    manifest = load_manifest(args.manifest)
    settings = manifest["settings"]

    run_name = settings["run_name"]
    run_mode = settings.get("run_mode", "multiple")

    if run_mode not in {"multiple", "pairwise"}:
        raise ValueError("run_mode must be 'multiple' or 'pairwise'")

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

    if run_mode == "multiple":
        run_multiple_mode(graphs, base_output, run_name, association_config, log)
    else:
        run_pairwise_mode(graphs, base_output, run_name, association_config, log)

    if tracker_residues:
        out_path = tracker_residues.dump_json()
        log.info(f"Residue tracking report saved to: {out_path}")

    if args.dashboard:
        log.info("Opening dashboard in the default web browser...")
        if run_mode == "multiple":
            dash_path = base_output / "MULTIPLE" / "Dashboard_Multiple.html"
            if dash_path.exists():
                webbrowser.open(f"file://{dash_path.resolve()}")
        else:
            dash_path = base_output / "PAIRWISE" / "Dashboard_Pairs.html"
            if dash_path.exists():
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

    ``multiple``
        Process all graphs together in a single association task.

    ``pairwise``
        Perform pairwise comparisons between all graph combinations.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If the configured ``run_mode`` is not ``"multiple"`` or ``"pairwise"``.
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
