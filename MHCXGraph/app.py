import base64
import json
import logging
import os
import webbrowser
from itertools import combinations
from pathlib import Path

from MHCXGraph.cli.cli_parser import parse_args
from MHCXGraph.core.residue_tracking import ResidueTracker
from MHCXGraph.core.tracking import init_tracker
from MHCXGraph.scripts.renumber_MHCI_imgt import load_mhci_templates, process_structure_file_mhci
from MHCXGraph.scripts.renumber_MHCII_imgt import load_mhcii_templates, process_structure_file_mhcii
from MHCXGraph.utils.logging_utils import setup_logging
from MHCXGraph.utils.preprocessing import create_graphs
from MHCXGraph.workflow.association import run_association_task
from MHCXGraph.workflow.manifest import build_association_config, load_manifest

def _assemble_dashboard(graph_data_injection, log):
    """Assemble the full self-contained dashboard HTML with every asset inlined.

    ``graph_data_injection`` is the raw string dropped where the dashboard reads
    ``masterData`` (i.e. ``const masterData = <graph_data_injection>;``). Pass
    ``json.dumps(export_data)`` for a data-baked dashboard, or a sentinel token
    for the standalone (data-free) shell. Mode-specific text patches and the
    output filename are the caller's concern — this only builds the HTML.
    """
    import base64
    import json  # noqa: F401  (kept for parity / callers that pre-dump)
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
    final_html = inject_js(final_html, "graph_data", graph_data_injection)
    final_html = inject_js(final_html, "init", init_js)
    final_html = inject_js(final_html, "viewer", viewer_js)
    final_html = inject_js(final_html, "theme", theme_js)
    final_html = inject_js(final_html, "grid", grid_js)
    final_html = inject_js(final_html, "analysis", analysis_js)
    final_html = inject_js(final_html, "graph", graph_js)
    final_html = inject_js(final_html, "structures", structures_js)
    final_html = inject_js(final_html, "modal", modal_js)

    return final_html


# Runtime-side screening patch, shared by the baked and standalone dashboards.
# For the baked dashboard it is a no-op wrapper the Python side already applied;
# the standalone applies the equivalent text swaps in the browser after the JSON
# is loaded (see the loader shell).
_SCREENING_META_OBSERVER = """
<script>
window.addEventListener('DOMContentLoaded', () => {
    if (typeof masterData !== 'undefined' && masterData.actual_mode === 'screening') {
        const observer = new MutationObserver(() => {
            const metaPanel = document.getElementById('metadata-panel');
            if (metaPanel && metaPanel.innerHTML.includes('pairwise')) {
                metaPanel.innerHTML = metaPanel.innerHTML.replace(/pairwise/g, 'screening');
            }
        });
        observer.observe(document.body, { childList: true, subtree: true });
    }
});
</script>
"""


def _apply_screening_patches(html):
    """Apply the screening-mode text swaps + metadata observer to assembled HTML."""
    html = html.replace("Pairwise View Mode", "Screening Mode (1 vs All)")
    html = html.replace("Global Pair Analysis", "Global Screening Analysis")
    html = html.replace("</body>", f"{_SCREENING_META_OBSERVER}\n</body>")
    return html


def create_master_dashboard(export_data, output_dir, log):
    """Build the data-baked interactive dashboard (unchanged behaviour)."""
    import json
    from pathlib import Path

    output_dir = Path(output_dir)

    final_html = _assemble_dashboard(json.dumps(export_data), log)

    actual_mode = export_data.get("actual_mode", export_data.get("mode"))
    if actual_mode == "screening":
        final_html = _apply_screening_patches(final_html)
        file_name = "Dashboard_Screening.html"
    else:
        mode = export_data.get("mode")
        file_name = "Dashboard_Pairwise.html" if mode == "pairwise" else "Dashboard_Multiple.html"

    full_path = output_dir / file_name
    with open(str(full_path), "w+", encoding="utf-8") as out:
        out.write(final_html)
    log.info(f"Interactive Dashboard saved to {full_path}")


# Sentinel injected where masterData is read. In the standalone shell it is
# replaced at RUNTIME (in the browser) with the user's exported JSON.
_STANDALONE_DATA_SENTINEL = '"__MHCX_STANDALONE_DATA__"'


# Loader shell. Placeholders __MHCX_TEMPLATE_B64__ / __MHCX_SENTINEL_LITERAL__
# are filled by create_standalone_dashboard(). No f-string / .format is used,
# so the {…} and ${…} below are passed through untouched.
_STANDALONE_SHELL = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>MHCXGraph — Standalone Dashboard Loader</title>
<style>
    :root { --accent: #2563eb; }
    * { box-sizing: border-box; }
    html, body { margin: 0; height: 100%; font-family: ui-sans-serif, system-ui, sans-serif; background: #0d0c11; color: #f4f4f5; }
    #mhcx-loader { position: fixed; inset: 0; display: flex; flex-direction: column; align-items: center; justify-content: center; gap: 18px; padding: 24px; text-align: center; z-index: 5; }
    #mhcx-loader h1 { font-size: 20px; font-weight: 800; margin: 0; letter-spacing: 0.3px; }
    #mhcx-loader p { margin: 0; color: #b0a9c0; font-size: 13px; max-width: 520px; line-height: 1.6; }
    #drop { width: min(560px, 92vw); border: 2px dashed #3a3550; border-radius: 14px; padding: 34px 24px; background: #14121a; transition: border-color .2s, background .2s; cursor: pointer; }
    #drop.drag { border-color: var(--accent); background: #171a2b; }
    #drop .big { font-size: 30px; opacity: .7; }
    #drop .hint { margin-top: 8px; color: #8c877c; font-size: 12px; }
    .btn { background: var(--accent); color: #fff; border: none; border-radius: 8px; padding: 9px 16px; font-weight: 600; cursor: pointer; font-size: 13px; }
    .btn:hover { background: #1d4ed8; }
    #err { color: #ff8a80; font-size: 13px; min-height: 16px; }
    #bar { display: none; position: fixed; top: 0; left: 0; right: 0; height: 40px; z-index: 10; background: #0d0c11; border-bottom: 1px solid #1f1b26; align-items: center; gap: 12px; padding: 0 12px; font-size: 12px; color: #b0a9c0; }
    #bar b { color: #f4f4f5; }
    #bar .spacer { flex: 1; }
    #frame { visibility: hidden; position: fixed; top: 40px; left: 0; right: 0; bottom: 0; width: 100%; height: calc(100% - 40px); border: 0; background: #050507; }
    a.ghost { color: #b0a9c0; text-decoration: none; border: 1px solid #1f1b26; padding: 6px 12px; border-radius: 8px; }
    a.ghost:hover { color: #fff; border-color: #3a3550; }
</style>
</head>
<body>

<div id="bar">
    <b>MHCXGraph</b>
    <span id="bar-file"></span>
    <span id="bar-mode"></span>
    <span class="spacer"></span>
    <button class="btn" id="reload-btn">Load different JSON</button>
</div>

<div id="mhcx-loader">
    <h1>MHCXGraph — Standalone Dashboard</h1>
    <p>This is a portable, self-contained dashboard. Load a <code>graph_data_*.json</code>
       exported by any MHCXGraph run and it renders exactly like the live dashboard.</p>
    <div id="drop">
        <div class="big">📂</div>
        <div>Drop your <b>graph_data_*.json</b> here, or click to choose a file</div>
        <div class="hint">Everything runs locally in your browser — nothing is uploaded.</div>
    </div>
    <div id="err"></div>
    <p style="font-size:11px;color:#6c6780;">Tip: to see 3D structures, load the PDB/CIF files from the sidebar once the dashboard opens.</p>
</div>

<iframe id="frame" title="MHCXGraph Dashboard"></iframe>

<input type="file" id="file-input" accept=".json,application/json" style="display:none">

<script>
"use strict";
// TEMPLATE ships as base64, not a JS string literal. A literal embed would
// need to escape every closing-script sequence to stay safe -- but the
// dashboard's own script tags (CDN includes, or vendored vis-network/3Dmol/
// plotly/fabric source if you keep local copies under assets/dashboard/js/)
// can legally contain an HTML comment-open marker followed later by another
// unescaped script-open tag. That combination flips the HTML tokenizer into
// a state where the *next* real closing tag doesn't end this element -- the
// rest spills onto the page as visible text instead of running as code.
// Base64's alphabet has no "<" in it at all, so that failure mode is
// categorically impossible here, whatever any bundled library contains.
const TEMPLATE_B64 = "__MHCX_TEMPLATE_B64__";
const SENTINEL = __MHCX_SENTINEL_LITERAL__;

function b64ToUtf8(b64) {
    const bin = atob(b64);
    const bytes = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
    return new TextDecoder('utf-8').decode(bytes);
}

let TEMPLATE;
try {
    TEMPLATE = b64ToUtf8(TEMPLATE_B64);
} catch (e) {
    document.body.innerHTML = '<pre style="padding:24px;color:#ff8a80;white-space:pre-wrap;">'
        + 'Failed to decode the embedded dashboard template: ' + e.message + '</pre>';
    throw e;
}

// Kept in sync with app.py::_apply_screening_patches (browser-side twin).
const SCREENING_OBSERVER = '<script>window.addEventListener("DOMContentLoaded",function(){if(typeof masterData!=="undefined"&&masterData.actual_mode==="screening"){var o=new MutationObserver(function(){var m=document.getElementById("metadata-panel");if(m&&m.innerHTML.indexOf("pairwise")>=0){m.innerHTML=m.innerHTML.replace(/pairwise/g,"screening");}});o.observe(document.body,{childList:true,subtree:true});}});<\/script>';

const loader = document.getElementById('mhcx-loader');
const drop = document.getElementById('drop');
const fileInput = document.getElementById('file-input');
const errEl = document.getElementById('err');
const frame = document.getElementById('frame');
const bar = document.getElementById('bar');
let currentUrl = null;

function fail(msg) { errEl.textContent = msg; }

function buildDoc(data) {
    // Re-serialize (canonical) and neutralize any stray closing-script tags in
    // string fields so the inlined data can't terminate this document's own
    // script element early.
    let jsonText = JSON.stringify(data).split('</').join('<\\/');
    let doc = TEMPLATE.split(SENTINEL).join(jsonText);
    if (data && data.actual_mode === 'screening') {
        doc = doc.split('Pairwise View Mode').join('Screening Mode (1 vs All)');
        doc = doc.split('Global Pair Analysis').join('Global Screening Analysis');
        doc = doc.replace('</body>', SCREENING_OBSERVER + '\n</body>');
    }
    return doc;
}

function render(data, fileName) {
    let doc;
    try { doc = buildDoc(data); }
    catch (e) { fail('Failed to assemble the dashboard: ' + e.message); return; }

    if (currentUrl) { URL.revokeObjectURL(currentUrl); currentUrl = null; }
    const blob = new Blob([doc], { type: 'text/html' });
    currentUrl = URL.createObjectURL(blob);

    // #frame stays visibility:hidden (never display:none) so it always has a
    // real, properly-sized box before its document loads -- otherwise the
    // loaded document's own "height: 100vh" can resolve against a stale or
    // degenerate viewport from when the iframe had no box at all, and the
    // whole flex layout collapses to content height instead of filling the
    // screen. Nudging a resize on load is a cheap extra safety net for the
    // same class of issue.
    frame.onload = () => {
        try { frame.contentWindow.dispatchEvent(new Event('resize')); } catch (e) {}
    };
    frame.src = currentUrl;

    loader.style.display = 'none';
    frame.style.visibility = 'visible';
    bar.style.display = 'flex';
    document.getElementById('bar-file').textContent = fileName ? ('· ' + fileName) : '';
    const mode = data.actual_mode || data.mode || 'unknown';
    document.getElementById('bar-mode').textContent = '· mode: ' + mode;
}

function ingest(file) {
    fail('');
    if (!file) return;
    const reader = new FileReader();
    reader.onload = () => {
        let data;
        try { data = JSON.parse(reader.result); }
        catch (e) { fail('That file is not valid JSON: ' + e.message); return; }
        if (!data || typeof data !== 'object' || (!data.pairs && !data.nodes)) {
            fail('This does not look like a MHCXGraph export (missing "pairs"/"nodes").');
            return;
        }
        render(data, file.name);
    };
    reader.onerror = () => fail('Could not read the file.');
    reader.readAsText(file);
}

drop.addEventListener('click', () => fileInput.click());
fileInput.addEventListener('change', e => ingest(e.target.files && e.target.files[0]));
['dragenter', 'dragover'].forEach(ev => drop.addEventListener(ev, e => { e.preventDefault(); drop.classList.add('drag'); }));
['dragleave', 'drop'].forEach(ev => drop.addEventListener(ev, e => { e.preventDefault(); drop.classList.remove('drag'); }));
drop.addEventListener('drop', e => { const f = e.dataTransfer && e.dataTransfer.files && e.dataTransfer.files[0]; ingest(f); });
document.getElementById('reload-btn').addEventListener('click', () => {
    frame.style.visibility = 'hidden'; bar.style.display = 'none'; loader.style.display = 'flex';
    if (currentUrl) { URL.revokeObjectURL(currentUrl); currentUrl = null; }
    frame.src = 'about:blank'; fileInput.value = '';
});
</script>
</body>
</html>
"""


def create_standalone_dashboard(output_dir, log, file_name="MHCXGraph_Standalone.html"):
    """Write a data-free dashboard that renders any exported ``graph_data_*.json``.

    The full dashboard (vis-network, 3DMol, Plotly, all app JS/CSS) is assembled
    once with a data sentinel, embedded as base64 (decoded client-side), and
    rendered in an <iframe> after the user picks a JSON file. Because the exact
    same template the baked dashboard uses is reused verbatim, the standalone
    stays in lockstep with the real dashboard — no second copy of the frontend
    to maintain.

    The template is base64, not a literal JS string: the assembled HTML is
    riddled with real "<script>...</script>" tags (CDN includes, or vendored
    vis-network/3Dmol/plotly/fabric source if you keep local copies under
    assets/dashboard/js/). Their content can legally contain "<!--" followed
    later by an unescaped "<script" -- the HTML tokenizer treats that as the
    start of the legacy script-hiding trick and won't close the outer <script>
    on the next "</script>" it sees, truncating this element early and
    spilling the remainder onto the page as visible text. A plain "</" ->
    "<\\/" escape only guards the closing tag, not that state transition.
    Base64's alphabet has no "<" in it at all, so the failure mode above is
    categorically impossible, at the cost of ~33% more bytes for this string.

    The JSON it consumes is precisely what ``write_master_json`` writes on every
    run (regardless of ``generate_dashboard``), so no change to the export
    payload is needed. Call this once, independent of any run — the resulting
    HTML has no data of its own and can load any graph_data_*.json export.

    Note: auto-loading 3D structures from ``protein_paths`` won't work from a
    portable file (blob origin can't read local paths); users upload PDBs via the
    sidebar's per-protein file inputs, exactly as in the baked dashboard offline.
    """
    from pathlib import Path

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    template = _assemble_dashboard(_STANDALONE_DATA_SENTINEL, log)
    template_b64 = base64.b64encode(template.encode("utf-8")).decode("ascii")
    sentinel_js_literal = json.dumps(_STANDALONE_DATA_SENTINEL)

    shell = _STANDALONE_SHELL.replace("__MHCX_TEMPLATE_B64__", template_b64)
    shell = shell.replace("__MHCX_SENTINEL_LITERAL__", sentinel_js_literal)

    full_path = output_dir / file_name
    with open(str(full_path), "w", encoding="utf-8") as out:
        out.write(shell)
    log.info(f"Standalone (data-free) dashboard saved to {full_path}")
    return full_path


def write_master_json(export_data, output_dir, log):
    """Write the raw dashboard payload as a standalone JSON file.

    This is the same ``master_export`` structure the dashboard consumes
    (``__GRAPH_DATA_JS_INJECTION__``), so a dashboard can be regenerated from
    it later. Written when dashboard generation is disabled — or always, if you
    want the data decoupled from the HTML. Cheap: a single json.dump, no asset
    loading and no template string assembly.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    actual_mode = export_data.get("actual_mode", export_data.get("mode"))
    if actual_mode == "screening":
        file_name = "graph_data_screening.json"
    elif export_data.get("mode") == "pairwise":
        file_name = "graph_data_pairwise.json"
    else:
        file_name = "graph_data_multiple.json"

    full_path = output_dir / file_name
    with open(str(full_path), "w", encoding="utf-8") as out:
        json.dump(export_data, out)
    log.info(f"Graph data JSON saved to {full_path}")
    return full_path


def emit_output(export_data, output_dir, generate_dashboard, log):
    """Route the assembled payload to HTML dashboard and/or raw JSON.

    generate_dashboard=True  -> build the self-contained (data-baked) HTML.
    generate_dashboard=False -> skip the HTML (fast, no assets).
    The JSON payload (graph_data_*.json) is always written, so it's there to
    load into MHCXGraph_Standalone.html later without re-running the pipeline
    (see: `mhcxgraph standalone-dashboard`).
    """
    if generate_dashboard:
        create_master_dashboard(export_data, output_dir, log)
    write_master_json(export_data, output_dir, log)


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


def run_multiple_mode(specs, base_output, run_name, config, log, generate_dashboard=True):
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

    graphs = [s.build() for s in specs] 

    G = run_association_task(
        graphs=graphs,
        output_path=target_dir,
        run_name=run_name,
        association_config=config,
        log=log,
    )

    if G and G.associated_graphs is not None:
        global_proteins = [clean_graph_name(g) for g in graphs]

        master_export = G.get_dashboard_data(global_proteins)

        master_export["mode"] = "multiple"
        master_export["run_name"] = run_name
        master_export["metadata"] = config

        emit_output(master_export, target_dir, generate_dashboard, log)


def clean_graph_name(graph):
    """Extract cleaned stem name from graph tuple."""
    name = Path(graph[1]).stem
    return name.replace("_nOH", "")


def run_pairwise_mode(specs, base_output, run_name, config, log, generate_dashboard=True):
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

    global_proteins = [clean_graph_name(s) for s in specs]

    master_export = {
        "mode": "pairwise",
        "run_name": run_name,
        "metadata": config,
        "proteins": global_proteins,
        "protein_paths": [str(Path(g[1]).resolve()) for g in specs],
        # Each protein's filtered graph is identical across every pair it
        # appears in. Store it ONCE here (keyed by protein name) instead of
        # duplicating it into every pair's payload — that duplication was the
        # main driver of the linear/superlinear RAM growth.
        "filtered_graphs": {},
        "pairs": {}
    }

    n = len(specs)
    for i in range(n):
        g1 = specs[i].build()          # pinned across the whole inner sweep
        name1 = clean_graph_name(g1)

        for j in range(i + 1, n):
            g2 = specs[j].build()      # built, used, dropped
            name2 = clean_graph_name(g2)

            pair_key = f"{name1}_vs_{name2}"
            G = run_association_task(
                graphs=[g1, g2],
                output_path=pair_base_dir / pair_key,
                run_name=f"{run_name}_{name1}_{name2}",
                association_config=config,
                log=log,
            )
            if G and G.associated_graphs is not None:
                # Pair payload WITHOUT the redundant per-protein filtered graphs.
                master_export["pairs"][pair_key] = G.get_dashboard_data(
                    global_proteins, include_filtered_graphs=False
                )
                # Store each protein's filtered graph exactly once.
                for prot_idx, gd in enumerate(G.graphs_data):
                    pname = gd["name"]
                    if pname not in master_export["filtered_graphs"]:
                        model_idx = global_proteins.index(pname)
                        master_export["filtered_graphs"][pname] = \
                            G.get_filtered_graph_data(prot_idx, model_idx)

            del g2

        del g1

    emit_output(master_export, pair_base_dir, generate_dashboard, log)

def run_screening_mode(ref_spec, target_specs, base_output, run_name, config, log, generate_dashboard=True):
    """
    Execute the association workflow in screening mode (1-vs-All).

    This mode compares a single reference graph against a collection of target 
    graphs. Each target is processed individually against the reference, and 
    the results are aggregated into a single interactive dashboard. To leverage 
    existing frontend logic, the dashboard payload mimics the "pairwise" mode 
    structure but includes an `actual_mode` flag to trigger specific UI text 
    replacements during HTML generation.

    Parameters
    ----------
    ref_graph : tuple
        A tuple containing the reference graph data produced by the preprocessing 
        stage. Typically structured as `(networkx.Graph, file_path, base_name)`.
    target_graphs : list of tuple
        A list of graph tuples to be compared against the reference graph.
    base_output : pathlib.Path
        The root directory where the screening results and the final HTML 
        dashboard will be saved.
    run_name : str
        A unique base identifier for the current execution run.
    config : dict[str, Any]
        The association configuration dictionary controlling the graph 
        association algorithm's parameters and thresholds.
    log : logging.Logger
        Logger instance used to record runtime progress, warnings, and errors.

    Returns
    -------
    None
    """
    if not target_specs:
        log.error("Screening mode requires at least 1 target graph alongside the reference.")
        return

    screening_base_dir = base_output / "SCREENING"

    ref_graph = ref_spec.build()
    ref_name = clean_graph_name(ref_graph)

    all_specs = [ref_spec, *list(target_specs)]
    global_proteins = [clean_graph_name(s) for s in all_specs]


    master_export = {
        "mode": "pairwise",
        "actual_mode": "screening",
        "reference_structure": ref_name,
        "run_name": run_name,
        "metadata": config,
        "proteins": global_proteins,
        "protein_paths": [str(Path(g[1]).resolve()) for g in all_specs],
        "filtered_graphs": {},
        "pairs": {}
    }


    for target_spec in target_specs:
        target_graph = target_spec.build()
        target_name = clean_graph_name(target_spec)

        pair_folder = f"{ref_name}_vs_{target_name}"
        pair_key = f"{ref_name}_vs_{target_name}"
        pair_run_name = f"{run_name}_{ref_name}_{target_name}"

        G = run_association_task(
            graphs=[ref_graph, target_graph],
            output_path=screening_base_dir / pair_folder,
            run_name=pair_run_name,
            association_config=config,
            log=log,
        )
        if G and G.associated_graphs is not None:
            master_export["pairs"][pair_key] = G.get_dashboard_data(
                global_proteins, include_filtered_graphs=False
            )
            for prot_idx, gd in enumerate(G.graphs_data):
                pname = gd["name"]
                if pname not in master_export["filtered_graphs"]:
                    model_idx = global_proteins.index(pname)
                    master_export["filtered_graphs"][pname] = \
                        G.get_filtered_graph_data(prot_idx, model_idx)

        del target_graph


    emit_output(master_export, screening_base_dir, generate_dashboard, log)


def run(args):
    manifest = load_manifest(args.manifest)
    settings = manifest["settings"]

    run_name = settings["run_name"]
    run_mode = settings.get("run_mode")

    if run_mode not in {"multiple", "pairwise", "screening"}:
        raise ValueError("run_mode must be 'multiple', 'pairwise' or 'screening'")

    base_output = Path(settings["output_path"])
    output_dir = base_output / run_name

    log = setup_logging(
        outdir=output_dir,
        debug=settings.get("debug_logs"),
        verbose=settings.get("verbose"),
    )

    tracker_residues = setup_trackers(output_dir=output_dir, settings=settings)
    association_config = build_association_config(settings, run_mode, tracker_residues)
    generate_dashboard = settings.get("generate_dashboard", True)

    specs = create_graphs(manifest)

    if run_mode == "multiple":
        run_multiple_mode(specs, base_output, run_name, association_config, log, generate_dashboard)
    elif run_mode == "pairwise":
        run_pairwise_mode(specs, base_output, run_name, association_config, log, generate_dashboard)
    elif run_mode == "screening":
        ref_name = settings.get("reference_structure")
        if not ref_name:
            raise ValueError("Screening mode requires 'reference_structure' to be defined in the manifest settings.")
        
        ref_spec = next((s for s in specs if clean_graph_name(s) == ref_name), None)

        if not ref_spec:
            raise ValueError(f"Reference structure '{ref_name}' not found among the input graphs.")

        target_specs = [s for s in specs if clean_graph_name(s) != ref_name]

        run_screening_mode(ref_spec, target_specs, base_output, run_name, association_config, log, generate_dashboard)

    if tracker_residues:
        out_path = tracker_residues.dump_json()
        log.info(f"Residue tracking report saved to: {out_path}")

    if args.dashboard and generate_dashboard:
        log.info("Opening dashboard in the default web browser...")
        dash_path = None
        if run_mode == "multiple":
            dash_path = base_output / "MULTIPLE" / "Dashboard_Multiple.html"
        elif run_mode == "pairwise":
            dash_path = base_output / "PAIRWISE" / "Dashboard_Pairwise.html"
        elif run_mode == "screening":
            dash_path = base_output / "SCREENING" / "Dashboard_Screening.html"

        if dash_path and dash_path.exists():
            webbrowser.open(f"file://{dash_path.resolve()}")
    elif args.dashboard and not generate_dashboard:
        log.info("Dashboard generation disabled (generate_dashboard=false); only JSON was written.")


def standalone_dashboard_command(args):
    """Generate the data-free MHCXGraph_Standalone.html loader on its own.

    No manifest, no pipeline run — the template it embeds carries no data
    (it's built around a sentinel), so this is a one-time download: reuse
    the same HTML file to load any graph_data_*.json export a run produced.
    """
    log = logging.getLogger("MHCXGraph")
    os.makedirs(args.output_dir, exist_ok=True)
    path = create_standalone_dashboard(args.output_dir, log)
    log.info(f"Open {path} and load a graph_data_*.json file exported by a previous run.")


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

    elif args.command == "standalone-dashboard":
        standalone_dashboard_command(args)

    elif args.command == "renumber":
        renumber(args)

    elif args.command == "heatmap":
        from MHCXGraph.scripts.create_heatmaps import create_heatmap

        create_heatmap(args)

if __name__ == "__main__":
    main()
