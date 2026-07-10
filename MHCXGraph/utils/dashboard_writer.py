"""Constant-memory writer for the self-contained MHCXGraph dashboard.

The original ``create_master_dashboard`` built the entire ``master_export``
dict, serialized it in one ``json.dumps``, and injected it into a template via
a chain of ~20 ``str.replace`` calls. At the peak that held, simultaneously:
the full dict, its JSON string (a second full copy as text), and the giant
HTML string being reallocated on every replace. Peak RAM therefore scaled with
the total data and with pair count.

This module keeps the output byte-for-byte equivalent and still fully
self-contained (everything embedded in one HTML), but never materializes the
whole graph payload at once:

* Each pair is streamed to its own temp file during the association loop and
  dropped from memory immediately (see app.run_*_mode).
* The HTML is assembled as an ordered sequence of segments written directly to
  the output file. The huge ``graph_data`` object is emitted incrementally:
  static header fields, then each pair read back from its temp file one at a
  time, then the deduplicated filtered_graphs. No single giant string, no
  repeated full-string copies, and at most one pair resident at a time.

The JS still receives exactly the same ``masterData`` object shape it did
before, so no frontend change is required.
"""

import json
from pathlib import Path


class DashboardStreamWriter:
    """Assemble the dashboard HTML with constant memory.

    Usage
    -----
    writer = DashboardStreamWriter(assets_dir, log)
    writer.write(
        output_path=...,
        file_name=...,
        header={...},              # top-level scalar/small fields of masterData
        pair_files=[(key, path)],  # each path holds one pair's JSON object
        filtered_graphs={...},     # deduped per-protein filtered graphs
        actual_mode=...,
    )
    """

    def __init__(self, assets_dir, log):
        self.assets_dir = Path(assets_dir)
        self.log = log

    # ------------------------------------------------------------------ assets
    def _load_asset(self, folder, filename, default=""):
        p = Path(folder) / filename
        return p.read_text(encoding="utf-8") if p.exists() else default

    def _b64_image(self, path):
        import base64
        with open(path, "rb") as fh:
            return base64.b64encode(fh.read()).decode("utf-8")

    def _build_template_segments(self, header, actual_mode):
        """Return (prefix, suffix) HTML strings split at the graph_data hole.

        Every injection EXCEPT the graph_data blob is done here, on the static
        template, exactly as the original did. We then split the result on the
        graph_data placeholder so the caller can stream the payload between the
        two halves.
        """
        import re

        assets_dir = self.assets_dir
        html_files = assets_dir / "dashboard"
        js_files = html_files / "js"

        def split_html(raw_html):
            css_match = re.search(r"<style>(.*?)</style>", raw_html, re.DOTALL)
            css = css_match.group(1) if css_match else ""
            clean = re.sub(r"<style>.*?</style>", "", raw_html, flags=re.DOTALL).strip()
            return css, clean

        def script_or_cdn(local_name, cdn_tag):
            p = js_files / local_name
            if p.exists():
                return f'<script>\n{p.read_text(encoding="utf-8")}\n</script>'
            return cdn_tag

        vis_injection = script_or_cdn(
            "vis-network.min.js",
            '<script type="text/javascript" src="https://unpkg.com/vis-network/standalone/umd/vis-network.min.js"></script>',
        )
        mol3d_injection = script_or_cdn(
            "3Dmol-min.js", '<script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>'
        )
        fabric_injection = script_or_cdn(
            "fabric.min.js",
            '<script src="https://cdnjs.cloudflare.com/ajax/libs/fabric.js/5.3.1/fabric.min.js"></script>',
        )
        plotly_injection = script_or_cdn(
            "plotly.min.js",
            '<script src="https://cdnjs.cloudflare.com/ajax/libs/plotly.js/2.32.0/plotly.min.js"></script>',
        )

        mhcx_logo_path = assets_dir / "images/MHCXGraph logo.png"
        mhcx_logo_injection = "<b>MHCXGraph</b>\n"
        favicon_injection = ""
        if mhcx_logo_path.exists():
            enc = self._b64_image(mhcx_logo_path)
            mhcx_logo_injection = f'<img src="data:image/png;base64,{enc}" alt="MHCXGraph Logo" style="width: 100%; height: 100%; margin-bottom: -10px;">'
            favicon_injection = f'<link rel="icon" type="image/png" href="data:image/png;base64,{enc}">'

        logo_injection = ""
        for name, cls in (("images/LNBio.png", "logo-light"), ("images/LNBio white.png", "logo-dark")):
            p = assets_dir / name
            if p.exists():
                enc = self._b64_image(p)
                logo_injection += f'<img src="data:image/png;base64,{enc}" alt="LNBio Logo" class="{cls}" style="height: 8rem; width: auto;">'
        if not logo_injection:
            self.log.debug("LNBio logos not found in assets/images. Skipping logo injection.")

        html_template = self._load_asset(html_files, "base.html")
        if not html_template:
            self.log.error(f"Template base.html not found at {html_files}.")
            return None, None

        sidebar_css, sidebar_dom = split_html(self._load_asset(html_files, "sidebar.html"))
        modal_css, modal_dom = split_html(self._load_asset(html_files, "export_modal.html"))

        # data.js contains the graph_data hole; inject data.js into main.js,
        # main.js into the template, so the hole survives into final_html.
        data_js = self._load_asset(js_files, "data.js")
        main_js = self._load_asset(js_files, "main.js")

        def inject(html, name, code):
            return html.replace(f"__{name.upper()}_JS_INJECTION__", code)

        final_html = html_template.replace("__FAVICON_INJECTION__", favicon_injection)
        final_html = inject(final_html, "vis", vis_injection)
        final_html = inject(final_html, "3Dmol", mol3d_injection)
        final_html = inject(final_html, "fabric", fabric_injection)
        final_html = inject(final_html, "plotly", plotly_injection)
        final_html = final_html.replace("__SIDEBAR_CSS_INJECTION__", sidebar_css)
        final_html = final_html.replace("__MODAL_CSS_INJECTION__", modal_css)
        final_html = final_html.replace("__SIDEBAR_HTML_INJECTION__", sidebar_dom)
        final_html = final_html.replace("__MHCXGRAPH_LOGO_INJECTION__", mhcx_logo_injection)
        final_html = final_html.replace("__LNBIO_LOGO_INJECTION__", logo_injection)
        final_html = final_html.replace("__MODAL_HTML_INJECTION__", modal_dom)

        # main.js still holds __DATA_JS_INJECTION__ = data.js (with the hole)
        main_js = inject(main_js, "data", data_js)
        final_html = inject(final_html, "main", main_js)

        # all remaining JS blocks (none contain the graph_data hole)
        for name, fname in (
            ("init", "init_functions.js"), ("viewer", "viewer.js"),
            ("theme", "theme.js"), ("grid", "grid.js"),
            ("analysis", "analysis.js"), ("graph", "graph.js"),
            ("structures", "structures.js"), ("modal", "export_modal.js"),
        ):
            final_html = inject(final_html, name, self._load_asset(js_files, fname))

        if actual_mode == "screening":
            final_html = final_html.replace("Pairwise View Mode", "Screening Mode (1 vs All)")
            final_html = final_html.replace("Global Pair Analysis", "Global Screening Analysis")
            patch_script = (
                "\n<script>\nwindow.addEventListener('DOMContentLoaded', () => {\n"
                "    if (typeof masterData !== 'undefined' && masterData.actual_mode === 'screening') {\n"
                "        const observer = new MutationObserver(() => {\n"
                "            const metaPanel = document.getElementById('metadata-panel');\n"
                "            if (metaPanel && metaPanel.innerHTML.includes('pairwise')) {\n"
                "                metaPanel.innerHTML = metaPanel.innerHTML.replace(/pairwise/g, 'screening');\n"
                "            }\n        });\n"
                "        observer.observe(document.body, { childList: true, subtree: true });\n"
                "    }\n});\n</script>\n"
            )
            final_html = final_html.replace("</body>", f"{patch_script}\n</body>")

        # Split at the graph_data hole -> stream payload between the halves.
        placeholder = "__GRAPH_DATA_JS_INJECTION__"
        if placeholder not in final_html:
            self.log.error("graph_data placeholder missing from template.")
            return None, None
        prefix, suffix = final_html.split(placeholder, 1)
        return prefix, suffix

    # ------------------------------------------------------------------- write
    def write(self, output_path, file_name, header, pair_files,
              filtered_graphs, actual_mode=None):
        """Stream the full dashboard to disk with constant memory.

        Parameters
        ----------
        header : dict
            Top-level masterData fields (mode, run_name, metadata, proteins,
            protein_paths, actual_mode, reference_structure...). Small.
        pair_files : list[tuple[str, str|Path]]
            (pair_key, path) where each path holds one pair's JSON object.
        filtered_graphs : dict
            Deduplicated per-protein filtered graph payloads. Modest.
        """
        output_path = Path(output_path)
        output_path.mkdir(parents=True, exist_ok=True)
        full_path = output_path / file_name

        prefix, suffix = self._build_template_segments(header, actual_mode)
        if prefix is None:
            return None

        with open(full_path, "w", encoding="utf-8") as out:
            out.write(prefix)

            # ---- stream the graph_data object literal, field by field ----
            out.write("{")

            # header scalars/small arrays: dump each value independently so we
            # never build a combined dict of the whole thing.
            first = True
            for key, value in header.items():
                if not first:
                    out.write(",")
                first = False
                out.write(json.dumps(key))
                out.write(":")
                out.write(json.dumps(value))

            # pairs: read each temp file and stream it in; one pair resident.
            out.write(',"pairs":{')
            for pi, (pair_key, pf) in enumerate(pair_files):
                if pi:
                    out.write(",")
                out.write(json.dumps(pair_key))
                out.write(":")
                # copy the pair JSON verbatim (it is already a JSON object)
                with open(pf, "r", encoding="utf-8") as pin:
                    for chunk in iter(lambda: pin.read(65536), ""):
                        out.write(chunk)
            out.write("}")

            # filtered_graphs: deduped, modest — stream per protein.
            out.write(',"filtered_graphs":{')
            for fi, (pname, fg) in enumerate(filtered_graphs.items()):
                if fi:
                    out.write(",")
                out.write(json.dumps(pname))
                out.write(":")
                out.write(json.dumps(fg))
            out.write("}")

            out.write("}")  # close graph_data object
            # ---------------------------------------------------------------

            out.write(suffix)

        self.log.info(f"Interactive Dashboard saved to {full_path}")
        return full_path
