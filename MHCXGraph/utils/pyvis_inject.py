from __future__ import annotations

import json
from importlib import resources
from typing import Any

import numpy as np


def _json_safe_matrix(mat: np.ndarray) -> str:
    arr = np.asarray(mat, dtype=float)
    safe = np.where(np.isfinite(arr), arr, None).tolist()
    return json.dumps(safe, separators=(",", ":"))


def _json_safe_index_map(d: dict[str, int]) -> str:
    return json.dumps(d, separators=(",", ":"))


def load_pyvis_std_js() -> str:
    # MHCXGraph/assets/pyvis_std_hover.js
    return resources.files("MHCXGraph").joinpath("assets/pyvis_std_hover.js").read_text(encoding="utf-8")


def inject_std_hover(
    html: str,
    *,
    std_matrix: np.ndarray,
    safe_node_index: dict[str, int],
) -> str:
    js_template = load_pyvis_std_js()

    js = js_template.replace("__STD_MATRIX__", _json_safe_matrix(std_matrix))
    js = js.replace("__NODE_INDEX__", _json_safe_index_map(safe_node_index))

    block = "<script type=\"text/javascript\">\n" + js + "\n</script>\n"

    if "</body>" in html:
        return html.replace("</body>", block + "</body>")
    return html + block
