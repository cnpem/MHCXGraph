(function () {
  if (typeof network === "undefined") return;

  const STD = __STD_MATRIX__;
  const IDX = __NODE_INDEX__;

  let tempEdgeIds = [];
  let pinnedNode = null;
  const origNodeColors = new Map();

  function clamp(x, a, b) { return Math.min(b, Math.max(a, x)); }
  function isNum(x) { return typeof x === "number" && isFinite(x); }

  function saveOrigColors() {
    const nodes = network.body.data.nodes.get();
    for (const n of nodes) {
      if (!origNodeColors.has(n.id)) origNodeColors.set(n.id, n.color);
    }
  }

  function restoreNodeColors() {
    if (!origNodeColors.size) return;
    const updates = [];
    for (const [id, color] of origNodeColors.entries()) updates.push({ id, color });
    network.body.data.nodes.update(updates);
  }

  function clearTemp() {
    if (tempEdgeIds.length) {
      network.body.data.edges.remove(tempEdgeIds);
      tempEdgeIds = [];
    }
    restoreNodeColors();
  }

  function rowMinMax(i) {
    let mn = Infinity, mx = -Infinity;
    const nodes = network.body.data.nodes.get();

    for (const n of nodes) {
      const j = IDX[n.id];
      if (j === undefined || j === i) continue;
      const s = STD[i][j];
      if (!isNum(s)) continue;
      if (s < mn) mn = s;
      if (s > mx) mx = s;
    }

    if (mn === Infinity || mx === -Infinity) { mn = 0.0; mx = 1.0; }
    if (mx - mn < 1e-12) mx = mn + 1e-12;
    return [mn, mx];
  }

  function tNorm(s, mn, mx) {
    if (!isNum(s)) s = mx;
    return clamp((s - mn) / (mx - mn), 0, 1);
  }

  function edgeRGBA(t) {
    const r = Math.floor(255 * t);
    const g = Math.floor(255 * (1 - t));
    const alphaStrong = 0.95;
    const alphaWeak = 0.12;
    const a = alphaStrong * (1 - t) + alphaWeak * t;
    return `rgba(${r},${g},0,${a.toFixed(3)})`;
  }

  function nodeAlpha(t) {
    const aStrong = 1.0;
    const aWeak = 0.15;
    return aStrong * (1 - t) + aWeak * t;
  }

  function setNodeOpacity(nodeId, alpha) {
    const orig = origNodeColors.get(nodeId);
    let border = "#222222";
    if (orig && typeof orig === "object" && orig.border) border = orig.border;

    const bg = `rgba(33,150,243,${alpha.toFixed(3)})`;
    network.body.data.nodes.update([{ id: nodeId, color: { background: bg, border: border } }]);
  }

  function drawFrom(nodeId) {
    clearTemp();
    saveOrigColors();

    const i = IDX[nodeId];
    if (i === undefined) return;

    const mm = rowMinMax(i);
    const mn = mm[0];
    const mx = mm[1];

    const nodes = network.body.data.nodes.get();
    const newEdges = [];

    setNodeOpacity(nodeId, 1.0);

    for (const n of nodes) {
      const j = IDX[n.id];
      if (j === undefined || j === i) continue;

      const s = STD[i][j];
      const t = tNorm(s, mn, mx);

      setNodeOpacity(n.id, nodeAlpha(t));

      const color = edgeRGBA(t);
      const label = isNum(s) ? s.toFixed(4) : "nan";

      newEdges.push({
        from: nodeId,
        to: n.id,
        color: { color: color },
        width: 1,
        physics: false,
        smooth: false,
        label: label,
        font: { align: "top", size: 10, strokeWidth: 4, strokeColor: "#ffffff" },
        title: "std_norm: " + label
      });
    }

    const added = network.body.data.edges.add(newEdges);
    tempEdgeIds = Array.isArray(added) ? added : [];
  }

  network.on("hoverNode", function (params) {
    if (pinnedNode !== null) return;
    drawFrom(params.node);
  });

  network.on("blurNode", function () {
    if (pinnedNode !== null) return;
    clearTemp();
  });

  network.on("click", function (params) {
    const clickedNode = (params.nodes && params.nodes.length) ? params.nodes[0] : null;

    if (!clickedNode) {
      pinnedNode = null;
      clearTemp();
      return;
    }

    if (pinnedNode === clickedNode) {
      pinnedNode = null;
      clearTemp();
    } else {
      pinnedNode = clickedNode;
      drawFrom(clickedNode);
    }
  });
})();
