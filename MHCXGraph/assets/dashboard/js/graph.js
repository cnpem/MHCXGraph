function getCSSVar(v) {
    const val = getComputedStyle(document.body)
        .getPropertyValue(v)
        .trim();

    return val;
}

function getStdColor(std, threshold) {
    if (std === null || std === undefined || isNaN(std)) return '#d1d5db'; 
    let ratio = std / (threshold + 0.0001);
    if (ratio > 1) ratio = 1; if (ratio < 0) ratio = 0;
    return `hsl(${120 * (1 - ratio)}, 100%, 45%)`;
}
function hexToHsl(hex) {
    hex = hex.replace('#', '');

    if (hex.length === 3) {
        hex = hex.split('').map(c => c + c).join('');
    }

    const r = parseInt(hex.slice(0, 2), 16) / 255;
    const g = parseInt(hex.slice(2, 4), 16) / 255;
    const b = parseInt(hex.slice(4, 6), 16) / 255;

    const max = Math.max(r, g, b);
    const min = Math.min(r, g, b);
    const d = max - min;

    let h = 0;
    let s = 0;
    let l = (max + min) / 2;

    if (d !== 0) {
        s = d / (1 - Math.abs(2 * l - 1));

        switch (max) {
            case r:
                h = 60 * (((g - b) / d) % 6);
                break;
            case g:
                h = 60 * ((b - r) / d + 2);
                break;
            case b:
                h = 60 * ((r - g) / d + 4);
                break;
        }
    }

    if (h < 0) h += 360;

    return { h, s: s * 100, l: l * 100 };
}

function hslToHex(h, s, l) {
    s /= 100;
    l /= 100;

    const c = (1 - Math.abs(2 * l - 1)) * s;
    const x = c * (1 - Math.abs((h / 60) % 2 - 1));
    const m = l - c / 2;

    let r = 0, g = 0, b = 0;

    if (h < 60) [r, g, b] = [c, x, 0];
    else if (h < 120) [r, g, b] = [x, c, 0];
    else if (h < 180) [r, g, b] = [0, c, x];
    else if (h < 240) [r, g, b] = [0, x, c];
    else if (h < 300) [r, g, b] = [x, 0, c];
    else [r, g, b] = [c, 0, x];

    const toHex = v => {
        const n = Math.round((v + m) * 255);
        return n.toString(16).padStart(2, '0');
    };

    return `#${toHex(r)}${toHex(g)}${toHex(b)}`;
}

function clamp(value, min, max) {
    return Math.max(min, Math.min(max, value));
}

function deriveVariantsFromBase(hex) {
    const { h, s, l } = hexToHsl(hex);

    const variants = [
        { dh: 0, ds: 0, dl: 0 },     
        { dh: 0, ds: 0, dl: -18 },   
        { dh: 0, ds: 0, dl: 18 },    
        { dh: 8, ds: -5, dl: -10 },  
        { dh: -8, ds: -5, dl: 10 }   
    ];

    return variants.map(v =>
        hslToHex(
            (h + v.dh + 360) % 360,
            clamp(s + v.ds, 35, 90),
            clamp(l + v.dl, 25, 80)
        )
    );
}

function buildExtendedPalette(basePalette, needed) {
    const extended = [];
    const seen = new Set();

    const families = basePalette.map(deriveVariantsFromBase);

    for (let layer = 0; extended.length < needed; layer++) {
        for (let i = 0; i < families.length && extended.length < needed; i++) {
            const family = families[i];
            if (layer < family.length) {
                const color = family[layer];
                if (!seen.has(color)) {
                    seen.add(color);
                    extended.push(color);
                }
            }
        }
    }

    return extended;
}

function clearDynamicStdEdges() {
    if (tempEdgeIds.length) { edgesDataset.remove(tempEdgeIds); tempEdgeIds = []; }
    if (currentGraphMode === 'associated') applyGraphFilters();
}

function drawDynamicStdEdges(nodeId) {
    if (currentGraphMode !== 'associated') return;
    clearDynamicStdEdges();
    
    const globalThresh = (graphData.metadata && graphData.metadata.global_distance_diff_threshold) || (masterData.metadata && masterData.metadata.global_distance_diff_threshold) || 2.0;
    let activeComp = null;
    for (let c of graphData.components) { if (c.id !== 0 && c.node_index_map && c.node_index_map[nodeId] !== undefined) { activeComp = c; break; } }
    if (!activeComp || !activeComp.std_matrix) return;

    const i = activeComp.node_index_map[nodeId];
    const newEdges = [], nodeUpdates = [];
    const isDark = document.body.getAttribute('data-theme') === 'dark';
    const dimBg = getCSSVar('--dim-bg');
    const dimBorder = getCSSVar('--dim-border');

    graphData.nodes.forEach(n => {
        if (n.id === nodeId) nodeUpdates.push({ id: n.id, color: { background: n.originalColor, border: themeBorder } });
        else nodeUpdates.push({ id: n.id, color: { background: dimBg, border: dimBorder } });
    });

    currentActiveNodes.forEach(targetId => {
        if (targetId === nodeId) return;
        const j = activeComp.node_index_map[targetId];
        if (j === undefined || activeComp.std_matrix[i][j] === null) return;
        
        const std = activeComp.std_matrix[i][j];
        const color = getStdColor(std, globalThresh);
        const targetNode = graphData.nodes.find(n => n.id === targetId);
        
        nodeUpdates.push({ id: targetId, color: { background: targetNode.originalColor, border: color } });
        newEdges.push({
            from: nodeId, to: targetId, color: { color: color, opacity: 0.9 }, width: optEdgeWidth * 2.5, physics: false, smooth: false,
            label: std.toFixed(3), 
            font: { align: "top", size: optLabelSize-2, color: getCSSVar('--text-main'), strokeWidth: 4, strokeColor: getCSSVar('--bg-panel'), face: "sans-serif" },
            title: "STD: " + std.toFixed(4)
        });
    });

    nodesDataset.update(
        nodeUpdates.map(n => ({
            id: n.id,
            color: {
                background: n.color.background,
                border: n.color.border
            }
        }))
    )
    const added = edgesDataset.add(newEdges);
    tempEdgeIds = Array.isArray(added) ? added : [added];
}

function assignAllPairGraphColors() {
    logDebug("Assigning global graph palettes for all pairs");
    try {
        // const pal = GRAPH_PALETTES[activePaletteName];
        let uniqueGroups = [...new Set(graphData.nodes.map(n => n.group))].sort();
        const basePalette = GRAPH_PALETTES[activePaletteName];
        const pal = buildExtendedPalette(basePalette, uniqueGroups.length);
        logDebug("Global palette selected", {
            activePaletteName,
            palette: pal
        });

        Object.entries(masterData.pairs).forEach(([pairKey, pData]) => {
            let uniqueGroups = [...new Set(pData.nodes.map(n => n.group))].sort();
            let cMap = {};
            uniqueGroups.forEach((g, i) => cMap[g] = pal[i % pal.length]);
            pData.nodes.forEach(n => n.originalColor = cMap[n.group]);

            logDebug("Pair palette assigned", {
                pairKey,
                uniqueGroups,
                cMap,
                sample: pData.nodes.slice(0, 5).map(n => ({
                    id: n.id,
                    group: n.group,
                    color: n.originalColor
                }))
            });
        });
    } catch(e) {
        logError("assignAllPairGraphColors failed", e);
    }
}
function assignGraphColorsToData() {
    try {
        if (!graphData) return;

        // const pal = GRAPH_PALETTES[activePaletteName];
        let uniqueGroups = [...new Set(graphData.nodes.map(n => n.group))].sort();
        const basePalette = GRAPH_PALETTES[activePaletteName];
        const pal = buildExtendedPalette(basePalette, uniqueGroups.length);

        graphNodeColors = {};
        uniqueGroups.forEach((g, i) => {
            graphNodeColors[g] = pal[i % pal.length];
        });

        graphData.nodes.forEach(n => {
            n.originalColor = graphNodeColors[n.group];
        });

        graphData.filtered_graphs.forEach(fg => {
            let fgGroups = [...new Set(fg.nodes.map(n => n.group))].sort();
            let fgColors = {};
            fgGroups.forEach((g, i) => fgColors[g] = pal[i % pal.length]);
            fg.nodes.forEach(n => n.originalColor = fgColors[n.group]);
        });

        logDebug("assignGraphColorsToData applied", {
            activePaletteName,
            uniqueGroups,
            graphNodeColors,
            sampleNodes: graphData.nodes.slice(0, 8).map(n => ({
                id: n.id,
                group: n.group,
                color: n.originalColor
            }))
        });
    } catch(e) {
        logError("assignGraphColorsToData crashed", e);
    }
}

function handleFilteredChange() {
    activeFilteredProtIdx = parseInt(document.querySelector('input[name="filtered_prot"]:checked').value);
    const thicknessMode = document.querySelector('input[name="edgeThickness"]:checked').value;
    const maxDist = (graphData.metadata && graphData.metadata.edge_threshold) || (masterData.metadata && masterData.metadata.edge_threshold) || 10.0;
    const fg = graphData.filtered_graphs[activeFilteredProtIdx];
    
    nodesDataset.clear(); edgesDataset.clear(); clearDynamicStdEdges();
    
    const coloredNodes = fg.nodes.map(n => ({...n, label: '<b>' + n.label + '</b>', color: { background: n.originalColor, border: themeBorder }, font: { color: themeText } }));
    const styledEdges = fg.edges.map(e => {
        let edgeWidth = optEdgeWidth;
        if (thicknessMode === 'distance' && e.raw_dist !== undefined && e.raw_dist !== null) {
            let dist = parseFloat(e.raw_dist);
            if (!isNaN(dist)) edgeWidth = optEdgeWidth * (4.0 - ((Math.max(0, Math.min(dist, maxDist))) / maxDist) * 3.8); 
        }
        return { ...e, color: { color: getCSSVar('--edge-default'), opacity: 1.0 }, width: Math.max(0.1, edgeWidth) };
    });

    nodesDataset.add(coloredNodes); edgesDataset.add(styledEdges);
    currentActiveNodes = new Set(fg.nodes.map(n => n.id));
    update3DViewer(); updateNodeListText();
    if (network) network.fit();
}

function handleTreeChange() {
    if (currentGraphMode !== 'associated') return;
    try {
        if (nodesDataset.length !== graphData.nodes.length) {
            nodesDataset.clear(); edgesDataset.clear();
            nodesDataset.add(graphData.nodes.map(n => ({...n, label: '<b>' + n.label + '</b>', font: { color: themeText } })));
            edgesDataset.add(graphData.edges);
        }

        let activeIds = new Set();
        document.querySelectorAll('.tree-cb:checked').forEach(cb => {
            let parts = cb.value.split('_');
            if (parts[0] === 'global') graphData.nodes.forEach(n => activeIds.add(n.id));
            else if (parts[0] === 'comp') { let c = graphData.components.find(c => c.id === parseInt(parts[1])); if(c) c.node_ids.forEach(id => activeIds.add(id)); } 
            else if (parts[0] === 'frame') { let c = graphData.components.find(c => c.id === parseInt(parts[1])); if(c) { let f = c.frames.find(f => f.id === parseInt(parts[2])); if(f) f.node_ids.forEach(id => activeIds.add(id)); } }
        });
        currentActiveNodes = activeIds;
        clearDynamicStdEdges(); applyGraphFilters(); update3DViewer(); updateNodeListText();
    } catch(e) { logError("handleTreeChange crashed", e); }
}

function handleGridTreeChange(pairKey) {
    let pData = masterData.pairs[pairKey];
    let dsNodes = gridNodesDatasets[pairKey]; let dsEdges = gridEdgesDatasets[pairKey];
    if (!dsNodes || !dsEdges) return;

    if (dsNodes.length !== pData.nodes.length) {
        dsNodes.clear(); dsEdges.clear();
        dsNodes.add(pData.nodes.map(n => ({...n, label: '<b>'+n.label+'</b>', font: { color: themeText }})));
        dsEdges.add(pData.edges.map(e => ({...e, color: {color: getCSSVar('--edge-default'), opacity: 0.8}, width: optEdgeWidth})));
    }

    let activeIds = new Set();
    let checkboxes = document.querySelectorAll(`.tree-cb-grid[data-pair="${pairKey}"]:checked`);
    checkboxes.forEach(cb => {
        let parts = cb.value.split('_');
        if (parts[0] === 'global') pData.nodes.forEach(n => activeIds.add(n.id));
        else if (parts[0] === 'comp') { let c = pData.components.find(c => c.id === parseInt(parts[1])); if(c) c.node_ids.forEach(id => activeIds.add(id)); } 
        else if (parts[0] === 'frame') { let c = pData.components.find(c => c.id === parseInt(parts[1])); if(c) { let f = c.frames.find(f => f.id === parseInt(parts[2])); if(f) f.node_ids.forEach(id => activeIds.add(id)); } }
    });
    
    if (!window.gridActiveNodes) window.gridActiveNodes = {};
    window.gridActiveNodes[pairKey] = activeIds;

    applyGraphFiltersGrid(pairKey, activeIds, dsNodes, dsEdges, pData);
    update3DViewerGridForPair(pairKey);
    document.getElementById('active-nodes-text').innerHTML = "<i>Active nodes highlighted in Grid.</i>";
}

function applyGraphFiltersGrid(pairKey, activeIds, dsNodes, dsEdges, pData) {
    if(!pData) pData = masterData.pairs[pairKey];
    if(!activeIds) activeIds = (window.gridActiveNodes && window.gridActiveNodes[pairKey]) ? window.gridActiveNodes[pairKey] : new Set(pData.nodes.map(n=>n.id));
    if(!dsNodes) dsNodes = gridNodesDatasets[pairKey];
    if(!dsEdges) dsEdges = gridEdgesDatasets[pairKey];
    
    const mode = document.querySelector('input[name="viewMode"]:checked').value;
    const isDark = document.body.getAttribute('data-theme') === 'dark';
    const dimBg = getCSSVar('--dim-bg');
    const dimBorder = getCSSVar('--dim-border');
    const dimText = getCSSVar('--dim-text');

    const nodeUpdates = [], edgeUpdates = [];
    pData.nodes.forEach(n => {
        const isActive = activeIds.has(n.id);
        if (mode === 'hide') nodeUpdates.push({ id: n.id, hidden: !isActive, color: { background: n.originalColor, border: themeBorder }, font: { color: themeText } });
        else nodeUpdates.push({ id: n.id, hidden: false, color: isActive ? { background: n.originalColor, border: themeBorder } : { background: dimBg, border: dimBorder }, font: { color: isActive ? themeText : dimText } });
    });

    pData.edges.forEach(e => {
        const isEdgeActive = activeIds.has(e.from) && activeIds.has(e.to);
        if (mode === 'hide') edgeUpdates.push({ id: e.id, hidden: !isEdgeActive, color: {color: getCSSVar('--edge-default')}, width: optEdgeWidth });
        else edgeUpdates.push({ id: e.id, hidden: false, color: isEdgeActive ? { color: getCSSVar('--edge-default'), opacity: 1.0 } : { color: getCSSVar('--edge-faded'), opacity: 0.15 }, width: isEdgeActive ? optEdgeWidth : Math.max(0.1, optEdgeWidth * 0.5) });
    });

    dsNodes.update(nodeUpdates);
    dsEdges.update(edgeUpdates);
}

function updateNodeListText() {
    const container = document.getElementById('active-nodes-text');
    if (currentActiveNodes.size === 0) { container.innerHTML = "<i>No nodes selected.</i>"; return; }

    let protNodes = {};
    graphData.proteins.forEach((p, i) => protNodes[i] = []);
    let nodesToIterate = currentGraphMode === 'associated' ? graphData.nodes : graphData.filtered_graphs[activeFilteredProtIdx].nodes;

    nodesToIterate.forEach(n => {
        if (currentActiveNodes.has(n.id)) {
            n.mapping.forEach(m => {
                let localIdx = graphData.proteins.indexOf(masterData.proteins[m.model_idx]);
                if (protNodes[localIdx]) protNodes[localIdx].push(`${m.chain}:${m.resn}:${m.resi}`);
            });
        }
    });

    let html = "";
    // const pal = GRAPH_PALETTES[activePaletteName];
    let uniqueGroups = [...new Set(graphData.nodes.map(n => n.group))].sort();
    const basePalette = GRAPH_PALETTES[activePaletteName];
    const pal = buildExtendedPalette(basePalette, uniqueGroups.length);
    for (let i = 0; i < graphData.proteins.length; i++) {
        if (protNodes[i] && protNodes[i].length > 0) {
            let uniqueNodes = [...new Set(protNodes[i])];
            let globalIdx = masterData.proteins.indexOf(graphData.proteins[i]);
            const color = pal[globalIdx % pal.length];
            html += `<div style="margin-bottom: 6px;"><b style="color:${color}">Prot ${i}:</b> ${uniqueNodes.join(', ')}</div>`;
        }
    }
    container.innerHTML = html;
}

function applyGraphFilters(force = false) {
    if (currentGraphMode !== 'associated') return;
    if (!force && tempEdgeIds.length > 0) return;

    const mode = document.querySelector('input[name="viewMode"]:checked').value;
    const isDark = document.body.getAttribute('data-theme') === 'dark';
    const dimBg = getCSSVar('--dim-bg');
    const dimBorder = getCSSVar('--dim-border');
    const dimText = getCSSVar('--dim-text');

    const nodeUpdates = [], edgeUpdates = [];

    graphData.nodes.forEach(n => {
        const isActive = currentActiveNodes.has(n.id);
        if (mode === 'hide') {
            nodeUpdates.push({
                id: n.id,
                hidden: !isActive,
                color: { background: n.originalColor, border: themeBorder },
                font: { color: themeText }
            });
        } else {
            nodeUpdates.push({
                id: n.id,
                hidden: false,
                color: isActive
                    ? { background: n.originalColor, border: themeBorder }
                    : { background: dimBg, border: dimBorder },
                font: { color: isActive ? themeText : dimText }
            });
        }
    });

    graphData.edges.forEach(e => {
        const isEdgeActive = currentActiveNodes.has(e.from) && currentActiveNodes.has(e.to);
        if (mode === 'hide') {
            edgeUpdates.push({
                id: e.id,
                hidden: !isEdgeActive,
                color: { color: getCSSVar('--edge-default') },
                width: optEdgeWidth
            });
        } else {
            edgeUpdates.push({
                id: e.id,
                hidden: false,
                color: isEdgeActive
                    ? { color: getCSSVar('--edge-default'), opacity: 1.0 }
                    : { color: getCSSVar('--edge-faded'), opacity: 0.15 },
                width: isEdgeActive ? optEdgeWidth : Math.max(0.1, optEdgeWidth * 0.5)
            });
        }
    });

    logDebug("applyGraphFilters updates", {
        force,
        activePaletteName,
        tempEdgeIdsLength: tempEdgeIds.length,
        nodeSample: nodeUpdates.slice(0, 8).map(n => ({
            id: n.id,
            color: n.color
        }))
    });
    nodesDataset.update(nodeUpdates);
    edgesDataset.update(edgeUpdates);
}
