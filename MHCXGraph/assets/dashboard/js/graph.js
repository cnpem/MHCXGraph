
function getStdColor(std, threshold) {
    if (std === null || std === undefined || isNaN(std)) return '#d1d5db'; 
    let ratio = std / (threshold + 0.0001);
    if (ratio > 1) ratio = 1; if (ratio < 0) ratio = 0;
    return `hsl(${120 * (1 - ratio)}, 100%, 45%)`;
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
    const dimBg = isDark ? 'rgba(55,65,81,0.3)' : 'rgba(229,231,235,0.3)';
    const dimBorder = isDark ? 'rgba(75,85,99,0.3)' : 'rgba(209,213,219,0.3)';

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
            font: { align: "top", size: optLabelSize-2, color: isDark ? '#ffffff' : '#000000', strokeWidth: 4, strokeColor: isDark ? '#111827' : '#ffffff', face: "sans-serif" },
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
        const pal = GRAPH_PALETTES[activePaletteName];
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

        const pal = GRAPH_PALETTES[activePaletteName];
        let uniqueGroups = [...new Set(graphData.nodes.map(n => n.group))].sort();

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
        return { ...e, color: { color: '#9ca3af', opacity: 1.0 }, width: Math.max(0.1, edgeWidth) };
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
        dsEdges.add(pData.edges.map(e => ({...e, color: {color: '#9ca3af', opacity: 0.8}, width: optEdgeWidth})));
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
    const dimBg = isDark ? 'rgba(55,65,81,0.3)' : 'rgba(229,231,235,0.3)';
    const dimBorder = isDark ? 'rgba(75,85,99,0.3)' : 'rgba(209,213,219,0.3)';
    const dimText = isDark ? 'rgba(156,163,175,0.5)' : 'rgba(107,114,128,0.5)';

    const nodeUpdates = [], edgeUpdates = [];
    pData.nodes.forEach(n => {
        const isActive = activeIds.has(n.id);
        if (mode === 'hide') nodeUpdates.push({ id: n.id, hidden: !isActive, color: { background: n.originalColor, border: themeBorder }, font: { color: themeText } });
        else nodeUpdates.push({ id: n.id, hidden: false, color: isActive ? { background: n.originalColor, border: themeBorder } : { background: dimBg, border: dimBorder }, font: { color: isActive ? themeText : dimText } });
    });

    pData.edges.forEach(e => {
        const isEdgeActive = activeIds.has(e.from) && activeIds.has(e.to);
        if (mode === 'hide') edgeUpdates.push({ id: e.id, hidden: !isEdgeActive, color: {color: '#9ca3af'}, width: optEdgeWidth });
        else edgeUpdates.push({ id: e.id, hidden: false, color: isEdgeActive ? { color: '#9ca3af', opacity: 1.0 } : { color: '#6b7280', opacity: 0.15 }, width: isEdgeActive ? optEdgeWidth : Math.max(0.1, optEdgeWidth * 0.5) });
    });

    dsNodes.update(nodeUpdates); dsEdges.update(edgeUpdates);
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
    const pal = GRAPH_PALETTES[activePaletteName];
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
    const dimBg = isDark ? 'rgba(55,65,81,0.3)' : 'rgba(229,231,235,0.3)';
    const dimBorder = isDark ? 'rgba(75,85,99,0.3)' : 'rgba(209,213,219,0.3)';
    const dimText = isDark ? 'rgba(156,163,175,0.5)' : 'rgba(107,114,128,0.5)';

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
                color: { color: '#9ca3af' },
                width: optEdgeWidth
            });
        } else {
            edgeUpdates.push({
                id: e.id,
                hidden: false,
                color: isEdgeActive
                    ? { color: '#9ca3af', opacity: 1.0 }
                    : { color: '#6b7280', opacity: 0.15 },
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

    nodesDataset.update(
        nodeUpdates.map(n => ({
            id: n.id,
            color: {
                background: n.color.background,
                border: n.color.border
            }
        }))
    )
    edgesDataset.update(edgeUpdates);
}
