function init() {
    logDebug("teste");
    logDebug(getCSSVar('--dim-bg'));
    logDebug(getCSSVar('--bg-mol'));
    logDebug("Mode Detected:", masterData.mode);
    updateThemeVars();
    initSplitter();        
    initMolColors();
    initUploads();
    initAdvancedOptions(); 

    document.querySelectorAll('input[name="viewMode"]').forEach(r => r.addEventListener('change', applyGraphFilters));
    document.getElementById('showLabels').addEventListener('change', update3DViewerOrGrid);
    document.getElementById('showStdEdges').addEventListener('change', function() {
        if (!this.checked) clearDynamicStdEdges();
        else if (pinnedNode && currentGraphMode === 'associated') drawDynamicStdEdges(pinnedNode);
    });

    const togglePhysics = (e) => { 
        if(network) network.setOptions({ physics: !e.target.checked }); 
        gridNetworks.forEach(n => n.network.setOptions({ physics: !e.target.checked }));
    };
    document.getElementById('freezePhysics').addEventListener('change', togglePhysics);
    document.getElementById('freezePhysicsFiltered').addEventListener('change', togglePhysics);
    
    if (masterData.mode === 'pairwise') {
        document.getElementById('pair-mode-panel').style.display = 'block';
        let sel = document.getElementById('pair-selector');
        Object.keys(masterData.pairs).forEach(k => { sel.innerHTML += `<option value="${k}">🔍 Focus: ${k.replace('_vs_', ' vs ')}</option>`; });
        
        assignAllPairGraphColors(); 
        autoLoadStructures().then(() => { switchPairView('GRID'); });
    } else {
        document.getElementById('pair-mode-panel').style.display = 'none';
        graphData = masterData;
        assignGraphColorsToData();
        initMetadata();
        initHierarchy();
        
        setTimeout(() => {
            const sn = document.getElementById('single-network');
            if(sn) sn.style.display = 'block';
            const sm = document.getElementById('single-mols');
            if(sm) sm.style.display = 'block';
            const vc = document.getElementById('view-controls');
            if(vc) vc.style.display = 'flex';
            const tb = document.getElementById('viewer-toolbar');
            if(tb) tb.style.display = 'flex';
            
            initNetwork();
            handleTreeChange();
            autoLoadStructures();
            switchGraphViewMode();
        }, 100);
    }
}

function initAdvancedOptions() {
    document.getElementById('optPalette').addEventListener('change', function(e) {
        activePaletteName = e.target.value;
        if (masterData.mode === 'pairwise' && document.getElementById('pair-selector').value === 'GRID') {
            assignAllPairGraphColors(); Object.keys(masterData.pairs).forEach(k => applyGraphFiltersGrid(k)); 
        } else {
            assignGraphColorsToData(); initMetadata(); 
            if (currentGraphMode === 'associated') handleTreeChange(); else handleFilteredChange();
        }
    });

    document.getElementById('optNodeSize').addEventListener('input', function(e) {
        optNodeSize = parseInt(e.target.value); document.getElementById('valNodeSize').innerText = optNodeSize;
        if(network) network.setOptions({ nodes: { size: optNodeSize } });
        gridNetworks.forEach(n => n.network.setOptions({ nodes: { size: optNodeSize } }));
    });
    
    document.getElementById('optLabelSize').addEventListener('input', function(e) {
        optLabelSize = parseInt(e.target.value); document.getElementById('valLabelSize').innerText = optLabelSize;
        if(network) network.setOptions({ nodes: { font: { size: optLabelSize } } });
        gridNetworks.forEach(n => n.network.setOptions({ nodes: { font: { size: optLabelSize-4 } } }));
    });

    document.getElementById('optEdgeWidth').addEventListener('input', function(e) {
        optEdgeWidth = parseFloat(e.target.value); document.getElementById('valEdgeWidth').innerText = optEdgeWidth.toFixed(1);
        if (masterData.mode === 'pairwise' && document.getElementById('pair-selector').value === 'GRID') {
            Object.keys(masterData.pairs).forEach(k => applyGraphFiltersGrid(k));
        } else {
            if (currentGraphMode === 'associated') applyGraphFilters(); else handleFilteredChange();
        }
    });

    document.getElementById('optSpringLen').addEventListener('input', function(e) {
        optSpringLen = parseInt(e.target.value); document.getElementById('valSpringLen').innerText = optSpringLen;
        if(network) network.setOptions({ physics: { barnesHut: { springLength: optSpringLen } } });
        gridNetworks.forEach(n => n.network.setOptions({ physics: { barnesHut: { springLength: optSpringLen } } }));
    });

    document.getElementById('optStraightEdges').addEventListener('change', function(e) {
        if(network) network.setOptions({ edges: { smooth: !e.target.checked } });
        gridNetworks.forEach(n => n.network.setOptions({ edges: { smooth: !e.target.checked } }));
    });
}


function initSplitter() {
    let isDragging = false;
    const splitter = document.getElementById('v-splitter');
    const topPanel = document.getElementById('top-panel');
    const bottomPanel = document.getElementById('bottom-panel');

    splitter.addEventListener('mousedown', function(e) { isDragging = true; document.body.style.cursor = 'row-resize'; e.preventDefault(); });
    document.addEventListener('mousemove', function(e) {
        if (!isDragging) return;
        const containerRect = document.getElementById('main-content').getBoundingClientRect();
        let newTopHeight = e.clientY - containerRect.top;
        if (newTopHeight < containerRect.height * 0.1) newTopHeight = containerRect.height * 0.1;
        if (newTopHeight > containerRect.height * 0.9) newTopHeight = containerRect.height * 0.9;
        topPanel.style.flex = `0 0 ${newTopHeight}px`; bottomPanel.style.flex = `1 1 0%`;
    });
    document.addEventListener('mouseup', function(e) { if (isDragging) { isDragging = false; document.body.style.cursor = 'default'; } });
}


function initMetadataGlobalFallback() {
    try {
        if (masterData.mode !== 'pairwise') return;
        const c = masterData.metadata || {}; 
        
        // Highlight the reference structure if in screening mode
        const refHtml = (masterData.actual_mode === 'screening' && masterData.reference_structure) 
            ? `<div style="margin-bottom: 10px; padding: 6px 10px; background: rgba(37, 99, 235, 0.1); border-left: 3px solid var(--btn-bg); border-radius: 4px; color: var(--text-main);"><b>Reference Structure:</b> ${masterData.reference_structure}</div>` 
            : '';

        document.getElementById('metadata-panel').innerHTML = `
            <h3>Execution Metadata</h3>
            <div style="margin-bottom: 10px;"><b>Run Name:</b> ${masterData.run_name || 'N/A'}</div>
            ${refHtml}
            <div style="margin-bottom: 12px; font-size: 12px; font-style: italic;">Select a specific pair focus or view colors in the grid to see node mappings.</div>
            <div><b>Parameters:</b>
                <ul style="margin: 5px 0 0 0; padding-left: 20px; line-height: 1.5;">
                    <li><b>Mode:</b> ${c.run_mode || 'N/A'}</li>
                    <li><b>Granularity:</b> ${c.node_granularity || 'N/A'}</li>
                    <li><b>Edge Thresh:</b> ${c.edge_threshold || 'N/A'} Å</li>
                    <li><b>Global Diff:</b> ${c.global_distance_diff_threshold || 'N/A'} Å</li>
                </ul>
            </div>
        `;
    } catch(e) { logError("initMetadataGlobalFallback failed", e); }
}

function initMetadata() {
    try {
        if (!graphData) return;
        const c = graphData.metadata || masterData.metadata || {};
        const pal = GRAPH_PALETTES[activePaletteName];
        
        let p_html = graphData.proteins.map((p, localIdx) => {
            let globalIdx = masterData.proteins.indexOf(p);
            const color = pal[globalIdx % pal.length];
            
            // Tag the reference structure in the tuple list
            const isRef = (masterData.actual_mode === 'screening' && p === masterData.reference_structure);
            const refIcon = isRef ? ' <span style="color: var(--text-muted); font-size: 11px;"><i>(Ref)</i></span>' : '';
            
            return `<div style="margin-bottom: 3px;"><span class="color-dot" style="background-color: ${color};"></span><b>Prot ${localIdx}:</b> ${p}${refIcon}</div>`;
        }).join('');

        // Highlight the reference structure at the top
        const refHtml = (masterData.actual_mode === 'screening' && masterData.reference_structure) 
            ? `<div style="margin-bottom: 10px; padding: 6px 10px; background: rgba(37, 99, 235, 0.1); border-left: 3px solid var(--btn-bg); border-radius: 4px; color: var(--text-main);"><b>Reference Structure:</b> ${masterData.reference_structure}</div>` 
            : '';

        document.getElementById('metadata-panel').innerHTML = `
            <h3>Execution Metadata</h3>
            <div style="margin-bottom: 10px;"><b>Run Name:</b> ${graphData.run_name || masterData.run_name || 'N/A'}</div>
            ${refHtml}
            <div style="margin-bottom: 12px;"><b>Protein Order (Tuple Layout):</b><div style="margin-top: 5px;">${p_html}</div></div>
            <div>
                <b>Parameters:</b>
                <ul style="margin: 5px 0 0 0; padding-left: 20px; line-height: 1.5;">
                    <li><b>Mode:</b> ${c.run_mode || 'N/A'}</li>
                    <li><b>Granularity:</b> ${c.node_granularity || 'N/A'}</li>
                    <li><b>Edge Thresh:</b> ${c.edge_threshold || 'N/A'} Å</li>
                    <li><b>Global Diff:</b> ${c.global_distance_diff_threshold || 'N/A'} Å</li>
                </ul>
            </div>
        `;
    } catch(e) { logError("initMetadata crashed", e); }
}

function initNetwork() {
    logDebug("Initializing Single Network");
    try {
        const container = document.getElementById('single-network');
        if (!container) return;
        if (network) { network.destroy(); network = null; }
        
        network = new vis.Network(container, { nodes: nodesDataset, edges: edgesDataset }, {
            nodes: { shape: 'dot', size: optNodeSize, font: { size: optLabelSize, multi: 'html', color: themeText, face: 'ui-sans-serif' }, borderWidth: 2 },
            edges: { smooth: !document.getElementById('optStraightEdges').checked },
            physics: { solver: 'repulsion', repulsion: { nodeDistance: 100, springLength: 100 }, stabilization: { iterations: 20 } },
            interaction: { hover: true, tooltipDelay: 150, dragNodes: true }
        });
        
        network.once("stabilizationIterationsDone", function () { if(document.getElementById('freezePhysics').checked) network.setOptions({ physics: false }); network.fit(); });
        network.on("hoverNode", function (params) { if (currentGraphMode !== 'associated' || pinnedNode !== null || !document.getElementById('showStdEdges').checked) return; drawDynamicStdEdges(params.node); });
        network.on("blurNode", function () { if (currentGraphMode !== 'associated' || pinnedNode !== null) return; clearDynamicStdEdges(); });
        network.on("click", function (params) {
            if (currentGraphMode !== 'associated') return;
            const clickedNode = (params.nodes && params.nodes.length) ? params.nodes[0] : null;
            if (!clickedNode) { pinnedNode = null; clearDynamicStdEdges(); return; }
            if (pinnedNode === clickedNode) { pinnedNode = null; clearDynamicStdEdges(); } 
            else { pinnedNode = clickedNode; if (document.getElementById('showStdEdges').checked) drawDynamicStdEdges(clickedNode); }
        });
    } catch(e) { logError("initNetwork crashed", e); }
}

function initUploads() {
    logDebug("Initializing Manual Upload Inputs globally.");
    const container = document.getElementById('upload-inputs');
    if (!container) return;
    container.innerHTML = '';
    masterData.proteins.forEach((protName, idx) => {
        let firstChainColor = '#999999';
        for (let key in molChainColors) { if (key.startsWith(`${idx}_`)) { firstChainColor = molChainColors[key]; break; } }
        let div = document.createElement('div');
        div.style.marginBottom = '8px'; div.style.fontSize = '12px';
        div.innerHTML = `<span class="color-dot" style="background-color: ${firstChainColor};"></span><b>${protName}:</b><br>
            <input type="file" id="file-prot-${idx}" accept=".pdb,.cif" style="width:100%; margin-top: 4px; color: var(--text-main);">`;
        container.appendChild(div);
    });
}

function initHierarchy() {
    const tree = document.getElementById('hierarchy-tree'); let html = '';
    if(!graphData || !graphData.components) return;
    graphData.components.forEach(comp => {
        if (comp.id === 0) {
            html += `<div class="tree-item tree-level-1"><span class="collapse-toggle" style="visibility:hidden;"></span><label><input type="checkbox" class="tree-cb" value="global" checked onchange="handleTreeChange()"> 🌐 Global Associated Graph</label></div>`;
        } else {
            const hasFrames = comp.frames.some(f => f.id !== 0);
            const toggleBtn = hasFrames ? `<span class="collapse-toggle" onclick="toggleComp('comp_${comp.id}_frames', this)">▼</span>` : `<span class="collapse-toggle" style="visibility:hidden;"></span>`;
            html += `<div class="tree-item tree-level-2">${toggleBtn}<label><input type="checkbox" class="tree-cb" value="comp_${comp.id}" onchange="handleTreeChange()"> 📦 Component ${comp.id}</label></div>`;
            if (hasFrames) {
                html += `<div id="comp_${comp.id}_frames" style="display:block;">`;
                comp.frames.forEach(frame => { if (frame.id !== 0) html += `<div class="tree-item tree-level-3"><span class="collapse-toggle" style="visibility:hidden;"></span><label><input type="checkbox" class="tree-cb" value="frame_${comp.id}_${frame.id}" onchange="handleTreeChange()"> ↳ Frame ${frame.id}</label></div>`; });
                html += `</div>`;
            }
        }
    });
    tree.innerHTML = html;

    const fTree = document.getElementById('filtered-tree'); let fHtml = '';
    graphData.filtered_graphs.forEach((fg, idx) => {
        fHtml += `<div class="tree-item tree-level-1"><span class="collapse-toggle" style="visibility:hidden;"></span><label><input type="radio" name="filtered_prot" value="${idx}" ${idx===0 ? 'checked' : ''} onchange="handleFilteredChange()"> Prot ${idx}: ${fg.name}</label></div>`;
    });
    fTree.innerHTML = fHtml;
}

function initHierarchyGrid() {
    const tree = document.getElementById('hierarchy-tree'); let html = '';
    if (!masterData.pairs) return;
    Object.keys(masterData.pairs).forEach(pairKey => {
        let pData = masterData.pairs[pairKey]; let pkSafe = pairKey.replace(/[^a-zA-Z0-9]/g, '_');
        html += `<div class="tree-item tree-level-1" style="background: var(--bg-control); border: 1px solid var(--border); margin-top: 10px;"><span class="collapse-toggle" onclick="toggleComp('grid_${pkSafe}', this)">▼</span><label style="color: var(--btn-bg);"><b>${pairKey.replace('_vs_', ' vs ')}</b></label></div>`;
        html += `<div id="grid_${pkSafe}" style="display:block; padding-left: 10px; border-left: 2px solid var(--border-light); margin-left: 5px;">`;
        
        pData.components.forEach(comp => {
            if (comp.id === 0) html += `<div class="tree-item tree-level-1"><span class="collapse-toggle" style="visibility:hidden;"></span><label><input type="checkbox" class="tree-cb-grid" data-pair="${pairKey}" value="global" checked onchange="handleGridTreeChange('${pairKey}')"> 🌐 Global Associated Graph</label></div>`;
            else {
                const hasFrames = comp.frames.some(f => f.id !== 0);
                const toggleBtn = hasFrames ? `<span class="collapse-toggle" onclick="toggleComp('grid_${pkSafe}_comp_${comp.id}_frames', this)">▼</span>` : `<span class="collapse-toggle" style="visibility:hidden;"></span>`;
                html += `<div class="tree-item tree-level-2">${toggleBtn}<label><input type="checkbox" class="tree-cb-grid" data-pair="${pairKey}" value="comp_${comp.id}" onchange="handleGridTreeChange('${pairKey}')"> 📦 Component ${comp.id}</label></div>`;
                if (hasFrames) {
                    html += `<div id="grid_${pkSafe}_comp_${comp.id}_frames" style="display:block;">`;
                    comp.frames.forEach(frame => { if (frame.id !== 0) html += `<div class="tree-item tree-level-3"><span class="collapse-toggle" style="visibility:hidden;"></span><label><input type="checkbox" class="tree-cb-grid" data-pair="${pairKey}" value="frame_${comp.id}_${frame.id}" onchange="handleGridTreeChange('${pairKey}')"> ↳ Frame ${frame.id}</label></div>`; });
                    html += `</div>`;
                }
            }
        });
        html += `</div>`;
    });
    tree.innerHTML = html;
    document.getElementById('filtered-tree').style.display = 'none';
    document.getElementById('edge-interactions-panel').style.display = 'block';
}
