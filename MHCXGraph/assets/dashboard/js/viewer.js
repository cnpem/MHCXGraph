function focusOnPair(pairKey) {
    logDebug(`Direct focus button clicked for: ${pairKey}`);
    document.getElementById('pair-selector').value = pairKey;
    switchPairView(pairKey);
}

function triggerRebuild() {
    if (masterData.mode === 'pairwise' && document.getElementById('pair-selector').value === 'GRID') buildGridBottom();
    else rebuildViewer();
}

function switchGraphViewMode() {
    currentGraphMode = document.getElementById('main-graph-selector').value;
    const singleNet = document.getElementById('single-network');
    const analysisWrapper = document.getElementById('analysis-wrapper');
    const bottomPanel = document.getElementById('bottom-panel');
    const vSplitter = document.getElementById('v-splitter');
    
    if (currentGraphMode === 'associated') {
        if(singleNet) singleNet.style.display = 'block'; 
        if(bottomPanel) bottomPanel.style.display = 'flex';
        if(vSplitter) vSplitter.style.display = 'flex';
        if(analysisWrapper) analysisWrapper.style.display = 'none';
        document.getElementById('hierarchy-tree').style.display = 'block'; document.getElementById('filtered-tree').style.display = 'none';
        document.getElementById('edge-interactions-panel').style.display = 'block'; document.getElementById('filtered-options-panel').style.display = 'none';
        document.getElementById('nav-title').parentElement.style.display = 'block';
        document.getElementById('nav-title').innerText = "Navigation (Multi-Select)";
        handleTreeChange();
    } else if (currentGraphMode === 'filtered') {
        if(singleNet) singleNet.style.display = 'block'; 
        if(bottomPanel) bottomPanel.style.display = 'flex';
        if(vSplitter) vSplitter.style.display = 'flex';
        if(analysisWrapper) analysisWrapper.style.display = 'none';
        document.getElementById('hierarchy-tree').style.display = 'none'; document.getElementById('filtered-tree').style.display = 'block';
        document.getElementById('edge-interactions-panel').style.display = 'none'; document.getElementById('filtered-options-panel').style.display = 'block';
        document.getElementById('nav-title').parentElement.style.display = 'block';
        document.getElementById('nav-title').innerText = "Select Protein";
        handleFilteredChange();
    } else if (currentGraphMode === 'analysis') {
        if(singleNet) singleNet.style.display = 'none'; 
        if(bottomPanel) bottomPanel.style.display = 'none'; // HIDES 3D VIEWER
        if(vSplitter) vSplitter.style.display = 'none';     // HIDES SPLITTER
        if(analysisWrapper) analysisWrapper.style.display = 'block';
        
        document.getElementById('nav-title').parentElement.style.display = 'none';
        document.getElementById('edge-interactions-panel').style.display = 'none'; document.getElementById('filtered-options-panel').style.display = 'none';
        
        let isGlobal = (masterData.mode === 'pairwise' && document.getElementById('pair-selector') && document.getElementById('pair-selector').value === 'GRID');
        renderAnalysis(isGlobal); 
    }
}

function switchPairView(val) {
    try {
        if (!val || val instanceof Event) val = document.getElementById('pair-selector').value;
        logDebug(`Routing View To: ${val}`);

        document.getElementById('analysis-wrapper').style.display = 'none';
        document.getElementById('top-panel').style.display = 'none';
        document.getElementById('bottom-panel').style.display = 'none';
        document.getElementById('v-splitter').style.display = 'none';
        document.getElementById('view-controls').style.display = 'none';
        
        const singleNet = document.getElementById('single-network');
        if (singleNet) singleNet.style.display = 'block';
        const singleMols = document.getElementById('single-mols');
        if (singleMols) singleMols.style.display = 'block';

        const dynGridNet = document.getElementById('dynamic-grid-networks');
        if (dynGridNet) dynGridNet.style.display = 'none';
        const dynGridMols = document.getElementById('dynamic-grid-mols');
        if (dynGridMols) dynGridMols.style.display = 'none';
        
        if (network) { network.destroy(); network = null; }
        gridNetworks.forEach(n => n.network.destroy());
        gridNetworks = []; gridNodesDatasets = {}; gridEdgesDatasets = {}; viewers = [];
        nodesDataset.clear(); edgesDataset.clear();
        
        if (val === 'GRID') {
            document.getElementById('main-graph-selector').value = 'associated';
            document.getElementById('nav-title').parentElement.style.display = 'block';
            document.getElementById('top-panel').style.display = 'block';
            document.getElementById('bottom-panel').style.display = 'flex';
            document.getElementById('v-splitter').style.display = 'flex';
            
            initMetadataGlobalFallback();
            initHierarchyGrid(); 
            
            setTimeout(() => { buildGridTop(); buildGridBottom(); }, 50);
            
        } else if (val === 'ANALYSIS') {
            document.getElementById('main-graph-selector').value = 'associated';
            document.getElementById('nav-title').parentElement.style.display = 'none';
            document.getElementById('analysis-wrapper').style.display = 'block';
            renderAnalysis(true); 
            
        } else {
            document.getElementById('main-graph-selector').value = 'associated';
            document.getElementById('nav-title').parentElement.style.display = 'block';
            document.getElementById('top-panel').style.display = 'block';
            document.getElementById('bottom-panel').style.display = 'flex';
            document.getElementById('v-splitter').style.display = 'flex';
            document.getElementById('view-controls').style.display = 'flex';
            
            graphData = masterData.pairs[val];
            
            assignGraphColorsToData();
            initMetadata();
            initHierarchy(); 
            
            const select = document.getElementById('view-selector');
            if (select) {
                select.innerHTML = '<option value="aligned">Aligned (Single View)</option><option value="separate">Separate (Grid View)</option>';
                loadedModels.forEach((m, idx) => {
                    if (m && graphData.proteins.includes(m.name)) select.innerHTML += `<option value="prot_${idx}">Prot ${graphData.proteins.indexOf(m.name)}: ${m.name}</option>`;
                });
            }
            
            setTimeout(() => {
                document.getElementById('viewer-toolbar').style.display = 'flex';
                initNetwork(); 
                handleTreeChange();
                rebuildViewer();
                switchGraphViewMode(); 
            }, 50);
        }
    } catch(e) { logError("switchPairView crashed", e); }
}


function hoverCallback(atom, viewer, event) {
    if (!atom.resn) return;
    const tt = document.getElementById("mol-tooltip");
    let pName = viewer.protName || "Unknown";
    if (viewer.layoutMode === 'aligned') {
        const loaded = loadedModels.find(lm => lm && lm.internalModel === atom.model);
        if (loaded) pName = loaded.name;
    }
    tt.innerHTML = `<b>Protein:</b> ${pName}<br><b>Residue:</b> ${atom.resn} ${atom.resi} (Chain ${atom.chain})`;
    tt.style.display = "block"; tt.style.left = (event.clientX + 15) + "px"; tt.style.top = (event.clientY + 15) + "px";
}

function unhoverCallback() { document.getElementById("mol-tooltip").style.display = "none"; }
    function rebuildViewer() {
        try {
            const container = document.getElementById('single-mols');
            container.innerHTML = ''; viewers = []; window.activeMolViewer = null; 
            
            const layoutMode = document.getElementById('view-selector').value;
            const isDark = document.body.getAttribute('data-theme') === 'dark';
            const bgMol = isDark ? "#1f2937" : "#ffffff";
            
            const validModels = loadedModels.filter(m => m !== null && graphData.proteins.includes(m.name));
            if (validModels.length === 0) return;

            if (layoutMode === 'separate') {
                container.style.display = 'flex';
                container.style.flexDirection = 'row';

                validModels.forEach((m, i) => {
                    const subDiv = document.createElement('div'); subDiv.className = 'mol-subviewer';
                    subDiv.style.flex = '1';
                    if (i === 0) subDiv.style.borderRight = '1px solid var(--border-light)';
                    
                    const titleDiv = document.createElement('div'); titleDiv.className = 'viewer-title';
                    const baseCol = get3DColor('A', m.protIdx);
                    titleDiv.innerHTML = `<span class="color-dot" style="background-color: ${baseCol};"></span>Prot ${graphData.proteins.indexOf(m.name)}: ${m.name}`;
                    subDiv.appendChild(titleDiv);
                    
                    const molDiv = document.createElement('div'); molDiv.style.width = '100%'; molDiv.style.height = '100%';
                    subDiv.appendChild(molDiv); container.appendChild(subDiv);
                    
                    const initViewer = () => {
                        if (molDiv.clientWidth === 0 || molDiv.clientHeight === 0) {
                            setTimeout(initViewer, 20); return;
                        }
                        let v = $3Dmol.createViewer(molDiv, { backgroundColor: bgMol });
                        let mObj = v.addModel(m.text, m.format); m.internalModel = mObj.getID();
                        v.setStyle({}, { cartoon: { colorfunc: function(atom) { return get3DColor(atom.chain, m.protIdx); } } });
                        
                        molDiv.addEventListener('mousedown', () => { window.activeMolViewer = v; });
                        molDiv.addEventListener('wheel', () => { window.activeMolViewer = v; }, {passive: true});
                        
                        v.protIdx = m.protIdx; v.protName = m.name; v.layoutMode = 'separate';
                        v.setHoverable({}, true, hoverCallback, unhoverCallback); v.zoomTo(); viewers.push(v);
                        update3DViewer();
                    };
                    initViewer();
                });
            } else if (layoutMode.startsWith('prot_')) {
                container.style.display = 'block';
                const targetIdx = parseInt(layoutMode.split('_')[1]);
                const targetModel = validModels.find(m => m.protIdx === targetIdx);
                
                if (targetModel) {
                    let molDiv = document.createElement('div'); molDiv.style.width = '100%'; molDiv.style.height = '100%';
                    container.appendChild(molDiv);

                    const initViewer = () => {
                        if (molDiv.clientWidth === 0 || molDiv.clientHeight === 0) {
                            setTimeout(initViewer, 20); return;
                        }
                        let v = $3Dmol.createViewer(molDiv, { backgroundColor: bgMol });
                        let mObj = v.addModel(targetModel.text, targetModel.format); targetModel.internalModel = mObj.getID();
                        v.setStyle({}, { cartoon: { colorfunc: function(atom) { return get3DColor(atom.chain, targetModel.protIdx); } } });
                        v.protIdx = targetModel.protIdx; v.protName = targetModel.name; v.layoutMode = 'single';
                        v.setHoverable({}, true, hoverCallback, unhoverCallback); v.zoomTo(); viewers.push(v);
                        update3DViewer();
                    };
                    initViewer();
                }
            } else {
                container.style.display = 'block';
                let molDiv = document.createElement('div'); molDiv.style.width = '100%'; molDiv.style.height = '100%';
                container.appendChild(molDiv);
                
                const initViewer = () => {
                    if (molDiv.clientWidth === 0 || molDiv.clientHeight === 0) {
                        setTimeout(initViewer, 20); return;
                    }
                    let v = $3Dmol.createViewer(molDiv, { backgroundColor: bgMol });
                    validModels.forEach(m => {
                        let mObj = v.addModel(m.text, m.format); m.internalModel = mObj.getID();
                        v.setStyle({model: m.internalModel}, { cartoon: { colorfunc: function(atom) { return get3DColor(atom.chain, m.protIdx); } } });
                    });
                    v.layoutMode = 'aligned'; v.setHoverable({}, true, hoverCallback, unhoverCallback); v.zoomTo(); viewers.push(v);
                    update3DViewer();
                };
                initViewer();
            }
        } catch(e) { logError("rebuildViewer crashed", e); }
    }

function update3DViewerGridForPair(pairKey) {
    if (viewers.length === 0 || loadedModels.length === 0) return;
    const showLabels = document.getElementById('showLabels').checked;
    const isDark = document.body.getAttribute('data-theme') === 'dark';
    const labelColor = window.exportLabelColor || (isDark ? "white" : "black");
    
    viewers.forEach(v => {
        if (v.pairKey !== pairKey) return;
        v.removeAllLabels();
        v.protIdxs.forEach(idx => {
            if(loadedModels[idx]) v.setStyle({}, { cartoon: { colorfunc: function(atom) { return get3DColor(atom.chain, idx); }, opacity: 0.5 } });
        });

        if(!window.gridActiveNodes || !window.gridActiveNodes[pairKey]) { v.render(); return; }
        
        let activeIds = window.gridActiveNodes[pairKey];
        let pData = masterData.pairs[pairKey];
        
        pData.nodes.forEach(n => {
            if (activeIds.has(n.id)) {
                n.mapping.forEach(m => {
                    let sel = { chain: m.chain, resi: parseInt(m.resi) }; 
                    const drawColor = get3DColor(m.chain, m.model_idx);
                    v.setStyle(sel, { cartoon: { color: drawColor }, stick: { colorscheme: 'whiteCarbon', radius: 0.25 } }, true);
                    if (showLabels) v.addLabel(`${m.resn}:${m.resi}`, { font: "sans-serif", fontSize: 13, fontColor: labelColor, showBackground: false, alignment: "center", position: sel }, sel, true);
                });
            }
        });
        v.render();
    });
}

function update3DViewer() {
    if (viewers.length === 0 || loadedModels.length === 0) return;
    const layoutMode = document.getElementById('view-selector').value;
    const showLabels = document.getElementById('showLabels').checked;
    const isDark = document.body.getAttribute('data-theme') === 'dark';
    const labelColor = window.exportLabelColor || (isDark ? "white" : "black");
    
    viewers.forEach(v => {
        v.removeAllLabels();
        if (layoutMode === 'aligned') {
            loadedModels.forEach(m => { if (m && graphData.proteins.includes(m.name)) v.setStyle({model: m.internalModel}, { cartoon: { colorfunc: function(atom) { return get3DColor(atom.chain, m.protIdx); }, opacity: 0.5 } }); });
        } else if (layoutMode === 'separate' || layoutMode.startsWith('prot_')) {
            v.setStyle({}, { cartoon: { colorfunc: function(atom) { return get3DColor(atom.chain, v.protIdx); }, opacity: 0.5 } });
        }
    });
    
    let nodesToIterate = currentGraphMode === 'associated' ? graphData.nodes : graphData.filtered_graphs[activeFilteredProtIdx].nodes;
    
    nodesToIterate.forEach(n => {
        if (currentActiveNodes.has(n.id)) {
            n.mapping.forEach(m => {
                const globalIdx = m.model_idx;
                const loaded = loadedModels.find(lm => lm && lm.protIdx === globalIdx);
                if (!loaded) return;

                let v = null, targetModel = 0;
                if (layoutMode === 'aligned') { v = viewers[0]; targetModel = loaded.internalModel; } 
                else if (layoutMode === 'separate') { v = viewers.find(vo => vo.protIdx === globalIdx); targetModel = 0; } 
                else if (layoutMode.startsWith('prot_')) { if (globalIdx !== parseInt(layoutMode.split('_')[1])) return; v = viewers[0]; targetModel = 0; }

                if (!v) return;

                let sel = { model: targetModel, chain: m.chain, resi: parseInt(m.resi) };
                const drawColor = get3DColor(m.chain, globalIdx);
                v.setStyle(sel, { cartoon: { color: drawColor }, stick: { colorscheme: 'whiteCarbon', radius: 0.25 } }, true);
                
                if (showLabels) v.addLabel(`${m.resn}:${m.resi}`, { font: "sans-serif", fontSize: 13, fontColor: labelColor, showBackground: false, alignment: "center", position: sel }, sel, true);
            });
        }
    });
    viewers.forEach(v => v.render());
}

function initMolColors(forceReset = false) {
    if (Object.keys(molChainColors).length > 0 && !forceReset) return;
    const isDark = document.body.getAttribute('data-theme') === 'dark';
    const pal = isDark ? MOL_DEFAULTS.dark : MOL_DEFAULTS.light;
    let i = 0;
    masterData.proteins.forEach((protName, pIdx) => {
        let chains = new Set();
        if (masterData.mode === 'pairwise') Object.values(masterData.pairs).forEach(pd => pd.nodes.forEach(n => n.mapping.forEach(m => { if (m.model_idx === pIdx) chains.add(m.chain); })));
        else masterData.nodes.forEach(n => n.mapping.forEach(m => { if (m.model_idx === pIdx) chains.add(m.chain); }));
        Array.from(chains).sort().forEach(c => { molChainColors[`${pIdx}_${c}`] = pal[i % pal.length]; i++; });
    });
    renderMolColorTree();
}

function renderMolColorTree() {
    const container = document.getElementById('color-tree');
    let html = '';
    masterData.proteins.forEach((protName, pIdx) => {
        if (masterData.mode === 'pairwise' && document.getElementById('pair-selector').value !== 'GRID' && document.getElementById('pair-selector').value !== 'ANALYSIS' && graphData && !graphData.proteins.includes(protName)) return;
        html += `<div class="tree-item tree-level-1" style="justify-content: space-between;"><div style="display: flex; align-items: center; gap: 6px;"><span class="collapse-toggle" onclick="toggleComp('colors_prot_${pIdx}', this)">▼</span><label>Prot ${pIdx}: ${protName}</label></div></div><div id="colors_prot_${pIdx}" style="display:block; padding-left: 10px;">`;
        let chains = new Set();
        if (masterData.mode === 'pairwise') Object.values(masterData.pairs).forEach(pd => pd.nodes.forEach(n => n.mapping.forEach(m => { if (m.model_idx === pIdx) chains.add(m.chain); })));
        else masterData.nodes.forEach(n => n.mapping.forEach(m => { if(m.model_idx === pIdx) chains.add(m.chain); }));
        Array.from(chains).sort().forEach(c => {
            const key = `${pIdx}_${c}`;
            html += `<div class="tree-item tree-level-2" style="justify-content: space-between; padding-right: 15px;"><label>Chain ${c}</label><input type="color" value="${molChainColors[key]}" onchange="updateMolChainColor(${pIdx}, '${c}', this.value)" style="width: 25px; height: 25px; padding: 0; border: none; background: none; cursor: pointer; border-radius: 4px;"></div>`;
        });
        html += `</div>`;
    });
    container.innerHTML = html;
}

function update3DViewerOrGrid() {
    if (masterData.mode === 'pairwise' && document.getElementById('pair-selector').value === 'GRID') {
        Object.keys(masterData.pairs).forEach(k => update3DViewerGridForPair(k));
    } else update3DViewer();
}


