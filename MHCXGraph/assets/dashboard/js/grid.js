function buildGridTop() {
    try {
        const topPanel = document.getElementById('top-panel');
        const singleNet = document.getElementById('single-network'); 
        if (singleNet) singleNet.style.display = 'none'; 

        let gridWrapper = document.getElementById('dynamic-grid-networks');
        if (!gridWrapper) {
            gridWrapper = document.createElement('div');
            gridWrapper.id = 'dynamic-grid-networks';
            gridWrapper.style.width = '100%'; gridWrapper.style.height = '100%';
            gridWrapper.style.padding = '15px'; gridWrapper.style.boxSizing = 'border-box';
            gridWrapper.style.overflowY = 'auto'; gridWrapper.style.background = 'var(--bg-body)';
            topPanel.appendChild(gridWrapper);
        }
        
        gridWrapper.innerHTML = ''; 
        gridWrapper.style.display = 'grid';
        
        let numPairs = Object.keys(masterData.pairs).length;
        let cols = numPairs > 4 ? 3 : (numPairs > 1 ? 2 : 1);
        gridWrapper.style.gridTemplateColumns = `repeat(${cols}, 1fr)`;
        gridWrapper.style.gridAutoRows = '350px';
        gridWrapper.style.gap = '15px';

        Object.keys(masterData.pairs).forEach(pairKey => {
            let pData = masterData.pairs[pairKey];
            
            let card = document.createElement('div');
            card.style.border = '1px solid var(--border)'; card.style.borderRadius = '8px';
            card.style.display = 'flex'; card.style.flexDirection = 'column';
            card.style.background = 'var(--bg-panel)'; card.style.overflow = 'hidden';
            
            let header = document.createElement('div');
            header.style.padding = '8px 12px'; header.style.fontWeight = 'bold'; header.style.background = 'var(--bg-control)';
            header.style.borderBottom = '1px solid var(--border)'; header.style.color = 'var(--text-main)';
            header.innerHTML = `<div style="display:flex; justify-content:space-between; align-items:center;">
                <span>${pairKey.replace('_vs_', ' x ')}</span>
                <button class="btn" style="padding: 4px 10px; font-size: 11px;" onclick="focusOnPair('${pairKey}')">🔍 Focus</button>
            </div>`;
            card.appendChild(header);
            
            let netDiv = document.createElement('div');
            netDiv.style.flex = '1'; netDiv.style.position = 'relative'; netDiv.style.minHeight = '0';
            card.appendChild(netDiv);
            gridWrapper.appendChild(card); 
            
            setTimeout(() => {
                let dsNodes = new vis.DataSet(pData.nodes.map(n => ({...n, label: '<b>'+n.label+'</b>', color: {background: n.originalColor || getCSSVar('--edge-default'), border: themeBorder}, font: {color: themeText}})));
                let dsEdges = new vis.DataSet(pData.edges.map(e => ({...e, color: {color: getCSSVar('--edge-faded'), opacity: 0.8}, width: optEdgeWidth})));
                gridNodesDatasets[pairKey] = dsNodes; gridEdgesDatasets[pairKey] = dsEdges;
                
                let net = new vis.Network(netDiv, {nodes: dsNodes, edges: dsEdges}, {
                    nodes: { shape: 'dot', size: optNodeSize, font: { size: optLabelSize-4, multi: 'html', face: 'ui-sans-serif' }, borderWidth: 2 },
                    edges: { smooth: !document.getElementById('optStraightEdges').checked },
                    physics: { solver: 'repulsion', repulsion: { nodeDistance: 100, springLength: 90 }, stabilization: { iterations: 40, fit: true } }, 
                    interaction: {zoomView: true, dragView: true}
                });
                
                net.once("stabilizationIterationsDone", function () { 
                    if (document.getElementById('freezePhysics').checked) net.setOptions({ physics: false }); 
                });
                gridNetworks.push({ pair: pairKey, network: net });
                handleGridTreeChange(pairKey); 
            }, 10);
        });
    } catch(e) { logError("buildGridTop crashed", e); }
}

function buildGridBottom() {
    try {
        viewers = []; 
        window.activeMolViewer = null;

        const bottomPanel = document.getElementById('bottom-panel');
        const singleMols = document.getElementById('single-mols'); 
        if (singleMols) singleMols.style.display = 'none'; 

        let gridWrapper = document.getElementById('dynamic-grid-mols');
        if (!gridWrapper) {
            gridWrapper = document.createElement('div');
            gridWrapper.id = 'dynamic-grid-mols';
            gridWrapper.style.width = '100%'; 
            gridWrapper.style.flex = '1'; 
            gridWrapper.style.padding = '15px'; 
            gridWrapper.style.boxSizing = 'border-box';
            gridWrapper.style.overflowY = 'auto'; 
            gridWrapper.style.background = 'var(--bg-body)';
            bottomPanel.appendChild(gridWrapper);
        }
        
        gridWrapper.innerHTML = ''; 
        gridWrapper.style.display = 'grid';
        
        const layoutMode = document.getElementById('view-selector').value;
        const isDark = document.body.getAttribute('data-theme') === 'dark';
        const bgMol = getCSSVar('--bg-mol');
        
        let numPairs = Object.keys(masterData.pairs).length;
        let cols = numPairs > 4 ? 3 : (numPairs > 1 ? 2 : 1);
        gridWrapper.style.gridTemplateColumns = `repeat(${cols}, 1fr)`;
        gridWrapper.style.gridAutoRows = '350px';
        gridWrapper.style.gap = '15px';

        Object.keys(masterData.pairs).forEach(pairKey => {
            let pData = masterData.pairs[pairKey];
            let p1_idx = masterData.proteins.indexOf(pData.proteins[0]);
            let p2_idx = masterData.proteins.indexOf(pData.proteins[1]);
            
            let card = document.createElement('div');
            card.style.border = '1px solid var(--border)'; card.style.borderRadius = '8px';
            card.style.display = 'flex'; card.style.flexDirection = 'column';
            card.style.background = 'var(--bg-panel)'; card.style.overflow = 'hidden';
            
            let header = document.createElement('div');
            header.style.padding = '8px 12px'; header.style.fontWeight = 'bold'; header.style.background = 'var(--bg-control)';
            header.style.borderBottom = '1px solid var(--border)'; header.style.textAlign = 'center'; header.style.color = 'var(--text-main)';
            header.innerText = pairKey.replace('_vs_', ' x ');
            card.appendChild(header);

            let body = document.createElement('div');
            body.style.flex = '1'; body.style.display = 'flex'; body.style.position = 'relative'; body.style.minHeight = '0';
            card.appendChild(body);
            gridWrapper.appendChild(card);

            if (layoutMode === 'separate') {
                [p1_idx, p2_idx].forEach(idx => {
                    let molDiv = document.createElement('div');
                    molDiv.style.flex = '1'; molDiv.style.height = '100%'; molDiv.style.position = 'relative';
                    if (idx === p1_idx) molDiv.style.borderRight = '1px solid var(--border-light)';
                    body.appendChild(molDiv);
                    
                    let titleDiv = document.createElement('div');
                    titleDiv.className = 'viewer-title';
                    titleDiv.style.fontSize = '11px';
                    titleDiv.style.padding = '2px 6px';
                    let pName = masterData.proteins[idx];
                    let baseCol = get3DColor('A', idx);
                    titleDiv.innerHTML = `<span class="color-dot" style="background-color: ${baseCol}; width: 8px; height: 8px;"></span>${pName}`;
                    molDiv.appendChild(titleDiv);

                    let canvasWrapper = document.createElement('div');
                    canvasWrapper.style.width = '100%'; canvasWrapper.style.height = '100%'; canvasWrapper.style.position = 'absolute';
                    molDiv.appendChild(canvasWrapper);

                    const initViewer = () => {
                        if (canvasWrapper.clientWidth === 0 || canvasWrapper.clientHeight === 0) {
                            setTimeout(initViewer, 20);
                            return;
                        }
                        let v = $3Dmol.createViewer(canvasWrapper, { backgroundColor: bgMol });
                        canvasWrapper.addEventListener('mousedown', () => { window.activeMolViewer = v; });
                        canvasWrapper.addEventListener('wheel', () => { window.activeMolViewer = v; }, {passive: true});
                        
                        if (loadedModels[idx]) {
                            v.addModel(loadedModels[idx].text, loadedModels[idx].format);
                            v.setStyle({}, {cartoon: {colorfunc: (atom) => get3DColor(atom.chain, idx)}});
                            v.zoomTo();
                        }
                        v.pairKey = pairKey; v.protIdxs = [idx]; viewers.push(v);
                        update3DViewerGridForPair(pairKey);
                    };
                    initViewer();
                });
            } else {
                let molDiv = document.createElement('div');
                molDiv.style.width = '100%'; molDiv.style.height = '100%'; molDiv.style.position = 'absolute';
                body.appendChild(molDiv);
                
                const initViewer = () => {
                    if (molDiv.clientWidth === 0 || molDiv.clientHeight === 0) {
                        setTimeout(initViewer, 20);
                        return;
                    }
                    let v = $3Dmol.createViewer(molDiv, { backgroundColor: bgMol });
                    molDiv.addEventListener('mousedown', () => { window.activeMolViewer = v; });
                    molDiv.addEventListener('wheel', () => { window.activeMolViewer = v; }, {passive: true});
                    
                    let hasMol = false;
                    [p1_idx, p2_idx].forEach(idx => {
                        if (loadedModels[idx]) {
                            let mObj = v.addModel(loadedModels[idx].text, loadedModels[idx].format);
                            v.setStyle({model: mObj.getID()}, {cartoon: {colorfunc: (atom) => get3DColor(atom.chain, idx)}});
                            hasMol = true;
                        }
                    });
                    v.pairKey = pairKey; v.protIdxs = [p1_idx, p2_idx]; 
                    if(hasMol) v.zoomTo(); viewers.push(v);
                    update3DViewerGridForPair(pairKey);
                };
                initViewer();
            }
        });
    } catch(e) { logError("buildGridBottom crashed", e); }
}
