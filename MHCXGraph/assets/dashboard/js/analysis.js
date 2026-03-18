function renderAnalysis(isGlobalPairMode = false) {
    logDebug(`Rendering Analysis Tab. Global: ${isGlobalPairMode}`);
    try {
        const mData = isGlobalPairMode ? masterData : graphData;
        if (!mData) {
            logDebug("renderAnalysis aborted: Data is null.");
            return;
        }

        const displayProts = isGlobalPairMode ? masterData.proteins : mData.proteins;
        const numDisplay = displayProts.length; 
        logDebug(`Number of proteins to display: ${numDisplay}`, displayProts);
        
        let matrix = Array(numDisplay).fill(null).map(() => Array(numDisplay).fill(null));

        let origSizes = [];
        if (isGlobalPairMode) {
             origSizes = Array(numDisplay).fill(0);
             Object.values(mData.pairs).forEach(pData => {
                 let p1 = displayProts.indexOf(pData.proteins[0]);
                 let p2 = displayProts.indexOf(pData.proteins[1]);
                 origSizes[p1] = pData.filtered_graphs[0].nodes.length;
                 origSizes[p2] = pData.filtered_graphs[1].nodes.length;
             });
        } else {
             origSizes = mData.filtered_graphs.map(fg => fg.nodes.length);
        }
        logDebug("Original sizes calculated:", origSizes);

        let largestPatch = 0;

        if (isGlobalPairMode) {
            Object.keys(mData.pairs).forEach(pairKey => {
                let pData = mData.pairs[pairKey];
                let p1 = displayProts.indexOf(pData.proteins[0]);
                let p2 = displayProts.indexOf(pData.proteins[1]);
                
                let unique0 = new Set(); let unique1 = new Set();
                pData.components.forEach(comp => {
                    if (comp.id === 0) return;
                    comp.frames.forEach(frame => {
                        if (frame.id === 0) return;
                        if (frame.node_ids.length > largestPatch) largestPatch = frame.node_ids.length;
                        frame.node_ids.forEach(nId => {
                            let n = pData.nodes.find(node => node.id === nId);
                            if (n) {
                                n.mapping.forEach(m => {
                                    let localIdx = pData.proteins.indexOf(masterData.proteins[m.model_idx]);
                                    if (localIdx === 0) unique0.add(`${m.chain}:${m.resn}:${m.resi}`);
                                    if (localIdx === 1) unique1.add(`${m.chain}:${m.resn}:${m.resi}`);
                                });
                            }
                        });
                    });
                });
                
                let assocNodesCount = unique0.size + unique1.size; 
                let origSum = pData.filtered_graphs[0].nodes.length + pData.filtered_graphs[1].nodes.length;
                let ratio = origSum > 0 ? (assocNodesCount / origSum) : 0;
                let val = parseFloat((ratio * 100).toFixed(1));
                
                let cellObj = { pct: val, shared: assocNodesCount, total: origSum, u1: unique0.size, u2: unique1.size };
                matrix[p1][p2] = cellObj; matrix[p2][p1] = cellObj; 
            });
        } else {
            mData.components.forEach(c => {
                if (c.id !== 0) {
                    c.frames.forEach(f => {
                        if (f.id !== 0 && f.node_ids.length > largestPatch) largestPatch = f.node_ids.length;
                    });
                }
            });

            for (let i = 0; i < numDisplay; i++) {
                for (let j = 0; j < numDisplay; j++) {
                    let unique_i = new Set();
                    let unique_j = new Set();
                    
                    let global_i = masterData.proteins.indexOf(displayProts[i]);
                    let global_j = masterData.proteins.indexOf(displayProts[j]);
                    
                    mData.components.forEach(comp => {
                        if (comp.id === 0) return;
                        comp.frames.forEach(frame => {
                            if (frame.id === 0) return;
                            frame.node_ids.forEach(nId => {
                                let n = mData.nodes.find(node => node.id === nId);
                                if (n) {
                                    n.mapping.forEach(m => {
                                        if (m.model_idx === global_i) unique_i.add(`${m.chain}:${m.resn}:${m.resi}`);
                                        if (m.model_idx === global_j) unique_j.add(`${m.chain}:${m.resn}:${m.resi}`);
                                    });
                                }
                            });
                        });
                    });

                    if (i === j) {
                        let origSize = origSizes[i];
                        let ratio = origSize > 0 ? (unique_i.size / origSize) : 0;
                        let val = parseFloat((ratio * 100).toFixed(1));
                        matrix[i][j] = { pct: val, shared: unique_i.size, total: origSize, u1: unique_i.size, u2: unique_i.size };
                    } else {
                        let assocNodesCount = unique_i.size + unique_j.size; 
                        let origSum = origSizes[i] + origSizes[j];
                        let ratio = origSum > 0 ? (assocNodesCount / origSum) : 0;
                        let val = parseFloat((ratio * 100).toFixed(1));
                        matrix[i][j] = { pct: val, shared: assocNodesCount, total: origSum, u1: unique_i.size, u2: unique_j.size };
                    }
                }
            }
        }
        logDebug("Matrix populated successfully.");

        const wrapper = document.getElementById('analysis-wrapper');
        let framesCount = 0;
        let compsCount = isGlobalPairMode ? Object.keys(mData.pairs).length : (mData.components.length - 1);

        if (!isGlobalPairMode) {
            mData.components.forEach(c => { if(c.id !== 0) framesCount += c.frames.length - 1; });
        } else {
            Object.values(mData.pairs).forEach(pData => {
                pData.components.forEach(c => { if(c.id !== 0) framesCount += c.frames.length - 1; });
            });
        }

        const safeMetadata = mData.metadata || masterData.metadata || {};
        const granularity = safeMetadata.node_granularity || 'N/A';

        logDebug("Injecting wrapper HTML...");
        wrapper.innerHTML = `
            <div class="summary-cards">
                <div class="card"><h4>Proteins Analyzed</h4><div class="value">${numDisplay}</div></div>
                <div class="card"><h4>Conserved Components</h4><div class="value">${compsCount}</div></div>
                <div class="card"><h4>Structural Frames</h4><div class="value">${framesCount}</div></div>
                <div class="card"><h4>Largest Patch</h4><div class="value">${largestPatch} <span style="font-size:12px; color:var(--text-muted)">nodes</span></div></div>
                <div class="card"><h4>Granularity</h4><div class="value" style="font-size:18px; margin-top:10px; color:var(--text-main);">${granularity}</div></div>
            </div>

            <div class="analysis-section" id="hm-section">
                <h2 id="analysis-title-1">
                    ${isGlobalPairMode ? 'Pairwise Coverage Heatmap (%)' : 'Protein Co-occurrence Heatmap'}
                </h2>
                <div class="heatmap-wrapper" id="heatmap-container"></div>
            </div>

            <div class="analysis-section">
                <h2>${isGlobalPairMode ? 'Pairwise Association Metrics' : 'Per-Frame Coverage & Properties'}</h2>
                <div style="max-height: 500px; overflow-y: auto; border: 1px solid var(--border-light); border-radius: 6px;">
                    <table class="metrics-table">
                        <thead><tr id="metrics-table-head"></tr></thead>
                        <tbody id="metrics-table-body"></tbody>
                    </table>
                </div>
            </div>
        `;

        logDebug("Building Heatmap...");
        const hmContainer = document.getElementById('heatmap-container');
        const grid = document.createElement('div');
        grid.className = 'heatmap-grid';
        grid.style.gridTemplateColumns = `auto repeat(${numDisplay}, 1fr)`;
        
        grid.appendChild(document.createElement('div')); 
        
        for (let i = 0; i < numDisplay; i++) {
            let h = document.createElement('div'); h.className = 'heatmap-header';
            h.innerHTML = `${displayProts[i]}<br><span style="font-size:10px; color:var(--text-muted); font-weight:normal;">(${origSizes[i]} nodes)</span>`;
            grid.appendChild(h);
        }
        
        for (let i = 0; i < numDisplay; i++) {
            let rh = document.createElement('div'); rh.className = 'heatmap-header';
            rh.style.justifyContent = "flex-end";
            rh.style.paddingRight = "15px";
            rh.innerHTML = `${displayProts[i]} <span style="font-size:10px; color:var(--text-muted); font-weight:normal; margin-left:5px;">(${origSizes[i]})</span>`;
            grid.appendChild(rh);
            
            for (let j = 0; j < numDisplay; j++) {
                let cell = document.createElement('div'); 
                let cellData = matrix[i][j];
                cell.className = 'heatmap-cell'; 
                
                let isDark = document.body.getAttribute('data-theme') === 'dark';
                
                if (!cellData || cellData.pct === 0) {
                    cell.innerText = cellData ? '0%' : '-';
                    cell.style.background = getCSSVar('--bg-control');
                    cell.style.color = getCSSVar('--text-faint');
                    cell.style.textShadow = 'none';
                    cell.style.cursor = 'default';
                } else {
                    let val = cellData.pct;
                    cell.innerText = `${val}%`;
                    
                    let intensity = val / 100.0;
                    let hue = 280 - (intensity * 230); 
                    let light = 35 + (intensity * 35);
                    cell.style.background = `hsl(${hue}, 85%, ${light}%)`;
                    
                    const tt = document.createElement('div');
                    tt.className = 'hm-tooltip';
                    let pName_i = displayProts[i];
                    let pName_j = displayProts[j];
                    
                    if (i === j) {
                        tt.innerHTML = `
                            <div style="border-bottom: 1px solid rgba(255,255,255,0.2); padding-bottom: 4px; margin-bottom: 4px; font-weight: bold;">
                                ${pName_i} (Self)
                            </div>
                            Conserved Nodes: <b>${cellData.shared}</b> / ${origSizes[i]}<br>
                            Self-Coverage: <b style="color: #FFD54F;">${val}%</b>
                        `;
                    } else {
                        tt.innerHTML = `
                            <div style="border-bottom: 1px solid rgba(255,255,255,0.2); padding-bottom: 4px; margin-bottom: 4px; font-weight: bold;">
                                ${pName_i} x ${pName_j}
                            </div>
                            ${pName_i} Match: <b>${cellData.u1}</b> / ${origSizes[i]}<br>
                            ${pName_j} Match: <b>${cellData.u2}</b> / ${origSizes[j]}<br>
                            Total Combined Coverage: <b style="color: #FFD54F;">${val}%</b>
                        `;
                    }
                    cell.appendChild(tt);
                }
                grid.appendChild(cell);
            }
        }
        hmContainer.appendChild(grid);

        logDebug("Building Metrics Table...");
        const tHead = document.getElementById('metrics-table-head');
        const tBody = document.getElementById('metrics-table-body');
        
        if (isGlobalPairMode) {
            tHead.innerHTML = `<th>Pair Name</th><th>Original Nodes (P1 + P2)</th><th>Unique Associated Residues</th><th>Total Coverage</th>`;
            let bodyHtml = '';
            Object.keys(mData.pairs).forEach(pairKey => {
                let pData = mData.pairs[pairKey];
                let unique0 = new Set(); let unique1 = new Set();
                pData.components.forEach(comp => {
                    if (comp.id === 0) return;
                    comp.frames.forEach(frame => {
                        if (frame.id === 0) return;
                        frame.node_ids.forEach(nId => {
                            let n = pData.nodes.find(node => node.id === nId);
                            if (n) {
                                n.mapping.forEach(m => {
                                    let localIdx = pData.proteins.indexOf(masterData.proteins[m.model_idx]);
                                    if (localIdx === 0) unique0.add(`${m.chain}:${m.resn}:${m.resi}`);
                                    if (localIdx === 1) unique1.add(`${m.chain}:${m.resn}:${m.resi}`);
                                });
                            }
                        });
                    });
                });
                
                let assocNodesCount = unique0.size + unique1.size; 
                let origSum = pData.filtered_graphs[0].nodes.length + pData.filtered_graphs[1].nodes.length;
                let ratio = origSum > 0 ? ((assocNodesCount / origSum)*100).toFixed(1) : 0;
                bodyHtml += `<tr><td><b>${pairKey.replace('_vs_', ' x ')}</b></td><td>${origSum}</td><td>${assocNodesCount}</td><td><span style="color:var(--btn-bg); font-weight:bold;">${ratio}%</span></td></tr>`;
            });
            tBody.innerHTML = bodyHtml || `<tr><td colspan="4" style="text-align:center;">No pairs found.</td></tr>`;
        } else {
            let headHtml = `<th>Comp</th><th>Frame</th><th>Residues</th><th>Edge Density</th><th>Residue Composition</th>`;
            
            let hasRmsd = false;
            if(mData.components.length > 1 && mData.components[1].frames.length > 1) {
                if (mData.components[1].frames[1].rmsds && Object.keys(mData.components[1].frames[1].rmsds).length > 0) hasRmsd = true;
            }
            if (hasRmsd) headHtml += `<th>Align RMSD</th>`;
            
            for (let i = 0; i < numDisplay; i++) {
                let pName = displayProts[i];
                headHtml += `<th>${pName} Cov.<br><span style="font-size:10px; color:var(--text-muted); text-transform:none;">(of ${origSizes[i]})</span></th>`;
            }
            tHead.innerHTML = headHtml;
            
            let bodyHtml = '';
            mData.components.forEach(comp => {
                if (comp.id === 0) return; 
                comp.frames.forEach(frame => {
                    if (frame.id === 0) return; 
                    
                    let frameNodes = frame.node_ids.length; 
                    let protUniqueResidues = Array(numDisplay).fill(null).map(() => new Set());
                    let classCounts = { 'Hydrophobic':0, 'Polar':0, 'Charged (+)':0, 'Charged (-)':0, 'Special':0, 'Unknown':0 };
                    
                    frame.node_ids.forEach(nId => {
                        let n = mData.nodes.find(node => node.id === nId);
                        if (n) {
                            n.mapping.forEach(m => {
                                let localIdx = displayProts.indexOf(masterData.proteins[m.model_idx]);
                                if (localIdx !== -1) {
                                    let resKey = `${m.chain}:${m.resn}:${m.resi}`;
                                    if(!protUniqueResidues[localIdx].has(resKey)){
                                        protUniqueResidues[localIdx].add(resKey);
                                        let resClass = RESIDUE_CLASSES[m.resn] || 'Unknown';
                                        classCounts[resClass]++;
                                    }
                                }
                            });
                        }
                    });
                    
                    let protCounts = protUniqueResidues.map(s => s.size);
                    let totalUniqueResidues = protCounts.reduce((a, b) => a + b, 0);
                    
                    let activeNodeSet = new Set(frame.node_ids);
                    let internalEdges = 0;
                    mData.edges.forEach(e => {
                        if (activeNodeSet.has(e.from) && activeNodeSet.has(e.to)) internalEdges++;
                    });
                    let maxEdges = (frameNodes * (frameNodes - 1)) / 2;
                    let edgeDensity = maxEdges > 0 ? ((internalEdges / maxEdges) * 100).toFixed(1) : 0;

                    let classHtml = '';
                    Object.keys(classCounts).forEach(cls => {
                        if (classCounts[cls] > 0) {
                            let pct = Math.round((classCounts[cls] / totalUniqueResidues) * 100);
                            classHtml += `<span class="class-badge" style="background-color:${CLASS_COLORS[cls]};" title="${classCounts[cls]} residues">${cls} ${pct}%</span>`;
                        }
                    });
                    
                    let row = `<tr>
                        <td><b>C${comp.id}</b></td>
                        <td>F${frame.id}</td>
                        <td>${frameNodes}</td>
                        <td>${edgeDensity}%</td>
                        <td style="max-width: 200px; line-height: 1.8;">${classHtml}</td>`;
                    
                    if (hasRmsd) {
                        let rmsdHtml = '';
                        if (frame.rmsds) {
                            Object.entries(frame.rmsds).forEach(([k, v]) => {
                                rmsdHtml += `<div style="font-size:11px;">${k}: <b>${v}Å</b></div>`;
                            });
                        }
                        row += `<td>${rmsdHtml || '-'}</td>`;
                    }

                    for (let i = 0; i < numDisplay; i++) {
                        let origSize = origSizes[i];
                        let coverage = origSize > 0 ? ((protCounts[i] / origSize) * 100).toFixed(1) : 0;
                        row += `<td><span style="color:var(--btn-bg); font-weight:bold;">${coverage}%</span><br><span style="font-size:10px; color:var(--text-muted);">${protCounts[i]} nodes</span></td>`;
                    }
                    row += `</tr>`; bodyHtml += row;
                });
            });
            tBody.innerHTML = bodyHtml || `<tr><td colspan="${5 + (hasRmsd?1:0) + numDisplay}" style="text-align:center; padding: 20px;">No component frames found.</td></tr>`;
        }
        logDebug("Analysis Tab rendering complete.");
    } catch (e) { logError("renderAnalysis crashed", e); }
}
