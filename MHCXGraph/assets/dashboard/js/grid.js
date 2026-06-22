/* =========================================================================
 * grid.js — lazy-loaded grid view for pairwise / screening dashboards
 *
 * Why this file looks different from before:
 *   - The previous version instantiated one vis.Network AND one (or two)
 *     $3Dmol viewers per pair on grid-open. With N proteins → P = N(N-1)/2
 *     pairs, this saturates the WebGL context budget (~8-16 per browser)
 *     and floods the main thread with physics simulations.
 *   - This version renders lightweight placeholder cards immediately and
 *     instantiates the heavy widgets only when a card is in the viewport
 *     (IntersectionObserver) or the user clicks "Load". A strict budget
 *     evicts the oldest off-screen instances when the budget is exceeded.
 *
 * Public API kept the same so the rest of the dashboard does not change:
 *   - buildGridTop()         build the vis-network grid (top panel)
 *   - buildGridBottom()      build the 3Dmol grid (bottom panel)
 *   - update3DViewerGridForPair(pairKey)   used by graph.js
 *   - applyGraphFiltersGrid(...)           used by graph.js
 *
 * New API used internally:
 *   - lazyMode, gridStats           current lazy state + stats
 *   - loadVisForPair / loadMolForPair
 *   - unloadVisForPair / unloadMolForPair
 *   - rebuildGridCards()            re-create placeholders after theme change
 * ========================================================================= */

/* --- Lazy-load budget and thresholds ------------------------------------- */
/* Tunable; chosen to keep most desktop browsers happy:
 *   - WebGL contexts cap at ~8-16 per page across most browsers/GPUs.
 *   - vis-network instances are CPU/memory heavy but not WebGL.
 *   - Above LAZY_PAIRS_THRESHOLD or HEAVY_NODES_THRESHOLD we enter lazy
 *     mode automatically. The user can override with "Load all".
 */
const LAZY_PAIRS_THRESHOLD = 12;
const HEAVY_NODES_THRESHOLD = 5000;
const MAX_VIS_CONCURRENT = 30;
const MAX_MOL_CONCURRENT = 8;

/* --- Module-local state -------------------------------------------------- */
let lazyMode = false;
let gridStats = { pairCount: 0, totalNodes: 0, totalEdges: 0 };
let visLoadedKeys = [];   /* LRU order: oldest first */
let molLoadedKeys = [];
let visObserver = null;
let molObserver = null;
const cardRefs = {};      /* pairKey -> { topCard, topBody, botCard, botBody } */

/* ========================================================================= */
/* Stats and threshold detection                                             */
/* ========================================================================= */

function computeGridStats() {
    const pairs = (masterData && masterData.pairs) || {};
    let totalNodes = 0, totalEdges = 0;
    const keys = Object.keys(pairs);
    keys.forEach(k => {
        const p = pairs[k];
        totalNodes += (p.nodes || []).length;
        totalEdges += (p.edges || []).length;
    });
    return { pairCount: keys.length, totalNodes, totalEdges };
}

function decideLazyMode(stats) {
    return stats.pairCount >= LAZY_PAIRS_THRESHOLD
        || stats.totalNodes >= HEAVY_NODES_THRESHOLD;
}

/* Quick per-pair stats for the placeholder header */
function pairSummary(pData) {
    const nNodes = (pData.nodes || []).length;
    const nEdges = (pData.edges || []).length;
    let nAssoc = 0;
    if (pData.components) {
        /* Component 0 == "global associated" in this codebase */
        const c0 = pData.components.find(c => c.id === 0);
        if (c0 && c0.node_ids) nAssoc = c0.node_ids.length;
        else nAssoc = nNodes;
    }
    const nComp = (pData.components || []).filter(c => c.id !== 0).length;
    return { nNodes, nEdges, nAssoc, nComp };
}

/* ========================================================================= */
/* Warning banner                                                            */
/* ========================================================================= */

function buildHeavyWarningBanner(gridWrapper) {
    const banner = document.createElement('div');
    banner.id = 'grid-heavy-warning';
    banner.style.cssText = [
        'grid-column: 1 / -1',
        'background: var(--bg-control)',
        'border: 1px solid var(--border)',
        'border-left: 4px solid var(--btn-bg)',
        'border-radius: 6px',
        'padding: 10px 14px',
        'margin-bottom: 4px',
        'color: var(--text-main)',
        'font-size: 13px',
        'display: flex',
        'flex-wrap: wrap',
        'gap: 10px',
        'align-items: center',
        'justify-content: space-between'
    ].join(';');
    banner.innerHTML = `
        <div>
            <b>⚠ Heavy dataset — lazy mode on.</b>
            <span style="color: var(--text-muted);">
                ${gridStats.pairCount} pairs · ${gridStats.totalNodes.toLocaleString()} nodes total.
                Graphs render only when a card scrolls into view
                (budget: ${MAX_VIS_CONCURRENT} graphs, ${MAX_MOL_CONCURRENT} 3D viewers at once).
            </span>
            <span id="grid-loaded-counter" style="margin-left: 8px; color: var(--text-faint);"></span>
        </div>
        <div style="display: flex; gap: 6px;">
            <button class="btn" style="padding: 4px 10px; font-size: 12px;"
                onclick="loadAllGridPairs(this)">Load all anyway</button>
            <button class="btn" style="padding: 4px 10px; font-size: 12px; background: var(--bg-panel); color: var(--text-main); border: 1px solid var(--border);"
                onclick="unloadAllGridPairs()">Unload all</button>
        </div>
    `;
    gridWrapper.appendChild(banner);
}

function updateLoadedCounter() {
    const el = document.getElementById('grid-loaded-counter');
    if (!el) return;
    el.textContent = `(loaded: ${visLoadedKeys.length} graphs, ${molLoadedKeys.length} 3D viewers)`;
}

/* User-facing controls referenced by the banner */
window.loadAllGridPairs = function(btn) {
    if (gridStats.pairCount > LAZY_PAIRS_THRESHOLD * 4) {
        if (!confirm(`This will instantiate ${gridStats.pairCount} graphs and 3D viewers at once. ` +
                     `On a heavy dataset this can freeze the page. Continue?`)) return;
    }
    if (btn) btn.disabled = true;
    /* Stagger to keep the UI responsive */
    const keys = Object.keys(masterData.pairs);
    let i = 0;
    function step() {
        const slice = keys.slice(i, i + 4);
        slice.forEach(k => { loadVisForPair(k); loadMolForPair(k); });
        i += 4;
        if (i < keys.length) setTimeout(step, 30);
        else if (btn) btn.disabled = false;
    }
    step();
};

window.unloadAllGridPairs = function() {
    Object.keys(masterData.pairs).forEach(k => {
        unloadVisForPair(k);
        unloadMolForPair(k);
    });
};

/* ========================================================================= */
/* Card construction (placeholder)                                           */
/* ========================================================================= */

function makeTopCard(pairKey, pData) {
    const s = pairSummary(pData);
    const card = document.createElement('div');
    card.dataset.pair = pairKey;
    card.style.cssText = [
        'border: 1px solid var(--border)',
        'border-radius: 8px',
        'display: flex',
        'flex-direction: column',
        'background: var(--bg-panel)',
        'overflow: hidden',
        'min-height: 0'
    ].join(';');

    const header = document.createElement('div');
    header.style.cssText = [
        'padding: 8px 12px',
        'font-weight: bold',
        'background: var(--bg-control)',
        'border-bottom: 1px solid var(--border)',
        'color: var(--text-main)'
    ].join(';');
    header.innerHTML = `
        <div style="display:flex; justify-content:space-between; align-items:center; gap: 8px;">
            <span style="overflow: hidden; text-overflow: ellipsis; white-space: nowrap;">
                ${pairKey.replace('_vs_', ' × ')}
            </span>
            <span style="font-weight: normal; font-size: 11px; color: var(--text-muted); white-space: nowrap;">
                ${s.nNodes}n · ${s.nEdges}e · ${s.nAssoc} assoc
            </span>
            <span style="display: flex; gap: 4px;">
                <button class="btn grid-toggle-load" style="padding: 3px 8px; font-size: 11px;"
                    onclick="toggleVisForPair('${pairKey}', this)">Load</button>
                <button class="btn" style="padding: 3px 8px; font-size: 11px;"
                    onclick="focusOnPair('${pairKey}')">🔍 Focus</button>
            </span>
        </div>
    `;
    card.appendChild(header);

    const body = document.createElement('div');
    body.style.cssText = 'flex: 1; position: relative; min-height: 0; background: var(--bg-panel);';
    /* Placeholder visible until load */
    const placeholder = document.createElement('div');
    placeholder.className = 'grid-placeholder';
    placeholder.style.cssText = [
        'position: absolute',
        'inset: 0',
        'display: flex',
        'align-items: center',
        'justify-content: center',
        'color: var(--text-faint)',
        'font-size: 12px',
        'text-align: center',
        'padding: 10px',
        'box-sizing: border-box'
    ].join(';');
    placeholder.innerHTML = `
        <div>
            <div style="font-size: 22px; opacity: 0.5;">⏸</div>
            <div>Graph not loaded</div>
            <div style="font-size: 11px; margin-top: 4px;">
                ${s.nComp} components · ${s.nNodes} nodes
            </div>
        </div>
    `;
    body.appendChild(placeholder);
    card.appendChild(body);

    return { card, body };
}

function makeBotCard(pairKey, pData) {
    const s = pairSummary(pData);
    const card = document.createElement('div');
    card.dataset.pair = pairKey;
    card.style.cssText = [
        'border: 1px solid var(--border)',
        'border-radius: 8px',
        'display: flex',
        'flex-direction: column',
        'background: var(--bg-panel)',
        'overflow: hidden',
        'min-height: 0'
    ].join(';');

    const header = document.createElement('div');
    header.style.cssText = [
        'padding: 8px 12px',
        'font-weight: bold',
        'background: var(--bg-control)',
        'border-bottom: 1px solid var(--border)',
        'text-align: center',
        'color: var(--text-main)',
        'display: flex',
        'justify-content: space-between',
        'gap: 8px'
    ].join(';');
    header.innerHTML = `
        <span style="overflow: hidden; text-overflow: ellipsis; white-space: nowrap;">
            ${pairKey.replace('_vs_', ' × ')}
        </span>
        <span style="font-weight: normal; font-size: 11px; color: var(--text-muted); white-space: nowrap;">
            ${s.nNodes}n
        </span>
        <button class="btn grid-toggle-mol" style="padding: 3px 8px; font-size: 11px;"
            onclick="toggleMolForPair('${pairKey}', this)">Load 3D</button>
    `;
    card.appendChild(header);

    const body = document.createElement('div');
    body.style.cssText = 'flex: 1; display: flex; position: relative; min-height: 0;';

    const placeholder = document.createElement('div');
    placeholder.className = 'grid-placeholder';
    placeholder.style.cssText = [
        'position: absolute',
        'inset: 0',
        'display: flex',
        'align-items: center',
        'justify-content: center',
        'color: var(--text-faint)',
        'font-size: 12px',
        'text-align: center'
    ].join(';');
    placeholder.innerHTML = `<div><div style="font-size: 22px; opacity: 0.5;">🧬</div><div>3D not loaded</div></div>`;
    body.appendChild(placeholder);
    card.appendChild(body);

    return { card, body };
}

/* ========================================================================= */
/* Loaders / unloaders                                                       */
/* ========================================================================= */

function loadVisForPair(pairKey) {
    if (gridNodesDatasets[pairKey]) return; /* already loaded */
    const refs = cardRefs[pairKey];
    if (!refs || !refs.topBody) return;
    const pData = masterData.pairs[pairKey];
    if (!pData) return;

    /* Budget eviction — drop the oldest loaded that isn't this one */
    while (visLoadedKeys.length >= MAX_VIS_CONCURRENT) {
        const victim = visLoadedKeys.find(k => k !== pairKey);
        if (!victim) break;
        unloadVisForPair(victim);
    }

    try {
        const placeholder = refs.topBody.querySelector('.grid-placeholder');
        if (placeholder) placeholder.style.display = 'none';

        const dsNodes = new vis.DataSet(pData.nodes.map(n => ({
            ...n,
            label: '<b>' + n.label + '</b>',
            color: { background: n.originalColor || getCSSVar('--edge-default'), border: themeBorder },
            font: { color: themeText }
        })));
        const dsEdges = new vis.DataSet(pData.edges.map(e => ({
            ...e,
            color: { color: getCSSVar('--edge-faded'), opacity: 0.8 },
            width: optEdgeWidth
        })));
        gridNodesDatasets[pairKey] = dsNodes;
        gridEdgesDatasets[pairKey] = dsEdges;

        const net = new vis.Network(refs.topBody, { nodes: dsNodes, edges: dsEdges }, {
            nodes: { shape: 'dot', size: optNodeSize, font: { size: optLabelSize - 4, multi: 'html', face: 'ui-sans-serif' }, borderWidth: 2 },
            edges: { smooth: !document.getElementById('optStraightEdges').checked },
            physics: { solver: 'repulsion', repulsion: { nodeDistance: 100, springLength: 90 }, stabilization: { iterations: 40, fit: true } },
            interaction: { zoomView: true, dragView: true }
        });
        net.once('stabilizationIterationsDone', () => {
            if (document.getElementById('freezePhysics').checked) net.setOptions({ physics: false });
        });

        gridNetworks.push({ pair: pairKey, network: net });
        visLoadedKeys.push(pairKey);

        /* Apply current tree filter to the just-loaded pair */
        handleGridTreeChange(pairKey);

        /* Update button label */
        const btn = refs.topCard && refs.topCard.querySelector('.grid-toggle-load');
        if (btn) btn.textContent = 'Unload';

        updateLoadedCounter();
    } catch (e) {
        logError(`loadVisForPair(${pairKey}) crashed`, e);
    }
}

function unloadVisForPair(pairKey) {
    const refs = cardRefs[pairKey];
    const idx = gridNetworks.findIndex(g => g.pair === pairKey);
    if (idx >= 0) {
        try { gridNetworks[idx].network.destroy(); } catch (e) {}
        gridNetworks.splice(idx, 1);
    }
    delete gridNodesDatasets[pairKey];
    delete gridEdgesDatasets[pairKey];
    visLoadedKeys = visLoadedKeys.filter(k => k !== pairKey);

    if (refs && refs.topBody) {
        /* Drop everything except the placeholder */
        const placeholder = refs.topBody.querySelector('.grid-placeholder');
        refs.topBody.innerHTML = '';
        if (placeholder) {
            placeholder.style.display = 'flex';
            refs.topBody.appendChild(placeholder);
        }
    }
    if (refs && refs.topCard) {
        const btn = refs.topCard.querySelector('.grid-toggle-load');
        if (btn) btn.textContent = 'Load';
    }
    updateLoadedCounter();
}

window.toggleVisForPair = function(pairKey, btn) {
    if (gridNodesDatasets[pairKey]) unloadVisForPair(pairKey);
    else loadVisForPair(pairKey);
};

function loadMolForPair(pairKey) {
    const refs = cardRefs[pairKey];
    if (!refs || !refs.botBody) return;
    if (molLoadedKeys.includes(pairKey)) return;
    const pData = masterData.pairs[pairKey];
    if (!pData) return;

    /* Budget eviction — drop the oldest loaded that isn't this one */
    while (molLoadedKeys.length >= MAX_MOL_CONCURRENT) {
        const victim = molLoadedKeys.find(k => k !== pairKey);
        if (!victim) break;
        unloadMolForPair(victim);
    }

    try {
        const layoutMode = document.getElementById('view-selector').value;
        const bgMol = getCSSVar('--bg-mol');
        const p1_idx = masterData.proteins.indexOf(pData.proteins[0]);
        const p2_idx = masterData.proteins.indexOf(pData.proteins[1]);

        const placeholder = refs.botBody.querySelector('.grid-placeholder');
        if (placeholder) placeholder.style.display = 'none';

        if (layoutMode === 'separate') {
            [p1_idx, p2_idx].forEach(idx => {
                const molDiv = document.createElement('div');
                molDiv.className = 'grid-mol-container';
                molDiv.style.cssText = 'flex: 1; height: 100%; position: relative;';
                if (idx === p1_idx) molDiv.style.borderRight = '1px solid var(--border-light)';
                refs.botBody.appendChild(molDiv);

                const titleDiv = document.createElement('div');
                titleDiv.className = 'viewer-title';
                titleDiv.style.cssText = 'font-size: 11px; padding: 2px 6px;';
                const pName = masterData.proteins[idx];
                const baseCol = get3DColor('A', idx);
                titleDiv.innerHTML = `<span class="color-dot" style="background-color: ${baseCol}; width: 8px; height: 8px;"></span>${pName}`;
                molDiv.appendChild(titleDiv);

                const canvasWrapper = document.createElement('div');
                canvasWrapper.style.cssText = 'width: 100%; height: 100%; position: absolute;';
                molDiv.appendChild(canvasWrapper);

                const initViewer = () => {
                    if (canvasWrapper.clientWidth === 0 || canvasWrapper.clientHeight === 0) {
                        setTimeout(initViewer, 20);
                        return;
                    }
                    const v = $3Dmol.createViewer(canvasWrapper, { backgroundColor: bgMol });
                    canvasWrapper.addEventListener('mousedown', () => { window.activeMolViewer = v; });
                    canvasWrapper.addEventListener('wheel', () => { window.activeMolViewer = v; }, { passive: true });
                    if (loadedModels[idx]) {
                        v.addModel(loadedModels[idx].text, loadedModels[idx].format);
                        v.setStyle({}, { cartoon: { colorfunc: (atom) => get3DColor(atom.chain, idx) } });
                        v.zoomTo();
                    }
                    v.pairKey = pairKey; v.protIdxs = [idx]; viewers.push(v);
                    update3DViewerGridForPair(pairKey);
                };
                initViewer();
            });
        } else {
            const molDiv = document.createElement('div');
            molDiv.className = 'grid-mol-container';
            molDiv.style.cssText = 'width: 100%; height: 100%; position: absolute;';
            refs.botBody.appendChild(molDiv);

            const initViewer = () => {
                if (molDiv.clientWidth === 0 || molDiv.clientHeight === 0) {
                    setTimeout(initViewer, 20);
                    return;
                }
                const v = $3Dmol.createViewer(molDiv, { backgroundColor: bgMol });
                molDiv.addEventListener('mousedown', () => { window.activeMolViewer = v; });
                molDiv.addEventListener('wheel', () => { window.activeMolViewer = v; }, { passive: true });
                let hasMol = false;
                [p1_idx, p2_idx].forEach(idx => {
                    if (loadedModels[idx]) {
                        const mObj = v.addModel(loadedModels[idx].text, loadedModels[idx].format);
                        v.setStyle({ model: mObj.getID() }, { cartoon: { colorfunc: (atom) => get3DColor(atom.chain, idx) } });
                        hasMol = true;
                    }
                });
                v.pairKey = pairKey; v.protIdxs = [p1_idx, p2_idx];
                if (hasMol) v.zoomTo();
                viewers.push(v);
                update3DViewerGridForPair(pairKey);
            };
            initViewer();
        }

        molLoadedKeys.push(pairKey);
        const btn = refs.botCard && refs.botCard.querySelector('.grid-toggle-mol');
        if (btn) btn.textContent = 'Unload 3D';
        updateLoadedCounter();
    } catch (e) {
        logError(`loadMolForPair(${pairKey}) crashed`, e);
    }
}

function unloadMolForPair(pairKey) {
    const refs = cardRefs[pairKey];
    /* Remove viewers associated with this pair */
    const remaining = [];
    viewers.forEach(v => {
        if (v.pairKey === pairKey) {
            try { if (v.removeAllModels) v.removeAllModels(); } catch (e) {}
            /* No clean .destroy() on 3Dmol viewers; we drop the canvas
             * via innerHTML below and rely on GC to release the WebGL ctx. */
        } else {
            remaining.push(v);
        }
    });
    viewers = remaining;
    if (window.activeMolViewer && window.activeMolViewer.pairKey === pairKey) {
        window.activeMolViewer = null;
    }
    molLoadedKeys = molLoadedKeys.filter(k => k !== pairKey);

    if (refs && refs.botBody) {
        const placeholder = refs.botBody.querySelector('.grid-placeholder');
        refs.botBody.innerHTML = '';
        if (placeholder) {
            placeholder.style.display = 'flex';
            refs.botBody.appendChild(placeholder);
        }
    }
    if (refs && refs.botCard) {
        const btn = refs.botCard.querySelector('.grid-toggle-mol');
        if (btn) btn.textContent = 'Load 3D';
    }
    updateLoadedCounter();
}

window.toggleMolForPair = function(pairKey, btn) {
    if (molLoadedKeys.includes(pairKey)) unloadMolForPair(pairKey);
    else loadMolForPair(pairKey);
};

/* ========================================================================= */
/* IntersectionObservers (auto-load on scroll into view)                     */
/* ========================================================================= */

function makeObserver(loader) {
    if (typeof IntersectionObserver === 'undefined') return null;
    return new IntersectionObserver(entries => {
        entries.forEach(e => {
            if (!e.isIntersecting) return;
            const k = e.target.dataset.pair;
            if (k) loader(k);
        });
    }, { rootMargin: '200px' });   /* preload slightly before visible */
}

/* ========================================================================= */
/* Public entry points (signatures preserved)                                */
/* ========================================================================= */

function buildGridTop() {
    try {
        const topPanel = document.getElementById('top-panel');
        const singleNet = document.getElementById('single-network');
        if (singleNet) singleNet.style.display = 'none';

        let gridWrapper = document.getElementById('dynamic-grid-networks');
        if (!gridWrapper) {
            gridWrapper = document.createElement('div');
            gridWrapper.id = 'dynamic-grid-networks';
            gridWrapper.style.width = '100%';
            gridWrapper.style.height = '100%';
            gridWrapper.style.padding = '15px';
            gridWrapper.style.boxSizing = 'border-box';
            gridWrapper.style.overflowY = 'auto';
            gridWrapper.style.background = 'var(--bg-body)';
            topPanel.appendChild(gridWrapper);
        }

        /* Tear down previous state */
        gridNetworks.forEach(g => { try { g.network.destroy(); } catch (e) {} });
        gridNetworks = []; gridNodesDatasets = {}; gridEdgesDatasets = {};
        visLoadedKeys = [];
        if (visObserver) visObserver.disconnect();

        /* Reset card refs (the bottom builder also writes here) */
        Object.keys(cardRefs).forEach(k => { cardRefs[k].topCard = null; cardRefs[k].topBody = null; });

        gridWrapper.innerHTML = '';
        gridWrapper.style.display = 'grid';

        gridStats = computeGridStats();
        lazyMode = decideLazyMode(gridStats);

        const numPairs = gridStats.pairCount;
        const cols = numPairs > 4 ? 3 : (numPairs > 1 ? 2 : 1);
        gridWrapper.style.gridTemplateColumns = `repeat(${cols}, 1fr)`;
        gridWrapper.style.gridAutoRows = '350px';
        gridWrapper.style.gap = '15px';

        if (lazyMode) buildHeavyWarningBanner(gridWrapper);

        visObserver = lazyMode ? makeObserver(loadVisForPair) : null;

        Object.keys(masterData.pairs).forEach(pairKey => {
            const pData = masterData.pairs[pairKey];
            const { card, body } = makeTopCard(pairKey, pData);
            gridWrapper.appendChild(card);
            cardRefs[pairKey] = Object.assign(cardRefs[pairKey] || {}, { topCard: card, topBody: body });

            if (!lazyMode) {
                /* Old behavior for small datasets: just load it */
                setTimeout(() => loadVisForPair(pairKey), 10);
            } else if (visObserver) {
                visObserver.observe(card);
            }
        });

        updateLoadedCounter();
    } catch (e) { logError('buildGridTop crashed', e); }
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

        molLoadedKeys = [];
        if (molObserver) molObserver.disconnect();

        Object.keys(cardRefs).forEach(k => { cardRefs[k].botCard = null; cardRefs[k].botBody = null; });

        gridWrapper.innerHTML = '';
        gridWrapper.style.display = 'grid';

        /* Make sure stats reflect the current state (may differ from top if top wasn't built yet) */
        if (gridStats.pairCount === 0) gridStats = computeGridStats();
        lazyMode = decideLazyMode(gridStats);

        const numPairs = gridStats.pairCount;
        const cols = numPairs > 4 ? 3 : (numPairs > 1 ? 2 : 1);
        gridWrapper.style.gridTemplateColumns = `repeat(${cols}, 1fr)`;
        gridWrapper.style.gridAutoRows = '350px';
        gridWrapper.style.gap = '15px';

        molObserver = lazyMode ? makeObserver(loadMolForPair) : null;

        Object.keys(masterData.pairs).forEach(pairKey => {
            const pData = masterData.pairs[pairKey];
            const { card, body } = makeBotCard(pairKey, pData);
            gridWrapper.appendChild(card);
            cardRefs[pairKey] = Object.assign(cardRefs[pairKey] || {}, { botCard: card, botBody: body });

            if (!lazyMode) {
                setTimeout(() => loadMolForPair(pairKey), 10);
            } else if (molObserver) {
                molObserver.observe(card);
            }
        });

        updateLoadedCounter();
    } catch (e) { logError('buildGridBottom crashed', e); }
}
