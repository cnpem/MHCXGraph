let fabCanvas = null;
let canvasW = 1920;
let canvasH = 1080;

// --- UNDO/REDO HISTORY ENGINE ---
let undoStack = [];
let redoStack = [];
let isHistoryAction = false;

// --- PREVIEW PAN / PALETTE STATE ---
let isPreviewPanning = false;
let lastPreviewPanX = 0;
let lastPreviewPanY = 0;
let isApplyingExportPalette = false;
let isSpaceDown = false;

let exportOriginalPalette = null;

function initHistory() {
    undoStack = [];
    redoStack = [];
    isHistoryAction = false;
    saveHistory();
}

function saveHistory() {
    if (isHistoryAction || !fabCanvas) return;
    redoStack = [];
    if (undoStack.length > 30) undoStack.shift();
    undoStack.push(JSON.stringify(fabCanvas.toJSON(['customId', 'customType', 'isSystemAsset', 'customLabel'])));
}

function undo() {
    if (undoStack.length > 1) {
        isHistoryAction = true;
        redoStack.push(undoStack.pop());
        let state = undoStack[undoStack.length - 1];
        fabCanvas.loadFromJSON(state, function() {
            fabCanvas.renderAll();
            renderElementsList();
            updateUIForSelection();
            isHistoryAction = false;
        });
    }
}

function redo() {
    if (redoStack.length > 0) {
        isHistoryAction = true;
        let state = redoStack.pop();
        undoStack.push(state);
        fabCanvas.loadFromJSON(state, function() {
            fabCanvas.renderAll();
            renderElementsList();
            updateUIForSelection();
            isHistoryAction = false;
        });
    }
}

// --- HELPERS ---
function getCurrentThemeBgColor() {
    const css = getComputedStyle(document.body);
    return css.getPropertyValue('--bg-panel').trim() || '#111827';
}

function getCurrentThemeTextColor() {
    const css = getComputedStyle(document.body);
    return css.getPropertyValue('--text-main').trim() || '#f9fafb';
}

function logExportDebug(msg, data) {
    if (typeof logDebug === 'function') logDebug(`[EXPORT] ${msg}`, data);
    else console.log(`[EXPORT] ${msg}`, data);
}

function getObjectDebugInfo(obj) {
    if (!obj) return null;
    return {
        customId: obj.customId,
        type: obj.type,
        width: obj.width,
        height: obj.height,
        scaleX: obj.scaleX,
        scaleY: obj.scaleY,
        scaledW: typeof obj.getScaledWidth === 'function' ? obj.getScaledWidth() : null,
        scaledH: typeof obj.getScaledHeight === 'function' ? obj.getScaledHeight() : null,
        left: obj.left,
        top: obj.top,
        visible: obj.visible
    };
}

function centerExportWorkspace() {
    const preview = document.getElementById('export-preview-area');
    const wrapper = document.getElementById('export-canvas-wrapper');
    if (!preview || !wrapper) return;

    const maxLeft = Math.max(0, wrapper.scrollWidth - preview.clientWidth);
    const maxTop = Math.max(0, wrapper.scrollHeight - preview.clientHeight);

    preview.scrollLeft = maxLeft / 2;
    preview.scrollTop = maxTop / 2;
}

function stopPreviewPan() {
    isPreviewPanning = false;
    const preview = document.getElementById('export-preview-area');
    if (preview) preview.style.cursor = 'default';
    if (fabCanvas) {
        fabCanvas.defaultCursor = 'default';
        fabCanvas.hoverCursor = 'move';
    }
}

// --- DYNAMIC IMAGE GENERATORS ---

async function cropTransparent(image) {
    try {
        const canvas = document.createElement('canvas');
        canvas.width = image.width || 500;
        canvas.height = image.height || 500;
        const ctx = canvas.getContext('2d');
        ctx.drawImage(image, 0, 0);

        const pixels = ctx.getImageData(0, 0, canvas.width, canvas.height);
        let bounds = { top: null, left: null, right: null, bottom: null };

        for (let i = 0; i < pixels.data.length; i += 4) {
            if (pixels.data[i + 3] !== 0) {
                const x = (i / 4) % canvas.width;
                const y = Math.floor((i / 4) / canvas.width);
                if (bounds.top === null) bounds.top = y;
                if (bounds.left === null || x < bounds.left) bounds.left = x;
                if (bounds.right === null || x > bounds.right) bounds.right = x;
                if (bounds.bottom === null || y > bounds.bottom) bounds.bottom = y;
            }
        }
        if (bounds.top === null) return image.src;

        const pad = 20;
        const cropX = Math.max(0, bounds.left - pad);
        const cropY = Math.max(0, bounds.top - pad);
        const cropW = Math.max(10, (bounds.right - bounds.left) + pad * 2);
        const cropH = Math.max(10, (bounds.bottom - bounds.top) + pad * 2);

        const croppedCanvas = document.createElement('canvas');
        croppedCanvas.width = cropW;
        croppedCanvas.height = cropH;
        const croppedCtx = croppedCanvas.getContext('2d');

        croppedCtx.drawImage(
            canvas,
            cropX, cropY, cropW, cropH,
            0, 0, cropW, cropH
        );

        return croppedCanvas.toDataURL('image/png');
    } catch(e) {
        logError("Crop failure", e);
        return image.src;
    }
}

function loadFabricImage(url, options) {
    return new Promise(resolve => {
        fabric.Image.fromURL(url, function(img) {
            if (img) {
                img.set(options);
                resolve(img);
            } else {
                resolve(null);
            }
        });
    });
}

function normalizeFabricImage(img) {

    if (!img) return;

    const w = img.width || img._element?.naturalWidth;
    const h = img.height || img._element?.naturalHeight;

    if (!w || !h) return;

    img.scaleX = 1;
    img.scaleY = 1;

    img.width = w;
    img.height = h;

    img.setCoords();
}
async function generateMetadataCard(textColor, bgColor, borderColor) {
    const c = document.createElement('canvas');
    const ctx = c.getContext('2d');
    c.width = 450;
    c.height = 250;

    ctx.fillStyle = bgColor;
    ctx.fillRect(0, 0, c.width, c.height);

    ctx.strokeStyle = borderColor;
    ctx.lineWidth = 4;
    ctx.strokeRect(2, 2, c.width - 4, c.height - 4);

    ctx.fillStyle = textColor;
    ctx.font = 'bold 24px sans-serif';
    ctx.textAlign = 'center';
    ctx.fillText("Dataset Metadata", c.width / 2, 40);

    ctx.font = '20px sans-serif';
    ctx.textAlign = 'left';

    let y = 90;
    let md = (typeof masterData !== 'undefined' && masterData.metadata) ? masterData.metadata : {};
    let lines = [
        `• Mode: ${typeof masterData !== 'undefined' ? masterData.mode : 'N/A'}`,
        `• Granularity: ${md.node_granularity || 'N/A'}`,
        `• Edge Threshold: ${md.edge_threshold || 'N/A'} Å`,
        `• Distance Diff: ${md.global_distance_diff_threshold || 'N/A'} Å`
    ];

    lines.forEach(l => {
        ctx.fillText(l, 30, y);
        y += 35;
    });

    return c.toDataURL('image/png');
}

// --- 1. INITIALIZE STUDIO ---
async function openExportStudio() {
    document.getElementById('export-modal').style.display = 'flex';

    const isDark = document.body.getAttribute('data-theme') === 'dark';
    const themeBg = getCurrentThemeBgColor();
    const themeText = getCurrentThemeTextColor();

    if (typeof activePaletteName !== 'undefined') window.originalPalette = activePaletteName;
    window.originalThemeText = themeText;

    document.getElementById('exp-bg-color').value = themeBg;
    document.getElementById('exp-data-color').value = themeText;
    document.getElementById('exp-transparent-bg').checked = false;
    document.getElementById('exp-w').value = 1920;
    document.getElementById('exp-h').value = 1080;

    canvasW = 1920;
    canvasH = 1080;

    if (fabCanvas) {
        fabCanvas.dispose();
    }

    fabCanvas = new fabric.Canvas('export-canvas', {
        width: canvasW,
        height: canvasH,
        preserveObjectStacking: true,
        backgroundColor: themeBg
    });

    if (fabCanvas.wrapperEl) {
        fabCanvas.wrapperEl.style.transformOrigin = 'top left';
        fabCanvas.wrapperEl.style.boxShadow = '0 10px 25px rgba(0,0,0,0.3)';
    }

    setupFabricInteractions();
    setupGlobalKeyboardShortcuts();

    if (typeof activePaletteName !== 'undefined') {
        let palDropdown = document.getElementById('exp-palette');
        if (palDropdown) palDropdown.value = activePaletteName;
    }

    let labelColor = themeText;
    await new Promise(r => setTimeout(r, 100));

    let loadPromises = [];

    // Capture Graph + Subtitle
    const graphCanvasUI = document.querySelector('#top-panel canvas');

    let gLabelStr = 'Associated Graph';

    if (typeof masterData !== 'undefined' && masterData.mode === 'pair') {
        let pSel = document.getElementById('pair-selector');
        if (pSel && pSel.value !== 'GRID' && pSel.value !== 'ANALYSIS') {
            gLabelStr = `Graph: ${pSel.value}`;
        }
    } else if (typeof graphData !== 'undefined' && graphData && graphData.proteins) {
        gLabelStr = `Associated Graph: ${graphData.proteins.join(' x ')}`;
    }
    if (graphCanvasUI) {
        let img = new Image();
        img.src = graphCanvasUI.toDataURL('image/png');
        await new Promise(r => img.onload = r);
        let croppedUrl = await cropTransparent(img);

        let p = loadFabricImage(croppedUrl, {
            customId: 'graph',
            customType: 'image',
            customLabel: gLabelStr,
            isSystemAsset: true,
            originX: 'center',
            originY: 'center'
        });
        loadPromises.push(p);

        let gLabel = new fabric.Textbox(gLabelStr, {
            customId: 'graph_label',
            customType: 'subtitle',
            isSystemAsset: true,
            fontSize: 28,
            fontWeight: 'bold',
            fill: labelColor,
            originX: 'center',
            originY: 'center',
            textAlign: 'center',
            width: 600
        });
        loadPromises.push(Promise.resolve(gLabel));
    }

    // Capture Plotly Heatmap + Subtitle
    let heatmapUrl = await generatePlotlyHeatmap();
    if (heatmapUrl) {
        let p = loadFabricImage(heatmapUrl, {
            customId: 'heatmap',
            customType: 'image',
            isSystemAsset: true,
            customLabel: 'Co-occurrence Heatmap',
            originX: 'center',
            originY: 'center'
        });
        loadPromises.push(p);

        let hmLabel = new fabric.Textbox("Co-occurrence Heatmap (%)", {
            customId: 'heatmap_label',
            customType: 'subtitle',
            isSystemAsset: true,
            fontSize: 28,
            fontWeight: 'bold',
            fill: labelColor,
            originX: 'center',
            originY: 'center',
            textAlign: 'center',
            width: 600
        });
        loadPromises.push(Promise.resolve(hmLabel));
    }

    // Capture 3DMol Viewers + Subtitles
    if (typeof viewers !== 'undefined' && viewers.length > 0) {
        for (let i = 0; i < viewers.length; i++) {
            let v = viewers[i];
            v.setBackgroundColor(0x000000, 0.0);
            v.render();

            let img = new Image();
            img.src = v.pngURI();
            await new Promise(r => img.onload = r);
            let croppedUrl = await cropTransparent(img);

            v.setBackgroundColor(isDark ? "#1f2937" : "#ffffff", 1.0);
            v.render();

            let pName = (v.layoutMode === 'aligned' && typeof graphData !== 'undefined')
                ? `Aligned Mols`
                : (v.protName || v.pairKey || `3D View ${i + 1}`);


            let p = loadFabricImage(croppedUrl, {
                customId: `mol_${i}`,
                customType: 'image',
                customLabel: pName,
                isSystemAsset: true,
                originX: 'center',
                originY: 'center'
            });

            loadPromises.push(p);
            let molLabel = new fabric.Textbox(pName, {
                customId: `mol_${i}_label`,
                customType: 'subtitle',
                isSystemAsset: true,
                fontSize: 28,
                fontWeight: 'bold',
                fill: labelColor,
                originX: 'center',
                originY: 'center',
                textAlign: 'center',
                width: 500
            });
            loadPromises.push(Promise.resolve(molLabel));
        }
    }

    // Generate Metadata Card
    let metaUrl = await generateMetadataCard(
        labelColor,
        isDark ? '#1f2937' : '#f3f4f6',
        isDark ? '#374151' : '#d1d5db'
    );

    let metaP = loadFabricImage(metaUrl, {
        customId: 'metadata',
        customType: 'image',
        customLabel: 'Dataset Metadata',
        isSystemAsset: true,
        originX: 'center',
        originY: 'center'
    });
    loadPromises.push(metaP);

    let loadedAssets = await Promise.all(loadPromises);
    logExportDebug('openExportStudio loaded assets', loadedAssets.map(getObjectDebugInfo));
    loadedAssets.forEach(asset => {

        if (!asset) return;

        if (asset.type === "image") {
            normalizeFabricImage(asset);
        }

        fabCanvas.add(asset);

    });
    // Main Title
    let runName = (typeof graphData !== 'undefined' && graphData && graphData.run_name)
        ? graphData.run_name
        : 'Output';

    let mainTitle = new fabric.Textbox(`MHCXGraph Analysis: ${runName}`, {
        customId: 'main_title',
        customType: 'title',
        isSystemAsset: false,
        left: fabCanvas.width / 2,
        top: 60,
        fontSize: 42,
        fontWeight: 'bold',
        fill: labelColor,
        originX: 'center',
        originY: 'center',
        textAlign: 'center',
        width: 1200
    });
    fabCanvas.add(mainTitle);

    document.getElementById('exp-layout').value = 'article_1';
    applyFabricSnapLayout();
    logExportDebug('layout applied after studio open');
    renderElementsList();

    const previewArea = document.getElementById('export-preview-area');
    previewArea.scrollTop = 0;
    previewArea.scrollLeft = 0;

    let zoom = Math.min((previewArea.clientWidth - 80) / 1920, 0.6);
    zoom = Math.max(zoom, 0.35);
    document.getElementById('exp-preview-zoom').value = zoom;
    updatePreviewZoom();

    requestAnimationFrame(() => {
        centerExportWorkspace();
    });
    if (typeof nodesDataset !== 'undefined' && nodesDataset) {
        const sampleIds = [0, 1, 2, 3, 4, 7];
        logExportDebug('dataset colors before graph capture', {
            sample: sampleIds.map(id => {
                const n = nodesDataset.get(id);
                return n ? {
                    id: n.id,
                    color: n.color
                } : null;
            }).filter(Boolean)
        });
    }
    initHistory();
}

// --- 2. PLOTLY HEATMAP ---
async function generatePlotlyHeatmap() {
    let isGlobalPairMode = (
        typeof masterData !== 'undefined' &&
        masterData.mode === 'pair' &&
        document.getElementById('pair-selector') &&
        document.getElementById('pair-selector').value === 'GRID'
    );

    const mData = isGlobalPairMode ? masterData : (typeof graphData !== 'undefined' ? graphData : null);
    if (!mData || !mData.proteins || mData.proteins.length === 0) return null;

    const labels = mData.proteins;
    const numDisplay = labels.length;
    let matrix = Array(numDisplay).fill(null).map(() => Array(numDisplay).fill(null));

    let origSizes = [];
    if (isGlobalPairMode) {
        origSizes = Array(numDisplay).fill(0);
        Object.values(mData.pairs).forEach(pData => {
            let p1 = labels.indexOf(pData.proteins[0]);
            let p2 = labels.indexOf(pData.proteins[1]);
            origSizes[p1] = pData.filtered_graphs[0].nodes.length;
            origSizes[p2] = pData.filtered_graphs[1].nodes.length;
        });
    } else {
        if (mData.filtered_graphs) origSizes = mData.filtered_graphs.map(fg => fg.nodes.length);
        else origSizes = Array(numDisplay).fill(0);
    }

    if (isGlobalPairMode) {
        Object.keys(mData.pairs).forEach(pairKey => {
            let pData = mData.pairs[pairKey];
            let p1 = labels.indexOf(pData.proteins[0]);
            let p2 = labels.indexOf(pData.proteins[1]);
            let unique0 = new Set();
            let unique1 = new Set();

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
            let ratio = origSum > 0 ? (assocNodesCount / origSum) : 0;
            matrix[p1][p2] = parseFloat((ratio * 100).toFixed(1));
            matrix[p2][p1] = matrix[p1][p2];
        });
    } else {
        for (let i = 0; i < numDisplay; i++) {
            for (let j = 0; j < numDisplay; j++) {
                let unique_i = new Set();
                let unique_j = new Set();
                let global_i = masterData.proteins.indexOf(labels[i]);
                let global_j = masterData.proteins.indexOf(labels[j]);

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

                let assocNodesCount = unique_i.size + unique_j.size;
                let origSum = (i === j) ? origSizes[i] : (origSizes[i] + origSizes[j]);
                let ratio = origSum > 0 ? (assocNodesCount / origSum) : 0;
                matrix[i][j] = parseFloat((ratio * 100).toFixed(1));
            }
        }
    }

    const isDark = document.body.getAttribute('data-theme') === 'dark';
    const textColor = document.getElementById('exp-data-color') ? document.getElementById('exp-data-color').value : '#f9fafb';
    const colorscaleElem = document.getElementById('exp-hm-colorscale');
    const colorscale = colorscaleElem ? colorscaleElem.value : 'Viridis';
    const annotElem = document.getElementById('exp-hm-annot');
    const showAnnot = annotElem ? annotElem.checked : true;

    const axisSize = parseInt(document.getElementById('exp-hm-axis-size')?.value || 18);
    const valSize = parseInt(document.getElementById('exp-hm-val-size')?.value || 16);

    const axisBold = document.getElementById('exp-hm-axis-bold') ? document.getElementById('exp-hm-axis-bold').checked : true;
    const valBold = document.getElementById('exp-hm-val-bold') ? document.getElementById('exp-hm-val-bold').checked : false;
    
    const axisWeight = axisBold ? 'bold' : 'normal';

    let maxVal = 0;
    for (let i = 0; i < numDisplay; i++) {
        for (let j = 0; j < numDisplay; j++) {
            if (matrix[i][j] > maxVal) maxVal = matrix[i][j];
        }
    }
    let midPoint = maxVal > 0 ? maxVal / 2 : 50;

    let textColorsMatrix = [];
    const darkHigh = ['Blues', 'Reds', 'YlGnBu'].includes(colorscale);

    for (let i = 0; i < numDisplay; i++) {
        let rowColors = [];
        for (let j = 0; j < numDisplay; j++) {
            let val = matrix[i][j];
            if (val === null || val === undefined) {
                rowColors.push(textColor);
            } else {
                if (darkHigh) {
                    rowColors.push(val > midPoint ? 'white' : 'black');
                } else {
                    // Viridis, Plasma, Inferno, Magma: HIGH values are bright
                    rowColors.push(val > midPoint ? 'black' : 'white');
                }
            }
        }
        textColorsMatrix.push(rowColors);
    }

    let cellTextTemplate = "";
    if (showAnnot) {
        cellTextTemplate = valBold ? "<b>%{z}%</b>" : "%{z}%";
    }
    var data = [{
        z: matrix,
        x: labels,
        y: labels,
        type: 'heatmap',
        colorscale: colorscale,
        showscale: true,
        texttemplate: cellTextTemplate,
        textfont: { color: textColorsMatrix, size: valSize },
        colorbar: { 
            tickfont: { color: textColor, size: axisSize, weight: axisWeight }
        }
    }];

    var layout = {
        plot_bgcolor: 'transparent',
        paper_bgcolor: 'transparent',
        xaxis: { 
            tickfont: { color: textColor, size: axisSize, weight: axisWeight },
            tickangle: -45,
            automargin: true
        },
        yaxis: { 
            tickfont: { color: textColor, size: axisSize, weight: axisWeight }, 
            autorange: 'reversed',
            automargin: true
        },
        width: 800,
        height: 800,
        // MARGENS AUMENTADAS: Para garantir que os nomes das proteínas caibam na imagem
        margin: { t: 20, b: 40, l: 40, r: 20 } 
    };
    const hiddenDiv = document.getElementById('plotly-hidden-container');
    await Plotly.newPlot(hiddenDiv, data, layout);
    let dataUrl = await Plotly.toImage(hiddenDiv, { format: 'png', width: 800, height: 800 });
    Plotly.purge(hiddenDiv);
    return dataUrl;
}

// --- 3. SCROLLING, ZOOMING & DRAWING ---
let currentDrawMode = 'none';
let drawStartX, drawStartY, tempShape;

function setDrawMode(mode) {
    currentDrawMode = mode;
    fabCanvas.isDrawingMode = (mode === 'brush');
    document.querySelectorAll('.tool-btn').forEach(b => b.classList.remove('active'));

    let color = document.getElementById('exp-element-color').value || '#000000';

    if (mode === 'brush') {
        document.getElementById('btn-draw-brush').classList.add('active');
        fabCanvas.freeDrawingBrush.color = color;
        fabCanvas.freeDrawingBrush.width = 5;
        fabCanvas.selection = false;
        fabCanvas.forEachObject(o => o.selectable = false);
        fabCanvas.defaultCursor = 'crosshair';
        fabCanvas.hoverCursor = 'crosshair';
    } else if (mode === 'rect' || mode === 'circle') {
        document.getElementById('btn-draw-' + mode).classList.add('active');
        fabCanvas.selection = false;
        fabCanvas.forEachObject(o => o.selectable = false);
        fabCanvas.defaultCursor = 'crosshair';
        fabCanvas.hoverCursor = 'crosshair';
    } else {
        fabCanvas.selection = true;
        fabCanvas.forEachObject(o => o.selectable = true);
        fabCanvas.defaultCursor = 'default';
        fabCanvas.hoverCursor = 'move';
    }
}

function updatePreviewZoom() {
    const zoom = parseFloat(document.getElementById('exp-preview-zoom').value);
    document.getElementById('exp-preview-zoom-val').innerText = Math.round(zoom * 100) + '%';

    const wrapper = document.getElementById('export-canvas-wrapper');
    if (!wrapper || !fabCanvas) return;

    wrapper.style.width = (canvasW * zoom) + 'px';
    wrapper.style.height = (canvasH * zoom) + 'px';

    if (fabCanvas.wrapperEl) {
        fabCanvas.wrapperEl.style.transformOrigin = 'top left';
        fabCanvas.wrapperEl.style.transform = `scale(${zoom})`;
    }

    fabCanvas.renderAll();
}

function setupGlobalKeyboardShortcuts() {
    if (!window.hasExportShortcutsBound) {
        window.addEventListener('keydown', function(e) {
            if (document.getElementById('export-modal').style.display === 'none') return;

            if (e.key === 'Delete' || e.key === 'Backspace') {
                if (e.target.tagName !== 'INPUT' && e.target.tagName !== 'TEXTAREA') {
                    const active = fabCanvas.getActiveObject();
                    if (active && !active.isEditing) deleteSelectedFabricAsset();
                }
            }
            if (e.code === 'Space') {
                if (e.target.tagName === 'INPUT' || e.target.tagName === 'TEXTAREA') return;
                e.preventDefault(); 
                if (!isSpaceDown) {
                    isSpaceDown = true;
                    const previewArea = document.getElementById('export-preview-area');
                    previewArea.style.cursor = 'grab';
                    if (fabCanvas) {
                        fabCanvas.defaultCursor = 'grab';
                        fabCanvas.selection = false;
                        // Temporarily disable hover/selection for all objects
                        fabCanvas.getObjects().forEach(o => {
                            o.tempSelectable = o.selectable;
                            o.selectable = false;
                            o.evented = false; 
                        });
                    }
                }
                return;
            }
            if (e.ctrlKey && e.key.toLowerCase() === 'a') {
                            e.preventDefault();
                            if (fabCanvas) {
                                fabCanvas.discardActiveObject();
                                // Grab everything that isn't locked or invisible
                                const allObjs = fabCanvas.getObjects().filter(o => o.visible && !o.locked);
                                if (allObjs.length) {
                                    const sel = new fabric.ActiveSelection(allObjs, { canvas: fabCanvas });
                                    fabCanvas.setActiveObject(sel);
                                    fabCanvas.requestRenderAll();
                                    if (typeof updateUIForSelection === 'function') updateUIForSelection();
                                }
                            }
                            return;
            }
            if (e.ctrlKey && e.shiftKey && e.key.toLowerCase() === 'z') {
                e.preventDefault();
                redo();
            } else if (e.ctrlKey && !e.shiftKey && e.key.toLowerCase() === 'z') {
                e.preventDefault();
                undo();
            }
        });
        window.addEventListener('keyup', function(e) {
            if (e.code === 'Space') {
                isSpaceDown = false;
                const previewArea = document.getElementById('export-preview-area');
                previewArea.style.cursor = 'default';
                if (fabCanvas) {
                    fabCanvas.defaultCursor = 'default';
                    fabCanvas.selection = true;
                    // Restore objects to their original selectable state
                    fabCanvas.getObjects().forEach(o => {
                        if (o.tempSelectable !== undefined) {
                            o.selectable = o.tempSelectable;
                            o.evented = true;
                        }
                    });
                }
                if (isPreviewPanning) stopPreviewPan();
            }
        });
        window.addEventListener('mouseup', stopPreviewPan);
        window.hasExportShortcutsBound = true;
    }
}

function setupFabricInteractions() {
    const previewArea = document.getElementById('export-preview-area');

    if (!window.hasExportWheelBound) {
        previewArea.addEventListener('wheel', (e) => {
            if (document.getElementById('export-modal').style.display === 'none') return;

            if (e.ctrlKey) {
                e.preventDefault();
                let zoomSlider = document.getElementById('exp-preview-zoom');
                let currentZoom = parseFloat(zoomSlider.value);
                let zoomChange = e.deltaY > 0 ? -0.05 : 0.05;
                zoomSlider.value = Math.min(Math.max(currentZoom + zoomChange, 0.1), 1.5);
                updatePreviewZoom();
                return;
            }

            e.preventDefault();

            if (e.shiftKey) {
                previewArea.scrollLeft += (Math.abs(e.deltaY) > Math.abs(e.deltaX) ? e.deltaY : e.deltaX);
            } else {
                previewArea.scrollTop += e.deltaY;
                if (e.deltaX) previewArea.scrollLeft += e.deltaX;
            }
        }, { passive: false });

        window.hasExportWheelBound = true;
    }

    fabCanvas.on('mouse:down', function(opt) {
        var evt = opt.e;

        
        if (isSpaceDown) {
            isPreviewPanning = true;
            lastPreviewPanX = evt.clientX;
            lastPreviewPanY = evt.clientY;
            
            const previewArea = document.getElementById('export-preview-area');
            previewArea.style.cursor = 'grabbing';
            fabCanvas.defaultCursor = 'grabbing';
            return;
        }

        // 🔴 FIX: Allow normal selection boxes! If clicking blank canvas, do nothing and let Fabric handle it.
        if (currentDrawMode === 'none' && !opt.target) {
            return; 
        }

        if (currentDrawMode === 'none' && !opt.target) {
            isPreviewPanning = true;
            lastPreviewPanX = evt.clientX;
            lastPreviewPanY = evt.clientY;
            previewArea.style.cursor = 'grabbing';
            fabCanvas.defaultCursor = 'grabbing';
            fabCanvas.hoverCursor = 'grabbing';
            return;
        }

        if (currentDrawMode === 'rect' || currentDrawMode === 'circle') {
            var pointer = fabCanvas.getPointer(evt);
            drawStartX = pointer.x;
            drawStartY = pointer.y;
            let color = document.getElementById('exp-element-color').value || '#000000';

            if (currentDrawMode === 'rect') {
                tempShape = new fabric.Rect({
                    left: drawStartX,
                    top: drawStartY,
                    originX: 'left',
                    originY: 'top',
                    width: 0,
                    height: 0,
                    fill: 'transparent',
                    stroke: color,
                    strokeWidth: 4
                });
            } else {
                tempShape = new fabric.Circle({
                    left: drawStartX,
                    top: drawStartY,
                    originX: 'left',
                    originY: 'top',
                    radius: 0,
                    fill: 'transparent',
                    stroke: color,
                    strokeWidth: 4
                });
            }

            fabCanvas.add(tempShape);
        }
    });

    fabCanvas.on('mouse:move', function(opt) {
        if (isPreviewPanning) {
            const evt = opt.e;
            const dx = evt.clientX - lastPreviewPanX;
            const dy = evt.clientY - lastPreviewPanY;

            previewArea.scrollLeft -= dx;
            previewArea.scrollTop -= dy;

            lastPreviewPanX = evt.clientX;
            lastPreviewPanY = evt.clientY;
            return;
        }

        if ((currentDrawMode === 'rect' || currentDrawMode === 'circle') && tempShape) {
            var pointer = fabCanvas.getPointer(opt.e);

            if (currentDrawMode === 'rect') {
                tempShape.set({
                    width: Math.abs(pointer.x - drawStartX),
                    height: Math.abs(pointer.y - drawStartY)
                });
                if (pointer.x < drawStartX) tempShape.set({ left: pointer.x });
                if (pointer.y < drawStartY) tempShape.set({ top: pointer.y });
            } else {
                let radius = Math.max(Math.abs(pointer.x - drawStartX), Math.abs(pointer.y - drawStartY)) / 2;
                tempShape.set({ radius: radius });
                if (pointer.x < drawStartX) tempShape.set({ left: pointer.x });
                if (pointer.y < drawStartY) tempShape.set({ top: pointer.y });
            }

            fabCanvas.renderAll();
        }
    });

    fabCanvas.on('mouse:up', function() {
        if (isPreviewPanning) {
            stopPreviewPan();
            return;
        }

        if ((currentDrawMode === 'rect' || currentDrawMode === 'circle') && tempShape) {
            tempShape.set({
                customId: currentDrawMode + '_' + Date.now(),
                customLabel: currentDrawMode.toUpperCase() + ' Shape'
            });
            tempShape.setCoords();
            tempShape = null;
            setDrawMode('none');
            saveHistory();
            renderElementsList();
        }
    });

    fabCanvas.on('object:modified', () => {
        document.getElementById('exp-layout').value = 'none';
        saveHistory();
    });

    fabCanvas.on('object:added', () => {
        if (!isHistoryAction) saveHistory();
    });

    fabCanvas.on('object:removed', () => {
        if (!isHistoryAction) saveHistory();
    });

    fabCanvas.on('path:created', function(opt) {
        opt.path.set({
            customId: 'path_' + Date.now(),
            customLabel: 'Brush Stroke'
        });
        setDrawMode('none');
        renderElementsList();
        saveHistory();
    });

    fabCanvas.on('selection:created', updateUIForSelection);
    fabCanvas.on('selection:updated', updateUIForSelection);
    fabCanvas.on('selection:cleared', () => {
        document.getElementById('exp-element-color').disabled = true;
        document.getElementById('btn-layer-up').disabled = true;
        document.getElementById('btn-layer-down').disabled = true;
        document.getElementById('exp-delete-btn').disabled = true;
    });
}

// --- 4. EDITOR TOOLS ---
function addFabricText() {
    let color = document.getElementById('exp-element-color').value || '#000000';
    let text = new fabric.Textbox('Double Click to Edit', {
        customId: 'text_' + Date.now(),
        customLabel: 'Text Block',
        isSystemAsset: false,
        left: fabCanvas.width / 2,
        top: fabCanvas.height / 2,
        fontSize: 48,
        fontWeight: 'bold',
        fill: color,
        originX: 'center',
        originY: 'center',
        width: 400
    });
    fabCanvas.add(text);
    fabCanvas.setActiveObject(text);
    renderElementsList();
}

// --- 5. ELEMENTS LIST & LAYERS ---
function renderElementsList() {
    if (!fabCanvas) return;
    const listDiv = document.getElementById('exp-elements-list');
    listDiv.innerHTML = '';

    let objs = [...fabCanvas.getObjects()].reverse();
    objs.forEach((obj) => {
        if (obj.customId === 'main_title') return;

        let lbl = document.createElement('label');
        lbl.className = 'export-cb-label';

        let displayTxt = obj.customLabel || (obj.type === 'path' ? 'Brush Stroke' : obj.type);
        if (obj.type === 'textbox' || obj.type === 'text') {
            displayTxt = `Text: "${(obj.text || '').substring(0, 20)}..."`;
        }

        lbl.innerHTML = `<input type="checkbox" ${obj.visible ? 'checked' : ''} onchange="toggleFabricAssetVisibility('${obj.customId || obj.type}', this.checked, ${fabCanvas.getObjects().indexOf(obj)})"> ${displayTxt}`;
        listDiv.appendChild(lbl);
    });
}

function toggleFabricAssetVisibility(id, isVisible, index) {
    let obj = fabCanvas.item(index);
    if (obj) {
        obj.set('visible', isVisible);
        if (obj.customType === 'image') {
            let sub = fabCanvas.getObjects().find(o => o.customId === obj.customId + '_label');
            if (sub) sub.set('visible', isVisible);
        }
        fabCanvas.renderAll();
        saveHistory();
    }
}

function moveLayer(dir) {
    const active = fabCanvas.getActiveObject();
    if (!active) return;

    if (dir === 'up') fabCanvas.bringForward(active);
    else if (dir === 'down') fabCanvas.sendBackwards(active);

    renderElementsList();
    saveHistory();
}

function updateUIForSelection() {
    const active = fabCanvas.getActiveObject();
    if (!active) return;

    let colorPicker = document.getElementById('exp-element-color');
    colorPicker.disabled = false;

    if (active.type === 'textbox' || active.type === 'text') colorPicker.value = active.fill;
    else if (active.type === 'path') colorPicker.value = active.stroke;
    else colorPicker.value = active.stroke || active.fill || '#000000';

    document.getElementById('exp-delete-btn').disabled = false;
    document.getElementById('btn-layer-up').disabled = false;
    document.getElementById('btn-layer-down').disabled = false;
}

function updateSelectedFabricColor() {
    const active = fabCanvas.getActiveObject();
    const newCol = document.getElementById('exp-element-color').value;

    if (fabCanvas.isDrawingMode) fabCanvas.freeDrawingBrush.color = newCol;
    if (!active) return;

    if (active.type === 'textbox' || active.type === 'text') active.set('fill', newCol);
    else if (active.type === 'path') active.set('stroke', newCol);
    else active.set('stroke', newCol);

    fabCanvas.renderAll();
    saveHistory();
}

function deleteSelectedFabricAsset() {
    const activeObjects = fabCanvas.getActiveObjects();
    if (activeObjects.length) {
        fabCanvas.discardActiveObject();
        activeObjects.forEach(function(object) {
            fabCanvas.remove(object);
        });
        renderElementsList();
        saveHistory();
    }
}

// --- 6. AUTO-LAYOUTS WITH SUBTITLES ---
function applyFabricSnapLayout() {
    const layout = document.getElementById('exp-layout').value;
    if (layout === 'none' || !fabCanvas) return;

    let objs = fabCanvas.getObjects();
    let sysImgs = objs.filter(o => o.isSystemAsset && o.visible && o.customType === 'image');
    if (sysImgs.length === 0) {
        logExportDebug('applyFabricSnapLayout aborted: no visible system images');
        return;
    }

    let graph = sysImgs.find(o => o.customId === 'graph');
    let heatmap = sysImgs.find(o => o.customId === 'heatmap');
    let meta = sysImgs.find(o => o.customId === 'metadata');
    let mols = sysImgs.filter(o => o.customId && o.customId.startsWith('mol_'));

    const gap = 60;
    const padding = 80;
    const headerHeight = 160;
    const subtitleGap = 15;
    const subtitleReserve = 60;

    const usableW = canvasW - padding * 2;
    const usableH = canvasH - padding * 2 - headerHeight;

    logExportDebug('applyFabricSnapLayout start', {
        layout,
        canvasW,
        canvasH,
        usableW,
        usableH,
        graph: getObjectDebugInfo(graph),
        heatmap: getObjectDebugInfo(heatmap),
        meta: getObjectDebugInfo(meta),
        mols: mols.map(getObjectDebugInfo)
    });

    const snapSubtitle = (imgObj) => {
        let sub = objs.find(o => o.customId === imgObj.customId + '_label');
        if (sub && sub.visible) {
            sub.set({
                left: imgObj.left,
                top: imgObj.top - (imgObj.getScaledHeight() / 2) - (sub.getScaledHeight() / 2) - subtitleGap
            });
            sub.setCoords();
        }
    };

    const resetImageScale = (img) => {
        if (!img) return;
        img.scaleX = 1;
        img.scaleY = 1;
        img.setCoords();
    };

    sysImgs.forEach(resetImageScale);

    if (layout === 'article_1') {
        let topAssets = [];
        if (graph) topAssets.push(graph);
        if (heatmap) topAssets.push(heatmap);
        if (meta) topAssets.push(meta);

        const topCount = topAssets.length;
        const molCount = mols.length;

        const topAvailableH = molCount > 0 ? usableH * 0.48 : usableH * 0.82;
        const bottomAvailableH = molCount > 0 ? usableH * 0.36 : 0;

        logExportDebug('article_1 sections', {
            topCount,
            molCount,
            topAvailableH,
            bottomAvailableH
        });

        let currentY = padding + headerHeight;

        if (topAssets.length > 0) {
            const totalNaturalW = topAssets.reduce((sum, a) => sum + a.width, 0);
            const maxNaturalH = Math.max(...topAssets.map(a => a.height));
            const totalGapW = gap * Math.max(0, topAssets.length - 1);

            const widthScale = (usableW - totalGapW) / totalNaturalW;
            const heightScale = (topAvailableH - subtitleReserve) / maxNaturalH;
            const topScale = Math.min(widthScale, heightScale, 2.0);

            logExportDebug('article_1 top scale', {
                totalNaturalW,
                maxNaturalH,
                widthScale,
                heightScale,
                topScale
            });

            let scaledWidths = topAssets.map(a => a.width * topScale);
            let totalScaledW = scaledWidths.reduce((a, b) => a + b, 0) + totalGapW;
            let startX = (canvasW - totalScaledW) / 2;

            let maxPlacedH = 0;

            topAssets.forEach((a, idx) => {
                a.scale(topScale);
                const w = a.getScaledWidth();
                const h = a.getScaledHeight();

                a.set({
                    left: startX + w / 2,
                    top: currentY + subtitleReserve + h / 2
                });
                a.setCoords();
                snapSubtitle(a);

                startX += w + gap;
                maxPlacedH = Math.max(maxPlacedH, h + subtitleReserve);

                logExportDebug('placed top asset', {
                    idx,
                    asset: getObjectDebugInfo(a)
                });
            });

            currentY += maxPlacedH + gap;
        }

        if (mols.length > 0) {
            const totalNaturalW = mols.reduce((sum, a) => sum + a.width, 0);
            const maxNaturalH = Math.max(...mols.map(a => a.height));
            const totalGapW = gap * Math.max(0, mols.length - 1);

            const widthScale = (usableW - totalGapW) / totalNaturalW;
            const heightScale = (bottomAvailableH - subtitleReserve) / maxNaturalH;
            const molScale = Math.min(widthScale, heightScale, 1.8);

            logExportDebug('article_1 mol scale', {
                totalNaturalW,
                maxNaturalH,
                widthScale,
                heightScale,
                molScale
            });

            let scaledWidths = mols.map(a => a.width * molScale);
            let totalScaledW = scaledWidths.reduce((a, b) => a + b, 0) + totalGapW;
            let startX = (canvasW - totalScaledW) / 2;

            mols.forEach((a, idx) => {
                a.scale(molScale);
                const w = a.getScaledWidth();
                const h = a.getScaledHeight();

                a.set({
                    left: startX + w / 2,
                    top: currentY + subtitleReserve + h / 2
                });
                a.setCoords();
                snapSubtitle(a);

                startX += w + gap;

                logExportDebug('placed mol asset', {
                    idx,
                    asset: getObjectDebugInfo(a)
                });
            });
        }
    }

    else if (layout === 'article_2') {
        const leftAssets = [];
        if (graph) leftAssets.push(graph);
        if (heatmap) leftAssets.push(heatmap);
        if (meta) leftAssets.push(meta);

        const leftW = usableW * 0.48;
        const rightW = usableW * 0.48;
        let leftY = padding + headerHeight;
        let rightY = padding + headerHeight;

        leftAssets.forEach((a, idx) => {
            resetImageScale(a);
            const scale = Math.min(
                leftW / a.width,
                (usableH / Math.max(1, leftAssets.length) - subtitleReserve) / a.height,
                1.8
            );
            a.scale(scale);
            a.set({
                left: padding + leftW / 2,
                top: leftY + subtitleReserve + a.getScaledHeight() / 2
            });
            a.setCoords();
            snapSubtitle(a);
            leftY += a.getScaledHeight() + subtitleReserve + gap;

            logExportDebug('article_2 left asset', { idx, asset: getObjectDebugInfo(a), scale });
        });

        mols.forEach((a, idx) => {
            resetImageScale(a);
            const scale = Math.min(
                rightW / a.width,
                (usableH / Math.max(1, mols.length) - subtitleReserve) / a.height,
                1.8
            );
            a.scale(scale);
            a.set({
                left: canvasW - padding - rightW / 2,
                top: rightY + subtitleReserve + a.getScaledHeight() / 2
            });
            a.setCoords();
            snapSubtitle(a);
            rightY += a.getScaledHeight() + subtitleReserve + gap;

            logExportDebug('article_2 right asset', { idx, asset: getObjectDebugInfo(a), scale });
        });
    }

    else if (layout === 'grid') {
        let cols = Math.ceil(Math.sqrt(sysImgs.length));
        let rows = Math.ceil(sysImgs.length / cols);
        const cellW = (usableW - gap * (cols - 1)) / cols;
        const cellH = (usableH - gap * (rows - 1)) / rows;

        sysImgs.forEach((a, i) => {
            resetImageScale(a);

            const scale = Math.min(
                cellW / a.width,
                (cellH - subtitleReserve) / a.height,
                1.8
            );

            a.scale(scale);

            let r = Math.floor(i / cols);
            let cIdx = i % cols;

            const x = padding + cIdx * (cellW + gap) + cellW / 2;
            const y = padding + headerHeight + r * (cellH + gap) + subtitleReserve + a.getScaledHeight() / 2;

            a.set({ left: x, top: y });
            a.setCoords();
            snapSubtitle(a);

            logExportDebug('grid asset', {
                idx: i,
                row: r,
                col: cIdx,
                cellW,
                cellH,
                scale,
                asset: getObjectDebugInfo(a)
            });
        });
    }

    fabCanvas.renderAll();
    saveHistory();

    logExportDebug('applyFabricSnapLayout end', {
        objects: sysImgs.map(getObjectDebugInfo)
    });
}

// --- 7. UPDATING GRAPH/DATA COLORS NATIVELY ---
async function applyDataNodesColor() {
    const newColor = document.getElementById('exp-data-color').value;

    if (typeof nodesDataset !== 'undefined' && nodesDataset) {
        let updates = [];
        nodesDataset.forEach(n => updates.push({ id: n.id, font: { color: newColor } }));
        nodesDataset.update(updates);
    }

    if (typeof gridNodesDatasets !== 'undefined' && gridNodesDatasets) {
        Object.entries(gridNodesDatasets).forEach(([key, ds]) => {
            let updates = [];
            ds.forEach(n => updates.push({ id: n.id, font: { color: newColor } }));
            ds.update(updates);
        });
    }

    if (typeof network !== 'undefined' && network) {
        network.setOptions({ nodes: { font: { color: newColor } } });
    }

    if (typeof gridNetworks !== 'undefined') {
        gridNetworks.forEach((n, idx) => {
            n.network.setOptions({ nodes: { font: { color: newColor } } });
        });
    }

    window.exportLabelColor = newColor;
    if (typeof update3DViewerOrGrid === 'function') {
        update3DViewerOrGrid();
    }

    await new Promise(r => setTimeout(r, 400));

    let graphObj = fabCanvas.getObjects().find(o => o.customId === 'graph');
    const graphCanvasUI = document.querySelector('#top-panel canvas');

    if (graphObj && graphCanvasUI) {
        let img = new Image();
        img.src = graphCanvasUI.toDataURL('image/png');
        await new Promise(r => img.onload = r);
        let croppedUrl = await cropTransparent(img);

        await new Promise(resolve => {
            const savedW = graphObj.getScaledWidth();
            const savedProps = {
                customId: graphObj.customId,
                customType: graphObj.customType,
                customLabel: graphObj.customLabel,
                isSystemAsset: graphObj.isSystemAsset,
                originX: graphObj.originX,
                originY: graphObj.originY,
                left: graphObj.left,
                top: graphObj.top,
                visible: graphObj.visible
            };
            const savedIdx = fabCanvas.getObjects().indexOf(graphObj);

            fabric.Image.fromURL(croppedUrl, function(fImg) {
                if (!fImg || !fImg.width) { resolve(); return; }
                const scale = savedW / fImg.width;
                fImg.set({ ...savedProps, scaleX: scale, scaleY: scale });
                
                isHistoryAction = true;
                fabCanvas.remove(graphObj);
                fabCanvas.insertAt(fImg, savedIdx);
                isHistoryAction = false;
                resolve();
            });
        });
    }

    if (typeof viewers !== 'undefined' && viewers.length > 0) {
        const isDark = document.body.getAttribute('data-theme') === 'dark';
        for (let i = 0; i < viewers.length; i++) {
            let molObj = fabCanvas.getObjects().find(o => o.customId === `mol_${i}`);
            if (!molObj) continue;

            let v = viewers[i];
            v.setBackgroundColor(0x000000, 0.0);
            v.render();
            let molImg = new Image();
            molImg.src = v.pngURI();
            await new Promise(r => molImg.onload = r);
            let croppedMolUrl = await cropTransparent(molImg);
            v.setBackgroundColor(isDark ? "#1f2937" : "#ffffff", 1.0);
            v.render();

            await new Promise(resolve => {
                const savedW = molObj.getScaledWidth();
                const savedH = molObj.getScaledHeight();
                const savedProps = {
                    customId: molObj.customId,
                    customType: molObj.customType,
                    customLabel: molObj.customLabel,
                    isSystemAsset: molObj.isSystemAsset,
                    originX: molObj.originX,
                    originY: molObj.originY,
                    left: molObj.left,
                    top: molObj.top,
                    visible: molObj.visible
                };
                const savedIdx = fabCanvas.getObjects().indexOf(molObj);

                fabric.Image.fromURL(croppedMolUrl, function(fImg) {
                    if (!fImg || !fImg.width) { resolve(); return; }
                    fImg.set({
                        ...savedProps,
                        scaleX: savedW / fImg.width,
                        scaleY: savedH / fImg.height,
                    });
                    isHistoryAction = true;
                    fabCanvas.remove(molObj);
                    fabCanvas.insertAt(fImg, savedIdx);
                    isHistoryAction = false;
                    resolve();
                });
            });
        }
    }

    isHistoryAction = false; 
    fabCanvas.renderAll();
    if (typeof updateHeatmapAsset === 'function') updateHeatmapAsset();
    saveHistory();
}

let paletteChangeLock = false;

async function changeExportPalette() {
    if (paletteChangeLock) return;
    paletteChangeLock = true;
    document.getElementById('exp-palette').disabled = true;

    try {
        let exportPalette = document.getElementById('exp-palette').value;
        if (!GRAPH_PALETTES[exportPalette]) return;

        activePaletteName = exportPalette;

        // Apply colors to datasets
        if (typeof assignAllPairGraphColors === 'function') assignAllPairGraphColors();
        if (typeof assignGraphColorsToData === 'function') assignGraphColorsToData();

        let isGlobal = (
            typeof masterData !== 'undefined' && masterData.mode === 'pair' &&
            document.getElementById('pair-selector') && document.getElementById('pair-selector').value === 'GRID'
        );

        if (isGlobal) {
            Object.keys(masterData.pairs).forEach(k => {
                if (typeof applyGraphFiltersGrid === 'function') applyGraphFiltersGrid(k);
            });
        } else {
            if (typeof currentGraphMode !== 'undefined' && currentGraphMode === 'associated') {
                if (typeof applyGraphFilters === 'function') applyGraphFilters(true);
            } else if (typeof handleFilteredChange === 'function') {
                handleFilteredChange();
            }
        }

        await new Promise(r => setTimeout(r, 400));

        let graphObj = fabCanvas.getObjects().find(o => o.customId === 'graph');
        const graphCanvasUI = document.querySelector('#top-panel canvas');
        
        if (graphObj && graphCanvasUI) {
            let img = new Image();
            img.src = graphCanvasUI.toDataURL('image/png');
            await new Promise(r => img.onload = r);
            let croppedUrl = await cropTransparent(img);

            await new Promise(resolve => {
                const savedW = graphObj.getScaledWidth();
                const savedProps = {
                    customId: graphObj.customId,
                    customType: graphObj.customType,
                    customLabel: graphObj.customLabel,
                    isSystemAsset: graphObj.isSystemAsset,
                    originX: graphObj.originX,
                    originY: graphObj.originY,
                    left: graphObj.left,
                    top: graphObj.top,
                    visible: graphObj.visible
                };
                const savedIdx = fabCanvas.getObjects().indexOf(graphObj);

                fabric.Image.fromURL(croppedUrl, function(fImg) {
                    if (!fImg || !fImg.width) { resolve(); return; }
                    const scale = savedW / fImg.width;
                    fImg.set({ ...savedProps, scaleX: scale, scaleY: scale });
                    
                    isHistoryAction = true;
                    fabCanvas.remove(graphObj);
                    fabCanvas.insertAt(fImg, savedIdx);
                    isHistoryAction = false;
                    resolve();
                });
            });
        }

        // Update the label text natively
        const graphLabel = fabCanvas.getObjects().find(o => o.customId === 'graph_label');
        if (graphLabel) {
            let pairSel = document.getElementById('pair-selector');
            let labelText = pairSel ? pairSel.value : 'Associated Graph';
            graphLabel.text = `Graph: ${labelText}`;
            graphLabel.setCoords();
        }

        fabCanvas.renderAll();
        saveHistory();

    } catch (e) {
        logError('changeExportPalette crashed', e);
    } finally {
        paletteChangeLock = false;
        document.getElementById('exp-palette').disabled = false;
    }
}

// Global UI Updates
function updateCanvasBackground() {
    if (!fabCanvas) return;

    const isTrans = document.getElementById('exp-transparent-bg').checked;
    fabCanvas.backgroundColor = isTrans ? null : document.getElementById('exp-bg-color').value;
    fabCanvas.renderAll();
    saveHistory();
}

function resizeCanvasWorkspace() {
    canvasW = parseInt(document.getElementById('exp-w').value) || 1920;
    canvasH = parseInt(document.getElementById('exp-h').value) || 1080;

    if (fabCanvas) {
        fabCanvas.setWidth(canvasW);
        fabCanvas.setHeight(canvasH);
        fabCanvas.renderAll();
        saveHistory();
    }

    updatePreviewZoom();
    requestAnimationFrame(() => {
        centerExportWorkspace();
    });
}

function downloadFabricExport() {
    if (!fabCanvas) return;

    fabCanvas.discardActiveObject();
    fabCanvas.renderAll();

    const dataURL = fabCanvas.toDataURL({
        format: 'png',
        quality: 1,
        multiplier: 1
    });

    const link = document.createElement('a');
    let runName = (typeof graphData !== 'undefined' && graphData) ? graphData.run_name : 'Export';
    link.download = `MHCXGraph_${runName}_Figure.png`;
    link.href = dataURL;
    link.click();
}

async function closeExportStudio() {
    document.getElementById('export-modal').style.display = 'none';

    if (window.originalPalette) {
        activePaletteName = window.originalPalette;
        if (typeof assignAllPairGraphColors === 'function') assignAllPairGraphColors();
        if (typeof assignGraphColorsToData === 'function') assignGraphColorsToData();

        let isGlobal = (
            typeof masterData !== 'undefined' &&
            masterData.mode === 'pair' &&
            document.getElementById('pair-selector') &&
            document.getElementById('pair-selector').value === 'GRID'
        );

        if (isGlobal) {
            Object.keys(masterData.pairs).forEach(k => {
                if (typeof applyGraphFiltersGrid === 'function') applyGraphFiltersGrid(k);
            });
        } else {
            if (typeof currentGraphMode !== 'undefined' && currentGraphMode === 'associated') {
                if (typeof applyGraphFilters === 'function') applyGraphFilters();
            } else if (typeof handleFilteredChange === 'function') {
                handleFilteredChange();
            }
        }
    }

    if (typeof updateThemeVars === 'function') updateThemeVars();

    let defaultThemeText = window.originalThemeText || '#f9fafb';

    if (typeof nodesDataset !== 'undefined' && nodesDataset) {
        let updates = [];
        nodesDataset.forEach(n => updates.push({ id: n.id, font: { color: defaultThemeText } }));
        nodesDataset.update(updates);
    }

    if (typeof gridNodesDatasets !== 'undefined' && gridNodesDatasets) {
        Object.values(gridNodesDatasets).forEach(ds => {
            let updates = [];
            ds.forEach(n => updates.push({ id: n.id, font: { color: defaultThemeText } }));
            ds.update(updates);
        });
    }
    await new Promise(r => setTimeout(r, 200));
    if (typeof network !== 'undefined' && network) {
        network.setOptions({ nodes: { font: { color: defaultThemeText } } });
        network.redraw();
    }

    if (typeof gridNetworks !== 'undefined') {
        gridNetworks.forEach(n => {
            n.network.setOptions({ nodes: { font: { color: defaultThemeText } } });
            n.network.redraw();
        });
    }

    window.exportLabelColor = null;
    if (typeof update3DViewerOrGrid === 'function') update3DViewerOrGrid();

    stopPreviewPan();

    if (fabCanvas) {
        fabCanvas.dispose();
        fabCanvas = null;
    }
    exportOriginalPalette = null;
}


async function updateHeatmapAsset() {
    if (!fabCanvas) return;
    
    let heatmapObj = fabCanvas.getObjects().find(o => o.customId === 'heatmap');
    if (!heatmapObj) return;

    // Reduz a opacidade para o usuário saber que está carregando
    heatmapObj.set({ opacity: 0.5 });
    fabCanvas.renderAll();

    let newUrl = await generatePlotlyHeatmap();
    if (!newUrl) {
        heatmapObj.set({ opacity: 1.0 });
        fabCanvas.renderAll();
        return;
    }

    let img = new Image();
    img.src = newUrl;
    await new Promise(r => img.onload = r);
    let croppedUrl = await cropTransparent(img);

    // Guarda as propriedades exatas para não perder a posição
    const savedProps = {
        customId: heatmapObj.customId,
        customType: heatmapObj.customType,
        customLabel: heatmapObj.customLabel,
        isSystemAsset: heatmapObj.isSystemAsset,
        originX: heatmapObj.originX,
        originY: heatmapObj.originY,
        left: heatmapObj.left,
        top: heatmapObj.top,
        scaleX: heatmapObj.scaleX,
        scaleY: heatmapObj.scaleY,
        visible: heatmapObj.visible,
        opacity: 1.0
    };
    
    const savedIdx = fabCanvas.getObjects().indexOf(heatmapObj);

    fabric.Image.fromURL(croppedUrl, function(fImg) {
        if (!fImg) return;
        fImg.set(savedProps);
        
        isHistoryAction = true;
        fabCanvas.remove(heatmapObj);
        fabCanvas.insertAt(fImg, savedIdx);
        isHistoryAction = false;
        
        fabCanvas.renderAll();
        saveHistory();
    });
}
