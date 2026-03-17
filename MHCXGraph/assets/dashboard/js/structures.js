function updateMolChainColor(pIdx, chain, color) {
    molChainColors[`${pIdx}_${chain}`] = color;
    update3DViewerOrGrid(); initUploads(); 
}

function resetMolColors() {
    initMolColors(true); update3DViewerOrGrid(); initUploads(); 
}

function get3DColor(chain, globalIdx) {
    return molChainColors[`${globalIdx}_${chain}`] || '#cccccc';
}

async function autoLoadStructures() {
    let loadedAny = false;
    loadedModels = new Array(masterData.proteins.length).fill(null);
    for (let i = 0; i < masterData.proteins.length; i++) {
        const path = masterData.protein_paths ? masterData.protein_paths[i] : null;
        if (!path) continue;
        try {
            const filename = path.split(/[/\\]/).pop();
            let response = await fetch(path).catch(() => null);
            if (!response || !response.ok) response = await fetch(filename).catch(() => null);
            if (response && response.ok) {
                const text = await response.text();
                const format = path.toLowerCase().endsWith('.cif') ? 'cif' : 'pdb';
                loadedModels[i] = { text, format, name: masterData.proteins[i], protIdx: i };
                loadedAny = true;
            }
        } catch (e) {}
    }

    if (loadedAny) {
        const select = document.getElementById('view-selector');
        if (select && graphData) {
            select.innerHTML = '<option value="aligned">Aligned (Single View)</option><option value="separate">Separate (Grid View)</option>';
            loadedModels.forEach((m, idx) => {
                if (m && graphData.proteins.includes(m.name)) select.innerHTML += `<option value="prot_${idx}">Prot ${graphData.proteins.indexOf(m.name)}: ${m.name}</option>`;
            });
        }
        if (masterData.mode !== 'pairwise' || (document.getElementById('pair-selector').value !== 'GRID' && document.getElementById('pair-selector').value !== 'ANALYSIS')) {
            rebuildViewer();
        } else if (masterData.mode === 'pairwise' && document.getElementById('pair-selector').value === 'GRID') {
            document.getElementById('viewer-toolbar').style.display = 'flex';
            buildGridBottom(); 
        }
    }
}

async function loadStructures() {
    if (!loadedModels || loadedModels.length === 0) loadedModels = new Array(masterData.proteins.length).fill(null);
    for (let i = 0; i < masterData.proteins.length; i++) {
        const fileInput = document.getElementById(`file-prot-${i}`);
        if (fileInput && fileInput.files.length > 0) {
            const text = await fileInput.files[0].text();
            loadedModels[i] = { text, format: fileInput.files[0].name.endsWith('.cif') ? 'cif' : 'pdb', name: masterData.proteins[i], protIdx: i };
        }
    }
    
    if (masterData.mode === 'pairwise' && document.getElementById('pair-selector').value === 'GRID') {
        document.getElementById('viewer-toolbar').style.display = 'flex';
        buildGridBottom(); 
    } else {
        rebuildViewer();
        handleTreeChange();
    }
}
