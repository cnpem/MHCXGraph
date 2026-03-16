function logDebug(msg, data) {
    let out = `[MHCX-DEBUG] ${msg}`;
    if (data !== undefined) {
        try { out += " | " + JSON.stringify(data); } catch(e) { out += " | [Complex Object]"; }
    }
    console.log(out);
}
function logError(msg, err) {
    console.error(`[MHCX-ERROR] ${msg}`);
    console.error(err);
}

function toggleSidebar() { document.getElementById('sidebar').classList.toggle('collapsed'); setTimeout(() => { if (network) network.fit(); }, 350); }

function toggleComp(targetId, el) { const t = document.getElementById(targetId); if (t.style.display === 'none') { t.style.display = 'block'; el.innerText = '▼'; } else { t.style.display = 'none'; el.innerText = '▶'; } }

__DATA_JS_INJECTION__

__INIT_JS_INJECTION__

__VIEWER_JS_INJECTION__

__GRID_JS_INJECTION__

__THEME_JS_INJECTION__

__ANALYSIS_JS_INJECTION__

__GRAPH_JS_INJECTION__

__STRUCTURES_JS_INJECTION__

__MODAL_JS_INJECTION__

window.onload = init;
setInterval(() => {
    const syncCheckbox = document.getElementById('syncViews');
    if (syncCheckbox && syncCheckbox.checked && window.activeMolViewer && viewers.length > 1) {
        const currentView = window.activeMolViewer.getView();
        viewers.forEach(v => { if (v !== window.activeMolViewer) v.setView(currentView); });
    }
}, 1000 / 30);
