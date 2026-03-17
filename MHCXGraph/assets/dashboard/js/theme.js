function toggleTheme() {
    const body = document.body;
    const newIsDark = !(body.getAttribute('data-theme') === 'dark');
    body.setAttribute('data-theme', newIsDark ? 'dark' : 'light');
    updateThemeVars();
    
    if(network) network.setOptions({ nodes: { font: { color: themeText } } });
    gridNetworks.forEach(n => n.network.setOptions({ nodes: { font: { color: themeText } } }));
    
    if (masterData.mode === 'pairwise' && document.getElementById('pair-selector').value === 'GRID') {
        Object.keys(masterData.pairs).forEach(k => applyGraphFiltersGrid(k));
    } else {
        if (currentGraphMode === 'associated') applyGraphFilters();
        else if (currentGraphMode === 'filtered') handleFilteredChange();
    }

    const bgMol = newIsDark ? "#1f2937" : "#ffffff"; 
    viewers.forEach(v => { v.setBackgroundColor(bgMol); v.render(); });
    resetMolColors();
}

function updateThemeVars() {
    const isDark = document.body.getAttribute('data-theme') === 'dark';
    themeText = isDark ? '#f9fafb' : '#111827'; themeBorder = isDark ? '#f9fafb' : '#111827';
}
