const masterData = __GRAPH_DATA_JS_INJECTION__;

let graphData = masterData.mode === 'pairwise' ? null : masterData; 

const RESIDUE_CLASSES = {
    'ALA':'Hydrophobic', 'VAL':'Hydrophobic', 'ILE':'Hydrophobic', 'LEU':'Hydrophobic', 'MET':'Hydrophobic', 'PHE':'Hydrophobic', 'TYR':'Hydrophobic', 'TRP':'Hydrophobic',
    'SER':'Polar', 'THR':'Polar', 'ASN':'Polar', 'GLN':'Polar', 'CYS':'Polar',
    'ARG':'Charged (+)', 'HIS':'Charged (+)', 'LYS':'Charged (+)',
    'ASP':'Charged (-)', 'GLU':'Charged (-)',
    'GLY':'Special', 'PRO':'Special'
};

const CLASS_COLORS = {
    'Hydrophobic': '#4CAF50', 'Polar': '#2196F3', 'Charged (+)': '#3F51B5', 'Charged (-)': '#F44336', 'Special': '#FF9800', 'Unknown': '#9E9E9E'
};

const GRAPH_PALETTES = {
    standard: ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'],
    colorblind: ['#E69F00', '#56B4E9', '#009E73', '#0072B2', '#CC79A7', '#999999', '#000000'], 
    high_contrast: ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe']
};

const MOL_DEFAULTS = {
    dark: ['#FFC107', '#4FC3F7', '#69F0AE', '#FF8A65', '#E040FB', '#CDDC39', '#448AFF', '#18FFFF'],
    light: ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#17becf']
};

let activePaletteName = 'colorblind';
let graphNodeColors = {}; 
let molChainColors = {};  

let currentGraphMode = 'associated';
let activeFilteredProtIdx = 0;

let optNodeSize = 18, optLabelSize = 22, optEdgeWidth = 1.0, optSpringLen = 150;
let themeText = getCSSVar('--theme-text');
let themeBorder = getCSSVar('--theme-border');

let network = null;
let gridNetworks = [];
let gridNodesDatasets = {};
let gridEdgesDatasets = {};

let nodesDataset = new vis.DataSet([]);
let edgesDataset = new vis.DataSet([]);
let originalNodeColors = null;
let currentActiveNodes = new Set();

let loadedModels = []; 
let viewers = []; 
let tempEdgeIds = [];
let pinnedNode = null;

