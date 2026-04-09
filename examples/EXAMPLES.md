# MHCXGraph — Examples

This directory contains a curated set of examples to validate and explore MHCXGraph's main functionalities. Each case is driven by a manifest file (JSON) that defines inputs, processing mode, selectors, and output path.

## Directory Structure

```
examples/
├── EXAMPLES.md
├── input
│   ├── pre-renumbered               # Original MHC-I structures before renumbering
│   └── renumbered                   # Renumbered structures ready for processing
├── manifests
│   ├── manifest-minimal.json
│   ├── manifest-multiple.json
│   ├── manifest-pairwise.json
│   └── manifest-screening.json
└── results                          # Output is written here (defined per manifest)
```

---

## Overview

The core workflow of MHCXGraph is:

```
run  →  (optional) heatmap
```

Two usage paths depending on whether you want to restrict the analysis to functionally relevant hotspot residues:

| Path | When to use |
|---|---|
| **A** — `run` → `heatmap` | General structural comparison; no residue filtering |
| **B** — `renumber` → `run` (with selectors) → `heatmap` | Analysis restricted to MHC hotspot residues (requires IMGT-standardized numbering) |

`renumber` and `heatmap` are accessory modules. `renumber` is only required when using the `MHC1` or `MHC2` selectors (see [Section 3](#3-using-hotspot-residue-selectors-conditional)).

---

## 1. Core processing (`run` module)

The `run` command is the main execution engine. It is driven by a manifest file that defines the inputs, processing mode, and structural selectors. For a detailed breakdown of all available configuration parameters, please refer to the [Manifest] section.

> [!WARNING]
> The input PDB files should not contain other macromolecules (such as TCRs or antibodies) interacting with the MHC in a way that affects the exposure of surface residues. The presence of these molecules can alter Relative Solvent Accessibility (RSA) calculations and structural connectivity, leading to inaccurate comparison results.

### 1a. Pairwise mode (`pairwise`)

Performs an all-vs-all comparison across the input dataset. Its primary goal is to quantify structural similarities between every possible pair of MHC structures.

```bash
MHCXGraph run manifests/manifest-pairwise.json
```

### 1b. Multiple mode (`multiple`)

Used for simultaneous comparison of a group of structures (e.g., identifying shared structural features within a subset in a single run).

```bash
MHCXGraph run manifests/manifest-multiple.json
```

### 1c. Screening mode (`screening`)

Optimized for high-speed one-to-many comparisons. It compares a specific reference structure against a directory of target PDB files to identify structural cross-reactive candidates.

```bash
MHCXGraph run manifests/manifest-screening.json
```

---

## 2. Visualization module (`heatmap`)

After a pairwise run, MHCXGraph can generate a similarity heatmap with hierarchical clustering. The similarity metric is a coverage-based index defined as the fraction of unique nodes in association graph frames relative to the total number of nodes in both input graphs.

```bash
MHCXGraph heatmap \
  -i results/pairwise/PAIRWISE/ \
  -o results/pairwise/ \
  -n similarity_heatmap.png
```

| Flag | Description |
|---|---|
| `-i` | Input directory (the `PAIRWISE/` folder from a pairwise run) |
| `-o` | Output directory for the heatmap image |
| `-n` | Output filename |

> **Note:** This module requires a completed pairwise run as input. Run `MHCXGraph run` with a `pairwise` manifest first.

### Output files

Besides the `.png` plot, the heatmap module generates several data files essential for downstream statistical analysis and identifying specific structural features:

1. `distance_matrix.csv`: Contains the numerical representation of dissimilarity between every pair of proteins.
2. `component_count_matrix.csv`: Provides data on the number of connected components shared between every pair.
3. `unique_nodes_graph_*.csv`: Detailed lists of specific residues involved in the shared nodes for every pair.

---

## 3. Using hotspot residue selectors (conditional)

By default, MHCXGraph processes the full structure. If you want to restrict the analysis to functionally relevant residues — such as those lining the peptide-binding groove or forming TCR contact points — you can use the `MHC1` or `MHC2` selectors in your manifest.

**When using these selectors, you must first run `renumber`.** The residue indices in `MHC1` and `MHC2` follow the IMGT domain numbering scheme; applying them to structures with non-standardized numbering will produce incorrect or incomplete results.

### Step 1 — Standardize numbering (`renumber`)

```bash
MHCXGraph renumber \
  -i examples/input/pre-renumbered/ \
  -o examples/input/renumbered/ \
  -c mhci \
  --suffix _renumbered
```

| Flag | Description |
|---|---|
| `-i` | Input directory with original PDB files |
| `-o` | Output directory for renumbered files |
| `-c` | Chain convention (`mhci` for MHC class I, `mhcii` for MHC class II) |
| `--suffix` | Suffix appended to each output filename |

### Step 2 — Configure selectors in the manifest

Two predefined hotspot sets are available:

**`MHC1`** — MHC class I residues involved in peptide presentation and TCR recognition. Covers positions along the α1 and α2 helices and the peptide-binding groove floor (chain A, IMGT numbering).

```json
"MHC1": {
  "chains": ["C"],
  "residues": {
    "A": [18, 19, 42, 43, 44, 54, 55, 56, 58, 59, 61, 62, 63, 64, 65, 66,
          68, 69, 70, 71, 72, 73, 75, 76, 79, 80, 83, 84, 89, 108, 109,
          142, 143, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155,
          156, 157, 158, 159, 161, 162, 163, 165, 166, 167, 169, 170, 171]
  }
}
```

**`MHC2`** — MHC class II residues at the peptide-binding interface, distributed across the α chain (chain A) and β chain (chain B).

```json
"MHC2": {
  "chains": ["C"],
  "residues": {
    "A": [37, 51, 52, 53, 55, 56, 58, 59, 60, 62, 63, 65, 66, 67, 69],
    "B": [56, 57, 59, 60, 61, 62, 63, 65, 66, 67, 68, 69, 70, 71, 72,
          73, 74, 77, 78, 81]
  }
}
```

> **Note:** Use `MHC1` for class I structures and `MHC2` for class II. Do not mix selectors across incompatible chain configurations.

### Step 3 — Run with the renumbered input

Point the manifest to the renumbered directory and proceed normally:

```bash
MHCXGraph run manifests/manifest-pairwise.json
```
