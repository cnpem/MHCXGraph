# MHCXGraph — Examples

This directory contains a curated set of examples to validate and explore MHCXGraph's main functionalities. Each case is driven by a manifest file (JSON) that defines inputs, processing mode, selectors, and output path.

## Directory Structure

```
examples/
├── EXAMPLES.md
├── input
│   ├── pre-renumbered               # Original MHC-I structures before renumbering
│   └── renumbered                   # Renumbered structures ready for processing
├── manifests
│   ├── manifest-minimal.json
│   ├── manifest-multiple.json
│   ├── manifest-pairwise.json
│   └── manifest-screening.json
└── results                          # Output is written here (defined per manifest)

```

## 1. Preprocessing module (`renumber`)

Before starting any analysis that uses residue selectors (like `MHC1` or `MHC2`), you must standardize the structure numbering to the IMGT domain numbering scheme. This ensures full compatibility between hotspot residue lists and the input structures.

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
| `-c` | Chain convention (`mhci`for MHC class I and `mhcii` for MHC class II) |
| `--suffix` | Suffix appended to each output filename |


## 2. Core processing (`run` module)
The `run` command is the main execution engine. It is driven by a manifest file that defines the inputs, processing mode, and structural selectors. For a detailed breakdown of all available configuration parameters, please refer to the [Manifest] section.

>[!WARNING]
>The input PDB files should not contain other macromolecules (such as TCRs or antibodies) interacting with the MHC in a way that affects the exposure of surface residues. The presence of these molecules can alter Relative Solvent Accessibility (RSA) calculations and structural connectivity, leading to inaccurate comparison results.


### 2a. Pairwise mode (`pairwise`)
The `pairwise` mode performs an all-vs-all comparison across the input dataset. Its primary goal is to quantify structural similarities between every possible pair of MHC structures.

```bash
MHCXGraph run manifests/manifest-pairwise.json
```

### 2b. Multiple mode (`multiple`)
Used for simultaneous comparison of a group of structures (e.g., identifying shared structural features within a subset in a single run).

```bash
MHCXGraph run manifests/manifest-multiple.json
```
### 2c. Screening mode (`screening`)
Optimized for high-speed one-to-many comparisons. It compares a specific reference structure against a  directory of target PDB files to identify structural cross-reactive candidates.

```bash
MHCXGraph run manifests/manifest-screening.json
```

## 3. Visualization module (`heatmap`)

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
