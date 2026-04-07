# Structural Similarity Detection in pMHC Complexes (MHCXGraph)

A python package for graph-based detection of structurally similar surface regions in protein complexes, with a focus on peptide–MHC (pMHC) systems involved in T cell receptor (TCR) recognition.

The methodology supports the investigation of TCR cross-reactivity by identifying conserved surface patterns that may be recognized by the same TCR across different pMHC structures.

See also:

- [Documentation](https://cnpem.github.io/MHCXGraph/)
- [GitHub repository](https://github.com/cnpem/MHCXGraph)


## Installation

**From PyPI (stable):**
To install the latest release on
[PyPI](https://pypi.org/project/MHCXGraph), run:

```bash
pip install MHCXGraph
```

**From source (development):**
```bash
git clone https://github.com/cnpem/MHCXGraph.git
pip install -e MHCXGraph
```

## Quick Start

```bash
MHCXGraph run manifest.json
```

MHCXGraph is configured via a JSON manifest file. An example is in `examples/minimal/manifest.json`. A minimal manifest structure:

```json
{
  "settings": {
    "run_name": "my-run",
    "run_mode": "pairwise",
    "output_path": "results/pairwise",
    "edge_threshold": 10,
    "node_granularity": "ca_only",
    "triad_rsa": false,
    "rsa_filter": 0.1,
    "asa_filter": 5,
    "local_distance_diff_threshold":1.0,
    "global_distance_diff_threshold":2.0,
    "distance_bin_width": 2,
    "close_tolerance": 0.1
  },
  "inputs": [
    {
      "path": "data/structures",
      "extensions": [".pdb", ".cif"],
      "selectors": [{ "name": "MHC1" }]
    }
  ],
  "selectors": {
    "MHC1": {
      "chains": ["C"],
      "structures": {},
      "residues": {
        "A": [18,19,42,43,44,54,55,56,58,59,61,62,63,64,65,66,68,69,70,71,72,73,75,76,79,80,83,84,89,108,109,142,143,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,161,162,163,165,166,167,169,170,171]
      }
    }
  }
}
```
> [!NOTE]
> Input structures in this example have been pre-processed with **`MHCXGraph renumber`**.

## Basic Configuration Parameters¹

| Category | Parameter | Description | Default |
|----------|-----------|-------------|---------|
| **Run** | `run_name` | Name of the run | `test` |
|         | `run_mode` | Execution mode: `pairwise`, `multiple`, or `screening` | `multiple`|
|         | `output_path` | Path for results output | `./outputs` |
|         | `reference_structure` | Path to reference structure (required for `screening` mode) | `None`|
| **Graph** | `node_granularity` | Atomic representation for residue nodes: `ca_only`, `all_atoms`, `backbone`, or `sidechain` | `all_atoms` |
|           | `edge_threshold` | Distance cutoff (Å) for defining edges between nodes | `8.5` |
|           | `include_ligands` | Include ligands as graph nodes | `true` |
|           | `include_noncanonical_residues` | Include modified amino acids as nodes | `true` |
|           | `include_waters` | Include water molecules as nodes | `true` |
|           | `filter_triads_by_chain` | Restrict triads to those with at least one node from the specified chain | `None` |
|           | `max_gap_helix` | Max residue gap between helices to treat them as continuous | `0` |
| **Triad comparison** | `distance_bin_width` | Width of distance bins for discretization | `2.0` |
|                      | `local_distance_diff_threshold` | Max distance difference (d1/d2/d3) between triads for association | `1.0` |
|                      | `global_distance_diff_threshold` | Max distance difference between non-adjacent nodes across structures (frame generation step) | `2.0` |
|                      | `close_tolerance` | Tolerance for placing distances at bin center | `0.1` |
| **Surface representation** | `rsa_filter` | Min RSA for canonical residues to be included as nodes | `0.1` |
|                  | `asa_filter` | Min ASA for non-canonical residues, waters, and ligands | `5` |
| **RSA tokens** ² | `triad_rsa` | Use RSA values in triad token representation | `false` |
|                  | `rsa_bin_width` | Width of RSA bins for discretization | `0.3` |
|                  | `rsa_diff_threshold` | Max RSA difference between triad nodes for association | `0.3` |
|                  | `close_tolerance_rsa` | Tolerance for placing RSA values at bin center | `0.01` |

> ¹ A complete description of all parameters can be found in the [documentation](https://cnpem.github.io/MHCXGraph/)

> ² RSA token parameters are only active when `triad_rsa: true`.

### Selectors

Selectors define which residues and chains are included as graph nodes. They can be specified by chain ID, residue index, or secondary structure type (`helix`, etc.), and referenced by name in the `inputs` block.

## Citation

- (manuscript in preparation)

## License

The software is licensed under the terms of the GNU Affero General Public License 3 (AGPL3) and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.
