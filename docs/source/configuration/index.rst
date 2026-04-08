Execution configuration (Manifest)
**********************************

**MHCXGraph** is a Python packaged designed primarily to be used as a command-line interface (CLI). This section describes the most important parameters, along with their use cases and recommended configurations.

All parameters are provided through a JSON configuration file, referred to as the `manifest`. The file name is user-defined. The manifest is organized into four top-level sections: `settings`, `classes`, `inputs`, and `selectors`.

An example of a manifest is shown below:

.. code-block:: json

   {
     "settings": {
       "run_name": "multiple-run",
       "run_mode": "multiple",
       "max_chunks": 5,
       "output_path": "examples/results/multiple",
       "debug_logs": false,
       "debug_tracking": false,
       "track_steps": false,
       "edge_threshold": 10,
       "node_granularity": "ca_only",
       "include_ligands": true,
       "include_noncanonical_residues": true,
       "include_waters": false,
       "triad_rsa": false,
       "rsa_filter": 0.1,
       "local_distance_diff_threshold": 1.0,
       "global_distance_diff_threshold": 2.0,
       "distance_bin_width": 2,
       "close_tolerance": 0.1
     },

     "inputs": [
       {
         "path": "examples/input/renumbered",
         "enable_tui": false,
         "extensions": [".pdb", ".cif"],
         "selectors": [
           { "name": "MHC1" }
         ]
       }
     ],

     "selectors": {
       "MHC1": {
         "chains": ["C"],
         "structures": {},
         "residues": {
           "A": [18,19,42,43,44,54,55,56,58,59,61,62,63,64,65,66,68,69,70,71,72,73,75,76,79,80,83,84,89,108,109,142,143,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,161,162,163,165,166,167,169,170,171]
         }
       },
       "MHC2": {
         "chains": ["C"],
         "residues": {
           "A": [37,51,52,53,55,56,58,59,60,62,63,65,66,67,69],
           "B": [56,57,59,60,61,62,63,65,66,67,68,69,70,71,72,73,74,77,78,81]
         }
       },
       "general": {
         "chains": ["C"],
         "structures": ["helix"],
         "residues": {}
       }
     }
   }

Settings
========

This section defines the execution parameters that control graph construction, residue filtering, and output generation. The table below provides an overview, followed by a detailed description of each parameter.

.. list-table:: Key Configuration Parameters¹
   :header-rows: 1
   :widths: 15 20 45 20

   * - Category
     - Parameter
     - Description
     - Default

   * - **Execution**
     - :mhcx-param:`run_name`
     - Name of the run
     - ``test``
   * -
     - :mhcx-param:`run_mode`
     - Execution mode: ``pairwise``, ``multiple``, or ``screening``
     - ``multiple``
   * -
     - :mhcx-param:`output_path`
     - Path for results output
     - ``./outputs``
   * -
     - :mhcx-param:`reference_structure`
     - Path to reference structure (required for ``screening`` mode)
     - ``None``

   * - **Graph**
     - :mhcx-param:`node_granularity`
     - Atomic representation for residue nodes: ``ca_only``, ``all_atoms``, ``backbone``, or ``sidechain``
     - ``all_atoms``
   * -
     - :mhcx-param:`edge_threshold`
     - Distance cutoff (Å) for defining edges between nodes
     - ``8.5``
   * -
     - :mhcx-param:`include_ligands`
     - Include ligands as graph nodes
     - ``true``
   * -
     - :mhcx-param:`include_noncanonical_residues`
     - Include modified amino acids as nodes
     - ``true``
   * -
     - :mhcx-param:`include_waters`
     - Include water molecules as nodes
     - ``true``
   * -
     - :mhcx-param:`filter_triads_by_chain`
     - Restrict triads to those with at least one node from the specified chain
     - ``None``
   * -
     - :mhcx-param:`max_gap_helix`
     - Max residue gap between helices to treat them as continuous
     - ``0``

   * - **Triad comparison**
     - :mhcx-param:`distance_bin_width`
     - Width of distance bins for discretization
     - ``2.0``
   * -
     - :mhcx-param:`local_distance_diff_threshold`
     - Max distance difference (d1/d2/d3) between triads for association
     - ``1.0``
   * -
     - :mhcx-param:`global_distance_diff_threshold`
     - Max distance difference between non-adjacent nodes across structures (frame generation step)
     - ``2.0``
   * -
     - :mhcx-param:`close_tolerance`
     - Tolerance for placing distances at bin center
     - ``0.1``

   * - **Surface representation**
     - :mhcx-param:`rsa_filter`
     - Min RSA for canonical residues to be included as nodes
     - ``0.1``
   * -
     - :mhcx-param:`asa_filter`
     - Min ASA for non-canonical residues, waters, and ligands
     - ``5``

   * - **RSA tokens²**
     - :mhcx-param:`triad_rsa`
     - Use RSA values in triad token representation
     - ``false``
   * -
     - :mhcx-param:`rsa_bin_width`
     - Width of RSA bins for discretization
     - ``0.3``
   * -
     - :mhcx-param:`rsa_diff_threshold`
     - Max RSA difference between triad nodes for association
     - ``0.3``
   * -
     - :mhcx-param:`close_tolerance_rsa`
     - Tolerance for placing RSA values at bin center
     - ``0.01``

.. mhcx-param:: run_name
  :type: string
  :default: test

  The name assigned to the current execution. The raw JSON file containing the graph nodes is saved as ``graph_{run_name}.json``.

.. mhcx-param:: run_mode
  :type: string
  :default: multiple

  There are three possible `run_mode` options, which determine how the execution is performed: `multiple`, `pairwise`, and `screening`.

  - ``multiple``: compares all input pMHC structures simultaneously, identifying regions that are common across all of them.

  - ``pairwise``: performs pairwise comparisons between all input pMHC structures.

  - ``screening``: uses the :mhcx-param:`reference_structure` to compare a reference pMHC against all other input structures. This mode is particularly useful when you have a target pMHC and want to evaluate it against a larger set of structures.

.. mhcx-param:: output_path
  :type: string
  :default: ./outputs

  The folder destination of all results generated. All the parents folder will be created even though they doesn't exist.

.. mhcx-param:: reference_structure
  :type: string | None
  :default: None

  The target structure that will be compared against all the other structures passed in the input parameter when the program is executed in `screening` mode.

.. mhcx-param:: node_granularity
  :type: string
  :default: all_atoms

  The ``node_granularity`` parameter determines which atoms of each residue are used to compute centroid positions. This choice strongly affects graph connectivity: in combination with :mhcx-param:`edge_threshold`, it can determine whether two residues are considered adjacent (i.e., whether an edge is formed between their corresponding nodes).

  Four granularities are available: ``ca_only``, ``all_atoms``, ``backbone``, and ``sidechain``.

  - ``ca_only``: uses only the C\ :math:`\alpha` atom of each residue.

  - ``all_atoms``: uses all heavy atoms in the residue.

  - ``backbone``: uses only heavy atoms belonging to the backbone.

  - ``sidechain``: uses only heavy atoms belonging to the side chain.

  .. note::

     Hydrogen atoms are not considered in any granularity. All centroid computations are performed using heavy atoms only.

.. mhcx-param:: edge_threshold
  :type: float
  :default: 8.5  

  The ``edge_threshold`` parameter defines the maximum distance (in Å) between the centroids of two residues for an edge to be created between their corresponding nodes.

.. mhcx-param:: include_ligands
  :type: bool
  :default: True

  If ``True``, ligand molecules bound to the structure are included as nodes and can participate in edge formation.

  .. note::

     This option does not include water molecules. Use :mhcx-param:`include_waters` to control their inclusion.

.. mhcx-param:: include_noncanonical_residues
  :type: bool
  :default: True

  If ``True``, non-canonical residues (e.g., modified amino acids) are included as nodes and can participate in edge formation.

.. mhcx-param:: include_waters
  :type: bool
  :default: True

  If ``True``, water molecules are included as nodes and can participate in edge formation.
