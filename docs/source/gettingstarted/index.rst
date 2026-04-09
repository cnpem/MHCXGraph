Getting Started
***************

This example provides the fastest way to run **MHCXGraph** with a minimal manifest configuration. It uses the ``multiple`` mode, reads structures from ``examples/input/renumbered``, applies the predefined ``MHC1`` selector, and writes the results to ``examples/results/minimal``.

Before running this example, make sure the input structures have already been renumbered if you intend to use the ``MHC1`` selector. The predefined ``MHC1`` and ``MHC2`` selectors follow IMGT numbering, as described in :ref:`hotspot_selectors_section`.

Run the example with:

.. code-block:: bash

   MHCXGraph run examples/manifests/manifest-minimal.json

The corresponding manifest is shown below:

.. code-block:: json

   {
     "settings": {
       "run_mode": "multiple",
       "run_name": "minimal",
       "output_path": "examples/results/minimal",
       "edge_threshold": 8.5,
       "node_granularity": "all_atoms",
       "triad_rsa": false,
       "rsa_filter": 0.1,
       "global_distance_diff_threshold": 2.0,
       "local_distance_diff_threshold": 1.0,
       "distance_bin_width": 2
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
         "residues": {
           "A": [18, 19, 42, 43, 44, 54, 55, 56, 58, 59, 61, 62, 63, 64, 65, 66,
                 68, 69, 70, 71, 72, 73, 75, 76, 79, 80, 83, 84, 89, 108, 109,
                 142, 143, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155,
                 156, 157, 158, 159, 161, 162, 163, 165, 166, 167, 169, 170, 171]
         }
       },

       "MHC2": {
         "chains": ["C"],
         "residues": {
           "A": [37, 51, 52, 53, 55, 56, 58, 59, 60, 62, 63, 65, 66, 67, 69],
           "B": [56, 57, 59, 60, 61, 62, 63, 65, 66, 67, 68, 69, 70, 71, 72, 73,
                 74, 77, 78, 81]
         }
       },

       "general": {
         "chains": ["C"],
         "structures": ["helix"],
         "residues": {}
       }
     }
   }

This manifest is intentionally minimal. It includes only the parameters required for a basic execution, while leaving the remaining options at their default values.

The main settings used in this example are:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Field
     - Description

   * - :mhcx-param:`run_mode`
     - Uses the ``multiple`` execution mode.

   * - :mhcx-param:`run_name`
     - Sets the name of the run to ``minimal``.

   * - :mhcx-param:`output_path`
     - Writes the output to ``examples/results/minimal``.

   * - :mhcx-param:`edge_threshold`
     - Sets the distance cutoff used to define graph edges.

   * - :mhcx-param:`node_granularity`
     - Uses ``all_atoms`` to compute residue centroids.

   * - :mhcx-param:`triad_rsa`
     - Disables RSA-based features in the triad token.

   * - :mhcx-param:`rsa_filter`
     - Restricts the analysis to solvent-exposed canonical residues.

   * - :mhcx-param:`global_distance_diff_threshold`
     - Sets the global distance consistency threshold used during frame generation.

   * - :mhcx-param:`local_distance_diff_threshold`
     - Sets the local distance consistency threshold used during triad comparison.

   * - :mhcx-param:`distance_bin_width`
     - Defines the discretization width for distance values.

For a complete description of these fields, see :ref:`manifest_section`.
