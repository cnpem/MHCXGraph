.. _manifest_section:

Execution configuration (Manifest)
**********************************

**MHCXGraph** is a Python package designed primarily to be used as a command-line interface (CLI). This section describes the most important parameters, along with their use cases and recommended configurations.

All parameters are provided through a JSON configuration file, referred to as the `manifest`. The file name is user-defined. The manifest is organized into four top-level sections: `settings`, `classes`, `inputs`, and `selectors`.

Overview
========

The manifest is organized into four main sections, each responsible for a different stage of the execution pipeline:

- :ref:`settings_section`  
  Defines execution parameters such as graph construction, filtering thresholds, and output configuration.

- :ref:`inputs_section`  
  Specifies which structure files are loaded and which selectors are applied to them.

- :ref:`selectors_section`  
  Defines residue selection rules used to restrict the analysis to specific regions of the structures.

- :ref:`classes_section`  
  Defines grouping rules that map residues or continuous descriptors into shared categories during triad construction.

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

.. _settings_section:

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

  .. image:: ../_static/images/node_granularity_light.svg
     :align: center
     :class: only-light
     :width: 650px
     :alt: Node granularity illustration

  .. image:: ../_static/images/node_granularity_dark.svg
     :align: center
     :class: only-dark
     :width: 650px
     :alt: Node granularity illustration
     
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

.. mhcx-param:: filter_triads_by_chain
  :type: list | None
  :default: None

  Restricts triads to those that include at least one residue from the specified chains.

  This parameter is particularly useful when multiple pMHC structures share the same MHC but differ in their bound peptides. In such cases, it restricts the analysis to the peptide and the surrounding MHC residues that interact with it.

  For example, passing ``["C"]`` ensures that, for any triad :math:`(v_1, v_2, v_3)`, at least one residue belongs to chain ``C``. Triads that do not satisfy this condition are discarded.

.. mhcx-param:: max_gap_helix
  :type: int
  :default: 0

  Maximum allowed gap (in residue index) between consecutive helical segments for them to be treated as a single continuous helix.

.. mhcx-param:: distance_bin_width
  :type: float
  :default: 2.0

  Defines the width of distance bins used to discretize inter-residue distances in triads. 

.. mhcx-param:: local_distance_diff_threshold
  :type: float
  :default: 1.0

  Maximum allowed difference between corresponding distances (:math:`d_1`, :math:`d_2`, :math:`d_3`) when comparing two triads. This threshold is applied during triad association to enforce local geometric consistency.

.. mhcx-param:: global_distance_diff_threshold
  :type: float
  :default: 2.0

  Maximum allowed distance difference between non-adjacent nodes across structures during frame generation. This parameter enforces global geometric consistency.

.. mhcx-param:: close_tolerance
  :type: float
  :default: 0.1

  Tolerance used when assigning distances to discretization bins. Values within this tolerance of a bin center may be associated with that bin.

.. mhcx-param:: rsa_filter
  :type: float
  :default: 0.1

  Minimum relative solvent accessibility (RSA) required for canonical residues to be included as nodes in the graph.

.. mhcx-param:: asa_filter
  :type: float
  :default: 5

  Minimum absolute solvent accessibility (ASA) required for non-canonical residues and ligands to be included as nodes.

.. mhcx-param:: triad_rsa
  :type: bool
  :default: False

  If ``True``, RSA values are included as features in the triad token representation, allowing solvent exposure to influence`` matching.

  .. note::
    
    When using :mhcx-param:`node_granularity` = ``ca_only``, it is recommended to set ``triad_rsa = False``, since centroids computed from C\ :math:`\alpha` atoms are not well correlated with RSA values.

.. mhcx-param:: rsa_bin_width
  :type: float
  :default: 0.3

  Defines the width of RSA bins used to discretize solvent accessibility values in triads.

.. mhcx-param:: rsa_diff_threshold
  :type: float
  :default: 0.3

  Maximum allowed difference in RSA values between corresponding nodes of two triads during association.

.. mhcx-param:: close_tolerance_rsa
  :type: float
  :default: 0.01

  Tolerance used when assigning RSA values to discretization bins, analogous to :mhcx-param:`close_tolerance` for distances.

.. _inputs_section:

Inputs
======

The ``inputs`` section defines which structure files are loaded for the analysis and which selectors are applied to them. Each entry in ``inputs`` is an input rule describing where the files are located, how they should be collected, and which selector definitions should be used.

This section allows the user to process either individual files or entire directories. When a directory is provided, all files matching the allowed extensions are considered. Selectors referenced in an input rule are resolved for each matching file before graph construction.

Each input rule may contain the following fields:

.. list-table:: Input fields
   :header-rows: 1
   :widths: 20 20 60

   * - Field
     - Type
     - Description

   * - ``path``
     - str or list[str]
     - Path to a structure file or directory. Multiple paths may also be provided.

   * - ``extensions``
     - list[str]
     - File extensions allowed when scanning directories, such as ``[".pdb", ".cif"]``.

   * - ``enable_tui``
     - bool
     - If ``True``, an interactive terminal selection is shown when the input path is a directory.

   * - ``selectors``
     - list[dict]
     - List of :ref:`selectors <selectors_section>` references applied to the files matched by this input rule.

An example is shown below:

.. code-block:: json

   "inputs": [
     {
       "path": "examples/input/renumbered",
       "enable_tui": false,
       "extensions": [".pdb", ".cif"],
       "selectors": [
         { "name": "MHC1" }
       ]
     }
   ]

In this example, all ``.pdb`` and ``.cif`` files inside ``examples/input/renumbered`` are considered as inputs, and the selector ``MHC1`` is applied to each file.

Multiple input rules may be defined in the same manifest. This makes it possible to process structures from different folders, apply different selector sets, or combine file-specific and directory-based rules in a single execution.

.. note::

   If ``enable_tui`` is set to ``False``, all matching files are collected automatically.

.. _selectors_section:

Selectors
=========

The ``selectors`` section defines how residues are chosen from each input structure after graph construction. A selector can restrict the analysis by chain, residue number, secondary structure, or a logical combination of these criteria.

Selectors are not applied globally by default. Instead, they are referenced from the ``inputs`` section, where one or more selector names can be attached to a given input rule. During execution, the selected constraints are resolved for each file and then combined with the solvent-exposure filtering step.

In practice, selectors allow the user to focus the analysis on structurally or biologically relevant regions, such as residues from a given chain, residues belonging to helices, or specific residue positions known to participate in recognition.

The ``selectors`` section is a dictionary in which each key is a selector name and each value defines a set of constraints. A selector may contain the following fields:

.. list-table:: Selector fields
   :header-rows: 1
   :widths: 20 20 60

   * - Field
     - Type
     - Description

   * - ``chains``
     - list[str]
     - Selects all nodes belonging to the specified chains.

   * - ``residues``
     - dict[str, list[int]]
     - Selects specific residue numbers for each chain.

   * - ``structures``
     - list[str] or dict[str, list[str]]
     - Selects nodes by secondary-structure annotation. It may be defined globally as a list, or per chain as a dictionary.

   * - ``logic``
     - str
     - Boolean expression combining named sets such as ``chains``, ``residues``, ``structures``, and ``exposed``.

A simple example is shown below:

.. code-block:: json

   "selectors": {
     "MHC1": {
       "chains": ["C"],
       "residues": {
         "A": [18, 19, 42, 43]
       }
     },
     "general": {
       "chains": ["C"],
       "structures": ["helix"],
       "residues": {}
     }
   }

In this example, ``MHC1`` restricts the selection to chain ``C`` together with a predefined set of residues from chain ``A``. The selector ``general`` restricts the graph to residues in chain ``C`` that are also annotated as helices.

By default, when no ``logic`` expression is provided, the sets defined by ``chains``, ``residues``, and ``structures`` are first combined by union, and the result is then intersected with the set of exposed residues. If none of these fields is provided, the selector reduces to the exposed-residue filter alone.

When a ``logic`` expression is provided, it is evaluated explicitly. The supported operators are ``&`` for intersection, ``|`` for union, and ``!`` for negation. Parentheses may also be used. For example:

.. code-block:: json

   "selectors": {
     "example_selector": {
       "chains": ["A", "C"],
       "structures": ["helix"],
       "logic": "exposed & chains:C & structures"
     }
   }

This expression keeps only exposed residues that belong to chain ``C`` and are part of a helix.

Selectors may also define chain-specific structural filters using a dictionary:

.. code-block:: json

   "selectors": {
     "mhc2_like": {
       "structures": {
         "A": ["helix"],
         "B": ["sheet"]
       }
     }
   }

In this case, residues are selected according to the allowed secondary-structure types for each chain.

.. note::

   After selector evaluation, isolated ligand and water nodes are removed automatically if they have no edges within the selected subgraph.

.. note::

   If a logical expression is provided and does not explicitly include ``exposed``, the exposure filter is still applied afterward.

.. _classes_section:

Classes
=======

The ``classes`` section defines optional grouping rules used during triad generation. These classes allow different residues or continuous values to be mapped into the same category before token construction.

This mechanism is particularly useful when the analysis should emphasize broader physicochemical similarity rather than exact identity. For example, residues with similar chemical properties can be assigned to the same class, allowing them to contribute to the same triad token.

The ``classes`` section may define classes for residues, distances, and RSA values.

.. list-table:: Class categories
   :header-rows: 1
   :widths: 20 20 60

   * - Field
     - Type
     - Description

   * - ``residues``
     - dict[str, list[str]]
     - Groups amino acids into user-defined residue classes.

   * - ``distance``
     - dict[str, list[float]] or similar
     - Defines custom distance classes used instead of automatic distance discretization.

   * - ``rsa``
     - dict[str, list[float]] or similar
     - Defines custom RSA classes used instead of automatic RSA discretization.

When residue classes are provided, each residue is mapped to its corresponding class name before the triad token is created. If no residue classes are defined, the residue identity itself is used.

For example, the following configuration groups amino acids by general physicochemical type:

.. code-block:: json

   "classes": {
     "residues": {
       "hydrophobic": ["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"],
       "polar": ["SER", "THR", "ASN", "GLN", "TYR", "CYS", "GLY"],
       "positive": ["LYS", "ARG", "HIS"],
       "negative": ["ASP", "GLU"]
     }
   }

With this definition, a triad containing ``LEU``, ``ILE``, and ``VAL`` will be represented using the shared class ``hydrophobic`` rather than the individual residue names.

This can increase the tolerance of the comparison by allowing chemically similar residues to match even when their exact identities differ.

Custom distance and RSA classes may also be provided. When these are defined, they replace the default discretization based on parameters such as :mhcx-param:`distance_bin_width`, :mhcx-param:`rsa_bin_width`, and their associated tolerances.

.. note::

   The ``classes`` section is optional. If it is omitted, residues, distances, and RSA values are handled using their default representations.

.. note::

   Residue classes affect triad tokens directly. As a result, changing these definitions may significantly alter the number and type of associations detected.
