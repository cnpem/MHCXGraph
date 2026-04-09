######################################
Welcome to MHCXGraph's documentation!
######################################

Welcome to the **MHCXGraph** documentation. This page is designed to help you get started with our graph-based structural similarity detection framework for biomolecular complexes.

**MHCXGraph** is an open-source Python package for detecting structurally conserved surface regions across protein structures, with a particular focus on peptide–MHC (pMHC) complexes involved in T cell receptor (TCR) recognition. The method represents biomolecular surfaces as graphs of solvent-exposed residues and identifies similarities through a triad-based geometric and physicochemical comparison framework.

The approach combines discretized representations of local residue environments with a restricted Cartesian product across multiple structures, enabling the detection of conserved spatial patterns that may underlie TCR cross-reactivity. The package also provides interactive visualization interfaces for exploring the resulting association graphs and structural correspondences.

.. Headings:
..  # with overline, for parts
..  * with overline, for chapters
..  = for sections
..  - for subsections
..  ^ for subsubsections
..  " for paragraphs

.. toctree::
  :maxdepth: 1
  :caption: Python Package
  :hidden:

  Installation <package/installation/index>
  Getting Started <gettingstarted/index>
  Execution Configuration (Manifest) <configuration/index>
  Examples <package/examples/index>
  API Reference <package/api_reference/index>
  GitHub repository <https://github.com/CNPEM/MHCXGraph>

.. toctree::
   :maxdepth: 1
   :caption: About

   about/index


