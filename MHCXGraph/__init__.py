"""
*MHCXGraph* is a Python package for working with MHC-peptide binding data
and generating immunogenicity graphs.

Command Line Interface
----------------------
Usage: MHCXGraph [-h] [--version] [-v] <path>

A python package for graph-based detection of structurally similar surface regions in protein complexes, with
a focus on peptide-MHC (pMHC) systems involved in T cell receptor (TCR) recognition.

positional arguments:
<path>         Path to the JSON manifest with complete set of parameters and settings.

options:
-h, --help     show this help message and exit
--version      Show MHCXGraph version number and exit.
-v, --verbose  Print extra information to standard output.

Example
-------
$ MHCXGraph examples/minimal/manifest.json

See also
--------
* GitHub repository: https://github.com/cnpem/MHCXGraph
* Documentation: https://cnpem.github.io/MHCXGraph
"""

__version__ = "0.1.1"
