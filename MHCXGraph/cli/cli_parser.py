import argparse
import os

from .. import __version__


def none_or_float(value):
    if value == "None":
        return None
    return float(value)


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "0"):
        return False
    raise argparse.ArgumentTypeError("Boolean value expected.")


def parse_args():
    """Configura e retorna os argumentos CLI."""
    parser = argparse.ArgumentParser(
        description="A python package for graph-based detection of "
        "structurally similar surface regions in protein complexes, with a "
        "focus on peptide-MHC (pMHC) systems involved in T cell receptor "
        "(TCR) recognition."
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"MHCXGraph v{__version__}",
        help="Show MHCXGraph version number and exit.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        help="Print extra information to standard output.",
        action="store_true",
        default=False,
    )

    subparsers = parser.add_subparsers(
        title="commands",
        dest="command",
        help="Choose a command to run."
    )

    subparsers.required = True

    parser_run = subparsers.add_parser(
        "run",
        help="Run the main MHCXGraph analysis using a manifest file."
    )
    parser_run.add_argument(
        "manifest",
        type=str,
        metavar="<path>",
        help="Path to the JSON manifest with complete set of parameters and settings.",
    )

    parser_run.add_argument(
            "--dashboard",
            help="Automatically open the interactive dashboard in the default web browser after execution.",
            action="store_true",
            default=False,
        )


    parser_heatmap = subparsers.add_parser("heatmap")
    parser_heatmap.add_argument("-i", '--input-dir', required=True, nargs='?', default=os.getcwd())
    parser_heatmap.add_argument("-o", "--output-dir", required=True, help="Output directory")
    parser_heatmap.add_argument("-n", "--name", required=True)

    parser_renumber = subparsers.add_parser(
        "renumber",
        help="Renumber MHC structures using IMGT mapping."
    )

    parser_renumber.add_argument("-i", "--input-dir", required=True, help="Input directory with .pdb/.cif/.mmcif files")
    parser_renumber.add_argument("-o", "--output-dir", required=True, help="Output directory")
    parser_renumber.add_argument("-c", "--mhc-class", required=True, help="The class of MHC to renumber.")
    parser_renumber.add_argument("--warn-score", type=float, default=50.0, help="Warn if the best alignment score is below this value")
    parser_renumber.add_argument("--debug", action="store_true", help="Print debug output")
    parser_renumber.add_argument("--suffix", default="", help="Optional suffix added before file extension in output")

    return parser.parse_args()
