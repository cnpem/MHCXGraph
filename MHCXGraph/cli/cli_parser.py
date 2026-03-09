import argparse


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
        "manifest",
        type=str,
        metavar="<path>",
        help="Path to the JSON manifest with complete set of parameters and settings.",
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

    parser.add_argument(
            "--dashboard",
            help="Automatically open the interactive dashboard in the default web browser after execution.",
            action="store_true",
            default=False,
        )

    return parser.parse_args()
