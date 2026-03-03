from __future__ import annotations

import json
import math
import re
import tempfile
from dataclasses import dataclass, field
from typing import Literal, TypedDict, cast

import networkx as nx
import numpy as np
import pandas as pd
from Bio.PDB.Chain import Chain
from Bio.PDB.DSSP import DSSP, residue_max_acc
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.Residue import Residue
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Structure import Structure
from Bio.PDB.vectors import Vector, rotaxis

from MHCXGraph.core.config import DSSPConfig, GraphConfig
from MHCXGraph.core.contact_map import contact_map_from_graph
from MHCXGraph.core.metadata import secondary_structure
from MHCXGraph.utils.logging_utils import get_log

log = get_log()

"""
Structural graph builder for pMHC using Bio.PDB and NetworkX.

Features
--------
- Residue and water centroid computation.
- ASA via Shrake–Rupley and RSA normalization (Tien et al. tables via Bio.PDB).
- Optional RSA via DSSP with fallback for non-canonicals.
- Distance-based graph construction (residue–residue and residue–water).
- Distance matrix and label export for reproducibility/debugging.

Requirements
------------
- biopython >= 1.79
- networkx
- numpy
- pandas
"""

AA1_TO_3 = {
    "A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","Q":"GLN","E":"GLU","G":"GLY",
    "H":"HIS","I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE","P":"PRO","S":"SER",
    "T":"THR","W":"TRP","Y":"TYR","V":"VAL"
}
CANONICAL_AA3 = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"
}
NONCANONICAL_TO_CANONICAL: dict[str, str] = {}
BACKBONE_NAMES = {"N", "CA", "C", "O", "OXT"}
WATER_NAMES = {
    "HOH", "H2O", "WAT", "SOL",
    "TIP3", "TIP3P", "TIP4P", "TIP5P",
    "SPC", "SPCE",
    "DOD", "D2O",
    "OH2",
}
HYDROGEN_RE = re.compile("[123 ]*H.*")

ResidueKind = Literal["canonical_aminoacid", "noncanonical_aminoacid", "ligand", "water"]
ResidueList = list[tuple[str, Residue, ResidueKind, np.ndarray]]


class AtomBundle(TypedDict):
    raw_df: pd.DataFrame
    node_centroids: pd.DataFrame
    ca_cb_map: dict[str, NodeCoords]


class ResidueInfo(TypedDict):
    canonical_aminoacid_residues: ResidueList
    noncanonical_aminoacid_residues: ResidueList
    ligands: ResidueList
    waters: ResidueList


class StructureDict(TypedDict):
    chains: list[Chain]
    residues: ResidueInfo


class NodeCoords(TypedDict):
    ca_coord: tuple[float, float, float]
    cb_coord: tuple[float, float, float]
    cb_is_virtual: bool


def _canonical_of(resname: str) -> str | None:
    """Return canonical 3-letter code for a residue name when available."""
    rn = (resname or "").strip().upper()
    if rn in CANONICAL_AA3:
        return rn
    return NONCANONICAL_TO_CANONICAL.get(rn)


def _is_water(res: Residue) -> bool:
    """Return True for water residues."""
    name = res.get_resname().strip().upper()
    return name in WATER_NAMES


def _heavy_atom_coords(res: Residue) -> np.ndarray:
    """Return heavy-atom coordinates of a residue as (N, 3) array."""
    coords: list[np.ndarray] = []

    for atom in res.get_atoms():
        element = getattr(atom, "element", "").strip().upper()

        if element:
            is_hydrogen = element in {"H", "D", "T"}
        else:
            fullname = atom.fullname.ljust(4)
            is_hydrogen = (not fullname[0].isalpha()) and (fullname[1] in {"H", "D", "T"})

        if not is_hydrogen:
            coords.append(atom.coord)  # type: ignore[attr-defined]

    if not coords:
        coords = [atom.coord for atom in res.get_atoms()]  # type: ignore[attr-defined]

    return np.asarray(coords, dtype=float)


def _node_id(chain_id: str, res: Residue, kind: str = "residue") -> str:
    """Build a stable node identifier: 'A:GLY:42' or 'A:HOH:2001'."""
    _, resseq, icode = res.id
    resname = res.get_resname()
    if kind == "water":
        return f"{chain_id}:HOH:{resseq}{(icode.strip() or '')}"
    return f"{chain_id}:{resname}:{resseq}{(icode.strip() or '')}"


def check_res_inconsistencies(res_tuples):
    inconsistencies = []
    for id_str, residue, _, _ in res_tuples:
        _, resseq, _ = residue.get_id()
        residue_number = int(resseq)
        residue_name_obj = residue.get_resname().strip()
        parts = id_str.split(":")
        if len(parts) == 3:
            _, residue_name_id, resseq_id = parts
            if residue_name_id != residue_name_obj:
                inconsistencies.append(
                    f"{id_str}: residue_name id='{residue_name_id}' vs obj='{residue_name_obj}'"
                )
            if str(resseq_id) != str(residue_number):
                inconsistencies.append(
                    f"{id_str}: residue_number id='{resseq_id}' vs obj='{residue_number}'"
                )

    return inconsistencies


@dataclass
class BuiltGraph:
    """
    Container returned by :class:`PDBGraphBuilder`.

    Attributes
    ----------
    graph : networkx.Graph
        Constructed graph with node/edge attributes.
    residue_index : list of tuple[str, Bio.PDB.Residue]
        Node id to residue pairing for amino acid residues used as the main
        distance base.
    residue_centroids : numpy.ndarray
        Centroids for ``residue_index`` in the same order, shape (N, 3).
    water_index : list of tuple[str, Bio.PDB.Residue]
        Node id to residue pairing for water residues (if included).
    water_centroids : numpy.ndarray or None
        Water centroids in the same order as ``water_index``, shape (W, 3).
    distance_matrix : numpy.ndarray or None
        Pairwise centroid distance matrix for residues in ``residue_index``.
    raw_pdb_df : pandas.DataFrame or None
        Atom-level table used to derive centroids and CA/CB coordinates.
    node_centroids : pandas.DataFrame or None
        DataFrame indexed by ``node_id`` with centroid coordinates
        ``x_coord``, ``y_coord``, ``z_coord``.
    dssp_df : pandas.DataFrame or None
        DSSP summary with an added "rsa" column, aligned to graph nodes.
    """
    graph: nx.Graph
    residue_index: list[tuple[str, Residue]]
    residue_centroids: np.ndarray
    water_index: list[tuple[str, Residue]] = field(default_factory=list)
    water_centroids: np.ndarray | None = None
    distance_matrix: np.ndarray | None = None
    raw_pdb_df: pd.DataFrame | None = None
    node_centroids: pd.DataFrame | None = None
    dssp_df: pd.DataFrame | None = None


class PDBGraphBuilder:
    """
    Build a structural graph from a PDB/mmCIF file.

    Parameters
    ----------
    pdb_path : str
        Path to the structure file.
    config : GraphBuildConfig, optional
        Graph construction options.

    Notes
    -----
    The distance matrix and the node labels are exported to ``.pmhc_tmp/``
    with filenames ``<stem>_distmat.npy`` and ``<stem>_residue_labels.txt``.
    """

    def __init__(self, pdb_path: str, config: GraphConfig | None = None) -> None:
        self.pdb_path = pdb_path
        self.config = config or GraphConfig()
        self.structure: Structure | None = None

    def load(self) -> None:
        """
        Load the structure into memory.

        Raises
        ------
        Exception
            If parsing fails.
        """
        if self.pdb_path.lower().endswith((".cif", ".mmcif")):
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)
        self.structure = parser.get_structure("struct", self.pdb_path)
        log.info("Structure loaded from %s", self.pdb_path)

    def _ensure_pdb_for_dssp(self) -> str:
        """
        Return a path to a PDB file suitable for mkdssp.
    If input is already PDB, return original path.
        If input is mmCIF, write a temporary PDB from self.structure and return it.
        """

        p = self.pdb_path.lower()
        if p.endswith(".pdb"):
            return self.pdb_path

        tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
        tmp_path = tmp.name
        tmp.close()

        io = PDBIO()
        io.set_structure(self.structure)
        io.save(tmp_path)
        return tmp_path

    def _compute_asa_rsa(
            self, res_tuples: ResidueList
            ) -> tuple[dict[str, tuple[float, float | None]], pd.DataFrame | None]:
        """
        Compute per-residue ASA and RSA.

        Parameters
        ----------
        res_tuples : list of tuple
            (node_id, Residue, centroid) for protein residues.

        Returns
        -------
        out : dict
            Mapping node_id -> (ASA_abs, RSA_rel in [0, 1]).
            RSA is guaranteed for non-canonical residues using SR-based fallback.
        dssp_df : pandas.DataFrame or None
            DSSP summary (when ``rsa_method='dssp'``) with "rsa" column added.

        Notes
        -----
        - For ``rsa_method='dssp'``, DSSP is used where available. Residues not
          covered by DSSP or lacking max-ASA are assigned RSA via SR fallback.
        - For ``rsa_method='sr'``, RSA = SR(ASA)/max-ASA for canonical residues;
          non-canonicals are normalized by a structure-wise ASA reference.
        """
        max_acc_table = residue_max_acc[self.config.dssp_acc_array]

        sr = ShrakeRupley(probe_radius=self.config.probe_radius, n_points=self.config.n_points)
        sr.compute(self.structure, level="R")
        out: dict[str, tuple[float, float | None]] = {}

        if self.config.rsa_method == "dssp":
            model = self.structure[self.config.model_index]

            dssp_input_path = self._ensure_pdb_for_dssp()

            dssp = DSSP(
                    model,
                    dssp_input_path,
                    dssp=self.config.dssp_exec,
                    acc_array=self.config.dssp_acc_array,
                    )
            idx2nid = {(res.get_parent().id, res.id): nid for nid, res, _, _ in res_tuples}
            rows: list[dict[str, object]] = []

            for key in dssp.keys():
                nid = idx2nid.get(key)
                if nid is None:
                    continue
                chain_id, (_, resnum, icode) = key
                t = dssp[key]
                dssp_index = t[0]
                aa1 = t[1]
                ss = t[2]
                rsa_rel = float(t[3])
                phi = float(t[4]); psi = float(t[5])
                nh_o_1_relidx, nh_o_1_energy = t[6], t[7]
                o_nh_1_relidx, o_nh_1_energy = t[8], t[9]
                nh_o_2_relidx, nh_o_2_energy = t[10], t[11]
                o_nh_2_relidx, o_nh_2_energy = t[12], t[13]
                aa3 = AA1_TO_3.get(aa1, "UNK")
                max_acc = float(max_acc_table.get(aa3, 0.0))
                asa_abs = rsa_rel * max_acc if max_acc > 0 else 0.0
                out[nid] = (asa_abs, rsa_rel)
                rows.append({
                    "chain": chain_id,
                    "resnum": int(resnum),
                    "icode": (icode or "").strip(),
                    "aa": aa3,
                    "ss": ss if ss != " " else "C",
                    "asa": asa_abs,
                    "phi": phi, "psi": psi,
                    "dssp_index": dssp_index,
                    "NH_O_1_relidx": nh_o_1_relidx, "NH_O_1_energy": nh_o_1_energy,
                    "O_NH_1_relidx": o_nh_1_relidx, "O_NH_1_energy": o_nh_1_energy,
                    "NH_O_2_relidx": nh_o_2_relidx, "NH_O_2_energy": nh_o_2_energy,
                    "O_NH_2_relidx": o_nh_2_relidx, "O_NH_2_energy": o_nh_2_energy,
                    "node_id": nid,
                    })

            dssp_df = pd.DataFrame(rows).set_index("node_id") if rows else None
            if dssp_df is not None and not dssp_df.empty:
                dssp_df["max_acc"] = dssp_df["aa"].map(max_acc_table.get).astype(float)
                with np.errstate(divide="ignore", invalid="ignore"):
                    dssp_df["rsa"] = dssp_df["asa"] / dssp_df["max_acc"]
                    dssp_df.loc[dssp_df["max_acc"] <= 0, "rsa"] = np.nan

            for nid, res, _, _ in res_tuples:
                if nid in out and out[nid][1] is not None:
                    continue
                asa_abs = float(getattr(res, "sasa", 0.0))
                resname = res.get_resname().strip().upper()
                aa_can = _canonical_of(resname)

                if aa_can is not None:
                    max_acc = float(max_acc_table.get(aa_can, 0.0))
                    rsa_rel = asa_abs / max_acc if max_acc > 0 else None

                else:
                    rsa_rel = None

                out[nid] = (asa_abs, None if rsa_rel is None else float(rsa_rel))

            return out, dssp_df

        for nid, res, _, _ in res_tuples:
            asa = float(getattr(res, "sasa", 0.0))
            aa3 = res.get_resname().strip().upper()
            aa_can = _canonical_of(aa3)

            if aa_can is not None:
                max_acc = float(max_acc_table.get(aa_can, 0.0))
                rsa = asa / max_acc if max_acc > 0 else None

            else:
                rsa = None

            out[nid] = (asa, None if rsa is None else float(rsa))

        return out, None

    def _centroid_mask_for_group(self, g: pd.DataFrame) -> pd.Series:
        """
        Compute a boolean mask selecting the atoms that contribute to the centroid
        of a residue or node group.

        The mask is built according to the granularity configuration used in
        the graph construction pipeline. The returned mask is aligned with the input
        DataFrame index and can be used to extract the coordinate subset before
        computing the centroid.

        Parameters
        ----------
        df : pandas.DataFrame
            DataFrame containing atom-level information for a single residue or
            graph node. Must contain at least the column ``"element"`` and the
            coordinate fields ``"x_coord"``, ``"y_coord"``, ``"z_coord"``.

        Returns
        -------
        numpy.ndarray of bool
            Boolean array of shape (N,) where N is the number of atoms in `df`.
            A value of ``True`` indicates that the corresponding atom should be
            included when computing the centroid.

        """
        names = g["atom_name"].str.upper()
        elems = g["element_symbol"].str.upper().fillna("")
        resn = str(g["residue_name"].iloc[0]).upper()

        elem_is_h = elems.isin({"H", "D", "T"})
        name_is_h = names.str.match(HYDROGEN_RE)
        is_hydrogen = elem_is_h | name_is_h

        heavy = ~is_hydrogen

        if resn in WATER_NAMES:
            mask = names.isin({"O", "OW", "OH2"})
            if not mask.any():
                mask = heavy
            return mask

        gran = getattr(self.config, "granularity", "all_atoms")

        if gran == "all_atoms":
            return heavy

        if gran == "backbone":
            return names.isin(BACKBONE_NAMES)

        if gran == "side_chain":
            mask = (~names.isin(BACKBONE_NAMES)) & heavy
            # fallback for GLY / unusual cases with no side chain heavy atoms
            if not mask.any():
                mask = names.eq("CA") | names.eq("CB")
                if not mask.any():
                    mask = heavy
            return mask

        # gran == "ca_only"
        mask = names.eq("CA")
        if not mask.any():
            mask = names.eq("CB")
        if not mask.any():
            mask = heavy
        return mask

    def _centroids_from_raw_df(self, raw_df: pd.DataFrame) -> pd.DataFrame:
        """
        Compute one centroid per residue node from an atom-level table.

        The centroid for each ``node_id`` is calculated using the same atom
        selection rules used for edge construction, as defined by
        ``_centroid_mask_for_group``.

        If the mask for a given residue selects no atoms (e.g., a glycine
        with ``granularity='side_chain'``), all atoms of that residue are
        used as a fallback.

        Parameters
        ----------
        raw_df : pandas.DataFrame
            Atom-level table with at least the following columns:
            ``node_id``, ``x_coord``, ``y_coord``, ``z_coord``,
            ``atom_name``, ``element_symbol``, and ``residue_name``.
            Typically generated by ``_make_raw_df``.

        Returns
        -------
        pandas.DataFrame
            A DataFrame indexed by ``node_id`` with columns:
            ``x_coord``, ``y_coord``, ``z_coord``, containing the centroid
            coordinates for each residue.

            The index name is set to ``"node_id"``. If ``raw_df`` is empty,
            an empty DataFrame with the correct structure is returned.
        """

        if raw_df.empty:
            return pd.DataFrame(
                columns=["x_coord", "y_coord", "z_coord"],
                index=pd.Index([], name="node_id"),
            )

        df = raw_df.set_index("node_id", drop=True)

        def reducer(g: pd.DataFrame) -> pd.Series:
            mask = self._centroid_mask_for_group(g)
            sub = g.loc[mask, ["x_coord", "y_coord", "z_coord"]]
            if sub.empty:
                sub = g[["x_coord", "y_coord", "z_coord"]]

            cx, cy, cz = sub.mean(axis=0).astype(float)
            return pd.Series(
                {"x_coord": cx, "y_coord": cy, "z_coord": cz}
            )

        out = df.groupby(level=0, sort=False).apply(reducer)

        out.index.name = "node_id"
        return out

    def _extract_ca_cb(self, raw_df: pd.DataFrame) -> dict[str, NodeCoords]:
        """
        Extracts CA, CB, N, and C atom coordinates for each residue node.

        Parameters
        ----------
        raw_df : pandas.DataFrame
            Input table with atomic information. Must contain the columns
            'node_id', 'atom_name', 'alt_loc', 'occupancy', 'b_factor',
            'atom_number', 'x_coord', 'y_coord', and 'z_coord'.

        Returns
        -------
        dict of str to NodeCoords
            Mapping from node_id to a dictionary with:
                'ca_coord' : tuple of float
                    The CA atom coordinates or (nan, nan, nan) if missing.
                'cb_coord' : tuple of float
                    The CB atom coordinates or (nan, nan, nan) if missing.
                'cb_is_virtual' : bool
                    Whether the CB was computed as a virtual glycine CB.
        Notes
        -----
        Atom selection follows a deterministic priority:
        empty altloc is preferred over non-empty,
        higher occupancy is preferred,
        lower B-factor is preferred,
        and lower atom_number is preferred.
        """

        if raw_df.empty:
            return {}

        atom_names = cast(pd.Series, raw_df["atom_name"])
        df = raw_df[atom_names.str.upper().isin(["CA", "CB", "N", "C"])].copy()

        if df.empty:
            return {}

        df["AN"] = cast(pd.Series, df["atom_name"]).str.upper()

        df["_alt_rank"] = [1 if str(x).strip() else 0 for x in df["alt_loc"]]
        df = df.fillna({
            "occupancy": 0.0,
            "b_factor": float("inf"),
            "atom_number": float("inf")
        })

        sort_cols: list[str] = ["node_id", "AN", "_alt_rank", "occupancy", "b_factor", "atom_number"]
        sort_asc: list[bool] = [True, True, True, False, True, True]

        df = df.sort_values(by=sort_cols, ascending=sort_asc, kind="mergesort")
        best = df.groupby(["node_id", "AN"], as_index=False).first()

        coords_dict: dict[str, dict[str, tuple[float, float, float]]] = {
            str(nid): {} for nid in raw_df["node_id"].unique()
        }

        for row in best.to_dict("records"):
            coords_dict[str(row["node_id"])][str(row["AN"])] = (
                float(row["x_coord"]),
                float(row["y_coord"]),
                float(row["z_coord"])
            )

        out: dict[str, NodeCoords] = {}
        make_virt_gly = getattr(self.config, "make_virtual_cb_for_gly", True)

        resnames = cast(pd.Series, raw_df.drop_duplicates("node_id").set_index("node_id")["residue_name"])
        resname_map = cast(dict[str, str], resnames.str.upper().to_dict())

        for nid, atoms in coords_dict.items():
            out[nid] = self._process_node_coords(nid, atoms, resname_map, make_virt_gly)

        return out

    def _process_node_coords(self, nid: str,
                             atoms: dict[str, tuple[float, float, float]],
                             resname_map: dict[str, str],
                             make_virt_gly: bool) -> NodeCoords:
        """
        Selects CA and CB coordinates for a given residue and computes a virtual CB
        for glycine when enabled.

        Parameters
        ----------
        nid : str
            Node identifier corresponding to the residue.
        atoms : dict of str to tuple of float
            Dictionary mapping atom names (e.g. 'CA', 'CB', 'N', 'C') to 3D coordinates.
        resname_map : dict of str to str
            Mapping from node_id to residue name.
        make_virt_gly : bool
            Whether to compute a virtual CB for glycine if CB is missing.

        Returns
        -------
        NodeCoords
            Dictionary with keys:
                'ca_coord' : tuple of float
                'cb_coord' : tuple of float
                'cb_is_virtual' : bool
        Notes
        -----
        Virtual CB computation follows standard geometric construction for glycine
        using CA, N, and C coordinates.
        """

        nan3 = (float("nan"), float("nan"), float("nan"))
        ca = atoms.get("CA", nan3)
        cb = atoms.get("CB", nan3)
        is_virt = False

        if math.isnan(cb[0]) and make_virt_gly and resname_map.get(nid) == "GLY":
            n = atoms.get("N", nan3)
            c = atoms.get("C", nan3)
            cb_virt = self._calc_virtual_cb(ca, n, c)

            if not math.isnan(cb_virt[0]):
                cb = cb_virt
                is_virt = True

        return {
            "ca_coord": ca,
            "cb_coord": cb,
            "cb_is_virtual": is_virt
        }

    def _calc_virtual_cb(self,
                         ca: tuple[float, float, float],
                         n: tuple[float, float, float],
                         c: tuple[float, float, float]) -> tuple[float, float, float]:
        """
        Calculates a virtual CB using the Biopython method:
        rotate the N–CA vector -120 degrees around the CA–C axis.

        Parameters
        ----------
        ca : tuple of float
            Coordinates of CA.
        n : tuple of float
            Coordinates of N.
        c : tuple of float
            Coordinates of C.

        Returns
        -------
        tuple of float
            Virtual CB coordinate. NaN triplet if invalid.
        """
        # Check for NaNs
        if any(math.isnan(x) for v in (ca, n, c) for x in v):
            return (float("nan"), float("nan"), float("nan"))

        # Convert to Biopython Vector
        v_ca = Vector(*ca)
        v_n = Vector(*n)
        v_c = Vector(*c)

        # Center at CA
        v_n0 = v_n - v_ca
        v_c0 = v_c - v_ca

        # Rotation of -120 degrees (in radians)
        angle = -math.pi * 120.0 / 180.0

        # Compute rotation matrix around CA–C axis
        rot = rotaxis(angle, v_c0)

        # Rotate N–CA vector
        v_cb0 = v_n0.left_multiply(rot)

        # Translate back to CA position
        v_cb = v_cb0 + v_ca

        return (float(v_cb[0]), float(v_cb[1]), float(v_cb[2]))


    def _make_raw_df(self, s_dict: StructureDict) -> pd.DataFrame:
        """
        Build atom-level DataFrames.

        Parameters
        ----------
        s_dict : StructureDict

        Returns
        -------
        raw_df : pandas.DataFrame
            All atoms (including waters/ligands).
        """
        rows = []
        model_idx = self.config.model_index

        categories: list[tuple[str, str]] = [
                    ("canonical_aminoacid_residues", "ATOM"),
                    ("noncanonical_aminoacid_residues", "ATOM"),
                    ("ligands", "HETATM"),
                    ("waters", "HETATM")
        ]

        for kind, record_name in categories:

            res_list = s_dict["residues"][kind]

            for node_id, res, _, _ in res_list:
                _, resseq, icode = res.id
                resname = res.get_resname()
                chain_id = res.parent.id
                residue_id = f"{chain_id}:{int(resseq)}{(icode.strip() or '')}"

                for atom in res.get_atoms():
                    try:
                        serial = atom.get_serial_number()
                    except Exception:
                        serial = None

                    name = atom.get_name().strip()

                    try:
                        alt = atom.get_altloc()
                    except Exception:
                        alt = ""

                    element = getattr(atom, "element", None)
                    if not element:
                        element = name[0] if name else ""

                    x, y, z = map(float, atom.coord)  # type: ignore[attr-defined]
                    occ = atom.get_occupancy()
                    bfac = atom.get_bfactor()

                    rows.append({
                        "record_name": record_name,
                        "atom_number": serial,
                        "atom_name": name,
                        "alt_loc": (alt or "").strip(),
                        "residue_name": resname,
                        "chain_id": chain_id,
                        "residue_number": int(resseq),
                        "insertion": (icode.strip() or ""),
                        "x_coord": x, "y_coord": y, "z_coord": z,
                        "occupancy": occ, "b_factor": bfac,
                        "element_symbol": (str(element).upper() if element else ""),
                        "charge": None,
                        "model_idx": model_idx,
                        "node_id": node_id,
                        "residue_id": residue_id,
                        "kind": kind
                    })

        raw_df = pd.DataFrame(rows)
        cols = [
            "record_name", "atom_number", "atom_name", "alt_loc", "residue_name", "chain_id",
            "residue_number", "insertion", "x_coord", "y_coord", "z_coord", "occupancy",
            "b_factor", "element_symbol", "charge", "model_idx", "node_id", "residue_id", "kind"
        ]

        raw_df = raw_df.reindex(columns=cols)

        return raw_df

    def _build_structure_dict(self) -> StructureDict:
        if self.structure is None:
            raise Exception("Structure not loaded")

        s_dict: StructureDict = {
            "chains": [],
            "residues": {
                "canonical_aminoacid_residues": [],
                "noncanonical_aminoacid_residues": [],
                "ligands": [],
                "waters": []
            }
        }

        model = self.structure[self.config.model_index]

        s_dict["chains"] = self._get_requested_chains(model)

        for ch in s_dict["chains"]:
            for res in ch.get_residues():
                self._process_residue(ch, res, s_dict["residues"])

        return s_dict

    def _get_requested_chains(self, model) -> list:
        """Helper to filter and validate requested chains."""
        if self.config.chains is None:
            return list(model)

        wanted = set(self.config.chains)
        chains = [ch for ch in model if ch.id in wanted]

        if not chains and not self.config.allow_empty_chains:
            msg = f"None of the requested chains were found: {self.config.chains}"
            raise ValueError(msg)

        return chains

    def _process_residue(self, ch, res, residues_dict: ResidueInfo) -> None:
        """Helper to categorize and store a single residue."""
        coords_heavy = _heavy_atom_coords(res)

        if is_aa(res, standard=False):
            key = "noncanonical_aminoacid_residues"
            kind = "noncanonical_aminoacid"

            if is_aa(res, standard=True):
                key = "canonical_aminoacid_residues"
                kind = "canonical_aminoacid"

            id_str = _node_id(ch.id, res)
            self._validate_aminoacid_id(id_str, res)

            residues_dict[key].append((id_str, res, kind, coords_heavy))

        elif _is_water(res):
            id_str = _node_id(ch.id, res, kind="water")
            residues_dict["waters"].append((id_str, res, "water", coords_heavy))

        else:
            id_str = _node_id(ch.id, res)
            residues_dict["ligands"].append((id_str, res, "ligand", coords_heavy))

    def _validate_aminoacid_id(self, id_str: str, res) -> None:
        """Helper to check for residue ID mismatches."""
        _, resseq, _ = res.get_id()
        residue_number = int(resseq)
        residue_name_obj = res.get_resname().strip()
        parts = id_str.split(":")

        if len(parts) == 3:
            _, residue_name_id, resseq_id = parts
            if residue_name_id != residue_name_obj:
                log.warning(
                    f"{id_str}: residue_name id='{residue_name_id}' vs obj='{residue_name_obj}'"
                )
            if str(resseq_id) != str(residue_number):
                log.warning(
                    f"{id_str}: residue_number id='{resseq_id}' vs obj='{residue_number}'"
                )

    def build_graph(self) -> BuiltGraph:
        """
        Run the full pipeline: load → select chains → ASA/RSA → distances → graph.

        Returns
        -------
        BuiltGraph
            Graph object and associated tables/arrays.

        Notes
        -----
        - The residue–residue distance matrix and node labels are saved under
          ``.pmhc_tmp/<stem>_distmat.npy`` and ``.pmhc_tmp/<stem>_residue_labels.txt``.
        - Waters are added as nodes with ``rsa=1.0`` and connected to nearby residues.
        """
        if self.structure is None:
            self.load()

        structure_dict = self._build_structure_dict()
        self.structure_dict = structure_dict

        canonical_nodes = structure_dict["residues"]["canonical_aminoacid_residues"]
        noncanonical_nodes = structure_dict["residues"]["noncanonical_aminoacid_residues"]
        ligand_nodes_all = structure_dict["residues"]["ligands"]
        water_nodes_all = structure_dict["residues"]["waters"]

        aa_nodes = canonical_nodes.copy()
        if self.config.include_noncanonical_residues:
            aa_nodes += noncanonical_nodes.copy()

        ligand_nodes = ligand_nodes_all if self.config.include_ligands else []
        water_nodes = water_nodes_all if self.config.include_waters else []

        all_nodes = aa_nodes + ligand_nodes + water_nodes

        raw_pdb_df = self._make_raw_df(structure_dict)
        ca_cb_map = self._extract_ca_cb(raw_pdb_df)
        node_centroids = self._centroids_from_raw_df(raw_pdb_df)

        asa_rsa: dict[str, tuple[float, float | None]] = {}
        dssp_df: pd.DataFrame | None = None
        if self.config.compute_rsa and aa_nodes:
            asa_rsa, dssp_df = self._compute_asa_rsa(aa_nodes)

        res_ids = [nid for nid, _, _, _ in aa_nodes]
        res_objects = [res for _, res, _, _ in aa_nodes]

        cent_rows = node_centroids.reindex(res_ids)
        if cent_rows[["x_coord", "y_coord", "z_coord"]].isna().any(axis=None):
            missing = cent_rows[cent_rows.isna().any(axis=1)].index.tolist()
            raise ValueError(f"Missing centroid coords for node_ids: {missing[:10]}")

        res_centroids = cent_rows[["x_coord", "y_coord", "z_coord"]].to_numpy(float)

        diff = res_centroids[:, None, :] - res_centroids[None, :, :]
        dist_mat = np.sqrt(np.sum(diff * diff, axis=2))

        G = nx.Graph()
        nan3 = (float("nan"), float("nan"), float("nan"))

        for nid, res, kind, _ in all_nodes:
            asa, rsa = asa_rsa.get(nid, (None, None))

            parent = res.get_parent()
            chain = parent.id if parent is not None else None

            extra = ca_cb_map.get(nid, {})
            ca_coord = extra.get("ca_coord", nan3)
            cb_coord = extra.get("cb_coord", nan3)
            cb_is_virtual = bool(extra.get("cb_is_virtual", False))

            c = node_centroids.loc[nid]
            cent = np.array([c["x_coord"], c["y_coord"], c["z_coord"]], dtype=float)

            if kind == "water":
                resname = "HOH"
                asa_val = None
                rsa_val = 1.0
                ca_coord = nan3
                cb_coord = nan3
                cb_is_virtual = False
            else:
                resname = res.get_resname()
                asa_val = None if asa is None else float(asa)
                rsa_val = None if rsa is None else float(rsa)

            if kind != "canonical_aminoacid":
                log.info("[EXTRA] %s kind=%s rsa=%s asa=%s", nid, kind, rsa_val, asa_val)

            G.add_node(
                nid,
                kind=kind,
                chain=chain,
                resname=resname,
                resseq=int(res.id[1]),
                icode=(res.id[2].strip() or ""),
                centroid=tuple(float(x) for x in cent),
                coords=cent,
                ca_coord=ca_coord,
                cb_coord=cb_coord,
                cb_is_virtual=cb_is_virtual,
                asa=asa_val,
                rsa=rsa_val,
            )

        residue_map = {nid: res for nid, res, _, _ in aa_nodes}
        secondary_structure(
            G,
            dssp_config=DSSPConfig(executable="mkdssp"),
            structure=self.structure,
            residue_map=residue_map,
            pdb_path=str(self.pdb_path),
        )

        cut = float(self.config.residue_distance_cutoff)
        iu, ju = np.triu_indices_from(dist_mat, k=1)
        keep = dist_mat[iu, ju] <= cut
        for i, j in zip(iu[keep], ju[keep], strict=True):
            G.add_edge(res_ids[i], res_ids[j], distance=float(dist_mat[i, j]), kind="res-res")

        water_ids: list[str] = []
        water_objs: list[Residue] = []
        water_centroids: np.ndarray | None = None

        if water_nodes:
            water_ids = [nid for nid, _, _, _ in water_nodes]
            water_objs = [res for _, res, _, _ in water_nodes]

            w_rows = node_centroids.reindex(water_ids)
            if w_rows[["x_coord", "y_coord", "z_coord"]].isna().any(axis=None):
                missing = w_rows[w_rows.isna().any(axis=1)].index.tolist()
                raise ValueError(f"Missing water centroid coords for node_ids: {missing[:10]}")

            water_centroids = w_rows[["x_coord", "y_coord", "z_coord"]].to_numpy(float)

            wd = water_centroids[:, None, :] - res_centroids[None, :, :]
            wdist = np.sqrt(np.sum(wd * wd, axis=2))

            wcut = float(self.config.water_distance_cutoff)
            hits = np.where((wdist > 0.0) & (wdist <= wcut))
            for wi, rj in zip(hits[0], hits[1], strict=True):
                G.add_edge(water_ids[wi], res_ids[rj], distance=float(wdist[wi, rj]), kind="wat-res")

        if ligand_nodes:
            ligand_ids = [nid for nid, _, _, _ in ligand_nodes]

            l_rows = node_centroids.reindex(ligand_ids)
            if l_rows[["x_coord", "y_coord", "z_coord"]].isna().any(axis=None):
                missing = l_rows[l_rows.isna().any(axis=1)].index.tolist()
                raise ValueError(f"Missing ligand centroid coords for node_ids: {missing[:10]}")

            ligand_centroids = l_rows[["x_coord", "y_coord", "z_coord"]].to_numpy(float)

            ld = ligand_centroids[:, None, :] - res_centroids[None, :, :]
            ldist = np.sqrt(np.sum(ld * ld, axis=2))

            hits = np.where((ldist > 0.0) & (ldist <= cut))
            for li, rj in zip(hits[0], hits[1], strict=True):
                G.add_edge(ligand_ids[li], res_ids[rj], distance=float(ldist[li, rj]), kind="lig-res")

        (
            contact_map,
            contact_node_order,
            residue_map_dict,
            residue_map_dict_all,
        ) = contact_map_from_graph(
            G,
            granularity=getattr(self.config, "granularity", "all_atoms"),
            exclude_kinds=(),
            fallback_to_centroid=True,
        )

        G.graph["contact_map"] = contact_map
        G.graph["contact_map_node_order"] = contact_node_order
        G.graph["residue_map_dict"] = residue_map_dict
        G.graph["residue_map_dict_all"] = residue_map_dict_all

        rsa_series = pd.Series(
            {nid: (float(d.get("rsa")) if d.get("rsa") is not None else np.nan) for nid, d in G.nodes(data=True)},
            name="rsa",
        )

        if dssp_df is None:
            dssp_df_all = rsa_series.to_frame()
        else:
            dssp_df_all = dssp_df.copy()
            missing_nodes = [nid for nid in G.nodes if nid not in dssp_df_all.index]
            if missing_nodes:
                dssp_df_all = pd.concat([dssp_df_all, pd.DataFrame(index=missing_nodes)], axis=0)
            if "rsa" not in dssp_df_all.columns:
                dssp_df_all["rsa"] = np.nan
            rsa_aligned = rsa_series.reindex(dssp_df_all.index)
            dssp_df_all["rsa"] = dssp_df_all["rsa"].where(dssp_df_all["rsa"].notna(), rsa_aligned)

        dssp_df = dssp_df_all

        built = BuiltGraph(
            graph=G,
            residue_index=list(zip(res_ids, res_objects, strict=True)),
            residue_centroids=res_centroids,
            water_index=list(zip(water_ids, water_objs, strict=True)),
            water_centroids=water_centroids,
            distance_matrix=dist_mat if self.config.store_distance_matrix else None,
            raw_pdb_df=raw_pdb_df,
            node_centroids=node_centroids,
            dssp_df=dssp_df,
        )
        return built

    @staticmethod
    def to_graphml(G: nx.Graph, path: str) -> None:
        """
        Export the graph to GraphML.

        Parameters
        ----------
        G : networkx.Graph
            Graph to export.
        path : str
            Output path.
        """
        nx.write_graphml(G, path)
        log.info("GraphML saved to %s", path)

    @staticmethod
    def to_json(G: nx.Graph, path: str) -> None:
        """
        Export the graph to JSON (nodes/edges with attributes).

        Parameters
        ----------
        G : networkx.Graph
            Graph to export.
        path : str
            Output path.
        """
        def _convert(v):
            if isinstance(v, (np.floating, np.float32, np.float64)):
                return float(v)
            if isinstance(v, (np.integer, np.int32, np.int64)):
                return int(v)
            if isinstance(v, (np.ndarray,)):
                return v.tolist()
            return v
        data = {
            "nodes": [
                {"id": n, **{k: _convert(v) for k, v in d.items()}}
                for n, d in G.nodes(data=True)
            ],
            "edges": [
                {"u": u, "v": v, **{k: _convert(w) for k, w in d.items()}}
                for u, v, d in G.edges(data=True)
            ],
        }
        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
        log.info("JSON saved to %s", path)
