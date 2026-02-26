from __future__ import annotations

import json
import logging
import tempfile
from dataclasses import dataclass, field

import networkx as nx
import numpy as np
import pandas as pd
from Bio.PDB import PDBIO, MMCIFParser, PDBParser, ShrakeRupley
from Bio.PDB.Chain import Chain
from Bio.PDB.DSSP import DSSP, residue_max_acc
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

from MHCXGraph.core.config import DSSPConfig, GraphConfig
from MHCXGraph.core.contact_map import contact_map_from_graph
from MHCXGraph.core.metadata import secondary_structure

log = logging.getLogger(__name__)

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


def _canonical_of(resname: str) -> str | None:
    """Return canonical 3-letter code for a residue name when available."""
    rn = (resname or "").strip().upper()
    if rn in CANONICAL_AA3:
        return rn
    return NONCANONICAL_TO_CANONICAL.get(rn)


def _is_protein_residue(res: Residue) -> bool:
    """Return True if residue is a standard (ATOM) amino-acid-like entry."""
    hetflag, _, _ = res.id
    return hetflag == " " and is_aa(res, standard=False)


def _is_water(res: Residue) -> bool:
    """Return True for water residues (HOH)."""
    return res.get_resname() == "HOH"


def _heavy_atom_coords(res: Residue) -> np.ndarray:
    """Return heavy-atom coordinates of a residue as (N, 3) array."""
    coords: list[np.ndarray] = []
    for atom in res.get_atoms():
        element = getattr(atom, "element", None)
        if element is None:
            if atom.get_name().upper().startswith("H"):
                continue
        else:
            if element.upper() == "H":
                continue
        coords.append(atom.coord)  # type: ignore[attr-defined]
    if not coords:
        coords = [atom.coord for atom in res.get_atoms()]  # type: ignore[attr-defined]
    return np.asarray(coords, dtype=float)


def _centroid(coords: np.ndarray) -> np.ndarray:
    """Return centroid of an (N, 3) coordinate array."""
    return coords.mean(axis=0)

def _node_id(chain_id: str, res: Residue, kind: str = "residue") -> str:
    """Build a stable node identifier: 'A:GLY:42' or 'A:HOH:2001'."""
    _, resseq, icode = res.id
    resname = res.get_resname()
    if kind == "water":
        return f"{chain_id}:HOH:{resseq}{(icode.strip() or '')}"
    return f"{chain_id}:{resname}:{resseq}{(icode.strip() or '')}"

def res_tuples_to_df(res_tuples):
    """
    Convert (node_id, Residue, centroid) tuples to a DataFrame.

    Parameters
    ----------
    res_tuples : list of tuple
        Each tuple is (node_id: str, residue: Bio.PDB.Residue, centroid: (3,) array).

    Returns
    -------
    df : pandas.DataFrame
        Columns: chain_id, residue_number, residue_name, insertion, x_coord, y_coord, z_coord.
    inconsistencies : list of str
        Messages describing id vs. object mismatches (if any).
    """
    rows, inconsistencies = [], []
    for id_str, residue, centroid in res_tuples:
        try:
            chain_id = id_str.split(":")[0]
        except Exception:
            chain_id = None
        hetflag, resseq, icode = residue.get_id()
        residue_number = int(resseq)
        insertion = "" if icode == " " else str(icode)
        residue_name_obj = residue.get_resname().strip()
        x, y, z = map(float, centroid)
        rows.append({
            "chain_id": chain_id,
            "residue_number": residue_number,
            "residue_name": residue_name_obj,
            "insertion": insertion,
            "x_coord": x,
            "y_coord": y,
            "z_coord": z,
        })
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
    df = pd.DataFrame(rows, columns=[
        "chain_id","residue_number","residue_name","insertion",
        "x_coord","y_coord","z_coord"
    ])
    return df, inconsistencies




@dataclass
class BuiltGraph:
    """
    Container returned by :class:`PDBGraphBuilder`.

    Attributes
    ----------
    graph : networkx.Graph
        Constructed graph with node/edge attributes.
    residue_index : list of (str, Bio.PDB.Residue)
        Node id to residue pairing for protein residues.
    residue_centroids : ndarray, shape (N, 3)
        Residue centroids used for distance calculations.
    water_index : list of (str, Bio.PDB.Residue)
        Node id to residue pairing for waters (if included).
    water_centroids : ndarray or None
        Water centroids, if any.
    distance_matrix : ndarray or None
        Residue–residue centroid distance matrix when requested.
    raw_pdb_df, pdb_df, rgroup_df : pandas.DataFrame or None
        Atom-level tables.
    dssp_df : pandas.DataFrame or None
        DSSP summary with an added "rsa" column; includes waters as rows (rsa=1.0).
    """
    graph: nx.Graph
    residue_index: list[tuple[str, Residue]]
    residue_centroids: np.ndarray
    water_index: list[tuple[str, Residue]] = field(default_factory=list)
    water_centroids: np.ndarray | None = None
    distance_matrix: np.ndarray | None = None
    raw_pdb_df: pd.DataFrame | None = None
    pdb_df: pd.DataFrame | None = None
    rgroup_df: pd.DataFrame | None = None
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
        assert self.structure is not None

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
        self, res_tuples: list[tuple[str, Residue, np.ndarray]]
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
        assert self.structure is not None

        sr = ShrakeRupley(probe_radius=self.config.probe_radius, n_points=self.config.n_points)
        sr.compute(self.structure, level="R")

        def _asa_ref_from_structure() -> float:
            vals: list[float] = []
            for _, res, _ in res_tuples:
                aa3 = res.get_resname().strip().upper()
                if aa3 in CANONICAL_AA3:
                    vals.append(float(getattr(res, "sasa", 0.0)))
            return max(vals) if vals else 200.0

        asa_ref = float(_asa_ref_from_structure())

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
            idx2nid = {(res.get_parent().id, res.id): nid for nid, res, _ in res_tuples}
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

            # Fill missing/undefined RSA using SR fallback
            for nid, res, _ in res_tuples:
                if nid in out and out[nid][1] is not None:
                    continue
                asa_abs = float(getattr(res, "sasa", 0.0))
                resname = res.get_resname().strip().upper()
                aa_can = _canonical_of(resname)

                if aa_can is not None:
                    max_acc = float(max_acc_table.get(aa_can, 0.0))
                    if max_acc > 0:
                        rsa_rel = asa_abs / max_acc
                    else:
                        rsa_rel = None
                else:
                    rsa_rel = None

                out[nid] = (asa_abs, None if rsa_rel is None else float(rsa_rel))

            return out, dssp_df

        # SR-only path
        for nid, res, _ in res_tuples:
            asa = float(getattr(res, "sasa", 0.0))
            aa3 = res.get_resname().strip().upper()
            aa_can = _canonical_of(aa3)

            if aa_can is not None:
                max_acc = float(max_acc_table.get(aa_can, 0.0))
                if max_acc > 0:
                    rsa = asa / max_acc
                else:
                    rsa = None
            else:
                rsa = None

            out[nid] = (asa, None if rsa is None else float(rsa))

        return out, None

    def _centroid_mask_for_group(self, g: pd.DataFrame) -> pd.Series:
        """
        Return a boolean mask selecting which atoms of a residue group `g` are
        used to compute the centroid, according to `config.node_granularity`.

        Notes
        -----
        - Heavy atoms => element != 'H'
        - Water residues: prefer O/OW/OH2; fallback to any heavy atom.
        - side_chain: if empty (e.g., GLY), fallback to CA (or CB).
        - ca_only: use CA; if missing, CB; if still missing, first heavy atom.
        """
        names = g["atom_name"].str.upper()
        elems = g["element_symbol"].str.upper().fillna("")
        resn  = str(g["residue_name"].iloc[0]).upper()
        heavy = elems != "H"

        # handle waters up front
        if resn == "HOH":
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
            mask = heavy  # ultimate fallback
        return mask

    def _centroids_from_raw_df(self, raw_df: pd.DataFrame) -> pd.DataFrame:
        """
        Compute one centroid per node_id using the same granularity rules as edges.

        Returns a DataFrame indexed by node_id with columns x_coord,y_coord,z_coord.
        """
        if raw_df.empty:
            return pd.DataFrame(columns=["x_coord", "y_coord", "z_coord"]).set_index(pd.Index([], name="node_id"))

        def _centroid_for_group(g: pd.DataFrame) -> pd.Series:
            mask = self._centroid_mask_for_group(g)
            sub = g.loc[mask, ["x_coord", "y_coord", "z_coord"]]
            if sub.empty:
                sub = g[["x_coord", "y_coord", "z_coord"]]
            cx, cy, cz = sub.mean(axis=0).values.astype(float)
            return pd.Series({"x_coord": cx, "y_coord": cy, "z_coord": cz})

        out = raw_df.groupby("node_id", sort=False, group_keys=False).apply(_centroid_for_group)
        out.index.name = "node_id"
        return out

    def _extract_ca_cb(self, raw_df: pd.DataFrame) -> dict[str, dict[str, tuple[float, float, float]]]:
        """
        Para cada node_id, retorna tuplas (x,y,z) para CA e CB.
        Se não existir, retorna (nan, nan, nan).
        Critério: prioriza altloc vazio sobre não vazio; maior occupancy; menor b_factor; menor atom_number.
        """
        if raw_df.empty:
            return {}

        df = raw_df.copy()
        df["AN"] = df["atom_name"].str.upper()
        df = df[df["AN"].isin(["CA", "CB", "N", "C"])]

        if df.empty:
            return {}

        # rank de escolha estável de conformero/altloc
        df["_alt_rank"] = (df["alt_loc"].astype(str).str.strip() != "").astype(int)
        df["_occ"] = df["occupancy"].fillna(0.0)
        df["_bfac"] = df["b_factor"].fillna(np.inf)
        df["_anum"] = df["atom_number"].fillna(np.inf)

        df = df.sort_values(
            ["node_id", "AN", "_alt_rank", "_occ", "_bfac", "_anum"],
            ascending=[True, True, True, False, True, True],
            kind="mergesort",
        )

        best = df.groupby(["node_id", "AN"], as_index=False).first()

        pick = lambda an: best[best["AN"] == an].set_index("node_id")[["x_coord", "y_coord", "z_coord"]]
        ca_tbl = pick("CA")
        cb_tbl = pick("CB")
        n_tbl  = pick("N")
        c_tbl  = pick("C")

        all_nodes = pd.Index(raw_df["node_id"].unique())
        ca_tbl = ca_tbl.reindex(all_nodes)
        cb_tbl = cb_tbl.reindex(all_nodes)
        n_tbl  = n_tbl.reindex(all_nodes)
        c_tbl  = c_tbl.reindex(all_nodes)

        # também precisamos do nome do resíduo por node_id
        resname_by_node = (
            raw_df.drop_duplicates("node_id")
                  .set_index("node_id")["residue_name"].str.upper()
                  .reindex(all_nodes)
        )

        nan3 = (float("nan"), float("nan"), float("nan"))
        out: dict[str, dict[str, tuple[float, float, float]]] = {}

        def _to_tuple(row) -> tuple[float, float, float]:
            if row is None or not hasattr(row, "values"):
                return nan3
            vals = row.values.astype(float)
            if np.isnan(vals).any():
                return nan3
            return (float(vals[0]), float(vals[1]), float(vals[2]))

        def _normalize(v: np.ndarray) -> np.ndarray:
            n = np.linalg.norm(v)
            if n == 0 or not np.isfinite(n):
                return v
            return v / n


        for nid in all_nodes:
            ca_t = _to_tuple(ca_tbl.loc[nid] if nid in ca_tbl.index else None)
            cb_t = _to_tuple(cb_tbl.loc[nid] if nid in cb_tbl.index else None)

            # Se CB está ausente e resíduo é GLY e config permite, tenta CB virtual
            if (np.isnan(cb_t[0]) or np.isnan(cb_t[1]) or np.isnan(cb_t[2])):
                if getattr(self.config, "make_virtual_cb_for_gly", True) and resname_by_node.get(nid, "") == "GLY":
                    ca_row = ca_tbl.loc[nid] if nid in ca_tbl.index else None
                    n_row  = n_tbl.loc[nid]  if nid in n_tbl.index  else None
                    c_row  = c_tbl.loc[nid]  if nid in c_tbl.index  else None

                    if all(r is not None and not np.isnan(r.values.astype(float)).any() for r in [ca_row, n_row, c_row]):
                        rCA = ca_row.values.astype(float)
                        rN  = n_row.values.astype(float)
                        rC  = c_row.values.astype(float)
                        nvec = _normalize(rN - rCA)
                        cvec = _normalize(rC - rCA)
                        b = _normalize(nvec + cvec)
                        cb_virtual = rCA - 1.522 * b
                        cb_t = (float(cb_virtual[0]), float(cb_virtual[1]), float(cb_virtual[2]))
                        out[nid] = {"ca_coord": ca_t, "cb_coord": cb_t, "cb_is_virtual": True}
                    else:
                        out[nid] = {"ca_coord": ca_t, "cb_coord": nan3, "cb_is_virtual": False}
                else:
                    out[nid] = {"ca_coord": ca_t, "cb_coord": cb_t, "cb_is_virtual": False}
            else:
                out[nid] = {"ca_coord": ca_t, "cb_coord": cb_t, "cb_is_virtual": False}

        return out

    def _centroid_pdb_df_from_raw(self, raw_df: pd.DataFrame) -> pd.DataFrame:
        """
        Build a centroid-style pdb_df: one representative atom per residue `node_id`,
        with (x, y, z) overwritten by the residue centroid computed from a
        configurable subset of atoms.

        Granularity
        -----------
        Controlled by `self.config.node_granularity`:
          - "all_atoms": all heavy atoms (element != 'H').
          - "backbone": heavy backbone atoms only (N, CA, C, O, OXT).
          - "side_chain": heavy side-chain atoms only; fallback to CA/CB if empty.
          - "ca_only": use CA (or CB if CA missing); centroid is that atom.

        Representative row
        -------------------
        We keep exactly one row per `node_id` by preferring CA, then CB, then
        water O/OW, then first heavy atom (stable). The (x,y,z) of that row are
        overwritten with the centroid computed above.
        """
        if raw_df.empty:
            return raw_df.copy()

        df = raw_df.copy()

        # Compute centroid per node_id with a custom mask per residue group
        def _centroid_for_group(g: pd.DataFrame) -> pd.Series:
            mask = self._centroid_mask_for_group(g)
            sub = g.loc[mask, ["x_coord", "y_coord", "z_coord"]]
            if sub.empty:  # paranoid fallback
                sub = g[["x_coord", "y_coord", "z_coord"]]
            cx, cy, cz = sub.mean(axis=0).values.astype(float)
            return pd.Series({"cx": cx, "cy": cy, "cz": cz})

        centroids = df.groupby("node_id", sort=False, group_keys=False).apply(_centroid_for_group)

        # Representative atom selection (stable) — prefer CA, then CB, water O/OW, else first heavy
        def _score(row) -> int:
            an = str(row["atom_name"]).upper()
            resn = str(row["residue_name"]).upper()
            if an == "CA":
                return 0
            if an == "CB":
                return 1
            if resn == "HOH" and an in {"O", "OW", "OH2"}:
                return 2
            return 3

        df["_repr_score"] = df.apply(_score, axis=1)

        df_sorted = df.sort_values(
            by=["_repr_score", "model_idx", "chain_id", "residue_number", "atom_number"],
            kind="mergesort",
        )

        reps = (
            df_sorted.groupby("node_id", as_index=False)
            .first()
            .drop(columns=["_repr_score"])
        )

        reps = reps.merge(centroids, left_on="node_id", right_index=True, how="left")
        for c_src, c_dst in (("cx", "x_coord"), ("cy", "y_coord"), ("cz", "z_coord")):
            reps[c_dst] = reps[c_src].where(reps[c_src].notna(), reps[c_dst])
            reps.drop(columns=[c_src], inplace=True)

        # Keep protein and waters; drop other ligands unless you want them later
        reps = reps.copy()

        cols = [
            "record_name","atom_number","atom_name","alt_loc","residue_name","chain_id",
            "residue_number","insertion","x_coord","y_coord","z_coord","occupancy",
            "b_factor","element_symbol","charge","model_idx","node_id","residue_id"
        ]
        reps = reps.reindex(columns=cols).reset_index(drop=True)

        return reps.reset_index(drop=True)

    def _make_atom_tables(self, chains: list[Chain]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Build atom-level DataFrames.

        Parameters
        ----------
        chains : list of Chain
            Chains to traverse.

        Returns
        -------
        raw_df : pandas.DataFrame
            All atoms (including waters/ligands).
        pdb_df : pandas.DataFrame
            Protein heavy atoms (no H, no waters).
        rgroup_df : pandas.DataFrame
            Side-chain heavy atoms (backbone removed).
        """
        rows = []
        model_idx = self.config.model_index
        for ch in chains:
            for res in ch.get_residues():
                hetflag, resseq, icode = res.id
                resname = res.get_resname()
                chain_id = ch.id
                is_water = (resname == "HOH")
                is_std = _is_protein_residue(res)
                record_name = "ATOM" if is_std else "HETATM"
                node_id = _node_id(chain_id, res, kind="water" if is_water else "residue")
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
                    })
        raw_df = pd.DataFrame(rows)
        pdb_df = raw_df[
            (raw_df["record_name"] == "ATOM") &
            (raw_df["element_symbol"] != "H") &
            (raw_df["residue_name"] != "HOH")
        ].copy()
        rgroup_df = pdb_df[~pdb_df["atom_name"].str.upper().isin(BACKBONE_NAMES)].copy()
        cols = [
            "record_name","atom_number","atom_name","alt_loc","residue_name","chain_id",
            "residue_number","insertion","x_coord","y_coord","z_coord","occupancy",
            "b_factor","element_symbol","charge","model_idx","node_id","residue_id"
        ]
        raw_df = raw_df.reindex(columns=cols)
        pdb_df = pdb_df.reindex(columns=cols)
        rgroup_df = rgroup_df.reindex(columns=cols)
        return raw_df, pdb_df, rgroup_df

    def _select_chains(self) -> list[Chain]:
        """
        Select chains according to configuration.

        Returns
        -------
        list of Chain
            Chains to be used.

        Raises
        ------
        ValueError
            If requested chains are not found and ``allow_empty_chains=False``.
        """
        assert self.structure is not None, "Structure not loaded"
        model = self.structure[self.config.model_index]
        if self.config.chains is None:
            return [ch for ch in model]
        chains_found: list[Chain] = []
        wanted = set(self.config.chains)
        for ch in model:
            if ch.id in wanted:
                chains_found.append(ch)
        if not chains_found and not self.config.allow_empty_chains:
            raise ValueError(f"None of the requested chains were found: {self.config.chains}")
        return chains_found

    def _collect_residues(self, chains: list[Chain]) -> list[tuple[str, Residue, np.ndarray]]:
        """
        Collect protein residues and centroids.

        Returns
        -------
        list of tuple
            (node_id, Residue, centroid) for standard residues only.
        """
        out: list[tuple[str, Residue, np.ndarray]] = []
        for ch in chains:
            for res in ch.get_residues():
                if not _is_protein_residue(res):
                    continue
                coords = _heavy_atom_coords(res)
                cent = _centroid(coords)
                out.append((_node_id(ch.id, res), res, cent))
        return out

    def _collect_waters(self, chains: list[Chain]) -> list[tuple[str, Residue, np.ndarray]]:
        """
        Collect water residues and centroids when enabled.

        Returns
        -------
        list of tuple
            (node_id, Residue, centroid) for waters (HOH). Empty if disabled.
        """
        if not self.config.include_waters:
            return []
        out: list[tuple[str, Residue, np.ndarray]] = []
        for ch in chains:
            for res in ch.get_residues():
                if _is_water(res):
                    coords = _heavy_atom_coords(res)
                    cent = _centroid(coords)
                    out.append((_node_id(ch.id, res, kind="water"), res, cent))
        return out

    def _collect_nonprotein_residues(
        self,
        chains: list[Chain],
    ) -> list[tuple[str, Residue, np.ndarray, str]]:
        """
        Collects non-protein residues (excluding water) and classifies them as:

        - 'noncanonical_residue': HETATM entries that are still amino acids (is_aa with standard=False)
        and when config.include_noncanonical_residues = True
        - 'ligand': HETATM entries that are not amino acids and when config.include_ligands = True

        Returns a list of (node_id, Residue, centroid, kind).
        """
        out: list[tuple[str, Residue, np.ndarray, str]] = []

        for ch in chains:
            for res in ch.get_residues():
                if _is_protein_residue(res) or _is_water(res):
                    continue

                hetflag, resseq, icode = res.id
                resname = res.get_resname().strip().upper()

                aa_like = is_aa(res, standard=False)

                kind: str | None = None
                if aa_like and self.config.include_noncanonical_residues:
                    kind = "noncanonical_residue"
                elif (not aa_like) and self.config.include_ligands:
                    kind = "ligand"

                if kind is None:
                    continue

                coords = _heavy_atom_coords(res)
                cent = _centroid(coords)
                nid = _node_id(ch.id, res, kind="residue")
                out.append((nid, res, cent, kind))

        return out

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
        chains = self._select_chains()

        res_tuples = self._collect_residues(chains)
        _, inconsist = res_tuples_to_df(res_tuples)
        for msg in inconsist:
            log.warning(msg)

        extra_tuples = self._collect_nonprotein_residues(chains)
        node_kind = {}

        for nid, res, cent in res_tuples:
            node_kind[nid] = "residue"
        for nid, res, cent, kind in extra_tuples:
            res_tuples.append((nid, res, cent))
            node_kind[nid] = kind

        water_tuples = self._collect_waters(chains)
        raw_pdb_df, pdb_df, rgroup_df = self._make_atom_tables(chains)
        pdb_df = self._centroid_pdb_df_from_raw(raw_pdb_df)
        ca_cb_map = self._extract_ca_cb(raw_pdb_df)
        centroid_tbl = self._centroids_from_raw_df(raw_pdb_df)
        nan3 = (float("nan"), float("nan"), float("nan"))
        if not res_tuples:
            raise ValueError("No protein residues found for the selected chains.")

        asa_rsa: dict[str, tuple[float, float | None]] = {}
        dssp_df: pd.DataFrame | None = None
        if self.config.compute_rsa:
            if res_tuples:
                asa_rsa, dssp_df = self._compute_asa_rsa(res_tuples)

        res_ids = [t[0] for t in res_tuples]
        res_objects = [t[1] for t in res_tuples]

        cent_rows = centroid_tbl.reindex(res_ids)
        if cent_rows[["x_coord", "y_coord", "z_coord"]].isna().any(axis=None):
            missing = cent_rows[cent_rows.isna().any(axis=1)].index.tolist()
            raise ValueError(f"Missing centroid coords for node_ids: {missing[:10]}")

        res_centroids = cent_rows[["x_coord", "y_coord", "z_coord"]].to_numpy(float)

        diff = res_centroids[:, None, :] - res_centroids[None, :, :]
        dist_mat = np.sqrt(np.sum(diff * diff, axis=2))

        G = nx.Graph()
        for nid, res, cent in res_tuples:
            asa, rsa = asa_rsa.get(nid, (None, None))

            if node_kind[nid] != "residue":
                print(f"[NONCAN] {nid}  RSA = {rsa}  ASA = {asa}")

            kind = node_kind.get(nid, "residue")
            extra = ca_cb_map.get(nid, {})
            ca_coord = extra.get("ca_coord", nan3)
            cb_coord = extra.get("cb_coord", nan3)
            cb_is_virtual = bool(extra.get("cb_is_virtual", False))
            c = centroid_tbl.loc[nid]
            cent = np.array([c["x_coord"], c["y_coord"], c["z_coord"]], dtype=float)

            G.add_node(
                nid,
                kind=kind,
                chain=res.get_parent().id,
                resname=res.get_resname(),
                resseq=int(res.id[1]),
                icode=(res.id[2].strip() or ""),
                centroid=tuple(float(x) for x in cent),
                coords=np.array(cent, dtype=float),
                ca_coord=ca_coord,
                cb_coord=cb_coord,
                cb_is_virtual=cb_is_virtual,
                asa=None if asa is None else float(asa),
                rsa=None if rsa is None else float(rsa),
            )

        residue_map = {nid: res for nid, res, _ in res_tuples}

        secondary_structure(
            G,
            dssp_config=DSSPConfig(executable="mkdssp"),
            structure=self.structure,
            residue_map=residue_map,
            pdb_path=str(self.pdb_path)
        )

        interacting_nodes = np.where(dist_mat <= self.config.residue_distance_cutoff)
        interacting_nodes = list(zip(interacting_nodes[0], interacting_nodes[1]))

        for a1, a2 in interacting_nodes:
            n1 = res_ids[a1]
            n2 = res_ids[a2]

            # row1 = pdb_df.set_index("node_id").loc[n1]
            # row2 = pdb_df.set_index("node_id").loc[n2]

            # n1_chain = row1["chain_id"]
            # n2_chain = row2["chain_id"]
            # n1_position = row1["residue_number"]
            # n2_position = row2["residue_number"]

            # condition_1 = n1_chain == n2_chain
            # condition_2 = (
            #     abs(n1_position - n2_position) <= 1
            # )

            d = float(dist_mat[a1, a2])

            # if not (condition_1 and condition_2):
            G.add_edge(n1, n2, distance=d, kind="res-res")

        water_ids: list[str] = []
        water_centroids: np.ndarray | None = None
        if water_tuples:
            water_ids = [t[0] for t in water_tuples]
            water_centroids = np.vstack([t[2] for t in water_tuples])
            for nid, res, cent in water_tuples:
                _, resseq, icode = res.id
                chain_id = res.get_parent().id  # type: ignore[attr-defined]
                G.add_node(
                    nid,
                    kind="water",
                    chain=chain_id,
                    resname="HOH",
                    resseq=int(resseq),
                    icode=(icode.strip() or ""),
                    centroid=tuple(float(x) for x in cent),
                    ca_coord=nan3,
                    cb_coord=nan3,
                    cb_is_virtual=False,
                    rsa=1.0,
                    asa=None
                )
            if water_centroids is not None:
                wd = water_centroids[:, None, :] - res_centroids[None, :, :]
                wdist = np.sqrt(np.sum(wd * wd, axis=2))
                wcut = self.config.water_distance_cutoff
                W, R = wdist.shape
                for wi in range(W):
                    for rj in range(R):
                        d = float(wdist[wi, rj])
                        if 0.0 < d <= wcut:
                            G.add_edge(water_ids[wi], res_ids[rj], distance=d, kind="wat-res")

        (
            contact_map,
            contact_node_order,
            residue_map_dict,
            residue_map_dict_all,
        ) = contact_map_from_graph(
            G,
            granularity=getattr(self.config, "granularity", "all_atoms"),
            exclude_kinds=(),              # include waters if included in graph
            fallback_to_centroid=True,
        )

        G.graph["contact_map"] = contact_map
        G.graph["contact_map_node_order"] = contact_node_order
        G.graph["residue_map_dict"] = residue_map_dict
        G.graph["residue_map_dict_all"] = residue_map_dict_all

        rsa_series = pd.Series(
            {nid: (float(d.get("rsa")) if d.get("rsa") is not None else np.nan)
             for nid, d in G.nodes(data=True)},
            name="rsa",
        )
        if dssp_df is None:
            dssp_df_all = rsa_series.to_frame()
        else:
            dssp_df_all = dssp_df.copy()
            missing = [nid for nid in G.nodes if nid not in dssp_df_all.index]
            if missing:
                dssp_df_all = pd.concat([dssp_df_all, pd.DataFrame(index=missing)], axis=0)
            if "rsa" not in dssp_df_all.columns:
                dssp_df_all["rsa"] = np.nan
            rsa_aligned = rsa_series.reindex(dssp_df_all.index)
            dssp_df_all["rsa"] = dssp_df_all["rsa"].where(dssp_df_all["rsa"].notna(), rsa_aligned)
        dssp_df = dssp_df_all

        built = BuiltGraph(
            graph=G,
            residue_index=list(zip(res_ids, res_objects)),
            residue_centroids=res_centroids,
            water_index=[(t[0], t[1]) for t in water_tuples],
            water_centroids=water_centroids,
            distance_matrix=dist_mat if self.config.store_distance_matrix else None,
            raw_pdb_df=raw_pdb_df,
            pdb_df=pdb_df,
            rgroup_df=rgroup_df,
            dssp_df=dssp_df
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


if __name__ == "__main__":
    import argparse
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")

    ap = argparse.ArgumentParser(description="Build a pMHC structural graph from PDB/mmCIF.")
    ap.add_argument("pdb_path", help="Path to a PDB or mmCIF file.")
    ap.add_argument("--chains", type=str, default=None, help="Comma-separated chain IDs, e.g., A,C.")
    ap.add_argument("--include_waters", action="store_true", help="Include water molecules as nodes.")
    ap.add_argument("--res_cut", type=float, default=10.0, help="Residue–residue cutoff distance in Å.")
    ap.add_argument("--wat_cut", type=float, default=6.0, help="Water–residue cutoff distance in Å.")
    ap.add_argument("--no_rsa", action="store_true", help="Disable ASA/RSA computation.")
    ap.add_argument("--graphml", type=str, default=None, help="Output GraphML path.")
    ap.add_argument("--json", type=str, default=None, help="Output JSON path.")

    args = ap.parse_args()

    cfg = GraphConfig(
        chains=None if args.chains is None else [c.strip() for c in args.chains.split(",") if c.strip()],
        include_waters=args.include_waters,
        residue_distance_cutoff=args.res_cut,
        water_distance_cutoff=args.wat_cut,
        compute_rsa=not args.no_rsa,
    )
    builder = PDBGraphBuilder(args.pdb_path, cfg)
    built = builder.build_graph()

    if args.graphml:
        PDBGraphBuilder.to_graphml(built.graph, args.graphml)
    if args.json:
        PDBGraphBuilder.to_json(built.graph, args.json)

    print(f"Graph built with {built.graph.number_of_nodes()} nodes and {built.graph.number_of_edges()} edges.")
    print("Distance matrix path:", built.graph.graph.get("distmat_path"))
    print("Residue labels path:", built.graph.graph.get("residue_labels_path"))
