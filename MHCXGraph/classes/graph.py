from __future__ import annotations

import base64
import copy
import json
import re
from pathlib import Path
from typing import Any, TypedDict

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from Bio.PDB import Atom, Chain, Model, Residue, Structure
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Superimposer import Superimposer

from MHCXGraph.core.config import GraphConfig
from MHCXGraph.core.pipeline import build_graph_with_config
from MHCXGraph.core.subgraphs import extract_subgraph
from MHCXGraph.utils.logging_utils import get_log
from MHCXGraph.utils.pyvis_inject import inject_fullscreen_css, inject_std_hover
from MHCXGraph.utils.tools import association_product

log = get_log()


class GraphData(TypedDict):
    id: int
    name: str
    graph: nx.Graph
    sorted_nodes: list[str]
    contact_map: np.ndarray
    residue_map: dict
    residue_map_all: dict
    rsa: np.ndarray
    pdb_file: str


class Graph:
    """
    Graph representation of a protein structure.

    This class loads a structural file (PDB or mmCIF), constructs a
    residue interaction graph using the configured pipeline, and
    provides utilities for extracting and exporting filtered subgraphs.

    Attributes
    ----------
    graph_path : str
        Path to the structure file used to construct the graph.

    config : GraphConfig
        Configuration object defining graph construction parameters.

    graph : networkx.Graph
        Graph representation of the protein structure.

    subgraphs : dict[str, networkx.Graph | list[str]]
        Dictionary storing named subgraphs or node selections.

    pdb_df : pandas.DataFrame or None
        DataFrame containing PDB atom-level information.

    raw_pdb_df : pandas.DataFrame or None
        DataFrame containing raw PDB data before filtering.

    dssp_df : pandas.DataFrame or None
        DataFrame containing DSSP secondary structure annotations.
    """

    def __init__(self, graph_path: str, config: GraphConfig):
        """
        Initialize a protein structure graph.

        Parameters
        ----------
        graph_path : str
            Path to a PDB or mmCIF file.

        config : GraphConfig
            Graph configuration object specifying node granularity,
            edge construction rules, and inclusion criteria for waters,
            ligands, and noncanonical residues.
        """

        self.graph_path = graph_path
        self.config = config

        self.graph: nx.Graph = build_graph_with_config(pdb_path=graph_path, config=self.config)

        self.subgraphs: dict[str, nx.Graph | list[str]] = {}
        self.pdb_df: pd.DataFrame | None = self.graph.graph.get("pdb_df")
        self.raw_pdb_df: pd.DataFrame | None = self.graph.graph.get("raw_pdb_df")
        self.dssp_df: pd.DataFrame | None = self.graph.graph.get("dssp_df")

    def get_subgraph(self, name: str):
        """
        Retrieve a stored subgraph by name.

        Parameters
        ----------
        name : str
            Name of the stored subgraph.

        Returns
        -------
        subgraph : networkx.Graph or list[str] or None
            Requested subgraph or node list if it exists. Returns
            ``None`` if the subgraph is not found.
        """

        if name not in self.subgraphs:
            log.warning(f"Can't find {name} in subgraph")
            return None
        else:
            return self.subgraphs[name]

    def create_subgraph(self, name: str, node_list: list | None = None, return_node_list: bool = False, **args):
        """
        Create and store a subgraph.

        A subgraph may be generated either from an explicit list of nodes
        or by applying filtering arguments passed to the subgraph
        extraction pipeline.

        Parameters
        ----------
        name : str
            Name used to store the subgraph.

        node_list : list, optional
            Explicit list of node identifiers used to build the subgraph.

        return_node_list : bool, default=False
            If True, return the node list of the generated subgraph.

        **args : dict
            Additional keyword arguments forwarded to the subgraph
            extraction function.

        Returns
        -------
        nodes : list or networkx.classes.reportviews.NodeView, optional
            Node identifiers of the generated subgraph if
            ``return_node_list`` is True.
        """

        if node_list is None:
            node_list = []

        if name in self.subgraphs:
            log.warning(f"You already have this subgraph created. Use graph.delete_subraph({name}) before creating it again.")
        elif not node_list:
            self.subgraphs[name] = extract_subgraph(g=self.graph, **args)
            log.info(f"Subgraph {name} created with success!")
        elif node_list:
            self.subgraphs[name] = self.graph.subgraph(node_list)

        if return_node_list:
            subgraphs_name = self.subgraphs[name]
            if isinstance(subgraphs_name, list):
                return subgraphs_name
            elif isinstance(subgraphs_name, nx.Graph):
                return subgraphs_name.nodes

    def _nx_to_serializable(self, g: nx.Graph) -> dict[str, Any]:
        """
        Convert a NetworkX graph into a JSON-serializable dictionary.

        Parameters
        ----------
        g : networkx.Graph
            Graph to be serialized.

        Returns
        -------
        data : dict[str, Any]
            Dictionary containing node identifiers, edge lists, and
            adjacency relationships encoded as strings.
        """

        nodes = [str(n) for n in g.nodes()]
        edges = [[str(u), str(v)] for u, v in g.edges()]
        neighbors = {
            str(n): [str(nb) for nb in g.neighbors(n)]
            for n in g.nodes()
        }

        return {
            "graph_path": self.graph_path,
            "nodes": nodes,
            "edges": edges,
            "neighbors": neighbors,
        }

    def save_filtered_pdb(self,
        g: nx.Graph,
        output_path: str | Path,
        name: str,
        use_cif: bool = False):
        """
        Save a filtered structural file containing only selected residues.

        Parameters
        ----------
        g : networkx.Graph
            Graph containing the subset of residues to export.

        output_path : str or pathlib.Path
            Directory where the filtered structure file will be written.

        name : str
            Base filename used when saving the structure.

        use_cif : bool, default=False
            If True, save the output in mmCIF format. Otherwise, save
            the structure as a PDB file.

        Returns
        -------
        None
        """

        path = str(self.graph_path)
        if path.lower().endswith((".cif", ".mmcif")):
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)

        orig = parser.get_structure("orig", self.graph_path)

        node_labels = set(g.nodes())
        new_struct = Structure.Structure("filtered")
        model = Model.Model(0)
        new_struct.add(model)

        for chain in orig[0]:
            new_chain = Chain.Chain(chain.id)
            for res in chain:
                chain_id = chain.id
                resname = res.get_resname()
                seqid = res.id[1]
                icode = res.id[2] if res.id[2] != " " else ""
                label = f"{chain_id}:{resname}:{seqid}{icode}"
                if label in node_labels:
                    new_res = Residue.Residue(res.id, resname, res.segid)
                    for atom in res:
                        new_atom = Atom.Atom(
                            atom.get_name(),
                            atom.get_coord().copy(),
                            atom.get_bfactor(),
                            atom.get_occupancy(),
                            atom.get_altloc(),
                            atom.get_fullname(),
                            atom.get_serial_number(),
                            element=atom.element
                        )
                        new_res.add(new_atom)
                    new_chain.add(new_res)
            if len(new_chain):
                model.add(new_chain)

        out_dir = Path(output_path)
        out_dir.mkdir(parents=True, exist_ok=True)
        out_file = out_dir / f"{name}.cif" if use_cif else out_dir / f"{name}.pdb"

        if use_cif:
            io = MMCIFIO()
        else:
            io = PDBIO()

        io.set_structure(new_struct)
        io.save(str(out_file))
        log.info(f"Filtered structure saved to {out_file}")


class AssociatedGraph:
    """
    Association engine for multiple protein graphs.

    This class performs structural association between several
    residue graphs, computes alignment frames, and generates
    visualization outputs and structural frames.

    Attributes
    ----------
    graphs : list[tuple]
        List containing graph tuples generated by preprocessing.

    output_path : pathlib.Path
        Directory where association outputs are written.

    run_name : str
        Identifier for the current association run.

    association_config : dict[str, Any]
        Configuration dictionary controlling association behavior.

    graphs_data : list[GraphData]
        Processed metadata extracted from the input graphs.

    associated_graphs : list or None
        Resulting association graph structures.
    """

    def __init__(self,
                graphs: list[tuple],
                output_path: str,
                run_name: str,
                association_config: dict[str, Any]):
        """
        Initialize an association task between protein graphs.

        Parameters
        ----------
        graphs : list[tuple]
            List of tuples produced during preprocessing containing
            graph objects and metadata.

        output_path : str
            Directory where association outputs will be saved.

        run_name : str
            Unique identifier for this run.

        association_config : dict[str, Any], optional
            Configuration dictionary controlling association behavior.
        """

        self.association_config = association_config

        self.track_residues = self.association_config.get("track_residues")
        self.graphs = graphs
        self.output_path = Path(output_path)
        self.run_name = run_name

        self.graphs_data = self._prepare_graph_data()

        result = association_product(graphs_data=self.graphs_data,
                                    config=self.association_config)

        if result is not None:
            self.__dict__.update(result)
            self.associated_graphs = result["AssociatedGraph"]
        else:
            self.associated_graphs = None

    def _parse_resnum_and_icode(self, resnum_str: str) -> tuple[int, str]:
        """
        Parse a residue number string into residue number and insertion code.

        Parameters
        ----------
        resnum_str : str
            Residue number string such as ``"6"`` or ``"6A"``.

        Returns
        -------
        resnum : int
            Residue sequence number.

        icode : str
            Insertion code or a blank space if none is present.
        """

        resnum_str = str(resnum_str).strip()
        if resnum_str.isdigit() or (resnum_str.startswith('-') and resnum_str[1:].isdigit()):
            return int(resnum_str), ' '
        match = re.match(r"\s*(-?\d+)\s*([A-Za-z]?)\s*", resnum_str)
        if not match:
            raise ValueError(f"Invalid residue resnum: {resnum_str!r}")
        resnum = int(match.group(1))
        icode = match.group(2) or ' '
        return resnum, icode

    def _find_residue_by_label(self, model, label: str):
        """
        Locate a residue in a structural model using a node label.

        Parameters
        ----------
        model : Bio.PDB.Model.Model
            Structural model containing residues.

        label : str
            Node label formatted as ``"CHAIN:RESNAME:RESNUM"``.

        Returns
        -------
        residue : Bio.PDB.Residue.Residue
            Residue matching the specified label.

        Raises
        ------
        KeyError
            If the residue cannot be found.
        """

        try:
            chain_id, resname, resnum_str = label.split(':')
        except ValueError:
            raise ValueError(f"Invalid node label format (expected 'CHAIN:RESNAME:RESNUM'): {label!r}")

        chain = model[chain_id]
        resnum, icode = self._parse_resnum_and_icode(resnum_str)

        candidates = []
        for res in chain.get_residues():
            hetflag, seq, ic = res.id
            if seq != resnum:
                continue
            if icode != ' ' and ic != icode:
                continue
            if res.get_resname().strip().upper() == resname.strip().upper():
                return res
            candidates.append(res)

        for res in candidates:
            return res

        for res in chain.get_residues():
            hetflag, seq, ic = res.id
            if seq == resnum:
                return res

        raise KeyError(f"Residue not found for label {label!r} (chain={chain_id}, resnum={resnum}, icode={icode!r})")

    def _prepare_graph_data(self) -> list[GraphData]:
        """
        Construct metadata required for graph association.

        The method extracts graph attributes such as contact maps,
        residue mappings, and solvent accessibility values.

        Returns
        -------
        graphs_data : list[GraphData]
            Structured metadata describing each input graph.
        """

        graphs_data = []
        for i, (g, pdb_file, name) in enumerate(self.graphs):
            contact_map = g.graph["contact_map"]
            residue_map = g.graph["residue_map_dict"]
            residue_map_all = g.graph["residue_map_dict_all"]

            sorted_nodes = sorted(list(g.nodes()))

            data: GraphData = {
                "id": i,
                "graph": g,
                "name": name,
                "sorted_nodes": sorted_nodes,
                "contact_map": contact_map,
                "residue_map": residue_map,
                "residue_map_all": residue_map_all,
                "rsa": g.graph["dssp_df"]["rsa"],
                "pdb_file": pdb_file
            }
            graphs_data.append(data)

        return graphs_data

    def create_pdb_per_protein(self):
        """
        Generate frame structures for each associated graph component.

        The method reconstructs residue subsets for each frame and writes
        them as separate structural models.

        Returns
        -------
        None
        """

        if not isinstance(self.associated_graphs, list):
            return

        for i in range(len(self.graphs)):
            pdb_file = self.graphs[i][1]
            parser = PDBParser(QUIET=True)
            orig_struct = parser.get_structure('orig', pdb_file)

            new_struct = Structure.Structure('frames')
            model = Model.Model(0)
            new_struct.add(model)

            for comp_id, comps in enumerate(self.associated_graphs):
                chain_counter = 1
                for frame_id, frame_graph in enumerate(comps[0]):
                    nodes = set(node[i] for node in frame_graph.nodes)

                    chain_id = f"K{comp_id}{chain_counter:03d}"
                    chain_counter += 1
                    chain = Chain.Chain(chain_id)

                    for node_label in nodes:
                        try:
                            orig_res = self._find_residue_by_label(orig_struct[0], node_label)
                        except KeyError:
                            log.warning(f"Residue {node_label} not found in {pdb_file}, skipping.")
                            continue

                        chain_name, resname, resnum_str = node_label.split(':')
                        resnum, icode = self._parse_resnum_and_icode(resnum_str)

                        hetflag, _, _ = orig_res.id
                        res_id = (hetflag, f"{resnum}:{chain_name}", icode)

                        new_res = Residue.Residue(res_id, orig_res.resname, orig_res.segid)

                        for atom in orig_res:
                            new_atom = Atom.Atom(
                                atom.get_name(),
                                atom.get_coord().copy(),
                                atom.get_bfactor(),
                                atom.get_occupancy(),
                                atom.get_altloc(),
                                atom.get_fullname(),
                                atom.get_serial_number(),
                                element=atom.element
                            )
                            new_res.add(new_atom)
                        chain.add(new_res)
                    model.add(chain)

            out_dir = self.output_path / "frames"
            out_dir.mkdir(exist_ok=True)
            use_cif = True
            if use_cif:
                io = MMCIFIO()
                out_file = out_dir / f"{Path(pdb_file).stem}_frames.cif"
            else:
                io = PDBIO()
                out_file = out_dir / f"{Path(pdb_file).stem}_frames.pdb"

            io.set_structure(new_struct)
            io.save(str(out_file))

            log.info(f"Structure saved in {out_file}")

    def _write_frame_multichain(self, comp_idx: int, frame_idx: int,
                                models: list, output_dir: Path):
        """
        Write a multi-chain mmCIF file containing aligned protein models.

        Parameters
        ----------
        comp_idx : int
            Component index of the association result.

        frame_idx : int
            Frame index within the component.

        models : list
            List of structural models to be written.

        output_dir : pathlib.Path
            Directory where the mmCIF file will be saved.

        Returns
        -------
        None
        """

        combo = Structure.Structure(f"comp{comp_idx}_frame{frame_idx}")
        combo_model = Model.Model(0)
        combo.add(combo_model)

        for prot_idx, prot_model in enumerate(models):
            for ch in prot_model:
                new_chain = copy.deepcopy(ch)

                new_chain.id = f"{ch.id}{prot_idx}"
                combo_model.add(new_chain)

        output_dir.mkdir(exist_ok=True)
        out_path = output_dir / f"comp{comp_idx}_frame{frame_idx}_all.cif"
        io = MMCIFIO()
        io.set_structure(combo)
        io.save(str(out_path))
        log.debug(f"[comp{comp_idx}_frame{frame_idx}] wrote {len(models)} proteins as "
            f"{len(combo_model)} chains to {out_path}")

    def align_all_frames(self):
        """
        Align protein structures across all associated graph frames.

        The alignment is performed using C-alpha atoms of residues
        shared across proteins. Frames lacking sufficient common
        residues are skipped.

        Returns
        -------
        None
        """

        parser = PDBParser(QUIET=True)

        if self.associated_graphs is not None:
            for comp_idx, (frame_graphs, _) in enumerate(self.associated_graphs):
                for frame_idx, assoc_graph in enumerate(frame_graphs):
                    nodes = list(assoc_graph.nodes())
                    if not nodes:
                        log.info(f"[comp{comp_idx}_frame{frame_idx}] graph vazio, pulando.")
                        continue

                    models = []
                    for prot_idx, (_, pdb_path, _) in enumerate(self.graphs):
                        struct = parser.get_structure(f"p{prot_idx}", pdb_path)
                        models.append(struct[0])

                    anchor_labels = []
                    for label in nodes:
                        try:
                            ref_res = self._find_residue_by_label(models[0], label[0])
                        except KeyError:
                            continue

                        try:
                            ref_res['CA']
                        except KeyError:
                            continue

                        ok = True
                        for prot_idx in range(1, len(models)):
                            try:
                                mob_res = self._find_residue_by_label(models[prot_idx], label[prot_idx])
                            except KeyError:
                                ok = False
                                break
                            try:
                                mob_res['CA']
                            except KeyError:
                                ok = False
                                break

                        if ok:
                            anchor_labels.append(label)

                    if len(anchor_labels) < 3:
                        log.info(
                            f"[comp{comp_idx}_frame{frame_idx}] "
                            f"menos de 3 resíduos com CA comum ({len(anchor_labels)} encontrados), "
                            "não dá para alinhar este frame, pulando."
                        )
                        continue

                    ref_cas = []
                    for label in anchor_labels:
                        ref_res = self._find_residue_by_label(models[0], label[0])
                        ref_cas.append(ref_res['CA'])

                    for prot_idx in range(1, len(models)):
                        mob_cas = []
                        for label in anchor_labels:
                            mob_res = self._find_residue_by_label(models[prot_idx], label[prot_idx])
                            mob_cas.append(mob_res['CA'])

                        sup = Superimposer()
                        sup.set_atoms(ref_cas, mob_cas)
                        sup.apply(models[prot_idx].get_atoms())

                        assoc_graph.graph[f"rmsd_p0_p{prot_idx}"] = round(sup.rms, 3)

                        log.debug(
                            f"[comp{comp_idx}_frame{frame_idx}] "
                            f"prot{prot_idx} <- prot0  RMSD={sup.rms:.2f}"
                        )

                    output_dir = self.output_path / "frames"
                    self._write_frame_multichain(comp_idx, frame_idx, models, output_dir)

    def draw_graph_interactive(self, show=False, save=True):
        """
        Generate interactive visualizations of associated graphs.

        Graphs are rendered using PyVis and can be displayed in a
        browser or saved as HTML files.

        Parameters
        ----------
        show : bool, default=False
            If True, display the generated visualization.

        save : bool, default=True
            If True, save the visualization to disk.

        Returns
        -------
        None
        """

        if not show and not save:
            log.info("You are not saving or viewing the graph. Please leave at least one of the parameters as true.")
            return

        if not isinstance(self.associated_graphs, list):
            log.warning(f"I can't draw the graph because it's {self.associated_graphs!r}")
            return

        out_dir = Path(self.output_path)
        out_dir.mkdir(parents=True, exist_ok=True)

        export_data = {
            "run_name": self.run_name,
            "metadata": self.association_config,
            "proteins": [gd['name'] for gd in self.graphs_data],
            "protein_paths": [gd['pdb_file'] for gd in self.graphs_data],
            "nodes": [],
            "edges": [],
            "components": [],
            "filtered_graphs": []
        }

        global_nodes = {}
        global_edges = {}
        node_id_counter = 0
        edge_id_counter = 0
        all_chain_ids = set()

        # Pass 1: Gather Associated Graphs Nodes, edges, and std_matrices
        for j, comps in enumerate(self.associated_graphs):
            comp_data = {
                "id": j, 
                "node_ids": [], 
                "frames": [],
                "std_matrix": None,
                "node_index_map": {}
            }
            
            if len(comps[0]) > 0:
                comp_std_matrix = comps[0][0].graph.get("edge_std_matrix")
                if comp_std_matrix is not None:
                    safe_matrix = [
                        [None if np.isnan(val) else float(val) for val in row]
                        for row in comp_std_matrix
                    ]
                    comp_data["std_matrix"] = safe_matrix

            for i, graph in enumerate(comps[0]):
                if graph.number_of_nodes() == 0:
                    continue
                
                frame_data = {"id": i, "node_ids": []}
                node_index_map = graph.graph.get("node_index_map")

                # Nodes
                for n in graph.nodes():
                    node_label = str(n)
                    try:
                        residues = [r for r in n]
                    except TypeError:
                        residues = [str(n)]
                        
                    mapping = []
                    chains = []
                    for prot_idx, res in enumerate(residues):
                        parts = str(res).split(':')
                        if len(parts) >= 3: 
                            chains.append(parts[0])
                            mapping.append({
                                "model_idx": prot_idx,
                                "chain": parts[0],
                                "resn": parts[1],
                                "resi": parts[2]
                            })
                    
                    chain_id = ''.join(chains) or "?"
                    all_chain_ids.add(chain_id)

                    if node_label not in global_nodes:
                        global_nodes[node_label] = {
                            "id": node_id_counter,
                            "label": node_label,
                            "title": f"Chains: {chain_id}\n{node_label}",
                            "group": chain_id,
                            "chain_id": chain_id,
                            "mapping": mapping,
                            "originalColor": None, 
                        }
                        node_id_counter += 1
                    
                    numeric_id = global_nodes[node_label]["id"]
                    if numeric_id not in comp_data["node_ids"]:
                        comp_data["node_ids"].append(numeric_id)
                    frame_data["node_ids"].append(numeric_id)
                    
                    if node_index_map is not None and n in node_index_map:
                        comp_data["node_index_map"][numeric_id] = int(node_index_map[n])

                # Edges
                for u, v, data in graph.edges(data=True):
                    u_id = global_nodes[str(u)]["id"]
                    v_id = global_nodes[str(v)]["id"]
                    edge_key = tuple(sorted([u_id, v_id]))
                    
                    if edge_key not in global_edges:
                        dist_tuple = []
                        for prot_idx, gd in enumerate(self.graphs_data):
                            u_res = u[prot_idx] if isinstance(u, tuple) else str(u)
                            v_res = v[prot_idx] if isinstance(v, tuple) else str(v)
                            try:
                                u_parts = u_res.split(":")
                                v_parts = v_res.split(":")
                                u_tup = (u_parts[0], u_parts[2], u_parts[1])
                                v_tup = (v_parts[0], v_parts[2], v_parts[1])
                                distance = gd["contact_map"][gd["residue_map_all"][u_tup], gd["residue_map_all"][v_tup]]
                                dist_tuple.append(f"{distance:.2f}")
                            except:
                                dist_tuple.append("N/A")
                        
                        dists = "(" + ", ".join(dist_tuple) + ")"
                        w = data.get('weight')
                        
                        std_val = None
                        if node_index_map is not None and comp_data["std_matrix"] is not None:
                            idx_u = node_index_map.get(u)
                            idx_v = node_index_map.get(v)
                            if idx_u is not None and idx_v is not None:
                                std_val = comp_data["std_matrix"][idx_u][idx_v]
                        
                        title = f"Distances: {dists}\n"
                        if w is not None: title += f"Weight: {w}\n"
                        if std_val is not None: title += f"STD: {std_val:.3f}"
                        
                        global_edges[edge_key] = {
                            "id": edge_id_counter,
                            "from": u_id,
                            "to": v_id,
                            "std": std_val,
                            "title": title,
                            "raw_dist": float(w) if isinstance(w, (int, float)) else None
                        }
                        edge_id_counter += 1

                comp_data["frames"].append(frame_data)
            export_data["components"].append(comp_data)

        export_data["nodes"] = list(global_nodes.values())
        export_data["edges"] = list(global_edges.values())

        filtered_graphs_data = []
        for prot_idx, (g, pdb_file, name) in enumerate(self.graphs):
            f_nodes = []
            f_edges = []
            
            for n in g.nodes():
                node_label = str(n)
                parts = node_label.split(':')
                chain_id = parts[0] if len(parts) >= 1 else "?"
                mapping = []
                if len(parts) >= 3:
                    mapping.append({
                        "model_idx": prot_idx,
                        "chain": parts[0],
                        "resn": parts[1],
                        "resi": parts[2]
                    })
                
                f_nodes.append({
                    "id": node_label,
                    "label": node_label,
                    "title": f"Chain: {chain_id}\n{node_label}",
                    "group": chain_id,
                    "chain_id": chain_id,
                    "mapping": mapping
                })
                
            for u, v, data in g.edges(data=True):
                w = data.get('weight')
                dist_val = data.get('distance', w)
                
                title = f"Distance: {float(dist_val):.2f} Å" if dist_val is not None else ""
                
                f_edges.append({
                    "id": f"{u}-{v}",
                    "from": str(u),
                    "to": str(v),
                    "title": title,
                    "raw_dist": float(dist_val) if isinstance(dist_val, (int, float)) else None
                })
                
            filtered_graphs_data.append({
                "id": prot_idx,
                "name": name,
                "nodes": f_nodes,
                "edges": f_edges
            })
            
        export_data["filtered_graphs"] = filtered_graphs_data

        # Safely load local Frameworks or fallback to CDNs
        assets_dir = Path(__file__).resolve().parent.parent / "assets"
        
        vis_local = assets_dir / "vis-network.min.js"
        if vis_local.exists():
            vis_injection = f'<script>\n{vis_local.read_text(encoding="utf-8")}\n</script>'
        else:
            log.warning("Local vis-network.min.js not found in assets/. Falling back to CDN.")
            vis_injection = '<script type="text/javascript" src="https://unpkg.com/vis-network/standalone/umd/vis-network.min.js"></script>'

        mol3d_local = assets_dir / "3Dmol-min.js"
        if mol3d_local.exists():
            mol3d_injection = f'<script>\n{mol3d_local.read_text(encoding="utf-8")}\n</script>'
        else:
            log.warning("Local 3Dmol-min.js not found in assets/. Falling back to CDN.")
            mol3d_injection = '<script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>'

        logo_dark_path = assets_dir / "LNBio white.png"
        logo_light_path = assets_dir / "LNBio.png"
        logo_injection = ""
        if logo_light_path.exists():
            with open(logo_light_path, "rb") as image_file:
                encoded = base64.b64encode(image_file.read()).decode("utf-8")
                logo_injection += f'<img src="data:image/png;base64,{encoded}" alt="LNBio Logo" class="logo-light" style="height: 7rem; width: auto;">'
        if logo_dark_path.exists():
            with open(logo_dark_path, "rb") as image_file:
                encoded = base64.b64encode(image_file.read()).decode("utf-8")
                logo_injection += f'<img src="data:image/png;base64,{encoded}" alt="LNBio Logo" class="logo-dark" style="height: 7rem; width: auto;">'
        if not logo_injection:
            log.debug("LNBio logos not found in assets/. Skipping logo injection.")

        mhcx_logo_path = assets_dir / "MHCXGraph logo.png"
        mhcx_logo_injection = "<h2>MHCXGraph</h2>" # Fallback text if image is missing
        favicon_injection = ""

        if mhcx_logo_path.exists():
            with open(mhcx_logo_path, "rb") as image_file:
                encoded = base64.b64encode(image_file.read()).decode("utf-8")
                mhcx_logo_injection = f'<img src="data:image/png;base64,{encoded}" alt="MHCXGraph Logo" style="max-width: 80%; height: auto;">'
                favicon_injection = f'<link rel="icon" type="image/png" href="data:image/png;base64,{encoded}">'
        else:
            log.debug("MHCXGraph logo not found in assets/. Using text fallback.")

        # Load Template and Inject Data
        template_path = assets_dir / "dashboard_template.html"
        try:
            with open(template_path, "r", encoding="utf-8") as f:
                html_template = f.read()
        except FileNotFoundError:
            log.error(f"Template not found at {template_path}. Please create it.")
            return

        json_data = json.dumps(export_data)
        final_html = html_template.replace("__GRAPH_DATA_INJECTION__", json_data)
        final_html = final_html.replace("__FAVICON_INJECTION__", favicon_injection)
        final_html = final_html.replace("__VIS_JS_INJECTION__", vis_injection)
        final_html = final_html.replace("__3DMOL_JS_INJECTION__", mol3d_injection)
        final_html = final_html.replace("__MHCXGRAPH_LOGO_INJECTION__", mhcx_logo_injection)
        final_html = final_html.replace("__LNBIO_LOGO_INJECTION__", logo_injection)

        if save:
            full_path = out_dir / "Dashboard.html"
            with open(str(full_path), "w+", encoding="utf-8") as out:
                out.write(final_html)
            log.info(f"Interactive Dashboard saved to {full_path}")
            
        if show:
            tmpfile = out_dir / "__preview.html"
            with open(str(tmpfile), "w+", encoding="utf-8") as out:
                out.write(final_html)

    def get_dashboard_data(self, global_proteins: list[str]) -> dict:
        """Extracts JSON serializable data for the dashboard injection."""
        global_idx = [global_proteins.index(gd['name']) for gd in self.graphs_data]
        
        export_data = {
            "proteins": [gd['name'] for gd in self.graphs_data],
            "protein_paths": [gd['pdb_file'] for gd in self.graphs_data],
            "nodes": [],
            "edges": [],
            "components": [],
            "filtered_graphs": []
        }

        global_nodes, global_edges = {}, {}
        node_id_counter, edge_id_counter = 0, 0

        for j, comps in enumerate(self.associated_graphs or []):
            comp_data = {"id": j, "node_ids": [], "frames": [], "std_matrix": None, "node_index_map": {}}
            if len(comps[0]) > 0 and comps[0][0].graph.get("edge_std_matrix") is not None:
                comp_data["std_matrix"] = [[None if np.isnan(val) else float(val) for val in row] for row in comps[0][0].graph.get("edge_std_matrix")]

            for i, frame_graph in enumerate(comps[0]):
                if frame_graph.number_of_nodes() == 0: continue
                
                # ---> EXTRACT THE SAVED RMSDS <---
                rmsds = {k: v for k, v in frame_graph.graph.items() if k.startswith("rmsd")}
                frame_data = {"id": i, "node_ids": [], "rmsds": rmsds}
                
                node_index_map = frame_graph.graph.get("node_index_map")

                for n in frame_graph.nodes():
                    node_label = str(n)
                    residues = [r for r in n] if isinstance(n, tuple) else [str(n)]
                    mapping, chains = [], []
                    
                    for prot_idx, res in enumerate(residues):
                        parts = str(res).split(':')
                        if len(parts) >= 3:
                            chains.append(parts[0])
                            
                            # ---> EXTRACT RSA FROM THE ORIGINAL GRAPH <---
                            orig_g = self.graphs[prot_idx][0]
                            rsa_val = orig_g.nodes.get(str(res), {}).get("rsa", "N/A")
                            if isinstance(rsa_val, float) and not np.isnan(rsa_val):
                                rsa_val = round(rsa_val, 3)
                            else:
                                rsa_val = "N/A"
                                
                            mapping.append({
                                "model_idx": global_idx[prot_idx], 
                                "chain": parts[0], "resn": parts[1], "resi": parts[2], 
                                "rsa": rsa_val # Pass RSA to frontend
                            })

                    chain_id = ''.join(chains) or "?"
                    if node_label not in global_nodes:
                        # ---> ADD RSA TO HOVER TITLE <---
                        title = f"Chains: {chain_id}\n{node_label}"
                        for m in mapping:
                            title += f"\nP{m['model_idx']} RSA: {m['rsa']}"
                            
                        global_nodes[node_label] = {
                            "id": node_id_counter, "label": node_label, "title": title,
                            "group": chain_id, "mapping": mapping, "originalColor": None,
                        }
                        node_id_counter += 1                   
                    numeric_id = global_nodes[node_label]["id"]
                    if numeric_id not in comp_data["node_ids"]: comp_data["node_ids"].append(numeric_id)
                    frame_data["node_ids"].append(numeric_id)
                    if node_index_map and n in node_index_map: comp_data["node_index_map"][numeric_id] = int(node_index_map[n])

                for u, v, data in frame_graph.edges(data=True):
                    u_id, v_id = global_nodes[str(u)]["id"], global_nodes[str(v)]["id"]
                    edge_key = tuple(sorted([u_id, v_id]))

                    if edge_key not in global_edges:
                        std_val = None
                        if node_index_map and comp_data["std_matrix"]:
                            idx_u, idx_v = node_index_map.get(u), node_index_map.get(v)
                            if idx_u is not None and idx_v is not None: std_val = comp_data["std_matrix"][idx_u][idx_v]

                        dist_tuple = []

                        for prot_idx, gd in enumerate(self.graphs_data):
                            try:
                                u_res = u[prot_idx] if isinstance(u, tuple) else str(u)
                                v_res = v[prot_idx] if isinstance(v, tuple) else str(v)

                                u_parts = u_res.split(":")
                                v_parts = v_res.split(":")

                                u_tup = (u_parts[0], u_parts[2], u_parts[1])
                                v_tup = (v_parts[0], v_parts[2], v_parts[1])

                                distance = gd["contact_map"][
                                    gd["residue_map_all"][u_tup],
                                    gd["residue_map_all"][v_tup]
                                ]

                                dist_tuple.append(distance)

                            except Exception:
                                dist_tuple.append(None)

                        dist_str = "(" + ", ".join(
                            f"{d:.2f}" if d is not None else "N/A"
                            for d in dist_tuple
                        ) + ")"

                        title = f"Distances: {dist_str}\n"

                        global_edges[edge_key] = {
                            "id": edge_id_counter,
                            "from": u_id,
                            "to": v_id,
                            "std": std_val,
                            "title": title,
                            "raw_dist": dist_tuple[0] if dist_tuple and isinstance(dist_tuple[0], (int, float)) else None
                        }

                        edge_id_counter += 1
                comp_data["frames"].append(frame_data)
            export_data["components"].append(comp_data)

        export_data["nodes"] = list(global_nodes.values())
        export_data["edges"] = list(global_edges.values())

        for prot_idx, (g, pdb_file, name) in enumerate(self.graphs):
            f_nodes, f_edges = [], []
            for n in g.nodes():
                parts = str(n).split(':')
                mapping = [{"model_idx": global_idx[prot_idx], "chain": parts[0], "resn": parts[1], "resi": parts[2]}] if len(parts) >= 3 else []
                f_nodes.append({"id": str(n), "label": str(n), "title": f"Chain: {parts[0] if len(parts)>=1 else '?'}\n{n}", "group": parts[0] if len(parts)>=1 else "?", "mapping": mapping})
            for u, v, data in g.edges(data=True):
                dist_val = data.get('distance')
                f_edges.append({"id": f"{u}-{v}", "from": str(u), "to": str(v), "title": f"Distance: {float(dist_val):.2f} Å" if dist_val is not None else "", "raw_dist": float(dist_val) if isinstance(dist_val, (int, float)) else None})
            export_data["filtered_graphs"].append({"id": global_idx[prot_idx], "name": name, "nodes": f_nodes, "edges": f_edges})

        return export_data
