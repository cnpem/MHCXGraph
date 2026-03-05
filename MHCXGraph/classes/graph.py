from __future__ import annotations

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
from pyvis.network import Network

from MHCXGraph.core.config import GraphConfig
from MHCXGraph.core.pipeline import build_graph_with_config
from MHCXGraph.core.subgraphs import extract_subgraph
from MHCXGraph.utils.logging_utils import get_log
from MHCXGraph.utils.pyvis_inject import inject_fullscreen_css, inject_std_hover
from MHCXGraph.utils.tools import association_product

log = get_log()


def _rgba_to_hex(rgba):
    r, g, b, _ = rgba
    return f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}"


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
    """Represents a protein structure graph (no external framework assumptions)."""

    def __init__(self, graph_path: str, config: GraphConfig):
        """
        Parameters
        ----------
        graph_path : str
            Path to a PDB or mmCIF file.
        config : GraphConfig, optional
            Unified graph configuration. If not provided, a sensible default is used.
        """
        self.graph_path = graph_path
        self.config = config

        self.graph: nx.Graph = build_graph_with_config(pdb_path=graph_path, config=self.config)

        self.subgraphs: dict[str, nx.Graph | list[str]] = {}
        self.pdb_df: pd.DataFrame | None = self.graph.graph.get("pdb_df")
        self.raw_pdb_df: pd.DataFrame | None = self.graph.graph.get("raw_pdb_df")
        self.dssp_df: pd.DataFrame | None = self.graph.graph.get("dssp_df")

    def get_subgraph(self, name: str):
        if name not in self.subgraphs:
            log.warning(f"Can't find {name} in subgraph")
            return None
        else:
            return self.subgraphs[name]

    def create_subgraph(self, name: str, node_list: list | None = None, return_node_list: bool = False, **args):
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
        Converte um nx.Graph qualquer em um dicionário simples
        com nodes, edges e vizinhos, tudo como string.
        Serve tanto para o grafo bruto quanto para subgrafos filtrados.
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

    def save_subgraph_view(
        self,
        g: nx.Graph,
        output_dir: str | Path,
        name: str,
        with_html: bool = True,
    ) -> None:
        out_dir = Path(output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        base = name

        # JSON "cru" do grafo
        data = self._nx_to_serializable(g)
        json_path = out_dir / f"{base}.json"
        json_path.write_text(json.dumps(data, indent=4), encoding="utf-8")
        log.vinfo(f"{name} JSON subgraph saved in {json_path}.", f"{name} JSON subgraph saved.")

        if not with_html:
            return

        # trabalhar em uma cópia para o HTML
        graph = g.copy()

        # chain_id + attrs mínimos para PyVis
        for n in graph.nodes():
            label = str(n)
            chain_id = label.split(":")[0] if ":" in label else "?"
            graph.nodes[n]["chain_id"] = chain_id
            graph.nodes[n]["label"] = label
            graph.nodes[n]["title"] = f"chain: {chain_id}\n{label}"
            graph.nodes[n]["size"] = 12
            graph.nodes[n]["group"] = chain_id

        chain_ids = sorted({data.get("chain_id", "?") for _, data in graph.nodes(data=True)})
        cmap = plt.cm.get_cmap("tab10", max(1, len(chain_ids)))
        palette = {cid: _rgba_to_hex(cmap(idx)) for idx, cid in enumerate(chain_ids)}

        for n in graph.nodes():
            cid = graph.nodes[n].get("chain_id", "?")
            graph.nodes[n]["color"] = palette.get(cid, "#999999")

        # LIMPEZA: deixar só atributos serializáveis básicos
        allowed_node_keys = {"label", "title", "color", "size", "group", "chain_id"}
        for _, data_node in graph.nodes(data=True):
            for k in list(data_node.keys()):
                if k not in allowed_node_keys:
                    del data_node[k]

        for _, _, data_edge in graph.edges(data=True):
            w = data_edge.get("weight")
            for k in list(data_edge.keys()):
                if k != "weight":
                    del data_edge[k]
            if not isinstance(w, (int, float)) and "weight" in data_edge:
                del data_edge["weight"]

        # renomear nós para ids simples
        safe_map = {n: f"v{idx}" for idx, n in enumerate(graph.nodes())}
        H = nx.relabel_nodes(graph, safe_map, copy=True)

        net = Network(
            height="100%",
            width="100%",
            bgcolor="#ffffff",
            notebook=False,
            cdn_resources="in_line",
            directed=graph.is_directed(),
        )
        net.from_nx(H)

        for (u, v, data) in H.edges(data=True):
            w = data.get("weight")
            if w is not None:
                title = f"weight: {w}"
                for e in net.edges:
                    if (e["from"] == u and e["to"] == v) or (
                        not graph.is_directed() and e["from"] == v and e["to"] == u
                    ):
                        e["title"] = title
                        if isinstance(w, (int, float)):
                            e.setdefault("value", float(w))

        net.set_options("""
        {
          "nodes": { "shape": "dot" },
          "interaction": { "hover": true, "tooltipDelay": 150, "zoomView": true, "dragView": true },
          "physics": {
            "enabled": true,
            "solver": "barnesHut",
            "barnesHut": {
              "gravitationalConstant": -8000,
              "centralGravity": 0.3,
              "springLength": 95,
              "springConstant": 0.04,
              "damping": 0.09,
              "avoidOverlap": 0.1
            },
            "minVelocity": 0.75
          }
        }
        """)

        html_path = out_dir / f"{base}.html"
        html = net.generate_html(notebook=False, local=True)
        html_path.write_text(html, encoding="utf-8")
        log.vinfo(f"{name} HTML subgraph saved in {html_path}.", f"{name} subgraph HTML saved.")


class AssociatedGraph:
    """Handles the association of multiple protein graphs."""

    def __init__(self,
                graphs: list[tuple],
                output_path: str,
                run_name: str,
                association_config: dict[str, Any] | None = None):
        """
        Initialize an AssociatedGraph instance with a reduced configuration.
 
        :param graphs: List of tuples (Graph instance, raw_data) from preprocessing.
        :param output_path: Where to save results.
        :param run_name: Unique run identifier.
        :param association_config: Dictionary with keys.
        """
        if association_config is None:
            self.association_config = {}
        else:
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
        Converte uma string de resnum (tipo '6', '6A') em (6, ' ') ou (6, 'A').
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
        Dado um rótulo tipo 'C:SEP:6', tenta achar o Residue correspondente no modelo,
        independente de ser ATOM ou HETATM (SEP, MSE, ligantes, etc.).
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

        # Fallback: bate só por número (se não achou com nome)
        for res in candidates:
            return res

        # Último fallback: qualquer resíduo com mesmo número na cadeia
        for res in chain.get_residues():
            hetflag, seq, ic = res.id
            if seq == resnum:
                return res

        raise KeyError(f"Residue not found for label {label!r} (chain={chain_id}, resnum={resnum}, icode={icode!r})")

    def _prepare_graph_data(self) -> list[GraphData]:
        """
        For each (Graph, raw) tuple, build a dictionary with the necessary data:
            - "graph": The Graph instance.
            - "contact_map": Output of build_contact_map(raw)[0].
            - "residue_map_all": Output of build_contact_map(raw)[2].
            - "rsa": np.array(g.graph["dssp_df"]["rsa"]).
            - "depth": g.depth.
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
        if not isinstance(self.associated_graphs, list):
            return

        for i in range(len(self.graphs)):
            pdb_file = self.graphs[i][1]
            parser = PDBParser(QUIET=True)
            orig_struct = parser.get_structure('orig', pdb_file)

            chain_counter = 1
            new_struct = Structure.Structure('frames')
            model = Model.Model(0)
            new_struct.add(model)

            for comp_id, comps in enumerate(self.associated_graphs):
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
        Write one mmCIF in which *all* proteins are merged into Model 0
        and distinguished only by their (renamed) chains.

        The original chain IDs are suffixed with the protein index:
        chain A in protein 0  -> 'A0'
        chain C in protein 2  -> 'C2'
        This preserves one-letter labels when you view the file in PyMOL
        (PyMOL truncates after the first character) but keeps unique
        `_atom_site.label_asym_id` in mmCIF so nothing collides.

        mmCIF permits multi-character chain IDs, so MMCIFIO will write them
        without complaints.
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
        For each component (frame) in self.associated_graphs:
        - recreate fresh instances of the PDB models in self.graphs
        - use model 0 as the reference
        - align each mobile protein using a subset of nodes
            for which ALL involved residues have a CA atom
        - ignore water, ligands, and anything that lacks a CA atom
        - save a multichain mmCIF with all models already aligned

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

                        log.debug(
                            f"[comp{comp_idx}_frame{frame_idx}] "
                            f"prot{prot_idx} <- prot0  RMSD={sup.rms:.2f}"
                        )

                    output_dir = self.output_path / "frames"
                    self._write_frame_multichain(comp_idx, frame_idx, models, output_dir)

    def draw_graph_interactive(self, show=False, save=True):

        if not show and not save:
            log.info("You are not saving or viewing the graph. Please leave at least one of the parameters as true.")
            return

        if not isinstance(self.associated_graphs, list):
            log.warning(f"I can't draw the graph because it's {self.associated_graphs!r}")
            return

        out_dir = Path(self.output_path)
        out_dir.mkdir(parents=True, exist_ok=True)

        for j, comps in enumerate(self.associated_graphs):
            for i, graph in enumerate(comps[0]):
                if graph.number_of_nodes() == 0:
                    log.warning(f"{j}:{i} graph has no nodes. Skipping.")
                    continue

                # chain_id per node
                for node in graph.nodes():
                    try:
                        residues = [r for r in node]  # flatten tuples
                    except TypeError:
                        residues = [str(node)]
                    chains = [str(r).split(':')[0] for r in residues]
                    graph.nodes[node]['chain_id'] = ''.join(chains) or "?"

                chain_ids = sorted({data.get('chain_id', '?') for _, data in graph.nodes(data=True)})
                cmap = plt.cm.get_cmap('tab10', max(1, len(chain_ids)))
                palette = {cid: _rgba_to_hex(cmap(idx)) for idx, cid in enumerate(chain_ids)}

                node_labels = {}
                for n in graph.nodes():
                    if isinstance(n, tuple) and n and isinstance(n[0], tuple):
                        combo1, combo2 = n
                        lab1 = repr(combo1).replace(" ", "")
                        lab2 = repr(combo2).replace(" ", "")
                        node_labels[n] = f"({lab1})({lab2})"
                    else:
                        node_labels[n] = str(n)

                for n in graph.nodes():
                    cid = graph.nodes[n].get('chain_id', '?')
                    graph.nodes[n]['label'] = node_labels[n]
                    graph.nodes[n]['title'] = f"chain: {cid}\n{node_labels[n]}"
                    graph.nodes[n]['color'] = palette.get(cid, "#999999")
                    graph.nodes[n]['size'] = 12
                    graph.nodes[n]['group'] = cid

                safe_map = {n: f"v{idx}" for idx, n in enumerate(graph.nodes())}
                H = nx.relabel_nodes(graph, safe_map, copy=True)

                net = Network(
                    height="100%",
                    width="100%",
                    bgcolor="#ffffff",
                    notebook=False,
                    cdn_resources="in_line",
                    directed=graph.is_directed()
                )
                net.from_nx(H)

                std_matrix = graph.graph.get("edge_std_matrix")
                node_index_map = graph.graph.get("node_index_map")

                if std_matrix is not None and node_index_map is not None:
                    safe_node_index = {safe_map[n]: node_index_map[n] for n in graph.nodes() if n in node_index_map}
                else:
                    safe_node_index = None

                # optional edge hover with weight
                for (u, v, data) in H.edges(data=True):
                    w = data.get('weight')
                    if w is not None:
                        # PyVis stores edges in net.edges, update matching one
                        # create a small title for hover
                        title = f"weight: {w}"
                        # there can be multiple edges, so update all matches
                        for e in net.edges:
                            if (e['from'] == u and e['to'] == v) or (not graph.is_directed() and e['from'] == v and e['to'] == u):
                                e['title'] = title
                                e.setdefault('value', float(w) if isinstance(w, (int, float)) else 1)

                net.set_options("""
                {
                "nodes": { "shape": "dot" },
                "interaction": { "hover": true, "tooltipDelay": 150, "zoomView": true, "dragView": true },
                "physics": {
                    "enabled": true,
                    "solver": "barnesHut",
                    "barnesHut": {
                    "gravitationalConstant": -8000,
                    "centralGravity": 0.3,
                    "springLength": 95,
                    "springConstant": 0.04,
                    "damping": 0.09,
                    "avoidOverlap": 0.1
                    },
                    "minVelocity": 0.75
                }
                }
                """)

                if save:
                    if i == 0 and j == 0:
                        filename = "Associated Graph.html"
                    elif i == 0 and j != 0:
                        filename = f"Associated Graph Component {j}.html"
                    else:
                        filename = f"Associated Graph Component {j} Frame {i}.html"
                    full = out_dir / filename
                    html = net.generate_html(
                        notebook=False,
                        local=True
                    )

                    if std_matrix is not None and safe_node_index is not None:
                        html = inject_std_hover(html, std_matrix=std_matrix, safe_node_index=safe_node_index)
                        html = inject_fullscreen_css(html)
                    with open(str(full), "w+", encoding="utf-8") as out:
                        out.write(html)

                    log.info(f"{j}: saved graph {i} to {full}")

                elif show:
                    tmpfile = out_dir / f"__preview_{j}_{i}.html"
                    html = net.generate_html(
                        notebook=False,
                        local=True
                    )

                    if std_matrix is not None and safe_node_index is not None:
                        html = inject_std_hover(html, std_matrix=std_matrix, safe_node_index=safe_node_index)
                        html = inject_fullscreen_css(html)

                    with open(str(tmpfile), "w+", encoding="utf-8") as out:
                        out.write(html)

