import json
from pathlib import Path

from MHCXGraph.classes.graph import AssociatedGraph
from MHCXGraph.utils.analysis import (
    _make_json_from_associated_graph,
    _save_eval_tables,
    evaluate_all_frames_nodes_weighted,
)


def run_association_task(graphs, output_path, run_name, association_config, log):
    """
    Helper function to run the association logic for a specific set of graphs
    and a specific output directory.
    """
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    task_config = association_config.copy()
    task_config["output_path"] = str(output_path)

    log.info(
        f"Starting run '{run_name}' in mode: {task_config['run_mode']} with {len(graphs)} graphs."
    )
    log.info(f"Output directory: {output_path}")

    G = AssociatedGraph(
        graphs=graphs,
        output_path=str(output_path),
        run_name=run_name,
        association_config=task_config,
    )

    if G.associated_graphs is None:
        log.warning(f"No associated graphs found for {run_name}.")
        return

    log.debug(f"Drawing Graph for {run_name}")
    G.draw_graph_interactive(show=False, save=True)
    # G.draw_graph(show=False, save=True)

    G.create_pdb_per_protein()
    G.align_all_frames()

    # log.debug("Growing Subgraph")
    # try:
    #     G.grow_subgraph_bfs()
    # except Exception as e:
    #     log.error(f"Unable to grow subgraphs with BFS. Error: {e}")

    graph_data = dict()
    for j, comps in enumerate(G.associated_graphs):
        graph_data[j] = {"comp": j, "frames": {}}
        for i in range(len(comps[0])):
            nodes = list(comps[0][i].nodes)
            edges = list(comps[0][i].edges)
            neighbors = {
                str(node): [str(neighbor) for neighbor in comps[0][i].neighbors(node)]
                for node in nodes
            }
            graph_data[j]["frames"][i] = {
                "nodes": nodes,
                "edges": edges,
                "neighbors": neighbors,
            }

    output_json = output_path / f"graph_{run_name}.json"
    with output_json.open("w") as f:
        json.dump(graph_data, f, indent=4)

    log.info(f"Graph data saved to {output_json}")

    analysis_json_path = output_path / f"graph_{run_name}.json"
    _make_json_from_associated_graph(G, analysis_json_path)
    log.info(f"Analysis graph JSON saved to {analysis_json_path}")

    try:
        df_fp_nodes, df_frames_nodes_w = evaluate_all_frames_nodes_weighted(
            analysis_json_path
        )
        if (
            not df_fp_nodes.empty
            and "frame_nodes_unique" in df_fp_nodes.columns
            and "total_nodes_associated" not in df_fp_nodes.columns
        ):
            df_fp_nodes["total_nodes_associated"] = df_fp_nodes["frame_nodes_unique"]

        _save_eval_tables(Path(output_path), df_fp_nodes, df_frames_nodes_w)
        log.info(f"Evaluation tables saved in {output_path}")
    except Exception as e:
        log.error(f"Failed to compute evaluation tables: {e}")
