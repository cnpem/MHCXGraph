"""
Generate clustered heatmaps for MHCXGraph component-count results.

SciPy's experimental Array API backend support is explicitly disabled before
SciPy or Seaborn is imported. All arrays passed to SciPy are converted to
NumPy ndarrays, so this module does not use PyTorch.
"""

import os

# IMPORTANT: this must run before importing scipy or seaborn.
# SciPy only enables alternative Array API backends when this variable is set.
os.environ.pop("SCIPY_ARRAY_API", None)

import argparse
import glob
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

plt.switch_backend("agg")


def extract_unique_aminoacids(file_path):
    """Extract unique residues per component and consolidate totals per protein."""
    try:
        with open(file_path, "r", encoding="utf-8") as file:
            data = json.load(file)

        if not isinstance(data, dict):
            return None

        results = {"components": {}}
        global_unique_0 = set()
        global_unique_1 = set()
        individual_rows = []

        sorted_comps = sorted(
            [key for key in data if key.isdigit() and key != "0"],
            key=int,
        )

        for comp_id in sorted_comps:
            comp_data = data[comp_id]

            try:
                comp_value = comp_data.get("comp")
                if comp_value is None or comp_value <= 0:
                    continue
            except (KeyError, TypeError):
                continue

            comp_unique_0 = set()
            comp_unique_1 = set()
            frames = comp_data.get("frames", {})

            for frame_id, frame_data in frames.items():
                if str(frame_id) == "0":
                    continue

                nodes = frame_data.get("nodes", [])

                for pair in nodes:
                    if isinstance(pair, list) and len(pair) == 2:
                        res_0, res_1 = pair

                        comp_unique_0.add(res_0)
                        comp_unique_1.add(res_1)
                        global_unique_0.add(res_0)
                        global_unique_1.add(res_1)

            results["components"][str(comp_id)] = {
                "prot_0_count": len(comp_unique_0),
                "prot_1_count": len(comp_unique_1),
            }

            for aminoacid in sorted(comp_unique_0):
                individual_rows.append(
                    {
                        "Component": comp_id,
                        "Protein": "0",
                        "Aminoacid": aminoacid,
                    }
                )

            for aminoacid in sorted(comp_unique_1):
                individual_rows.append(
                    {
                        "Component": comp_id,
                        "Protein": "1",
                        "Aminoacid": aminoacid,
                    }
                )

        results["global_prot_0_count"] = len(global_unique_0)
        results["global_prot_1_count"] = len(global_unique_1)
        results["individual_csv_data"] = individual_rows

        return results

    except (OSError, json.JSONDecodeError, TypeError, ValueError) as error:
        print(f"Error while processing {file_path}: {error}")
        return None


def extract_original_graph_info(file_path):
    """Read original graph sizes and protein names from a graph JSON file."""
    try:
        with open(file_path, "r", encoding="utf-8") as file:
            data = json.load(file)

        original_graphs = data.get("original_graphs", {})
        node_counts = {}
        protein_names = {}

        for graph_id, graph_data in original_graphs.items():
            nodes = graph_data.get("nodes", [])
            node_counts[int(graph_id)] = len(nodes)
            protein_names[f"Protein_Name_{graph_id}"] = graph_data.get(
                "name",
                f"Graph_{graph_id}",
            )

        return node_counts, protein_names

    except (OSError, json.JSONDecodeError, TypeError, ValueError) as error:
        print(f"Error while reading graph metadata from {file_path}: {error}")
        return {}, {}


def process_directories(directory_path, output_path):
    """Process graph JSON files and create the component-count CSV files."""
    directory_path = Path(directory_path)
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    recursive_pattern = str(directory_path / "**" / "graph_*.json")
    root_pattern = str(directory_path / "graph_*.json")

    json_files = sorted(
        set(
            glob.glob(recursive_pattern, recursive=True)
            + glob.glob(root_pattern)
        )
    )

    summary_data = []
    all_comp_cols = set()

    for path_string in json_files:
        path = Path(path_string)
        file_name = path.name

        data_extracted = extract_unique_aminoacids(path)
        orig_counts, prot_names = extract_original_graph_info(path)

        if data_extracted is None:
            continue

        row = {"File": file_name}
        row.update(prot_names)

        sum_orig = sum(orig_counts.values())

        for graph_id, count in orig_counts.items():
            row[f"Original_Graph_{graph_id}"] = count

        unique_0 = data_extracted["global_prot_0_count"]
        unique_1 = data_extracted["global_prot_1_count"]

        row["Unique_Prot_0"] = unique_0
        row["Unique_Prot_1"] = unique_1

        for comp_id, counts in data_extracted["components"].items():
            prot_0_column = f"Comp_{comp_id}_Prot_0"
            prot_1_column = f"Comp_{comp_id}_Prot_1"

            row[prot_0_column] = counts["prot_0_count"]
            row[prot_1_column] = counts["prot_1_count"]

            all_comp_cols.update([prot_0_column, prot_1_column])

        row["total_prot_comp"] = unique_0 + unique_1
        row["ratio_total_prot_comp"] = (
            round((unique_0 + unique_1) / sum_orig, 4)
            if sum_orig > 0
            else 0
        )

        summary_data.append(row)

        if data_extracted["individual_csv_data"]:
            individual_df = pd.DataFrame(
                data_extracted["individual_csv_data"]
            )

            individual_csv_path = (
                output_path
                / f"unique_nodes_{path.stem}.csv"
            )

            individual_df.to_csv(individual_csv_path, index=False)
            print(f"Individual file generated: {individual_csv_path}")

    if not summary_data:
        return None

    summary_df = pd.DataFrame(summary_data)

    protein_name_columns = sorted(
        column
        for column in summary_df.columns
        if "Protein_Name_" in column
    )

    original_graph_columns = sorted(
        column
        for column in summary_df.columns
        if "Original_Graph_" in column
    )

    unique_columns = ["Unique_Prot_0", "Unique_Prot_1"]

    component_columns = sorted(
        all_comp_cols,
        key=lambda column: (
            int(column.split("_")[1]),
            column.split("_")[3],
        ),
    )

    final_order = (
        ["File"]
        + protein_name_columns
        + original_graph_columns
        + unique_columns
        + component_columns
        + ["total_prot_comp", "ratio_total_prot_comp"]
    )

    summary_df = summary_df.reindex(columns=final_order).fillna(0)

    integer_columns = (
        original_graph_columns
        + unique_columns
        + component_columns
        + ["total_prot_comp"]
    )

    for column in integer_columns:
        summary_df[column] = summary_df[column].astype(int)

    global_matrix_path = output_path / "component_count_matrix.csv"
    summary_df.to_csv(global_matrix_path, index=False)

    print(f"\nGlobal matrix saved: {global_matrix_path}")
    return summary_df


def create_distance_matrix(csv_path):
    """Create a symmetric dissimilarity matrix from pairwise similarity ratios."""
    dataframe = pd.read_csv(csv_path)

    proteins = sorted(
        set(dataframe["Protein_Name_0"])
        | set(dataframe["Protein_Name_1"])
    )

    matrix = pd.DataFrame(
        np.zeros((len(proteins), len(proteins)), dtype=np.float64),
        index=proteins,
        columns=proteins,
    )

    for _, row in dataframe.iterrows():
        protein_1 = row["Protein_Name_0"]
        protein_2 = row["Protein_Name_1"]
        similarity = float(row["ratio_total_prot_comp"])
        distance = 1.0 - similarity

        matrix.loc[protein_1, protein_2] = distance
        matrix.loc[protein_2, protein_1] = distance

    return matrix


def build_component_matrix(comp_df):
    """Build a symmetric matrix containing total component-node counts."""
    proteins = sorted(
        set(comp_df["Protein_Name_0"])
        | set(comp_df["Protein_Name_1"])
    )

    number_of_proteins = len(proteins)
    component_matrix = np.zeros(
        (number_of_proteins, number_of_proteins),
        dtype=np.float64,
    )

    protein_to_index = {
        protein: index
        for index, protein in enumerate(proteins)
    }

    protein_to_original = {}

    for _, row in comp_df.iterrows():
        protein_0 = row["Protein_Name_0"]
        protein_1 = row["Protein_Name_1"]
        value = float(row["total_prot_comp"])

        index_0 = protein_to_index[protein_0]
        index_1 = protein_to_index[protein_1]

        component_matrix[index_0, index_1] = value
        component_matrix[index_1, index_0] = value

        protein_to_original[protein_0] = row["Original_Graph_0"]
        protein_to_original[protein_1] = row["Original_Graph_1"]

    for protein in proteins:
        index = protein_to_index[protein]
        component_matrix[index, index] = (
            float(protein_to_original.get(protein, 0)) * 2
        )

    return pd.DataFrame(
        component_matrix,
        index=proteins,
        columns=proteins,
    )


def build_ratio_matrix(comp_df):
    """Build a symmetric similarity matrix with a diagonal of 1.0."""
    proteins = sorted(
        set(comp_df["Protein_Name_0"])
        | set(comp_df["Protein_Name_1"])
    )

    number_of_proteins = len(proteins)
    ratio_matrix = np.zeros(
        (number_of_proteins, number_of_proteins),
        dtype=np.float64,
    )

    protein_to_index = {
        protein: index
        for index, protein in enumerate(proteins)
    }

    for _, row in comp_df.iterrows():
        protein_0 = row["Protein_Name_0"]
        protein_1 = row["Protein_Name_1"]
        value = float(row["ratio_total_prot_comp"])

        index_0 = protein_to_index[protein_0]
        index_1 = protein_to_index[protein_1]

        ratio_matrix[index_0, index_1] = value
        ratio_matrix[index_1, index_0] = value

    np.fill_diagonal(ratio_matrix, 1.0)

    return pd.DataFrame(
        ratio_matrix,
        index=proteins,
        columns=proteins,
    )


def calculate_linkage(distance_dataframe):
    """
    Calculate average-linkage hierarchical clustering with SciPy.

    The explicit NumPy conversion prevents alternative array backends from
    being selected even if a caller supplied another array-like object.
    """
    distance_matrix = np.asarray(
        distance_dataframe.to_numpy(),
        dtype=np.float64,
        order="C",
    ).copy()

    if distance_matrix.ndim != 2:
        raise ValueError("The distance matrix must be two-dimensional.")

    rows, columns = distance_matrix.shape

    if rows != columns:
        raise ValueError("The distance matrix must be square.")

    if rows < 2:
        raise ValueError(
            "At least two proteins are required for hierarchical clustering."
        )

    if not np.isfinite(distance_matrix).all():
        raise ValueError(
            "The distance matrix contains NaN or infinite values."
        )

    if not np.allclose(
        distance_matrix,
        distance_matrix.T,
        rtol=1e-10,
        atol=1e-12,
    ):
        raise ValueError("The distance matrix must be symmetric.")

    np.fill_diagonal(distance_matrix, 0.0)

    condensed_distance = squareform(
        distance_matrix,
        checks=False,
    )

    return linkage(
        condensed_distance,
        method="average",
        optimal_ordering=False,
    )


def create_heatmap(args):
    """Generate the clustered heatmap and its supporting CSV matrices."""
    input_directory = Path(args.input_dir)
    output_directory = Path(args.output_dir)
    output_directory.mkdir(parents=True, exist_ok=True)

    print(f"Processing directories inside {input_directory}")

    component_dataframe = process_directories(
        input_directory,
        output_directory,
    )

    if component_dataframe is None:
        raise RuntimeError(
            f"No valid graph_*.json files were found inside "
            f"{input_directory}."
        )

    component_matrix_path = (
        output_directory / "component_count_matrix.csv"
    )

    print(
        "Creating distance matrix from: "
        f"{component_matrix_path}"
    )

    distance_dataframe = create_distance_matrix(
        component_matrix_path
    )

    distance_dataframe.to_csv(
        output_directory / "distance_matrix.csv"
    )

    labels = distance_dataframe.index.tolist()
    linkage_matrix = calculate_linkage(distance_dataframe)

    comp_df = pd.read_csv(component_matrix_path)
    comp_df_full = build_component_matrix(comp_df)
    ratio_df_full = build_ratio_matrix(comp_df)

    # Seaborn receives a precomputed SciPy linkage matrix, so it does not
    # recalculate pairwise distances. It still uses SciPy to interpret and
    # render the dendrogram layout.
    cluster_grid = sns.clustermap(
        ratio_df_full,
        cmap="viridis",
        row_linkage=linkage_matrix,
        col_linkage=linkage_matrix,
        figsize=(10, 8),
        dendrogram_ratio=0.15,
        cbar_kws={"label": "Similarity index"},
    )

    row_order = cluster_grid.dendrogram_row.reordered_ind
    ordered_labels = [labels[index] for index in row_order]

    component_ordered = comp_df_full.reindex(
        index=ordered_labels,
        columns=ordered_labels,
    )

    cluster_grid.ax_heatmap.set_xticklabels(
        ordered_labels,
        fontsize=11,
        rotation=90,
    )

    cluster_grid.ax_heatmap.set_yticklabels(
        ordered_labels,
        fontsize=11,
        rotation=0,
    )

    heatmap_axis = cluster_grid.ax_heatmap

    for row_index in range(len(ordered_labels)):
        for column_index in range(len(ordered_labels)):
            value = component_ordered.iloc[
                row_index,
                column_index,
            ]

            heatmap_axis.text(
                column_index + 0.5,
                row_index + 0.5,
                f"{int(value)}",
                horizontalalignment="center",
                verticalalignment="center",
                fontsize=11,
                color="black",
            )

    output_file = output_directory / args.name

    cluster_grid.fig.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(cluster_grid.fig)
    print(f"Heatmap saved: {output_file}")

    return output_file


def build_argument_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Create a hierarchically clustered protein-similarity heatmap."
        )
    )

    parser.add_argument(
        "-i",
        "--input-dir",
        required=True,
        help="Directory containing graph_*.json files.",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        help="Directory where CSV files and the heatmap will be written.",
    )

    parser.add_argument(
        "-n",
        "--name",
        default="heatmap.png",
        help="Output image filename. Default: heatmap.png",
    )

    return parser


def main():
    parser = build_argument_parser()
    arguments = parser.parse_args()
    create_heatmap(arguments)


if __name__ == "__main__":
    main()
