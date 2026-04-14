import argparse
import glob
import json
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path

plt.switch_backend('agg')


# ==============================
# Substituições de scipy
# ==============================

def squareform(X):
    """
    Converte entre forma condensada e quadrada.
    - 1D (condensada, n*(n-1)/2) → matriz quadrada (n x n)
    - 2D quadrada (n x n)        → vetor condensado 1D
    Equivalente a scipy.spatial.distance.squareform.
    """
    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        n = int((1 + np.sqrt(1 + 8 * X.shape[0])) / 2)
        D = np.zeros((n, n), dtype=float)
        idx = 0
        for i in range(n - 1):
            length = n - i - 1
            D[i, i + 1:] = X[idx: idx + length]
            D[i + 1:, i] = X[idx: idx + length]
            idx += length
        return D
    elif X.ndim == 2 and X.shape[0] == X.shape[1]:
        n = X.shape[0]
        idx = np.triu_indices(n, k=1)
        return X[idx]
    else:
        raise ValueError("X deve ser 1D (condensada) ou 2D quadrada.")


def linkage(condensed_dist, method="average"):
    """
    Clustering hierárquico average linkage (UPGMA).
    Aceita vetor condensado (saída de squareform sobre matriz quadrada).
    Equivalente a scipy.cluster.hierarchy.linkage(..., method='average').
    Retorna Z: array (n-1, 4) — [cluster_i, cluster_j, distância, tamanho].
    """
    if method != "average":
        raise NotImplementedError("Apenas method='average' está implementado.")

    condensed_dist = np.asarray(condensed_dist, dtype=float)
    n = int((1 + np.sqrt(1 + 8 * condensed_dist.shape[0])) / 2)

    # Monta matriz quadrada
    D = squareform(condensed_dist)
    np.fill_diagonal(D, np.inf)

    sizes = {i: 1 for i in range(n)}
    active = list(range(n))
    label_to_pos = {i: i for i in range(n)}
    next_id = n
    Z = np.zeros((n - 1, 4), dtype=float)

    for step in range(n - 1):
        pos_list = [label_to_pos[lbl] for lbl in active]
        sub = D[np.ix_(pos_list, pos_list)]
        flat_min = np.argmin(sub)
        row_local, col_local = divmod(int(flat_min), len(pos_list))

        lbl_i = active[row_local]
        lbl_j = active[col_local]
        pi, pj = label_to_pos[lbl_i], label_to_pos[lbl_j]

        dist_ij = D[pi, pj]
        sz_i, sz_j = sizes[lbl_i], sizes[lbl_j]
        sz_new = sz_i + sz_j
        new_lbl = next_id
        next_id += 1

        Z[step] = [lbl_i, lbl_j, dist_ij, sz_new]

        pos_new = pi
        for lbl_k in active:
            if lbl_k in (lbl_i, lbl_j):
                continue
            pk = label_to_pos[lbl_k]
            d_new_k = (sz_i * D[pi, pk] + sz_j * D[pj, pk]) / sz_new
            D[pos_new, pk] = d_new_k
            D[pk, pos_new] = d_new_k

        D[pos_new, pos_new] = np.inf
        active.remove(lbl_j)
        active.remove(lbl_i)
        active.append(new_lbl)
        label_to_pos[new_lbl] = pos_new
        sizes[new_lbl] = sz_new

    return Z


# ==============================
# Funções originais do script
# ==============================

def extract_unique_aminoacids(file_path):
    """Extrai resíduos únicos por componente e consolida o global por proteína."""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            data = json.load(f)

        if not isinstance(data, dict):
            return None

        results = {"components": {}}
        global_unique_0 = set()
        global_unique_1 = set()
        individual_rows = []

        sorted_comps = sorted(
            [k for k in data.keys() if k.isdigit() and k != "0"],
            key=int
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
                        res_0 = pair[0]
                        res_1 = pair[1]
                        comp_unique_0.add(res_0)
                        comp_unique_1.add(res_1)
                        global_unique_0.add(res_0)
                        global_unique_1.add(res_1)

            results["components"][str(comp_id)] = {
                "prot_0_count": len(comp_unique_0),
                "prot_1_count": len(comp_unique_1)
            }

            for aa in sorted(list(comp_unique_0)):
                individual_rows.append({"Component": comp_id, "Protein": "0", "Aminoacid": aa})
            for aa in sorted(list(comp_unique_1)):
                individual_rows.append({"Component": comp_id, "Protein": "1", "Aminoacid": aa})

        results["global_prot_0_count"] = len(global_unique_0)
        results["global_prot_1_count"] = len(global_unique_1)
        results["individual_csv_data"] = individual_rows

        return results

    except Exception as e:
        print(f"Erro ao processar {file_path}: {e}")
        return None


def extract_original_graph_info(file_path):
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        original_graphs = data.get("original_graphs", {})
        node_counts = {}
        protein_names = {}
        for graph_id, graph_data in original_graphs.items():
            nodes = graph_data.get("nodes", [])
            node_counts[int(graph_id)] = len(nodes)
            protein_names[f"Protein_Name_{graph_id}"] = graph_data.get("name", f"Graph_{graph_id}")
        return node_counts, protein_names
    except:
        return {}, {}


def process_directories(directory_path, output_path):
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    search_pattern = os.path.join(directory_path, "**/graph_*.json")
    root_pattern = os.path.join(directory_path, "graph_*.json")
    json_files = list(set(glob.glob(search_pattern, recursive=True) + glob.glob(root_pattern)))

    summary_data = []
    all_comp_cols = set()

    for path in json_files:
        file_name = os.path.basename(path)
        data_extracted = extract_unique_aminoacids(path)
        orig_counts, prot_names = extract_original_graph_info(path)

        if data_extracted is not None:
            row = {"File": file_name}
            row.update(prot_names)

            sum_orig = sum(orig_counts.values())
            for g_id, count in orig_counts.items():
                row[f"Original_Graph_{g_id}"] = count

            u0 = data_extracted["global_prot_0_count"]
            u1 = data_extracted["global_prot_1_count"]
            row["Unique_Prot_0"] = u0
            row["Unique_Prot_1"] = u1

            for comp_id, counts in data_extracted["components"].items():
                c0_name, c1_name = f"Comp_{comp_id}_Prot_0", f"Comp_{comp_id}_Prot_1"
                row[c0_name] = counts["prot_0_count"]
                row[c1_name] = counts["prot_1_count"]
                all_comp_cols.update([c0_name, c1_name])

            row["total_prot_comp"] = u0 + u1
            row["ratio_total_prot_comp"] = round((u0 + u1) / sum_orig, 4) if sum_orig > 0 else 0

            summary_data.append(row)

            if data_extracted["individual_csv_data"]:
                df_ind = pd.DataFrame(data_extracted["individual_csv_data"])
                csv_name = f"unique_nodes_{file_name.replace('.json', '.csv')}"
                individual_csv_path = os.path.join(output_path, csv_name)
                df_ind.to_csv(individual_csv_path, index=False)
                print(f"Individual file generated: {individual_csv_path}")

    if summary_data:
        df = pd.DataFrame(summary_data)
        prot_names = sorted([c for c in df.columns if "Protein_Name_" in c])
        orig_graphs = sorted([c for c in df.columns if "Original_Graph_" in c])
        uniques = ["Unique_Prot_0", "Unique_Prot_1"]
        comps = sorted(list(all_comp_cols), key=lambda x: (int(x.split('_')[1]), x.split('_')[3]))

        final_order = ["File"] + prot_names + orig_graphs + uniques + comps + ["total_prot_comp", "ratio_total_prot_comp"]
        df = df.reindex(columns=final_order).fillna(0)

        for col in orig_graphs + uniques + comps + ["total_prot_comp"]:
            df[col] = df[col].astype(int)

        global_matrix_path = os.path.join(output_path, "component_count_matrix.csv")
        df.to_csv(global_matrix_path, index=False)
        print(f"\nGlobal Matrix saved: {global_matrix_path}")
        return df

    return None


def create_distance_matrix(csv_path):
    df = pd.read_csv(csv_path)
    proteins = sorted(set(df['Protein_Name_0']).union(set(df['Protein_Name_1'])))
    n = len(proteins)
    matrix = pd.DataFrame(np.zeros((n, n)), index=proteins, columns=proteins)
    matrix_values = matrix.values.copy()
    np.fill_diagonal(matrix_values, 0.0)

    for _, row in df.iterrows():
        p1, p2, dist = row['Protein_Name_0'], row['Protein_Name_1'], row['ratio_total_prot_comp']
        matrix.loc[p1, p2] = 1 - dist
        matrix.loc[p2, p1] = 1 - dist

    return matrix


def build_component_matrix(comp_df):
    proteins = sorted(set(
        comp_df['Protein_Name_0'].tolist() +
        comp_df['Protein_Name_1'].tolist()
    ))

    n = len(proteins)
    comp_matrix = np.zeros((n, n))
    protein_to_idx = {protein: i for i, protein in enumerate(proteins)}
    protein_to_original = {}

    for _, row in comp_df.iterrows():
        p0 = row['Protein_Name_0']
        p1 = row['Protein_Name_1']
        value = row['total_prot_comp']
        i = protein_to_idx[p0]
        j = protein_to_idx[p1]
        comp_matrix[i, j] = value
        comp_matrix[j, i] = value
        protein_to_original[p0] = row['Original_Graph_0']
        protein_to_original[p1] = row['Original_Graph_1']

    for protein in proteins:
        idx = protein_to_idx[protein]
        comp_matrix[idx, idx] = protein_to_original.get(protein, 0) * 2

    return pd.DataFrame(comp_matrix, index=proteins, columns=proteins)


def build_ratio_matrix(comp_df):
    """Matriz de similaridade (ratio_total_prot_comp), diagonal = 1.0."""
    proteins = sorted(set(
        comp_df['Protein_Name_0'].tolist() +
        comp_df['Protein_Name_1'].tolist()
    ))
    n = len(proteins)
    matrix = np.zeros((n, n))
    protein_to_idx = {protein: i for i, protein in enumerate(proteins)}

    for _, row in comp_df.iterrows():
        p0, p1 = row['Protein_Name_0'], row['Protein_Name_1']
        val = row['ratio_total_prot_comp']
        i, j = protein_to_idx[p0], protein_to_idx[p1]
        matrix[i, j] = val
        matrix[j, i] = val

    np.fill_diagonal(matrix, 1.0)
    return pd.DataFrame(matrix, index=proteins, columns=proteins)


def create_heatmap(args):
    args.input_dir, args.output_dir = Path(args.input_dir), Path(args.output_dir)
    print(f"Processing directories inside {args.input_dir}")
    process_directories(args.input_dir, args.output_dir)
    print(f"Creating distance matrices given: {args.output_dir / 'distance_matrix.csv'}")
    matrix = create_distance_matrix(args.output_dir / "component_count_matrix.csv")
    matrix.to_csv(args.output_dir / "distance_matrix.csv")
    dist_df = matrix
    dist_matrix = dist_df.values
    labels = dist_df.index.tolist()

    condensed_dist = squareform(dist_matrix)
    linkage_matrix = linkage(condensed_dist, method="average")

    # ==============================
    # 2. MATRIZ Total_Comp_Sum + Ratio
    # ==============================
    comp_df = pd.read_csv(args.output_dir / "component_count_matrix.csv")
    comp_df_full = build_component_matrix(comp_df)
    ratio_df_full = build_ratio_matrix(comp_df)

    # ==============================
    # 3. CLUSTERMAP
    # ==============================
    g = sns.clustermap(
        ratio_df_full,
        cmap='viridis',
        row_linkage=linkage_matrix,
        col_linkage=linkage_matrix,
        figsize=(10, 8),
        dendrogram_ratio=0.15,
        cbar_kws={'label': 'Similarity index'}
    )

    row_order = g.dendrogram_row.reordered_ind
    ordered_labels = [labels[i] for i in row_order]

    comp_ordered = comp_df_full.reindex(
        index=ordered_labels,
        columns=ordered_labels
    )

    g.ax_heatmap.set_xticklabels(ordered_labels, fontsize=11, rotation=90)
    g.ax_heatmap.set_yticklabels(ordered_labels, fontsize=11, rotation=0)

    # ==============================
    # 6. Inserir valores nas células
    # ==============================
    ax = g.ax_heatmap

    for i in range(len(ordered_labels)):
        for j in range(len(ordered_labels)):
            value = comp_ordered.iloc[i, j]
            ax.text(
                j + 0.5,
                i + 0.5,
                f"{int(value)}",
                ha='center',
                va='center',
                fontsize=11,
                color='black'
            )

    plt.tight_layout()
    plt.savefig(args.output_dir / args.name, dpi=300, bbox_inches='tight')
