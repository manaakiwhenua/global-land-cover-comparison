import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from pathlib import Path
import sys

def compute_jaccard_index(intersection, union):
    """Compute Jaccard index."""
    return intersection / union if union != 0 else 0

def load_confusion_matrix(csv_path):
    """Load a confusion matrix from a CSV file, excluding 'Total' rows and columns."""
    df = pd.read_csv(csv_path, index_col=0)
    if 'Total' in df.columns:
        df = df.drop(columns='Total')
    if 'Total' in df.index:
        df = df.drop(index='Total')
    df.index = df.index.astype(int)
    df.columns = df.columns.astype(int)
    return df

def get_node_color(concept, source, sources_config):
    """Retrieve node color from configuration."""
    colourmap = sources_config[source]['colourmap']
    default_color = "#ffffff"  # Default/placeholder color
    return colourmap.get(list(colourmap.keys())[0], default_color)  # Use the first color in the map

def main():
    # Snakemake variables
    outputs = snakemake.output
    params = snakemake.params
    inputs = snakemake.input['confusion_matrices']
    region = snakemake.wildcards.region

    # Jaccard correspondence directly from parameters
    correspondence = params.jaccard_correspondence
    sources_config = params.sources_config

    G = nx.Graph()
    threshold = 0.1

    # Process each pair of sources for the given region
    for csv_path in inputs:
        print(csv_path)
        src1, src2 = Path(csv_path).stem.split('--')

        try:
            confusion_matrix = load_confusion_matrix(csv_path)
        except Exception as e:
            print(f"Error loading {csv_path}: {e}", file=sys.stderr)
            continue

        for class_info in correspondence:
            concept = list(class_info.keys())[0]
            datasets = class_info[concept]

            if src1 not in datasets or src2 not in datasets:
                continue

            values1 = datasets.get(src1, [])
            values2 = datasets.get(src2, [])

            values1 = [v for v in values1 if v in confusion_matrix.index]
            values2 = [v for v in values2 if v in confusion_matrix.columns]

            if not values1 or not values2:
                print(f"{concept} - Missing classes for {src1} or {src2} in confusion matrix.", file=sys.stderr)
                continue

            try:
                intersection_count = confusion_matrix.loc[values1, values2].sum().sum()
                union_count = (confusion_matrix.loc[values1].sum().sum() +
                               confusion_matrix.loc[:, values2].sum().sum() -
                               intersection_count)

                jaccard_index = compute_jaccard_index(intersection_count, union_count)

                if jaccard_index > threshold:
                    label1 = f"{concept}-{src1}"
                    label2 = f"{concept}-{src2}"
                    G.add_node(label1, color=get_node_color(concept, src1, sources_config))
                    G.add_node(label2, color=get_node_color(concept, src2, sources_config))
                    G.add_edge(label1, label2, weight=jaccard_index)
            except KeyError as ke:
                print(f"KeyError for {concept}: {ke}", file=sys.stderr)

    # Prepare node colors
    node_color_list = [G.nodes[node]['color'] for node in G.nodes]

    pos = nx.spring_layout(G)
    plt.figure(figsize=(12, 8))
    nx.draw_networkx_nodes(G, pos, node_size=700, node_color=node_color_list)
    edges = G.edges(data=True)
    edge_colors = [d['weight'] for _, _, d in edges]
    
    # Use a non-linear scaling for edge widths
    edge_widths = [(weight * 5)**0.5 for weight in edge_colors]

    edge_collection = nx.draw_networkx_edges(G, pos, edgelist=edges, edge_color=edge_colors,
                                             edge_cmap=plt.cm.Blues, edge_vmin=0, edge_vmax=1,
                                             width=edge_widths)
    nx.draw_networkx_labels(G, pos, font_size=10)

    # Correctly associate colorbar with the edge collection
    plt.colorbar(edge_collection, label='Jaccard Index')
    plt.title(f"Jaccard Index Network for {region}")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(outputs.png)
    plt.close()

if __name__ == "__main__":
    main()
