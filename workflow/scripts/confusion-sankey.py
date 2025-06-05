import numpy as np
import pandas as pd
import plotly.graph_objects as go
from pathlib import Path

def format_hectares(count):
    hectares = count / 1000  # since each cell is 10 m²
    if hectares >= 10_000:
        kha = hectares / 1000
        return f"{kha:.1f} kha"
    elif hectares >= 1:
        return f"{hectares:,.0f} ha"
    else:
        return f"{hectares:.2f} ha"

format_func = format_hectares

left_colors = snakemake.params['colourmap1']
right_colors = snakemake.params['colourmap2']
left_labels = snakemake.params['labels1']
right_labels = snakemake.params['labels2']
source1 = snakemake.wildcards['source1']
source2 = snakemake.wildcards['source2']
region = snakemake.wildcards['region']

def plot_sankey_from_confusion(csv_path):
    # Load full matrix including totals
    df = pd.read_csv(csv_path, index_col=0)

    # Separate out the "core" matrix (excluding last row/col)
    core = df.iloc[:-1, :-1]

    # Convert index and columns to integers for sorting
    core.index = core.index.astype(int)
    core.columns = core.columns.astype(int)

    # Sort left and right indices based on their numeric IDs
    left_order = sorted(core.index)
    right_order = sorted(core.columns)

    # Reorder the core matrix
    core = core.loc[left_order, right_order]

    # Nodes
    left_nodes = [f"[{left}] {left_labels[left]}" for left in left_order]
    right_nodes = [f"[{right}] {right_labels[right]}" for right in right_order]

    # Prepare colors
    node_colors = [left_colors.get(left, "#a6cee3") for left in left_order] + \
                  [right_colors.get(right, "#1f78b4") for right in right_order]

    # Setting up the data for the Sankey diagram
    source = []
    target = []
    value = []
    link_colors = []  # To store the link colors based on their magnitude

    # Threshold and scaling setup
    threshold = core.values.max() * 0.01  # Set threshold as 1% of the maximum flow
    min_opacity = 0.2
    max_opacity = 0.8

    for i, left in enumerate(left_order):
        for j, right in enumerate(right_order):
            flow_value = core.loc[left, right]
            if flow_value > 0:
                source.append(i)
                target.append(len(left_order) + j)
                value.append(flow_value)
                
                # Calculate opacity based on flow value
                if flow_value > threshold:
                    opacity = min_opacity + (max_opacity - min_opacity) * (flow_value / core.values.max()) ** 0.5
                else:
                    opacity = min_opacity
                link_colors.append(f"rgba(31, 120, 180, {opacity})")  # Adjust color based on calculated opacity

    # Data for the Sankey diagram
    sankey_data = go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=left_nodes + right_nodes,
            color=node_colors
        ),
        link=dict(
            source=source,
            target=target,
            value=value,
            color=link_colors
        )
    )

    # Create figure
    fig = go.Figure(sankey_data)
    fig.update_layout(title_text=f"Sankey Diagram: {source1} → {source2} ({region})", font_size=10)

    # Save to file
    fig.write_image(snakemake.output['png'], scale=2)
    fig.write_html(snakemake.output['html'])

# Make sure the output directory exists
Path(snakemake.output['png']).parent.mkdir(parents=True, exist_ok=True)

# Generate the diagram
plot_sankey_from_confusion(snakemake.input['csv'])
