import numpy as np
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, PowerNorm
import matplotlib.patches as mpatches
import pandas as pd
import seaborn as sns

def format_compact(num):
    """Format number compactly (e.g. 1500 -> 1.5k, 2_500_000 -> 2.5m)"""
    if num >= 1_000_000:
        return f"{num/1_000_000:.1f}m"
    elif num >= 1_000:
        return f"{num/1_000:.1f}k"
    else:
        return str(num)

def format_km2(num):
    area_km2 = num * 10 / 1_000_000
    return f"{area_km2:.1f} km²" if area_km2 >= 0.1 else f"{area_km2:.2f} km²"

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

row_colors = snakemake.params['colourmap1']
col_colors = snakemake.params['colourmap2']
row_labels = snakemake.params['labels1']
col_labels = snakemake.params['labels2']
source1 = snakemake.wildcards['source1']
source2 = snakemake.wildcards['source2']
region = snakemake.wildcards['region']

def plot_confusion_with_totals(csv_path):
    # Load full matrix including totals
    df = pd.read_csv(csv_path, index_col=0)

    # Separate out the "core" matrix for color scaling (excluding last row/col)
    core = df.iloc[:-1, :-1]
    row_totals = df.iloc[:-1, -1]
    col_totals = df.iloc[-1, :-1]
    grand_total = df.iloc[-1, -1]

    # Format annotations for the whole matrix
    full_annot = df.map(lambda x: format_func(int(x)))

    # Plot only the core heatmap
    if max(len(core.columns), len(core.index)) > 10:
        width, height = (24, 12) if len(core.columns) > len(core.index) else (12, 24)
    else:
        width, height = 12, 10
    fig, ax = plt.subplots(figsize=(width, height))  # width, height

    # Create the heatmap
    sns.heatmap(core,
                cmap='Purples',
                alpha=0.5,
                annot=full_annot.iloc[:-1, :-1],
                fmt='',
                cbar=False,
                linewidths=0.5,
                linecolor='gray',
                norm=PowerNorm(gamma=0.3, vmin=1, vmax=core.values.max()),
                ax=ax)

    # Set ticks and labels
    ax.set_xticks(np.arange(len(core.columns)) + 0.5)
    ax.set_yticks(np.arange(len(core.index)) + 0.5)
    ax.set_xticklabels(core.columns, rotation=0, ha='center', va='bottom')
    ax.set_yticklabels(core.index, rotation=0)
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()

    # Color row labels background
    for tick_label in ax.get_yticklabels():
        label_text = tick_label.get_text()
        if int(label_text) in row_colors:
            tick_label.set_bbox(dict(facecolor=row_colors[int(label_text)], edgecolor='none', boxstyle='round,pad=0.3'))

    # Color column labels background
    for tick_label in ax.get_xticklabels():
        label_text = tick_label.get_text()
        if int(label_text) in col_colors:
            tick_label.set_bbox(dict(facecolor=col_colors[int(label_text)], edgecolor='none', boxstyle='round,pad=0.3'))

    # Manually annotate the total column (right side)
    for i, val in enumerate(row_totals):
        ax.text(len(core.columns) + 0.5, i + 0.5, format_func(int(val)),
                ha='center', va='center', fontsize=9, color='black')

    # Manually annotate the total row (bottom)
    for j, val in enumerate(col_totals):
        ax.text(j + 0.5, len(core.index) + 0.5, format_func(int(val)),
                ha='center', va='center', fontsize=9, color='black')

    # Annotate bottom-right total
    ax.text(len(core.columns) + 0.5, len(core.index) + 0.5,
            format_func(int(grand_total)), ha='center', va='center',
            fontsize=9, color='black')

    # Extend limits to accommodate total row/col
    ax.set_xlim(0, len(core.columns) + 1)
    ax.set_ylim(len(core.index) + 1, 0)  # y-axis is inverted

    # Add labels
    ax.set_xlabel(source2)
    ax.set_ylabel(source1)
    plt.title(f"Confusion Matrix: {source1} — {source2} ({region})", fontsize=12)

    # Draw lines on the edges of cells
    for i in range(core.shape[0]):
        for j in range(core.shape[1]):
            row_prop = core.iloc[i, j] / row_totals[i]
            col_prop = core.iloc[i, j] / col_totals[j]
            line_length_col = col_prop  # adjust for aesthetics
            line_length_row = row_prop  # adjust for aesthetics

            # Draw vertical line for column proportion
            ax.add_line(plt.Line2D([j + 1, j + 1], [i + 0.5 - line_length_col / 2, i + 0.5 + line_length_col / 2], color='black', linewidth=2))
            
            # Draw horizontal line for row proportion
            ax.add_line(plt.Line2D([j + 0.5 - line_length_row / 2, j + 0.5 + line_length_row / 2], [i + 1, i + 1], color='black', linewidth=2))

    plt.savefig(snakemake.output['png'], dpi=300)
    plt.savefig(snakemake.output['svg'])
    # plt.show()

Path(snakemake.output['png']).parent.mkdir(parents=True, exist_ok=True)

plot_confusion_with_totals(snakemake.input['csv'])
