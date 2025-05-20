from collections import Counter
import csv
import numpy as np
import os
from pathlib import Path

import rasterio

def window_generator(dataset, window_size):
    width = dataset.width
    height = dataset.height

    for top in range(0, height, window_size):
        for left in range(0, width, window_size):
            bottom = min(top + window_size, height)
            right = min(left + window_size, width)
            window = rasterio.windows.Window(left, top, right - left, bottom - top)
            yield window

def compute_confusion(src1_path, src2_path, window_size=512):
    counter = Counter()
    with rasterio.open(src1_path) as ds1, rasterio.open(src2_path) as ds2:
        assert ds1.width == ds2.width and ds1.height == ds2.height, "Rasters must be aligned"

        for window in window_generator(ds1, window_size):
            a = ds1.read(1, window=window)
            b = ds2.read(1, window=window)
            mask = (a != ds1.nodata) & (b != ds2.nodata)
            pairs = zip(a[mask].flatten(), b[mask].flatten())
            counter.update(pairs)
    return counter

def compute_confusion_csv(src1_path, src2_path, output_csv_path, window_size=512):
    counter = compute_confusion(src1_path, src2_path, window_size)
    
    # Find all unique labels
    labels = sorted(set(k for k, _ in counter.keys()) | set(k for _, k in counter.keys()))
    label_indices = {k: i for i, k in enumerate(labels)}
    
    # Create confusion matrix
    matrix = np.zeros((len(labels), len(labels)), dtype=int)
    for (a_label, b_label), count in counter.items():
        i = label_indices[a_label]
        j = label_indices[b_label]
        matrix[i, j] = count
    
    nonzero_rows = np.any(matrix != 0, axis=1)
    nonzero_cols = np.any(matrix != 0, axis=0)

    matrix = matrix[np.ix_(nonzero_rows, nonzero_cols)]
    row_labels = [label for i, label in enumerate(labels) if nonzero_rows[i]]
    col_labels = [label for i, label in enumerate(labels) if nonzero_cols[i]]

    row_totals = matrix.sum(axis=1)
    col_totals = matrix.sum(axis=0)
    grand_total = matrix.sum()

    num_rows, num_cols = matrix.shape
    # Create a matrix with totals, size (num_rows+1, num_cols+1)
    matrix_with_totals = np.zeros((num_rows + 1, num_cols + 1), dtype=int)
    matrix_with_totals[:num_rows, :num_cols] = matrix
    matrix_with_totals[:num_rows, num_cols] = row_totals
    matrix_with_totals[num_rows, :num_cols] = col_totals
    matrix_with_totals[num_rows, num_cols] = grand_total

    # Write CSV
    with open(output_csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([""] + col_labels + ["Total"])
        for i, row_label in enumerate(row_labels):
            writer.writerow([row_label] + matrix_with_totals[i, :-1].tolist() + [matrix_with_totals[i, -1]])
        writer.writerow(["Total"] + matrix_with_totals[-1, :-1].tolist() + [matrix_with_totals[-1, -1]])


Path(snakemake.output[0]).parent.mkdir(parents=True, exist_ok=True)

compute_confusion_csv(
    snakemake.input.src1,
    snakemake.input.src2,
    snakemake.output[0],
    snakemake.params.window_size
)