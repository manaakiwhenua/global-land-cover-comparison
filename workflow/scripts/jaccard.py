import pandas as pd
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
def main():
    # Snakemake variables
    outputs = snakemake.output
    params = snakemake.params
    inputs = snakemake.input['confusion_matrices']
    region = snakemake.wildcards.region

    correspondence = params.jaccard_correspondence
    sources_config = params.sources_config

    jaccard_records = []

    for csv_path in inputs:
        # Verify selected sources and file handling
        src1, src2 = Path(csv_path).stem.split('--')
        print(f"Evaluating Sources: {src1}, {src2}")

        try:
            confusion_matrix = load_confusion_matrix(csv_path)
        except Exception as e:
            print(f"Error loading {csv_path}: {e}", file=sys.stderr)
            continue

        for class_info in correspondence:
            concept = list(class_info.keys())[0]
            datasets = class_info[concept]

            # Evaluate source availability in datasets for corresponding concept
            if src1 not in datasets or src2 not in datasets:
                print(f"Skipping: {concept} as {src1} or {src2} not in datasets.", file=sys.stderr)
                continue

            # Attempt retrieval of corresponding values
            values1 = datasets.get(src1, [])
            values2 = datasets.get(src2, [])

            # print(f"Using Classes - Source1 ({src1}): {values1}, Source2 ({src2}): {values2}")

            # Safely filter valid IDs only
            values1 = [v for v in values1 if v in confusion_matrix.index]
            values2 = [v for v in values2 if v in confusion_matrix.columns]
            area1 = confusion_matrix.loc[values1].sum().sum()
            area2 = confusion_matrix.loc[:, values2].sum().sum()

            if not values1 or not values2:
                print(f"{concept} - Missing classes in confusion matrix.", file=sys.stderr)
                continue

            try:
                intersection_count = confusion_matrix.loc[values1, values2].sum().sum()
                union_count = (area1 + area2 - intersection_count)

                jaccard_index = compute_jaccard_index(intersection_count, union_count)
                jaccard_records.append({
                    "Concept": concept,
                    "Source1": src1,
                    "Source2": src2,
                    "Jaccard Index": jaccard_index,
                    "Area1_ha": area1 / 100,
                    "Area2_ha": area2 / 100,
                    "Area1_kha": area1 / 100000,
                    "Area2_kha": area2 / 100000,
                    "intersection_count": intersection_count,
                    "union_count": union_count
                })
            except KeyError as ke:
                print(f"Key Error for {concept}: {ke}", file=sys.stderr)

    # Generate resulting DataFrame for CSV
    jaccard_df = pd.DataFrame(jaccard_records)
    jaccard_df.to_csv(outputs.csv, index=False)
    print(f"Jaccard index CSV saved to {outputs.csv}")

if __name__ == "__main__":
    main()

