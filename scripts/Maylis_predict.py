"""
Maylis Prediction Script

Usage:
    python3 Maylis_predict.py <input_tsv>

Arguments:
    input_tsv: Path to the TSV file containing TPM expression data (with gene ID in first column)

Example:
    python3 Maylis_predict.py TPM_test.tsv
    python3 Maylis_predict.py TPM_test.tsv results/

Arguments:
    input_tsv:  Path to the TSV file containing TPM expression data (with gene ID in first column)
    output_dir: (Optional) Directory where outputs will be saved. Defaults to current directory.
"""

import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sklearn
from pickle import load
import urllib.request
import os


def subgroup_classification(sample_data, features, model):
    ## We remove gene ID versions (e.g. ENSG00000000460.16) so that any genome version used for alignment with STAR+RSEM is recognized by the classifier.
    sample_data["gene_id_without_version"] = sample_data["gene_id"].apply(lambda x: x.split('.', 1)[0])
    # Genes with 2 transcripts (e.g. with "_PAR_Y" in the second transcript) are deleted.
    sample_data.drop_duplicates(subset="gene_id_without_version", keep='first', inplace=True)
    del sample_data['gene_id']
    sample_data = sample_data.rename(columns={"gene_id_without_version": "gene_id"})
    # Shift the gene_id column to the first position
    cols = list(sample_data)
    cols.insert(0, cols.pop(cols.index("gene_id")))
    sample_data = sample_data.loc[:, cols]

    ## We keep the genes (ENSG) needed for the classifier (protein-coding genes, excluding those on the sex chromosomes).
    # Gene order must remain the same as in features.tsv
    sample_data = features.merge(sample_data, on='gene_id', how='left')
    # For required genes not present in the sample, we replace their expression after merge (NaN) by 0.
    sample_data['TPM'] = sample_data['TPM'].fillna(0)

    ## Data formatting
    sample_data.set_index("gene_id", inplace=True)
    # Log transformation on gene expression data
    sample_data = sample_data.astype(float)
    sample_data = np.log2(sample_data + 1)
    sample_array = sample_data[["TPM"]].T.to_numpy()

    ## Applying the classification model to the sample
    prediction_proba = model.predict_proba(sample_array)
    # We add the equivalences with the names of AML subtypes
    subgroups = [["KMT2At", "KMT2A translocations"],
                 ["t(15;17)", "t(15;17)(q24.1;q21.2)/PML-RARA"],
                 ["IDH2m", "IDH2 R172 mutations"],
                 ["NPM1m", "NPM1 mutations"],
                 ["Chr.Spl.", "Mutations in genes encoding chromatin compaction and/or splicing regulation"],
                 ["Inter.", "Normal or intermediate karyotype"],
                 ["Complex", "Complex karyotype and/or TP53 mutations"],
                 ["Mono5/7", "Monosomy or deletion of 5(q) and/or 7(q)"],
                 ["NUP98t", "NUP98 translocations"],
                 ["t(8;21)", "t(8;21)(q22;q22.1)/RUNX1-RUNX1T1"],
                 ["inv(16)", "inv(16)(p13.1q22)/CBFB-MYH11"],
                 ["CEBPAm", "CEBPA b-Zip in-frame mutations"],
                 ["MECOMr", "MECOM rearrangements"],
                 ["t(6;9)", "t(6;9)(p22.3;q34.1)/DEK-NUP214"],
                 ["others", "AML others"],
                 ["ETV6t", "ETV6 translocations"],
                 ["CBFA2T3t", "inv(16)(p13q24)/CBFA2T3-GLIS2"]]
    # Create the pandas DataFrame
    prediction_df = pd.DataFrame(subgroups, columns=['Subgroup_def', 'Subgroup_fr'])
    prediction_df['Probability'] = prediction_proba.flatten()
    # Sort probabilities in descending order
    prediction_df = prediction_df.sort_values(by="Probability", axis=0, ascending=False)
    prediction_df.reset_index(drop=True, inplace=True)
    prediction_df["Probability"] = prediction_df["Probability"] * 100

    return prediction_df

def main():

    # Set-up classifier
    # Download files if they are not there.
    if not os.path.exists("/opt/models"):
        os.makedirs("/opt/models")
    # File containing genes (features) taken into account by the classifier
    if not os.path.exists("/opt/models/features.tsv"):
        urllib.request.urlretrieve("https://nextcloud.computecanada.ca/index.php/s/Ax6fakA2aCSrj3e/download",
                                   "/opt/models/features.tsv")
    # Classifier
    if not os.path.exists("/opt/models/classifier.pkl"):
        urllib.request.urlretrieve("https://nextcloud.computecanada.ca/index.php/s/58Eb4gMGtiBYBWN/download",
                                   "/opt/models/classifier.pkl")

    # Gene expression data (TSV) from sample BA3192R from patient 2475 of the BeatAML cohort
    if not os.path.exists("/opt/test_data"):
        os.makedirs("/opt/test_data")
    if not os.path.exists("/opt/test_data/sample_example.tsv"):
        urllib.request.urlretrieve("https://nextcloud.computecanada.ca/index.php/s/sP3kmHJ87cDgPXN/download",
                                   "/opt/test_data/sample_example.tsv")

    color_palette = {
        "KMT2At": "#AFDFEF",
        "t(15;17)": "#08306B",
        "IDH2m": "#E6AB02",
        "NPM1m": "#C70E7B",
        "Chr.Spl.": "#009999",
        "Inter.": "#0079BA",
        "Complex": "#78222E",
        "Mono5/7": "#EE7A1A",
        "NUP98t": "#FFC0CB",
        "t(8;21)": "#6C6C9D",
        "inv(16)": "#DB7676",
        "CEBPAm": "#CF5BE3",
        "MECOMr": "#00A4DB",
        "t(6;9)": "#1C7A3F",
        "others": "#FDAE6B",
        "ETV6t": "#7F38CF",
        "CBFA2T3t": "#A9CC51"
    }

    features = pd.read_csv('/opt/models/features.tsv', sep='\t')
    with open("/opt/models/classifier.pkl", "rb") as f:
        model = load(f)

    # Check if input file is provided
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Error: Wrong number of arguments")
        print("\nUsage:")
        print("    python3 Maylis_predict.py <input_tsv> [output_dir]")
        print("\nExample:")
        print("    python3 Maylis_predict.py TPM_test.tsv")
        print("    python3 Maylis_predict.py TPM_test.tsv results/")
        sys.exit(1)

    # Get input file path
    input_tsv = sys.argv[1]

    # Get output directory (default to current directory)
    output_dir = sys.argv[2] if len(sys.argv) == 3 else "."
    os.makedirs(output_dir, exist_ok=True)

    # Read the TSV file
    try:
        print(f"Loading data from: {input_tsv}")
        # Fix 2: removed index_col=0 so that the first column (gene_id) is kept as a regular column
        test = pd.read_csv(input_tsv, sep='\t')
        print(f"Data loaded successfully. Shape: {test.shape}")
    except FileNotFoundError:
        print(f"Error: File '{input_tsv}' not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading TSV file: {e}")
        sys.exit(1)

    # Fix 3: pass the DataFrame 'test' instead of the file path string 'input_tsv'
    # Fix 4: removed the duplicate/incorrectly-indented model loading block that was here
    sample_prediction = subgroup_classification(test, features, model)

    # Derive output base name from the input file name (strip directory and extension)
    base_name = os.path.splitext(os.path.basename(input_tsv))[0]
    output_table = os.path.join(output_dir, f"{base_name}_predictions.tsv")
    output_figure = os.path.join(output_dir, f"{base_name}_predictions.png")

    # Save the prediction table
    sample_prediction.to_csv(output_table, sep='\t', index=False)
    print(f"Prediction table saved to: {output_table}")

    # Display classification results
    plt.figure(figsize=(12, 8))
    bar = sns.barplot(sample_prediction, x="Subgroup_def", y="Probability", hue="Subgroup_def", palette=color_palette, legend=False)
    for i in range(len(sample_prediction)):
        bar.text(i, sample_prediction.loc[i, "Probability"]+0.3, str(round(sample_prediction.loc[i, "Probability"], 1)),
                 fontdict=dict(fontsize=8),
                 horizontalalignment='center')
    plt.text(7, 80,
             f"Most probable AML subtypes",
             fontsize=12,
             fontweight="bold",
             ha="left",
             va="baseline",
             bbox=dict(facecolor='none', edgecolor='black', boxstyle='round', pad=0.4))
    plt.text(7, 66,
             f"1. {sample_prediction.loc[0, 'Subgroup_fr']}: {(sample_prediction.loc[0, 'Probability']):.1f}%\n2. {sample_prediction.loc[1, 'Subgroup_fr']}: {(sample_prediction.loc[1, 'Probability']):.1f}%\n3. {sample_prediction.loc[2, 'Subgroup_fr']}: {(sample_prediction.loc[2, 'Probability']):.1f}%",
             fontsize=12,
             ha="left",
             va="baseline")
    sns.despine()
    plt.ylim(0, 102)
    plt.xticks(fontsize=8.5)
    plt.xlabel("AML subtypes", fontsize=10, fontdict=dict(weight='bold'))
    plt.ylabel("Classification probability", fontsize=10, fontdict=dict(weight='bold'))
    plt.tight_layout()
    plt.savefig(output_figure, dpi=150)
    print(f"Prediction figure saved to: {output_figure}")


if __name__ == "__main__":
    main()