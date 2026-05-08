"""
Script to format input matrix to match AttentionAML reference format.

This script:
1. Reads the reference format file to get the expected genes
2. Reads the submitted matrix
3. Automatically detects if genes are in rows or columns
4. Matches genes by ID or symbol
5. Keeps only genes that exist in the reference format
6. Adds missing genes with zero values
7. Reorders genes to match the reference format
8. Outputs in the same format (CSV/TSV) as the reference file

Usage:
    python3 format_matrix.py <input_matrix> <reference_format> <output_matrix> [--gene-id-col COLUMN] [--strip-versions]

Arguments:
    --gene-id-col:    Optional column name for gene IDs if genes are in rows (default: auto-detect)
    --strip-versions: Ignore Ensembl version suffixes when matching genes
                      (e.g. treats ENSG00000000003.17 and ENSG00000000003.15 as the same gene)

Example:
    python3 format_matrix.py my_data.csv TPM_test_format.csv formatted_output.csv
    python3 format_matrix.py my_data.tsv TPM_test_format.tsv formatted_output.tsv --gene-id-col gene_symbol
    python3 format_matrix.py my_data.csv TPM_test_format.csv formatted_output.csv --strip-versions
"""

import sys
import pandas as pd
import numpy as np
import argparse
import os


def detect_file_format(filepath):
    """
    Detect if file is CSV or TSV based on extension and content.

    Args:
        filepath: Path to the file

    Returns:
        'csv' or 'tsv'
    """
    ext = os.path.splitext(filepath)[1].lower()
    if ext in ['.tsv', '.txt']:
        return 'tsv'
    elif ext == '.csv':
        return 'csv'
    else:
        with open(filepath, 'r') as f:
            first_line = f.readline()
            return 'tsv' if '\t' in first_line else 'csv'


def read_matrix(filepath):
    """
    Read matrix file (CSV or TSV) with auto-detection.

    Args:
        filepath: Path to the file

    Returns:
        pandas DataFrame and file format ('csv' or 'tsv')
    """
    file_format = detect_file_format(filepath)
    sep = '\t' if file_format == 'tsv' else ','
    df = pd.read_csv(filepath, sep=sep, index_col=0)
    return df, file_format


def write_matrix(df, filepath, file_format):
    """
    Write matrix file in specified format (CSV or TSV).

    Args:
        df: pandas DataFrame
        filepath: Path to save the file
        file_format: 'csv' or 'tsv'
    """
    sep = '\t' if file_format == 'tsv' else ','
    df.to_csv(filepath, sep=sep)


def strip_version(gene_name):
    """Strip Ensembl version suffix (e.g. ENSG00000000003.17 → ENSG00000000003)."""
    return str(gene_name).rsplit('.', 1)[0]


def detect_gene_orientation(df, reference_genes):
    """
    Detect if genes are in rows or columns.

    Args:
        df: Input dataframe
        reference_genes: List of reference gene names

    Returns:
        'columns' if genes are columns, 'rows' if genes are rows
    """
    col_overlap = len(set(df.columns) & set(reference_genes))
    row_overlap = len(set(df.index) & set(reference_genes))
    if col_overlap > row_overlap:
        return 'columns'
    elif row_overlap > col_overlap:
        return 'rows'
    return 'columns'


def find_gene_column(df, reference_genes):
    """
    Find which column in the dataframe contains gene identifiers.

    Args:
        df: Input dataframe
        reference_genes: List of reference gene names

    Returns:
        Column name with best gene match, or None
    """
    best_col, best_overlap = None, 0
    for col in df.columns:
        overlap = len(set(df[col].astype(str)) & set(reference_genes))
        if overlap > best_overlap:
            best_overlap, best_col = overlap, col
    return best_col if best_overlap > 0 else None


def normalize_gene_names(genes, strip_versions=False):
    """
    Normalize gene names for matching (case-insensitive, strip whitespace).

    Args:
        genes: List of gene names
        strip_versions: If True, also strip Ensembl version suffixes

    Returns:
        Dictionary mapping normalized names to original names
    """
    normalized = {}
    for gene in genes:
        key = str(gene).strip().upper()
        if strip_versions:
            key = strip_version(key)
        normalized[key] = gene
    return normalized


def format_matrix(input_file, reference_file, output_file,
                  gene_id_col=None, strip_versions=False):
    """
    Format input matrix to match reference format.

    Args:
        input_file: Path to input matrix
        reference_file: Path to reference format
        output_file: Path to save formatted output
        gene_id_col: Optional column name for gene IDs (if genes are in rows)
        strip_versions: If True, ignore Ensembl version suffixes when matching genes
    """
    print(f"Reading reference format from: {reference_file}")
    reference_df, ref_format = read_matrix(reference_file)
    print(f"Reference file format: {ref_format.upper()}")

    ref_orientation = 'columns' if reference_df.shape[1] > reference_df.shape[0] else 'rows'

    if ref_orientation == 'columns':
        reference_genes = reference_df.columns.tolist()
        print(f"Reference format has {len(reference_genes)} genes as columns")
    else:
        reference_genes = reference_df.index.tolist()
        print(f"Reference format has {len(reference_genes)} genes as rows")
        reference_df = reference_df.T
        reference_genes = reference_df.columns.tolist()

    print(f"\nReading input matrix from: {input_file}")
    input_df, input_format = read_matrix(input_file)
    print(f"Input file format: {input_format.upper()}")
    print(f"Input matrix has {input_df.shape[0]} rows and {input_df.shape[1]} columns")

    if strip_versions:
        print("Version stripping enabled: gene version suffixes will be ignored during matching")

    orientation = detect_gene_orientation(input_df, reference_genes)
    print(f"Detected genes as: {orientation}")

    if orientation == 'rows':
        if gene_id_col:
            if gene_id_col in input_df.columns:
                input_df = input_df.set_index(gene_id_col)
                print(f"Using '{gene_id_col}' column as gene identifiers")
            else:
                print(f"Warning: Specified column '{gene_id_col}' not found. Using index.")
        else:
            gene_col = find_gene_column(input_df, reference_genes)
            if gene_col:
                print(f"Auto-detected gene column: '{gene_col}'")
                input_df = input_df.set_index(gene_col)
            else:
                print("Using existing index as gene identifiers")

        input_df = input_df.T
        print(f"Transposed input matrix: now {input_df.shape[0]} rows × {input_df.shape[1]} columns")

    # Normalize and map gene names
    ref_normalized = normalize_gene_names(reference_genes, strip_versions)
    input_normalized = normalize_gene_names(input_df.columns, strip_versions)

    gene_mapping = {}
    for norm_input, input_gene in input_normalized.items():
        if norm_input in ref_normalized:
            gene_mapping[input_gene] = ref_normalized[norm_input]

    input_df = input_df.rename(columns=gene_mapping)

    matching_genes = [gene for gene in reference_genes if gene in input_df.columns]
    missing_genes  = [gene for gene in reference_genes if gene not in input_df.columns]
    extra_genes    = [gene for gene in input_df.columns  if gene not in reference_genes]

    print(f"\nGene analysis:")
    print(f"  - Matching genes: {len(matching_genes)}")
    print(f"  - Missing genes (will be filled with zeros): {len(missing_genes)}")
    print(f"  - Extra genes (will be removed): {len(extra_genes)}")

    # Keep only matching genes first
    output_df = input_df[matching_genes].copy()

    # Add all missing-gene columns in one concat instead of a loop (avoids PerformanceWarning)
    if missing_genes:
        zeros_df = pd.DataFrame(
            0.0,
            index=output_df.index,
            columns=missing_genes
        )
        output_df = pd.concat([output_df, zeros_df], axis=1)

    # Reorder to match reference
    output_df = output_df[reference_genes]

    if ref_orientation == 'rows':
        output_df = output_df.T
        print(f"\nTransposing output to match reference format (genes as rows)")

    print(f"\nSaving formatted matrix to: {output_file}")
    print(f"Output file format: {ref_format.upper()} (matching reference)")
    write_matrix(output_df, output_file, ref_format)
    print(f"✓ Output matrix saved with {output_df.shape[0]} rows and {output_df.shape[1]} columns")

    if missing_genes:
        print(f"\nNote: {len(missing_genes)} genes were added with zero values")
        preview = missing_genes[:10]
        suffix = '...' if len(missing_genes) > 10 else ''
        print(f"{'First 10 m' if suffix else 'M'}issing genes: {', '.join(preview)}{suffix}")

    if extra_genes:
        print(f"\nNote: {len(extra_genes)} extra genes were removed")
        preview = extra_genes[:10]
        suffix = '...' if len(extra_genes) > 10 else ''
        print(f"First 10 removed genes: {', '.join(preview)}{suffix}")


def main():
    parser = argparse.ArgumentParser(
        description='Format input matrix to match AttentionAML reference format.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 format_matrix.py my_data.csv TPM_test_format.csv formatted_output.csv
  python3 format_matrix.py my_data.tsv TPM_test_format.tsv formatted_output.tsv --gene-id-col gene_symbol
  python3 format_matrix.py my_data.csv ref.csv out.csv --strip-versions
        """
    )
    parser.add_argument('input_file',      help='Path to input matrix (CSV or TSV)')
    parser.add_argument('reference_file',  help='Path to reference format (CSV or TSV)')
    parser.add_argument('output_file',     help='Path to save formatted output')
    parser.add_argument('--gene-id-col',   help='Column name for gene IDs (if genes are in rows)', default=None)
    parser.add_argument('--strip-versions',
                        action='store_true',
                        help='Ignore Ensembl version suffixes when matching genes '
                             '(e.g. treats ENSG00000000003.17 and ENSG00000000003.15 as the same gene)')

    args = parser.parse_args()

    try:
        format_matrix(args.input_file, args.reference_file, args.output_file,
                      args.gene_id_col, args.strip_versions)
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()