#!/usr/bin/env python3
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
    python3 format_matrix.py <input_matrix> <reference_format> <output_matrix> [--gene-id-col COLUMN]

Arguments:
    --gene-id-col: Optional column name for gene IDs if genes are in rows (default: auto-detect)

Example:
    python3 format_matrix.py my_data.csv TPM_test_format.csv formatted_output.csv
    python3 format_matrix.py my_data.tsv TPM_test_format.tsv formatted_output.tsv --gene-id-col gene_symbol
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
        # Try to detect from content
        with open(filepath, 'r') as f:
            first_line = f.readline()
            if '\t' in first_line:
                return 'tsv'
            else:
                return 'csv'


def read_matrix(filepath):
    """
    Read matrix file (CSV or TSV) with auto-detection.
    
    Args:
        filepath: Path to the file
        
    Returns:
        pandas DataFrame and file format ('csv' or 'tsv')
    """
    file_format = detect_file_format(filepath)
    
    if file_format == 'tsv':
        df = pd.read_csv(filepath, sep='\t', index_col=0)
    else:
        df = pd.read_csv(filepath, index_col=0)
    
    return df, file_format


def write_matrix(df, filepath, file_format):
    """
    Write matrix file in specified format (CSV or TSV).
    
    Args:
        df: pandas DataFrame
        filepath: Path to save the file
        file_format: 'csv' or 'tsv'
    """
    if file_format == 'tsv':
        df.to_csv(filepath, sep='\t')
    else:
        df.to_csv(filepath)


def detect_gene_orientation(df, reference_genes):
    """
    Detect if genes are in rows or columns.
    
    Args:
        df: Input dataframe
        reference_genes: List of reference gene names
        
    Returns:
        'columns' if genes are columns, 'rows' if genes are rows
    """
    # Check overlap with columns
    col_overlap = len(set(df.columns) & set(reference_genes))
    
    # Check overlap with index
    row_overlap = len(set(df.index) & set(reference_genes))
    
    if col_overlap > row_overlap:
        return 'columns'
    elif row_overlap > col_overlap:
        return 'rows'
    else:
        # If equal or both zero, default to columns
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
    best_col = None
    best_overlap = 0
    
    for col in df.columns:
        overlap = len(set(df[col].astype(str)) & set(reference_genes))
        if overlap > best_overlap:
            best_overlap = overlap
            best_col = col
    
    return best_col if best_overlap > 0 else None


def normalize_gene_names(genes):
    """
    Normalize gene names for matching (case-insensitive, strip whitespace).
    
    Args:
        genes: List of gene names
        
    Returns:
        Dictionary mapping normalized names to original names
    """
    normalized = {}
    for gene in genes:
        normalized[str(gene).strip().upper()] = gene
    return normalized


def format_matrix(input_file, reference_file, output_file, gene_id_col=None):
    """
    Format input matrix to match reference format.
    
    Args:
        input_file: Path to input matrix
        reference_file: Path to reference format
        output_file: Path to save formatted output
        gene_id_col: Optional column name for gene IDs (if genes are in rows)
    """
    print(f"Reading reference format from: {reference_file}")
    reference_df, ref_format = read_matrix(reference_file)
    print(f"Reference file format: {ref_format.upper()}")
    
    # Detect if reference has genes as rows or columns
    ref_orientation = 'columns' if reference_df.shape[1] > reference_df.shape[0] else 'rows'
    
    if ref_orientation == 'columns':
        reference_genes = reference_df.columns.tolist()
        print(f"Reference format has {len(reference_genes)} genes as columns")
    else:
        reference_genes = reference_df.index.tolist()
        print(f"Reference format has {len(reference_genes)} genes as rows")
        # Transpose for uniform processing
        reference_df = reference_df.T
        reference_genes = reference_df.columns.tolist()
    
    print(f"\nReading input matrix from: {input_file}")
    input_df, input_format = read_matrix(input_file)
    print(f"Input file format: {input_format.upper()}")
    print(f"Input matrix has {input_df.shape[0]} rows and {input_df.shape[1]} columns")
    
    # Detect orientation
    orientation = detect_gene_orientation(input_df, reference_genes)
    print(f"Detected genes as: {orientation}")
    
    # Handle genes in rows
    if orientation == 'rows':
        # If gene_id_col specified, use it
        if gene_id_col:
            if gene_id_col in input_df.columns:
                input_df = input_df.set_index(gene_id_col)
                print(f"Using '{gene_id_col}' column as gene identifiers")
            else:
                print(f"Warning: Specified column '{gene_id_col}' not found. Using index.")
        else:
            # Try to find gene column automatically
            gene_col = find_gene_column(input_df, reference_genes)
            if gene_col:
                print(f"Auto-detected gene column: '{gene_col}'")
                input_df = input_df.set_index(gene_col)
            else:
                print("Using existing index as gene identifiers")
        
        # Transpose so genes are columns
        input_df = input_df.T
        print(f"Transposed input matrix: now {input_df.shape[0]} rows × {input_df.shape[1]} columns")
    
    # Normalize gene names for matching
    ref_normalized = normalize_gene_names(reference_genes)
    input_normalized = normalize_gene_names(input_df.columns)
    
    # Create mapping from input genes to reference genes
    gene_mapping = {}
    for norm_input, input_gene in input_normalized.items():
        if norm_input in ref_normalized:
            gene_mapping[input_gene] = ref_normalized[norm_input]
    
    # Rename input columns to match reference
    input_df = input_df.rename(columns=gene_mapping)
    
    # Find matching, missing, and extra genes
    matching_genes = [gene for gene in reference_genes if gene in input_df.columns]
    missing_genes = [gene for gene in reference_genes if gene not in input_df.columns]
    extra_genes = [gene for gene in input_df.columns if gene not in reference_genes]
    
    print(f"\nGene analysis:")
    print(f"  - Matching genes: {len(matching_genes)}")
    print(f"  - Missing genes (will be filled with zeros): {len(missing_genes)}")
    print(f"  - Extra genes (will be removed): {len(extra_genes)}")
    
    # Create output dataframe with matching genes
    output_df = input_df[matching_genes].copy()
    
    # Add missing genes with zeros
    for gene in missing_genes:
        output_df[gene] = 0.0
    
    # Reorder genes to match reference format
    output_df = output_df[reference_genes]
    
    # Transpose back if reference was originally in rows
    if ref_orientation == 'rows':
        output_df = output_df.T
        print(f"\nTransposing output to match reference format (genes as rows)")
    
    # Save to output file in the same format as reference
    print(f"\nSaving formatted matrix to: {output_file}")
    print(f"Output file format: {ref_format.upper()} (matching reference)")
    write_matrix(output_df, output_file, ref_format)
    print(f"✓ Output matrix saved with {output_df.shape[0]} rows and {output_df.shape[1]} columns")
    
    # Print summary statistics
    if missing_genes:
        print(f"\nNote: {len(missing_genes)} genes were added with zero values")
        if len(missing_genes) <= 10:
            print(f"Missing genes: {', '.join(missing_genes)}")
        else:
            print(f"First 10 missing genes: {', '.join(missing_genes[:10])}...")
    
    if extra_genes:
        print(f"\nNote: {len(extra_genes)} extra genes were removed")
        if len(extra_genes) <= 10:
            print(f"Removed genes: {', '.join(extra_genes)}")
        else:
            print(f"First 10 removed genes: {', '.join(extra_genes[:10])}...")


def main():
    parser = argparse.ArgumentParser(
        description='Format input matrix to match AttentionAML reference format.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 format_matrix.py my_data.csv TPM_test_format.csv formatted_output.csv
  python3 format_matrix.py my_data.tsv TPM_test_format.tsv formatted_output.tsv --gene-id-col gene_symbol
  
Note: Output format (CSV/TSV) will automatically match the reference file format.
        """
    )
    
    parser.add_argument('input_file', help='Path to input matrix (CSV or TSV)')
    parser.add_argument('reference_file', help='Path to reference format (CSV or TSV)')
    parser.add_argument('output_file', help='Path to save formatted output')
    parser.add_argument('--gene-id-col', help='Column name for gene IDs (if genes are in rows)', default=None)
    
    args = parser.parse_args()
    
    try:
        format_matrix(args.input_file, args.reference_file, args.output_file, args.gene_id_col)
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