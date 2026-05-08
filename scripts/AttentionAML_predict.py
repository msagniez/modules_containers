"""
un AttentionAML Prediction Script

Usage:
    python3 AttentionAML_predict.py <input_csv>

Arguments:
    input_csv: Path to the CSV file containing TPM expression data (with gene index in first column)

Example:
    python3 AttentionAML_predict.py TPM_test.csv
"""

import sys
import pandas as pd

# Add AttentionAML to path
sys.path.append('AttentionAML')
from AttentionAML import AttentionAML


def main():
    # Check if input file is provided
    if len(sys.argv) != 2:
        print("Error: Missing input file argument")
        print("\nUsage:")
        print("    python3 AttentionAML_predict.py <input_csv>")
        print("\nExample:")
        print("    python3 AttentionAML_predict.py TPM_test.csv")
        sys.exit(1)
    
    # Get input file path
    input_csv = sys.argv[1]
    
    # Read the CSV file
    try:
        print(f"Loading data from: {input_csv}")
        test = pd.read_csv(input_csv, index_col=0)
        print(f"Data loaded successfully. Shape: {test.shape}")
    except FileNotFoundError:
        print(f"Error: File '{input_csv}' not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        sys.exit(1)
    
    # Run prediction
    try:
        print("\nRunning AttentionAML prediction...")
        results = AttentionAML.Predict(Exp=test, exp_type='TPM')
        print("\nPrediction completed successfully!")
        return results
    except Exception as e:
        print(f"Error during prediction: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()