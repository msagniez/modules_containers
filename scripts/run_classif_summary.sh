#!/bin/bash
set -euo pipefail
 
# First (and only) argument is the output filename stem.
# The stem is used for both the CSV and the PDF.
# Examples:
#   Classification-summary          → /data/Classification-summary.csv / .pdf
#   MPXXXX_Classification-summary   → /data/MPXXXX_Classification-summary.csv / .pdf
STEM="${1:-Classification-summary}"
 
CSV="/data/${STEM}.csv"
PDF="/data/${STEM}.pdf"
 
echo "=== Step 1: build classification summary CSV ==="
python /opt/scripts/build_classification_summary.py \
    --classifiers-dir /data \
    --output "$CSV"
 
echo ""
echo "=== Step 2: plot classification summary PDF ==="
Rscript /opt/scripts/plot_classification_summary.R \
    --input  "$CSV" \
    --output "$PDF"
 
echo ""
echo "All done. Output files:"
echo "  $CSV"
echo "  $PDF"