#!/usr/bin/env bash
# Reproduce every result, table and figure. Run from the repository root after placing raw files in data/raw/ (see README).
set -euo pipefail
cd "$(dirname "$0")/scripts"
python 01_parse_metadata.py
python 02_load_expression.py
python 03_quality_control.py
python 03b_sex_verification.py
python 04_extract_markers.py
python 05_build_signatures.py
python 06_primary_analysis.py          # primary + sensitivity analyses, 1000-pair random gene-set test (~1 min)
python 07_gse63060_sensitivity.py      # alternative handling of GSE63060 controls
python 08_random_geneset_null.py 10000 20261123   # extended 10,000-pair test (~4 min)
python 09_figures.py
echo "Done. Outputs in results/intermediate and results/figures."
