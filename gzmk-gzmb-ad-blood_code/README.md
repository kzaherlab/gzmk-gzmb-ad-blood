# gzmk-gzmb-ad-blood

Analysis code for:

> Zaher, K.A.; [co-authors]. Granzyme B–Skewed Effector-Memory CD8⁺ T-Cell Transcriptional Signatures in Whole Blood from Individuals with Mild Cognitive Impairment and Alzheimer's Disease: A Three-Cohort Analysis. *[Journal]* **2026/2027**, [volume], [article]. https://doi.org/[DOI]

The code scores GZMK⁺ and GZMB⁺ effector-memory CD8⁺ T-cell gene signatures (derived from the Terekhova et al. 2023 lifespan blood atlas) in three public whole-blood microarray cohorts, computes the K/B index (standardized GZMK-Tem score minus GZMB-Tem score), and compares it between cognitively normal, MCI and AD participants with per-cohort regression and random-effects meta-analysis.

**Interpretation note.** The K/B index is a relative bulk-transcriptional measure. It does not count GZMK⁺ or GZMB⁺ CD8⁺ T cells.

## 1. Get the input data (all public)

Place these files in `data/raw/` with exactly these names (or edit `scripts/config.py`):

| File name in `data/raw/` | Source |
| --- | --- |
| `GSE63060_series_matrix.txt.gz` | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63060/matrix/ |
| `GSE63061_series_matrix.txt.gz` | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63061/matrix/ |
| `GSE140829_series_matrix.txt.gz` | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140829/matrix/ |
| `GPL6947_full_table.txt` | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6947 → "Download full table" |
| `GPL10558_full_table.txt` | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10558 → "Download full table" |
| `GPL15988_full_table.txt` | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL15988 → "Download full table" |
| `Terekhova2023_TableS5_cluster_markers.xlsx` | Supplementary cluster-marker table (file `mmc6.xlsx`) of Terekhova et al., *Immunity* 2023, https://doi.org/10.1016/j.immuni.2023.10.013 |

Keep the `.gz` files compressed. The derived signature gene lists are also provided in `results/derived/` so the signatures can be inspected without the atlas file.

## 2. Install and run

```bash
python -m venv .venv && source .venv/bin/activate
pip install numpy scipy pandas matplotlib openpyxl   # versions in requirements.txt
./run_all.sh
```

Runtime is about 6 minutes on a laptop. Random seeds are fixed (20261116 for the prespecified analyses; 20261123 for the extended random-gene-set test).

## 3. Scripts and what they produce

| Script | Produces | Manuscript item |
| --- | --- | --- |
| `01_parse_metadata.py` | `meta_raw.csv` (diagnosis, age, sex, APOE, batch) | Table 1 |
| `02_load_expression.py` | per-cohort gene-level expression (probe → symbol; highest-mean probe per gene) | — |
| `03_quality_control.py` | `keep_samples.json`, `qc_log.csv` (diagnosis, sex check, PCA outliers, duplicates) | Table S2 |
| `03b_sex_verification.py` | `sex_verification.json`, `GSE63060_sexcheck.csv` | Figure S2A–B, Table S6 |
| `04_extract_markers.py` | `terekhova_markers.json` | — |
| `05_build_signatures.py` | `signatures.json`, `TableS1_signature_filtering.csv` | Figure 1A–B, Table S1 |
| `06_primary_analysis.py` | `results.json`, `scores_all.csv`, `gene_level.csv` | Figures 1C, 2, 3, 4A, 5; Tables 2, S3–S5 |
| `07_gse63060_sensitivity.py` | `qc_scenarios.json` | Figure S2C, Table S5 |
| `08_random_geneset_null.py N SEED` | `null_N.json` | Figure 4B–C, Table S5 |
| `09_figures.py` | Figures 1–5, S1, S2 (PNG + PDF) | all figures |

Outputs go to `results/intermediate/` and `results/figures/`.

## 4. Methods in brief

- **Signatures.** The 100 top-ranked markers of the atlas clusters "Tem GZMK+" and "Tem GZMB+"; shared genes, immediate-early/dissociation genes, Ig/HLA class II/ribosomal/mitochondrial genes and top myeloid/B-cell markers removed; genes on all three platforms kept if, in every cohort, the leave-one-out correlation with the rest of the signature was ≥ 0.2 and exceeded the correlation with a neutrophil proxy (no diagnosis labels used).
- **Scoring.** ssGSEA (GSVA formulation, α = 0.25, no normalization), standardized to controls within cohort.
- **Models.** OLS per cohort with age, sex and T-cell, CD8, NK and neutrophil marker proxies (plus batch in GSE140829); SMD = coefficient / residual SD; REML random-effects pooling.
- **Specificity.** Each signature gene replaced by a random gene from the same mean-expression decile; empirical p = (r + 1)/(N + 1).

## 5. Known data issue

In GSE63060, recorded sex disagrees with expression-inferred sex for 46 of 104 controls (none of 225 cases), and recorded sex among controls is no better than chance. The prespecified analysis excludes the 46 discordant controls; `07_gse63060_sensitivity.py` reruns the analysis with all controls retained (inferred sex) and with all GSE63060 controls excluded.

## License

MIT (see `LICENSE`). Data remain subject to the terms of their original sources.
