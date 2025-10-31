# Ischemia Multi-omics Prediction Model

This repository contains a predictive model and associated gene sets for classifying ischemic conditions into four risk scores (RS1-RS4). The model is implemented using the `singscore` methodology for sample scoring.

## Contents

- Prediction model
- Gene sets:
  - RS1 gene signature
  - RS2 gene signature
  - RS3 gene signature
  - RS4 gene signature
- Usage vignette

## Installation

```R
# Install required packages
install.packages("singscore")
```

## Scoring

This repository provides two scripts for scoring new samples:

### RNA Subtype Scoring (`src/rna_subtype_scoring.R`)

**Inputs:**
- Expression matrix CSV (`genes = rows`, `samples = columns`, first column = gene symbols)

**Outputs:**
- `predictions.csv` (sample IDs and predicted RS)
- `rs_distribution.txt` (distribution of RS labels)

**Run:**
Edit the config variables at the top of the script, then run:
```bash
Rscript src/rna_subtype_scoring.R
```

### Methylation Subtype Scoring (`src/methylation_subtype_scoring.R`)

**Inputs:**
- New data CSV (`samples = rows`, `features = columns`)
- Feature names CSV (column `feature`)
- Trained model RDS

**Outputs:**
- `methylation_subtype_predictions.csv` (sample predictions and probabilities)
- `class_distribution.txt` (class counts)

**Run:**
Edit the config variables at the top of the script, then run:
```bash
Rscript src/methylation_subtype_scoring.R
```

## Scoring scripts

This repository includes two R scripts in `src/` to score new samples and produce predictions. Edit the configuration variables at the top of each script to point to your input files and output directory, then run them with `Rscript` or from an R session.

### RNA scoring — `src/rna_subtype_scoring.R`

Purpose

- Compute RNA-based RS subtype scores and predict RS class for each sample.

Inputs

- `expression_matrix_file`: path to a CSV expression matrix (genes = rows, samples = columns). The first column should contain gene symbols; the remaining columns are sample expression values.
- `outdir`: output directory where results will be saved.

Dependencies

- `singscore`, `vroom`, `tidyverse` (install via `install.packages()` or `pacman::p_load`).

What the script does

- Downloads gene sets (`genesets/RS*_genes.csv`) and the RNA model (`models/rs_model.rds`) from the repository (main branch) by default.
- Ranks genes per sample with `rankGenes()` (from `singscore`), computes `simpleScore()` for each RS gene set, scales the scores, and applies the `rs_model` to predict RS class.

Outputs

- `predictions.csv` — table of sample IDs and predicted RS.
- `rs_distribution.txt` — distribution of predicted RS labels.
- `session.log` and `session.RData` for reproducibility/debugging.

Run

Edit the variables at the top of `src/rna_subtype_scoring.R` then run:

```bash
Rscript src/rna_subtype_scoring.R
```

### Methylation scoring — `src/methylation_subtype_scoring.R`

Purpose

- Score new methylation samples using a trained classifier and output predicted subtype/probabilities.

Inputs

- `new_data_file`: CSV path to new data. Format expected by the script: samples = rows, features = columns.
- `feature_names_file`: CSV containing the feature names used during training (column named `feature`). The script checks that all required features are present.
- `model_file`: path to the trained model RDS used for prediction (for example, a caret model saved with `saveRDS()`).
- `scoring_outdir`: directory where scoring outputs will be written.

Dependencies

- `tidyverse`, `data.table`, `caret` (install via `install.packages()`).

What the script does

- Validates that the new data contains all required features, applies an inverse-normal (van der Waerden) transform per feature, orders features to match the training set, and runs the model to produce class predictions and class probabilities.

Outputs

- `methylation_subtype_predictions.csv` — sample-level predictions and per-class probabilities.
- `class_distribution.txt` — count of predicted classes.

Run

Edit the variables at the top of `src/methylation_subtype_scoring.R` then run:

```bash
Rscript src/methylation_subtype_scoring.R
```

Notes

- The scripts are intentionally simple: they expect you to edit the top-of-file configuration variables with absolute or repo-relative paths. If you prefer programmatic invocation, you can refactor them to accept command-line arguments (e.g., with `argparse` or `optparse`).
- Model and gene set files in this repository:
    - `models/rs_model.rds` — RNA subtype prediction model used by `src/rna_subtype_scoring.R`.
    - `models/ms_knn_model.rds` — example methylation model (if present).
    - `genesets/RS*_genes.csv` — gene sets used for scoring.

## Citation

If you use this model in your research, please cite the following publication: [TBD]

## Contact

For questions or issues, please open an issue in this repository or contact [mm12865@nyu.edu]
