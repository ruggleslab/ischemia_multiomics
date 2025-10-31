# Ischemia Multi-omics Prediction Model

This repository contains a predictive model and associated gene sets for classifying ischemic conditions into four risk scores (RS1-RS4). The model is implemented using the `singscore` methodology for sample scoring.

## Contents

- Prediction model
- Gene sets
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

## Citation

If you use this model in your research, please cite the following publication: [TBD]

## Contact

For questions or issues, please open an issue in this repository or contact [mm12865@nyu.edu]
