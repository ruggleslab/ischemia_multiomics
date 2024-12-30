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

## Usage

### Basic workflow

1. Load your expression data
2. Score samples using singscore
3. Apply the prediction model

```R
library(singscore)

# Load gene sets
rs1_genes <- read.csv("genesets/RS1_genes.csv", header = FALSE)$X
rs2_genes <- read.csv("genesets/RS2_genes.csv", header = FALSE)$X
rs3_genes <- read.csv("genesets/RS3_genes.csv", header = FALSE)$X
rs4_genes <- read.csv("genesets/RS4_genes.csv", header = FALSE)$X
rs_model <- readRDS("models/rs_model.rds")

# Calculate rankings for samples
rankings <- rankGenes(expression_matrix)

# Score samples using all gene sets
scores <- list(rs1_genes, rs2_genes, rs3_genes, rs4_genes) %>%
    purrr::map(~simpleScore(rankings, .)) %>%
    setNames(c("RS1", "RS2", "RS3", "RS4"))

# Apply model
predictions <- predict(rs_model, scores)
```

## Citation

If you use this model in your research, please cite the following publication: [TBD]

## Contact

For questions or issues, please open an issue in this repository or contact [mm12865@nyu.edu]
