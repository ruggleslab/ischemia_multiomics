# Author: Matthew Muller
# Date: 2025-07-21

# ======================== USER CONFIGURATION ========================
# Path to the output directory where results will be saved
outdir <- "path/to/output"

# Path to expression matrix file (CSV format)
# Expression Matrix format: samples = columns, genes = rows
# First column should contain gene symbols, remaining columns should be sample data
expression_matrix_file <- "path/to/expression_matrix.csv"

# ======================== LIBRARIES ========================
library(singscore)
library(vroom)
library(tidyverse)

# Function to calculate RS scores and predictions
calculate_rs_scores <- function(expression_matrix, rs1_genes, rs2_genes, rs3_genes, rs4_genes, rs_model) {
    # Calculate rankings for samples
    rankings <- rankGenes(expression_matrix)

    # Score samples using all gene sets
    scores <- list(rs1_genes, rs2_genes, rs3_genes, rs4_genes) %>%
        purrr::map(~ simpleScore(rankings, .)$TotalScore) %>%
        # scale the scores
        purrr::map(~ scale(.)) %>% # There may be batch variation, so we scale the scores in most cases
        setNames(c("UpScore__rna_subtype1", "UpScore__rna_subtype2", "UpScore__rna_subtype3", "UpScore__rna_subtype4")) %>%
        bind_rows()

    # Apply model
    predictions <- predict(rs_model, scores)
    predictions_df <- data.frame(
        sample_id = colnames(expression_matrix),
        RS = predictions,
        row.names = colnames(expression_matrix)
    )

    return(predictions_df)
}

# ======================== CODE ========================
# Create output directory if it doesn't exist
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Load expression matrix
cat("Loading expression matrix from:", expression_matrix_file, "\n")
expression_data <- vroom(expression_matrix_file) %>% column_to_rownames(names(.)[1])

# Download gene sets from GitHub repository
cat("Downloading gene sets and model...\n")
base_url <- "https://raw.githubusercontent.com/ruggleslab/ischemia_multiomics/main"
rs1_genes <- read.csv(url(file.path(base_url, "genesets/RS1_genes.csv")), header = FALSE)$V1
rs2_genes <- read.csv(url(file.path(base_url, "genesets/RS2_genes.csv")), header = FALSE)$V1
rs3_genes <- read.csv(url(file.path(base_url, "genesets/RS3_genes.csv")), header = FALSE)$V1
rs4_genes <- read.csv(url(file.path(base_url, "genesets/RS4_genes.csv")), header = FALSE)$V1
rs_model <- readRDS(url(file.path(base_url, "models/rs_model.rds")))

# Score the samples
cat("Calculating RNA subtype scores...\n")
predictions <- calculate_rs_scores(
    expression_data,
    rs1_genes,
    rs2_genes,
    rs3_genes,
    rs4_genes,
    rs_model
)

# Print and save RS distribution
tab <- with(predictions, table(RS))
cat("RNA Subtype (RS) distribution:\n")
print(tab)
writeLines(paste0("RS distribution:\n", paste(capture.output(print(tab)), collapse = "\n")), file.path(outdir, "rs_distribution.txt"))

# Save predictions to output directory
cat("Saving results to:", outdir, "\n")
write_csv(predictions, file.path(outdir, "predictions.csv"))

# ======================== END ========================
writeLines(capture.output(sessionInfo()), file.path(outdir, "session.log"))
save.image(file.path(outdir, "session.RData"))
