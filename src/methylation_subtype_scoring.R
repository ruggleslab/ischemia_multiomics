# Author: Matthew Muller
# Date: 2025-07-21

# ======================== USER CONFIGURATION ========================

# Path to your new data file (samples = rows, features = columns)
new_data_file <- "path/to/your/new_data.csv"

# Output directory for scoring results
scoring_outdir <- "output/methylation_scoring_results"
dir.create(scoring_outdir, recursive = TRUE, showWarnings = FALSE)

# Load the model and feature names
feature_names_file <- "path/to/your/feature_names.csv"
feature_names <- read_csv(feature_names_file)$feature

# Load the model to use for predictions
model_file <- "path/to/your/trained_model.rds"
trained_model <- readRDS(model_file)

# ======================== SETUP ========================
# Load required libraries
library(tidyverse)
library(data.table)
library(caret)

# ======================== LOAD AND PREPARE NEW DATA ========================
# Read the new data
new_data <- fread(new_data_file)

# Check if all required features are present
missing_features <- setdiff(feature_names, colnames(new_data))
if (length(missing_features) > 0) {
    stop("ERROR: The following required features are missing from the new data: ", paste(missing_features, collapse = ", "))
}

# Normalize each feature using the inverse normal transformation (vdw)
normalize_feature <- function(x) {
  qnorm(rank(x, na.last = "keep") / (sum(!is.na(x)) + 1))
}
new_data[, (feature_names) := lapply(.SD, normalize_feature), .SDcols = feature_names]

# Select and order features to match training data
new_data_processed <- new_data[, ..feature_names]

# ======================== GENERATE PREDICTIONS ========================
# Generate class predictions
pred_class <- predict(trained_model, new_data_processed)
pred_prob <- predict(trained_model, new_data_processed, type = "prob")

# Create results dataframe
results_df <- data.frame(
  sample_id = rownames(new_data_processed),
  prediction = as.character(pred_class),
  stringsAsFactors = FALSE
)

# Add probabilities for each class
for (class_name in colnames(pred_prob)) {
  results_df[[paste0("prob_", class_name)]] <- pred_prob[[class_name]]
}

# ======================== SAVE RESULTS ========================
# Save main results file
results_file <- file.path(scoring_outdir, "methylation_subtype_predictions.csv")
write_csv(results_df, results_file)

# Save summary statistics
table_summary <- table(results_df$prediction)
writeLines(
  paste("Class distribution in predictions:\n", 
        paste(names(table_summary), table_summary, sep = ": ", collapse = "\n")),
  file = file.path(scoring_outdir, "class_distribution.txt")
)
