###########################################################################
#
#                            Configuration File
#
###########################################################################
# Author: Matthew Muller
# Date: 2025-08-08
# Script Name: config.R
# Description: Configuration file for ISCHEMIA multiomics analysis
#
# This file contains all the configurable paths and parameters used across
# the analysis scripts. This allows for easy modification without changing
# individual scripts.

# ======================== PROJECT PATHS ========================
# Base project directory
PROJECT_ROOT <- here::here()

# Data directories
DATA_DIR <- file.path(PROJECT_ROOT, "data")
GENESETS_DIR <- file.path(PROJECT_ROOT, "genesets")
MODELS_DIR <- file.path(PROJECT_ROOT, "models")

# Output directories
OUTPUT_DIR <- file.path(PROJECT_ROOT, "output")
FIGURES_DIR <- file.path(OUTPUT_DIR, "figures")
TABLES_DIR <- file.path(OUTPUT_DIR, "tables")

# Ensure output directories exist
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)

# ======================== ANALYSIS PARAMETERS ========================
# Random seed for reproducibility
# Use this seed in all scripts to ensure consistent results
RANDOM_SEED <- 42

# QC parameters for expression data
MIN_READ_COUNT <- 2000000  # Minimum reads per sample
MIN_GENE_COUNT <- 16       # Minimum counts for gene inclusion
MIN_SAMPLES_FRACTION <- 0.5 # Fraction of samples gene must be expressed in

# Analysis parameters
WGCNA_POWER <- 16
WGCNA_MIN_MODULE_SIZE <- 40
NMF_RANK <- 4

# Variance gene selection percentiles (for PCA and heatmaps)
UPPER_PERCENTILE <- 0.9999
LOWER_PERCENTILE <- 0.95

# ======================== FILE PATHS ========================
# These paths should be updated to point to your actual data files
# or set as environment variables

# Count data file
COUNT_FILE <- file.path(DATA_DIR, "quant.featurecounts.counts.rev.txt")

# Metadata file
METADATA_FILE <- file.path(DATA_DIR, "metasheet_8_20220418.csv")

# Biomarker data
BIOMARKER_FILE <- file.path(DATA_DIR, "biomarkers_20250906", "20241025_biomarker_data_cleaned.csv")

# Color guide file
COLOR_GUIDE_FILE <- file.path(DATA_DIR, "ischemia2021_colorguide.csv")

# ======================== HELPER FUNCTIONS ========================
# Function to create output subdirectories
create_output_dir <- function(subdir) {
    full_path <- file.path(OUTPUT_DIR, subdir)
    dir.create(full_path, recursive = TRUE, showWarnings = FALSE)
    return(full_path)
}

# Function to generate timestamped filenames
timestamp_filename <- function(base_name, extension = ".csv") {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    paste0(base_name, "_", timestamp, extension)
}

# Set seed for reproducibility
set.seed(RANDOM_SEED)
