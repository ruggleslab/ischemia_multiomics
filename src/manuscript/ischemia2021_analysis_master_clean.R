###########################################################################
#
#                            ISCHEMIA RNA-seq Analysis Master Script
#
###########################################################################
# Author: ISCHEMIA Study Team
# Date: 2025-08-08
# Script Name: ischemia2021_analysis_master_clean.R
# Description: Master script for ISCHEMIA RNA-seq differential expression analysis
#
# This script performs differential expression analysis for:
# - Inducible ischemia severity (DEGRISCH, IMGDEGIS)
# - Coronary artery disease severity (CTMULT50, CTANYD50, DUKESCORE)
#
# Analysis groups:
# Inducible ischemia severity:
# 1. DEGRISCH (all testing modalities): 0-1 = none-mild; 2 = moderate; 3 = severe
# 2. IMGDEGIS (ischemia by imaging): 0-1 = none-mild; 2 = moderate; 3 = severe
#
# Coronary artery disease severity:
# 1. CTMULT50 (multivessel disease by CTA): 0,1, 2 non-evaluable
# 2. CTANYD50 (any obstructive disease): 0, 1, 2 non-evaluable
# 3. DUKESCORE: groups - 1-2, 3, 4-5, >6

# ======================== SETUP ========================
# Load configuration and utilities
source("config.R")
source("utils.R")

# Create output directories
analysis_dir <- setup_analysis_dirs("rna_differential_expression")
process_dir <- file.path(analysis_dir, "processed_data")
log_file <- file.path(analysis_dir, "logs", "analysis_log.txt")

log_analysis_step("Starting ISCHEMIA RNA-seq analysis", log_file)

# ======================== DATA LOADING ========================
log_analysis_step("Loading count data", log_file)

# Load count data (update COUNT_FILE path in config.R)
if (!file.exists(COUNT_FILE)) {
    stop("Count file not found. Please update COUNT_FILE path in config.R")
}

count_table <- read.table(COUNT_FILE, header = TRUE, row.names = 1, 
                         sep = "\t", check.names = FALSE, comment = "")

log_analysis_step(paste("Loaded count data:", nrow(count_table), "genes,", 
                       ncol(count_table), "samples"), log_file)

# ======================== METADATA PROCESSING ========================
log_analysis_step("Loading and processing metadata", log_file)

# Load metadata (update METADATA_FILE path in config.R)
if (!file.exists(METADATA_FILE)) {
    stop("Metadata file not found. Please update METADATA_FILE path in config.R")
}

metadata <- read.table(METADATA_FILE, sep = ",", header = TRUE, 
                      stringsAsFactors = FALSE, na.strings = c(NA, "NA", ""))

# Standardize sample IDs
metadata[,1] <- as.character(metadata[,1])
colnames(metadata)[1] <- "sample"

# Find overlapping samples between metadata and count data
sample_ids <- sort(intersect(metadata[,1], colnames(count_table)))
sample_ids_ordered <- sample_ids[match(metadata[,1][metadata[,1] %in% colnames(count_table)], sample_ids)]

# Check for samples only in one dataset
samples_only_meta <- setdiff(metadata[,1], colnames(count_table))
samples_only_counts <- setdiff(colnames(count_table), metadata[,1])

if (length(samples_only_meta) > 0) {
    log_analysis_step(paste("Warning:", length(samples_only_meta), 
                           "samples in metadata not in count table"), log_file)
}
if (length(samples_only_counts) > 0) {
    log_analysis_step(paste("Warning:", length(samples_only_counts), 
                           "samples in count table not in metadata"), log_file)
}

# Filter to overlapping samples
metadata_filtered <- metadata[match(sample_ids_ordered, metadata[,1]),]
count_table_filtered <- count_table[,sample_ids_ordered]

# Verify alignment
if (!identical(colnames(count_table_filtered), metadata_filtered[,1])) {
    stop("Sample ordering mismatch between count data and metadata")
}

log_analysis_step(paste("Final dataset:", nrow(metadata_filtered), "samples"), log_file)

# ======================== QUALITY CONTROL ========================
log_analysis_step("Performing quality control", log_file)

# Remove problematic genes (hemoglobin genes)
remove_genes <- c("HBB", "HBA1", "HBA2")
hemoglobin_counts <- t(count_table_filtered[rownames(count_table_filtered) %in% remove_genes,])

if (nrow(hemoglobin_counts) > 0) {
    hemoglobin_proportions <- cbind.data.frame(
        sample = rownames(hemoglobin_counts),
        hemoglobin_counts / colSums(count_table_filtered)
    )
    
    # Save hemoglobin proportions for QC
    write.csv(hemoglobin_proportions, 
             file.path(process_dir, "hemoglobin_proportions.csv"),
             row.names = FALSE)
    
    log_analysis_step(paste("Removed", length(remove_genes), "hemoglobin genes"), log_file)
}

# Filter count table
count_table_clean <- count_table_filtered[!rownames(count_table_filtered) %in% remove_genes,]

# Perform expression-based filtering
count_table_qc <- qc_expression_data(count_table_clean, min_counts = 10, min_samples = 0.1)

# ======================== OUTLIER DETECTION ========================
log_analysis_step("Detecting potential outliers", log_file)

# Calculate basic QC metrics
qc_metrics <- data.frame(
    sample = colnames(count_table_qc),
    total_counts = colSums(count_table_qc),
    detected_genes = colSums(count_table_qc > 0),
    stringsAsFactors = FALSE
)

# Detect outliers based on total counts
count_outliers <- detect_outliers(qc_metrics$total_counts)
gene_outliers <- detect_outliers(qc_metrics$detected_genes)

if (length(count_outliers) > 0) {
    log_analysis_step(paste("Potential count outliers:", 
                           paste(qc_metrics$sample[count_outliers], collapse = ", ")), log_file)
}

if (length(gene_outliers) > 0) {
    log_analysis_step(paste("Potential gene detection outliers:", 
                           paste(qc_metrics$sample[gene_outliers], collapse = ", ")), log_file)
}

# Save QC metrics
write.csv(qc_metrics, file.path(process_dir, "qc_metrics.csv"), row.names = FALSE)

# ======================== SAVE PROCESSED DATA ========================
log_analysis_step("Saving processed data", log_file)

# Save cleaned count matrix
write.table(count_table_qc, 
           file = file.path(process_dir, "count_table_cleaned.txt"),
           sep = "\t", quote = FALSE, col.names = NA)

# Save filtered metadata
write.csv(metadata_filtered, 
         file = file.path(process_dir, "metadata_filtered.csv"),
         row.names = FALSE)

# Save session info for reproducibility
writeLines(capture.output(sessionInfo()), 
          file.path(analysis_dir, "logs", "session_info.txt"))

log_analysis_step("Data preprocessing completed", log_file)

# ======================== DESEQ2 ANALYSIS SETUP ========================
log_analysis_step("Setting up DESeq2 analysis", log_file)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
    countData = count_table_qc,
    colData = metadata_filtered,
    design = ~ 1  # Will be updated for specific comparisons
)

# Perform variance stabilizing transformation for visualization
vst_data <- vst(dds, blind = TRUE)

# Save transformed data
write.table(assay(vst_data), 
           file = file.path(process_dir, "vst_transformed_data.txt"),
           sep = "\t", quote = FALSE, col.names = NA)

# ======================== DIFFERENTIAL EXPRESSION COMPARISONS ========================
# This section would contain the specific DE comparisons
# The original file continues with specific analyses that would need 
# individual cleaning based on the experimental design

log_analysis_step("Master analysis script completed. Run specific comparison scripts for DE analysis.", log_file)

# Print summary
cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Input samples:", ncol(count_table), "\n")
cat("Input genes:", nrow(count_table), "\n")
cat("Final samples:", ncol(count_table_qc), "\n")
cat("Final genes:", nrow(count_table_qc), "\n")
cat("Potential outliers detected:", length(unique(c(count_outliers, gene_outliers))), "\n")
cat("Output directory:", analysis_dir, "\n")
cat("========================\n")
