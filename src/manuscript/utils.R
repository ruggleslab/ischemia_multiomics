###########################################################################
#
#                            Utility Functions
#
###########################################################################
# Author: Matthew Muller
# Date: 2025-08-08
# Script Name: utils.R
# Description: Common utility functions used across ISCHEMIA analysis scripts

# ======================== SETUP ========================
# Load required packages
required_packages <- c(
    "tidyverse", "ggplot2", "reshape2", "DESeq2", "grid", "gridExtra", 
    "scales", "ggrepel", "tools", "Hmisc", "NMF", "caret", "blacksheepr", 
    "mlbench", "pROC", "WGCNA", "multiROC", "scattermore", "dendextend",
    "psych", "glue", "rmatt", "broom", "survival", "purrr", 
    "ComplexHeatmap", "circlize", "here"
)

# Function to install and load packages
load_packages <- function(packages) {
    for (pkg in packages) {
        if (!require(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)) {
            message(paste("Installing package:", pkg))
            install.packages(pkg)
            require(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
        }
    }
}

# Load all required packages
load_packages(required_packages)

# ======================== COMMON FUNCTIONS ========================

# Function to standardize ggplot theme
get_publication_theme <- function() {
    theme_bw() +
    theme(
        text = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        strip.text = element_text(size = 11),
        plot.title = element_text(size = 14, hjust = 0.5)
    )
}

# Function to save plots with consistent formatting
save_publication_plot <- function(plot, filename, width = 8, height = 6, dpi = 300) {
    ggsave(
        filename = filename,
        plot = plot,
        width = width,
        height = height,
        dpi = dpi,
        units = "in"
    )
}

# Function to create standardized output directory structure
setup_analysis_dirs <- function(analysis_name) {
    base_dir <- file.path("output", analysis_name)
    dirs <- c("figures", "tables", "processed_data", "logs")
    
    for (dir in dirs) {
        dir.create(file.path(base_dir, dir), recursive = TRUE, showWarnings = FALSE)
    }
    
    return(base_dir)
}

# Function to log analysis steps
log_analysis_step <- function(message, log_file = NULL) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    log_message <- paste(timestamp, "-", message)
    
    # Print to console
    message(log_message)
    
    # Write to log file if specified
    if (!is.null(log_file)) {
        write(log_message, file = log_file, append = TRUE)
    }
}

# Function to create summary statistics table
create_summary_stats <- function(data, group_var, numeric_vars) {
    summary_stats <- data %>%
        group_by(!!sym(group_var)) %>%
        summarise(
            across(all_of(numeric_vars), 
                   list(
                       n = ~sum(!is.na(.)),
                       mean = ~mean(., na.rm = TRUE),
                       sd = ~sd(., na.rm = TRUE),
                       median = ~median(., na.rm = TRUE),
                       q25 = ~quantile(., 0.25, na.rm = TRUE),
                       q75 = ~quantile(., 0.75, na.rm = TRUE)
                   ),
                   .names = "{.col}_{.fn}"),
            .groups = "drop"
        )
    
    return(summary_stats)
}

# Function to perform quality control on expression data
qc_expression_data <- function(count_matrix, min_counts = 10, min_samples = 0.1) {
    # Calculate the minimum number of samples
    min_sample_count <- ceiling(ncol(count_matrix) * min_samples)
    
    # Filter genes with low expression
    keep_genes <- rowSums(count_matrix >= min_counts) >= min_sample_count
    
    filtered_matrix <- count_matrix[keep_genes, ]
    
    log_analysis_step(paste(
        "QC filtering: Retained", nrow(filtered_matrix), "of", nrow(count_matrix), 
        "genes (",
        round(nrow(filtered_matrix)/nrow(count_matrix)*100, 1), "%)"
    ))
    
    return(filtered_matrix)
}

# Function to detect outliers using IQR method
detect_outliers <- function(x, multiplier = 1.5) {
    q25 <- quantile(x, 0.25, na.rm = TRUE)
    q75 <- quantile(x, 0.75, na.rm = TRUE)
    iqr <- q75 - q25
    
    lower_bound <- q25 - multiplier * iqr
    upper_bound <- q75 + multiplier * iqr
    
    outliers <- which(x < lower_bound | x > upper_bound)
    return(outliers)
}

# Function to create volcano plot
create_volcano_plot <- function(results_df, title = "Volcano Plot", 
                               p_threshold = 0.05, fc_threshold = 1) {
    
    # Prepare data
    plot_data <- results_df %>%
        mutate(
            significant = ifelse(padj < p_threshold & abs(log2FoldChange) > fc_threshold, 
                                "Significant", "Not Significant"),
            log10_padj = -log10(padj)
        )
    
    # Create plot
    volcano_plot <- ggplot(plot_data, aes(x = log2FoldChange, y = log10_padj, 
                                         color = significant)) +
        geom_point(alpha = 0.6) +
        geom_vline(xintercept = c(-fc_threshold, fc_threshold), 
                   linetype = "dashed", color = "gray") +
        geom_hline(yintercept = -log10(p_threshold), 
                   linetype = "dashed", color = "gray") +
        scale_color_manual(values = c("Not Significant" = "gray", 
                                     "Significant" = "red")) +
        labs(x = "log2 Fold Change", y = "-log10 Adjusted P-value", 
             title = title) +
        get_publication_theme() +
        theme(legend.title = element_blank())
    
    return(volcano_plot)
}

# Function to perform pathway enrichment analysis placeholder
# (This would need actual pathway databases to implement fully)
perform_pathway_analysis <- function(gene_list, background_genes = NULL) {
    # Placeholder function - would integrate with pathway databases
    # like KEGG, GO, Reactome, etc.
    
    log_analysis_step("Pathway analysis function called - implement with actual pathway database")
    
    return(data.frame(
        pathway = "Placeholder pathway",
        p_value = 0.05,
        genes = paste(head(gene_list, 5), collapse = ", ")
    ))
}

# Set default ggplot theme
theme_set(get_publication_theme())
