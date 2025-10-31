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

# Function to process tabulation table with standardized recoding
#' Process Tabulation Table
#'
#' Standardizes the tabulation table by recoding variables consistently
#' across all analysis scripts
#'
#' @param data Raw tabulation table
#' @return Processed tabulation table with standardized variables
process_tabulation_table <- function(data) {
    data %>%
        # Recode censoring variables (2 -> 0 for consistency)
        mutate_at(vars(starts_with("C_")), ~ ifelse(.x == 2, 0, .x)) %>%
        # Standardize demographic variables
        mutate(
            # Recode missing ethnicity
            ETHNIC = case_when(
                ETHNIC == "99" ~ NA_character_,
                TRUE ~ ETHNIC
            ),
            # Standardize race categories
            RACE = case_when(
                RACE == "White" ~ "White",
                RACE == "Black or African American" ~ "Black",
                RACE == "Asian" ~ "Asian",
                RACE %in% c("American Indian or Alaska Native",
                           "Native Hawaiian or Other Pacific Islander",
                           "Multiple Races") ~ "Other",
                TRUE ~ NA_character_
            ),
            # Create dialysis indicator
            NO_DIALYSIS = case_when(
                DIALYSIS == "Yes" ~ 0,
                DIALYSIS == "No" ~ 1,
                TRUE ~ NA_real_
            ),
            # Create composite subtypes if RNA and methylation clusters exist
            composite_subtype = if ("rna_4cluster" %in% colnames(.) &
                                   "meth_3cluster" %in% colnames(.)) {
                interaction(
                    gsub("nmf_cluster_", "RS", rna_4cluster),
                    gsub("meth_3cluster_", "MS", meth_3cluster),
                    drop = TRUE, sep = "_"
                )
            } else {
                NA_character_
            },
            # Create composites of interest
            composites_of_interest = if (!is.na(composite_subtype[1])) {
                factor(
                    composite_subtype,
                    levels = c("RS2_MS3", "RS4_MS3"),
                    labels = c("CS1", "CS2")
                )
            } else {
                NA_character_
            },
            # Create one-vs-rest comparisons
            cs1_v_rest = case_when(
                composites_of_interest == "CS1" ~ "CS1",
                TRUE ~ "Rest"
            ),
            cs2_v_rest = case_when(
                composites_of_interest == "CS2" ~ "CS2",
                TRUE ~ "Rest"
            )
        )
}

# Function to safely load data files with error checking
#' Load Data File Safely
#'
#' Loads a data file with existence and format checking
#'
#' @param file_path Path to the data file
#' @param file_type Type of file ("csv", "tsv", "txt")
#' @param required_columns Optional vector of required column names
#' @return Data frame or matrix
load_data_safe <- function(file_path, file_type = "csv", required_columns = NULL) {
    # Check file exists
    if (!file.exists(file_path)) {
        stop("Data file not found: ", file_path,
             "\nPlease check the file path in config.R or your data directory")
    }

    # Load based on file type
    data <- switch(file_type,
        "csv" = read.csv(file_path, row.names = 1, stringsAsFactors = FALSE),
        "tsv" = read.table(file_path, sep = "\t", header = TRUE, row.names = 1,
                          stringsAsFactors = FALSE),
        "txt" = read.table(file_path, sep = "\t", header = TRUE, row.names = 1,
                          stringsAsFactors = FALSE),
        stop("Unsupported file type: ", file_type)
    )

    # Check for required columns
    if (!is.null(required_columns)) {
        missing_cols <- setdiff(required_columns, colnames(data))
        if (length(missing_cols) > 0) {
            stop("Missing required columns in ", basename(file_path), ": ",
                 paste(missing_cols, collapse = ", "))
        }
    }

    log_analysis_step(paste("Loaded data from", basename(file_path),
                           "- Dimensions:", nrow(data), "x", ncol(data)))

    return(data)
}

# Function to calculate odds ratios for one-vs-rest comparisons
#' Calculate One-vs-Rest Odds Ratios
#'
#' Performs logistic regression to calculate odds ratios for each level
#' of a grouping variable compared to all other levels
#'
#' @param data Data frame containing the variables
#' @param group_var Character string naming the grouping variable
#' @param test_vars Character vector of test variable names
#' @param control_vars Character vector of control variable names (default: demographics)
#' @return Data frame with odds ratios, confidence intervals, and p-values
calculate_odds_ratios_one_vs_rest <- function(data, group_var, test_vars,
                                             control_vars = c("AGE_RAND", "SEX", "RACE", "ETHNIC")) {
    # Drop missing data
    d <- drop_na(data, all_of(c(group_var, test_vars, control_vars)))

    # One versus rest approach for each level of the group variable
    unique_levels <- unique(d[[group_var]])

    one_vs_rest_results <- map(unique_levels, ~ {
        level <- .x
        # Create binary variable for current level vs rest
        d_ovr <- d %>%
            mutate(!!glue("{group_var}_ovr") := ifelse(!!sym(group_var) == level, 1, 0))

        ovr_results <- map(test_vars, ~ {
            test_var <- .x
            # Create a logistic regression model for each test variable
            controls_formula <- paste(control_vars, collapse = " + ")
            fmla <- as.formula(glue("{group_var}_ovr ~ {test_var} + {controls_formula}"))

            tryCatch({
                model <- glm(fmla, data = d_ovr, family = binomial)
                model_tidy <- tidy(model)
                model_tidy <- model_tidy %>%
                    filter(term == test_var) %>%
                    mutate(
                        variable = test_var,
                        odds_ratio = exp(estimate),
                        lower_ci = exp(estimate - 1.96 * std.error),
                        upper_ci = exp(estimate + 1.96 * std.error),
                        level = level
                    ) %>%
                    select(variable, level, odds_ratio, lower_ci, upper_ci, p.value)
                return(model_tidy)
            }, error = function(e) {
                warning(glue("Model failed for {test_var} in {group_var} = {level}: {e$message}"))
                return(NULL)
            })
        })
        ovr_results <- bind_rows(ovr_results)
        return(ovr_results)
    })

    one_vs_rest_results <- bind_rows(one_vs_rest_results)
    one_vs_rest_results$group <- group_var

    return(one_vs_rest_results)
}

# Set default ggplot theme
theme_set(get_publication_theme())
