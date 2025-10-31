###########################################################################
#
#              Multiomics Subtypes and Drug Interactions
#
###########################################################################
# Purpose: Identify potential drug targets based on differential expression
#          in RNA/methylation subtypes using DGiDB
#
# Inputs:  - DESeq2 results for RNA subtypes
#          - Differential methylation results
#          - Composite subtype outlier analyses
#
# Outputs: - Drug-gene interaction tables
#          - Filtered high-confidence drug targets
#          - Visualization of drug interactions by subtype
###########################################################################

experiment <- "muller_multiomic_subtypes_and_drugs"
outdir <- file.path("output", experiment)
dir.create(outdir, showWarnings = FALSE)

# ======================== LIBRARIES ========================
library(tidyverse)
library(glue)
library(rmatt)
library(broom)
library(survival)
library(pROC)
library(purrr)
library(ComplexHeatmap)
library(circlize)
library(vroom)
library(patchwork)

# TODO: Use RANDOM_SEED from config.R instead of hardcoded value
set.seed(420)
theme_set(theme_bw())

source("code/muller_helpers.R")

# ======================== Load DGE Results ========================
dir.create("output/dgidb", showWarnings = FALSE)

rgx <- "deseq_results_comp_(RNAsubtype_RNAtype[1-4]_vs_NOT)"
dge_files <- list.files(
    "output/run4_rmoutliers2_asr_control/deseq",
    pattern = rgx,
    full.names = TRUE,
    recursive = TRUE
)
rna_dges <- map(dge_files, ~ {
    dge <- vroom(.x, col_types = cols())
    dge <- dge %>% select(gene = `...1`, log2FoldChange, pvalue, padj)
    return(dge)
})
names(rna_dges) <- c("RS1", "RS2", "RS3", "RS4")

rgx <- "DMP_methsubtype__methtype[1-3]_v_NOTmethtype[1-3]_DMP_table.csv"
meth_dge_files <- list.files(
    "output/run4_rmoutliers2_asr_control/meth_processing/run3/asr_corrected_DMP",
    pattern = rgx,
    full.names = TRUE,
    recursive = TRUE
)
meth_dges <- map(meth_dge_files, ~ {
    dge <- vroom(.x, col_types = cols())
    dge <- dge %>% select(gene, log2FoldChange = deltaBeta, pvalue, padj)
    return(dge)
})
names(meth_dges) <- c("MS1", "MS2", "MS3")

params <- read.csv("output/muller_composite_subtype_outlier_analysis/params.csv")

cs1_outliers_merged <- read_and_merge_cs_outliers(
    "output/muller_composite_subtype_outlier_analysis/meth_deva_ms3/negative/CS1/CS1_outliers.csv",
    "output/muller_composite_subtype_outlier_analysis/meth_deva_ms3/positive/CS1/CS1_outliers.csv",
    "output/muller_composite_subtype_outlier_analysis/meth_heatmaps/CS1/CS1_genes_fraction.csv",
    cs_name = "CS1"
)

cs2_outliers_merged <- read_and_merge_cs_outliers(
    "output/muller_composite_subtype_outlier_analysis/rna_deva_rs4/negative/CS2/CS2_outliers.csv",
    "output/muller_composite_subtype_outlier_analysis/rna_deva_rs4/positive/CS2/CS2_outliers.csv",
    "output/muller_composite_subtype_outlier_analysis/rna_heatmaps/CS2/CS2_genes.csv",
    cs_name = "CS2"
)

# Convert composite subtype outlier results to DGE format
cs1_dge <- cs1_outliers_merged %>%
    select(gene, log2FoldChange = odds_ratio, pvalue, padj) %>%
    mutate(
        log2FoldChange = log2FoldChange - 1,
    )
cs2_dge <- cs2_outliers_merged %>%
    select(gene, log2FoldChange = odds_ratio, pvalue, padj) %>%
    mutate(
        log2FoldChange = log2FoldChange - 1,
    )
css <- list(
    CS1 = cs1_dge,
    CS2 = cs2_dge
)

dges <- c(rna_dges, meth_dges, css)

# ======================== DGiDb Drug Interactions ========================
dir.create(file.path(outdir, "dgidb"), showWarnings = FALSE)

rs_drugs <- map(dges, ~ {
    probe_for_drug_interactions(
        dge = .x,
        dgidb = NULL,
    )
})
rs_drugs <- bind_rows(rs_drugs, .id = "group")
write_csv(rs_drugs, file.path(outdir, "dgidb", "rs_drugs.csv"))

# Filter for significant results and summarize top targets per subtype
rs_drugs_filtered <- rs_drugs %>%
    filter(padj < 0.05) %>%
    mutate(
        slogp = -log10(padj) * sign(log2FoldChange),
        interaction_score = as.numeric(interaction_score),
    ) %>%
    group_by(group, gene_name) %>%
    summarise(
        slogp = mean(slogp, na.rm = TRUE),
        log2FoldChange = mean(log2FoldChange, na.rm = TRUE),
        pvalue = mean(pvalue, na.rm = TRUE),
        padj = mean(padj, na.rm = TRUE), interaction_score = mean(interaction_score, na.rm = TRUE),
        drugs = paste(unique(drug_name), collapse = ", "),
        interaction_type = paste(unique(interaction_type), collapse = ", "),
        approved = any(approved),
        immunotherapy = any(immunotherapy),
        anti_neoplastic = any(anti_neoplastic),
    ) %>%
    slice_max(order_by = abs(slogp), n = 10)
write_csv(rs_drugs_filtered, file.path(outdir, "dgidb", "rs_drugs_filtered.csv"))

# Filter specifically for composite subtype signature genes
rs_drugs_filtered_cs_genes <- rs_drugs %>%
    filter(gene_name %in% unlist(map(css, ~ .x$gene))) %>%
    filter(padj < 0.05) %>%
    mutate(
        slogp = -log10(padj) * sign(log2FoldChange),
        interaction_score = as.numeric(interaction_score),
    ) %>%
    group_by(group, gene_name) %>%
    summarise(
        slogp = mean(slogp, na.rm = TRUE),
        log2FoldChange = mean(log2FoldChange, na.rm = TRUE),
        pvalue = mean(pvalue, na.rm = TRUE),
        padj = mean(padj, na.rm = TRUE),
        interaction_score = mean(interaction_score, na.rm = TRUE),
        drugs = paste(unique(drug_name), collapse = ", "),
        interaction_type = paste(unique(interaction_type), collapse = ", "),
        approved = any(approved),
        immunotherapy = any(immunotherapy),
        anti_neoplastic = any(anti_neoplastic),
    )
write_csv(rs_drugs_filtered_cs_genes, file.path(outdir, "dgidb", "rs_drugs_filtered_cs_genes.csv"))

# Create visualization of drug interactions by subtype group
create_drug_interaction_plots <- function(rs_drugs_filtered) {
    plot_groups <- list(
        RS = "RNA Subtypes (RS)",
        MS = "Methylation Subtypes (MS)",
        CS = "Composite Subtypes (CS)"
    )

    rs_plots <- map(names(plot_groups), ~ {
        rs_drugs_filtered %>%
            filter(str_starts(group, .x)) %>%
            plot_drug_interactions(
                x = "gene_name",
                y = "group",
                color = "slogp",
                size = "interaction_score"
            ) +
            ggtitle(plot_groups[[.x]])
    })
    names(rs_plots) <- names(plot_groups)

    rs_plot <- rs_plots$RS / rs_plots$MS / rs_plots$CS +
        plot_layout(heights = c(4, 3, 2))

    return(rs_plot)
}

rs_plot <- create_drug_interaction_plots(rs_drugs_filtered)
ggsave(filename = file.path(outdir, "dgidb", "rs_drugs_plot.pdf"), plot = rs_plot, width = 8, height = 8)

rs_plot_cs_genes <- create_drug_interaction_plots(rs_drugs_filtered_cs_genes)
ggsave(filename = file.path(outdir, "dgidb", "rs_drugs_plot_cs_genes.pdf"), plot = rs_plot_cs_genes, width = 8, height = 8)

# ======================== Gene Expression Validation ========================
dir.create(file.path(outdir, "spot_checks"), showWarnings = FALSE)

counts_rna <- vroom("output/run4_rmoutliers2_asr_control/rna_processing/normcounttab.txt", col_types = cols())
colnames(counts_rna)[1] <- "gene"

# TODO: Move this path to config.R
tabulation_table <- read.table(paste0("output/run4_rmoutliers2_asr_control/integrative_analyses/run9-figs/cohort_tabulations/tabulation_table.csv"), sep = ",", row.names = 1, header = TRUE)

tabulation_table <- tabulation_table %>%
    mutate_at(vars(starts_with("C_")), ~ ifelse(.x == 2, 0, .x)) %>%
    mutate(
        ETHNIC = case_when(
            ETHNIC == "99" ~ NA_character_,
            TRUE ~ ETHNIC
        ),
        RACE = case_when(
            RACE == "White" ~ "White",
            RACE == "Black or African American" ~ "Black",
            RACE == "Asian" ~ "Asian",
            RACE == "American Indian or Alaska Native" ~ "Other",
            RACE == "Native Hawaiian or Other Pacific Islander" ~ "Other",
            RACE == "Multiple Races" ~ "Other",
            TRUE ~ NA_character_
        ),
        NO_DIALYSIS = case_when(
            DIALYSIS == "Yes" ~ 0,
            DIALYSIS == "No" ~ 1,
            TRUE ~ NA_real_
        ),
        composite_subtype = interaction(
            gsub("nmf_cluster_", "RS", rna_4cluster),
            gsub("meth_3cluster_", "MS", meth_3cluster),
            drop = TRUE, sep = "_"
        ),
        composites_of_interest = factor(
            composite_subtype,
            levels = c("RS2_MS3", "RS4_MS3"),
            labels = c("CS1", "CS2")
        ),
        cs1_v_rest = case_when(
            composites_of_interest == "CS1" ~ "CS1",
            TRUE ~ "Rest"
        ),
        cs2_v_rest = case_when(
            composites_of_interest == "CS2" ~ "CS2",
            TRUE ~ "Rest"
        )
    )

add_genes_of_interest <- function(counts_rna, tabulation_table, genes_of_interest) {
    counts_long <- counts_rna %>%
        filter(gene %in% genes_of_interest) %>%
        pivot_longer(-gene, names_to = "PATNUM", values_to = "expression")

    merged_data <- tabulation_table %>%
        inner_join(counts_long, by = "PATNUM")

    return(merged_data)
}

genes_of_interest <- c("PDCD1")
merged_data <- add_genes_of_interest(counts_rna, tabulation_table, genes_of_interest)
write_csv(merged_data, file.path(outdir, "spot_checks", "genes_of_interest_expression_by_subtype.csv"))

create_gene_expression_boxplot <- function(merged_data, gene_name, subtype_column, subtype_labels = NULL, subtype_order = NULL) {
    plot_data <- merged_data %>%
        filter(gene == gene_name) %>%
        drop_na(.data[[subtype_column]], expression)

    if (is.null(subtype_labels)) {
        subtype_labels <- unique(plot_data[[subtype_column]])
    }

    if (is.null(subtype_order)) {
        subtype_order <- unique(plot_data[[subtype_column]])
    }

    plot <- plot_data %>%
        ggplot(aes(x = factor(.data[[subtype_column]], levels = subtype_order, labels = subtype_labels), y = expression)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.05) +
        labs(x = str_to_title(str_replace(subtype_column, "_", " ")), y = glue("{gene_name} Expression")) +
        theme_minimal()

    return(plot)
}

pdcd1_plot <- create_gene_expression_boxplot(
    merged_data,
    "PDCD1",
    "rna_4cluster",
    c("RS1", "RS2", "RS3", "RS4")
)
ggsave(filename = file.path(outdir, "spot_checks", "PDCD1_expression_by_RNA_subtype.pdf"), plot = pdcd1_plot, width = 6, height = 4)

pdcd1_meth_plot <- create_gene_expression_boxplot(
    merged_data,
    "PDCD1",
    "meth_3cluster",
    c("MS1", "MS2", "MS3")
)
ggsave(filename = file.path(outdir, "spot_checks", "PDCD1_expression_by_Methylation_subtype.pdf"), plot = pdcd1_meth_plot, width = 6, height = 4)

pdcd1_composite_plot <- create_gene_expression_boxplot(
    merged_data,
    "PDCD1",
    "composite_subtype",
    subtype_labels = c("RS1_MS1", "RS1_MS2", "RS1_MS3", "RS2_MS1", "RS2_MS2", "RS2_MS3", "RS3_MS1", "RS3_MS2", "RS3_MS3", "RS4_MS1", "RS4_MS2", "RS4_MS3"),
    subtype_order = c("RS1_MS1", "RS1_MS2", "RS1_MS3", "RS2_MS1", "RS2_MS2", "RS2_MS3", "RS3_MS1", "RS3_MS2", "RS3_MS3", "RS4_MS1", "RS4_MS2", "RS4_MS3")
) +
    geom_boxplot(outlier.shape = NA, aes(fill = composites_of_interest), alpha = 0.7) +
    theme(legend.position = "none")
ggsave(filename = file.path(outdir, "spot_checks", "PDCD1_expression_by_Composite_subtype.pdf"), plot = pdcd1_composite_plot, width = 10, height = 4)
