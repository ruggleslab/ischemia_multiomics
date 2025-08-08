###########################################################################
#
#                            muller_multiomic_subtypes_and_drugs
#
###########################################################################
# Author: Matthew Muller
# Date: 2024-09-22
# Script Name: muller_multiomic_subtypes_and_drugs
# Output directory:
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

set.seed(420)
theme_set(theme_bw())

source("code/muller_helpers.R")

# ======================== DGiDb ========================
dir.create("output/dgidb", showWarnings = FALSE)

# get the dges we want to probe
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

# get the methylation dges
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

# now get the composite subtype outlier analyses
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

# Combine the composite subtype outlier analyses with the DGE lists
cs1_dge <- cs1_outliers_merged %>%
    select(gene, log2FoldChange = odds_ratio, pvalue, padj) %>%
    mutate(
        log2FoldChange = log2FoldChange - 1,  # Adjusting the odds ratio to log2 scale
    )
cs2_dge <- cs2_outliers_merged %>%
    select(gene, log2FoldChange = odds_ratio, pvalue, padj) %>%
    mutate(
        log2FoldChange = log2FoldChange - 1,  # Adjusting the odds ratio to log2 scale
    )
css <- list(
    CS1 = cs1_dge,
    CS2 = cs2_dge
)

# Combine RNA and methylation DGE lists
dges <- c(rna_dges, meth_dges, css)

#======================== DGiDb ========================
dir.create(file.path(outdir, "dgidb"), showWarnings = FALSE)

# Probe for drug interactions
rs_drugs <- map(dges, ~ {
    probe_for_drug_interactions(
        dge = .x,
        dgidb = NULL,
    )
})
rs_drugs <- bind_rows(rs_drugs, .id = "group")
write_csv(rs_drugs, file.path(outdir, "dgidb", "rs_drugs.csv"))

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
        padj = mean(padj, na.rm = TRUE),        interaction_score = mean(interaction_score, na.rm = TRUE),
        drugs = paste(unique(drug_name), collapse = ", "),
        interaction_type = paste(unique(interaction_type), collapse = ", "),
        approved = any(approved),
        immunotherapy = any(immunotherapy),
        anti_neoplastic = any(anti_neoplastic),
    ) %>%
    slice_max(order_by = abs(slogp), n = 10)
write_csv(rs_drugs_filtered, file.path(outdir, "dgidb", "rs_drugs_filtered.csv"))

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
    
# Create separate plots for MS, RS, and CS using map
plot_groups <- list(
    RS = "RNA Subtypes (RS)",
    MS = "Methylation Subtypes (MS)",
    CS = "Composite Subtypes (CS)"
)
create_drug_interaction_plots <- function(rs_drugs_filtered) {
    # Create separate plots for MS, RS, and CS using map
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
    
    # Combine plots
    rs_plot <- rs_plots$RS / rs_plots$MS / rs_plots$CS +
        plot_layout(heights = c(4, 3, 2))
    
    return(rs_plot)
}

# Use the function
rs_plot <- create_drug_interaction_plots(rs_drugs_filtered)
ggsave(filename = file.path(outdir, "dgidb", "rs_drugs_plot.pdf"), plot = rs_plot, width = 8, height = 8)

# Use on the filtered CS genes
rs_plot_cs_genes <- create_drug_interaction_plots(rs_drugs_filtered_cs_genes)
ggsave(filename = file.path(outdir, "dgidb", "rs_drugs_plot_cs_genes.pdf"), plot = rs_plot_cs_genes, width = 8, height = 8)
