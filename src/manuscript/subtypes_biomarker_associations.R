###########################################################################
#
#           Multiomics Subtypes and Biomarker Associations
#
###########################################################################
# Purpose: Analyze associations between RNA/methylation subtypes and
#          circulating biomarkers, compare risk prediction models
#
# Inputs:  - Cohort tabulation table with subtype assignments
#          - Biomarker measurements
#          - Outcome data
#
# Outputs: - Biomarker summary statistics by subtype
#          - Odds ratio analyses
#          - ROC curve comparisons
###########################################################################

# Setup
experiment <- "muller_multiomic_subtypes_and_biomarkers"
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

# TODO: Use RANDOM_SEED from config.R instead of hardcoded value
set.seed(420)
theme_set(theme_bw())

# ======================== LOAD DATA ========================
# TODO: Move these paths to config.R
colorguidefile <- "data/ischemia2021_colorguide.csv"
colorguide <- read.table(colorguidefile, sep = ",", header = TRUE, comment.char = "",
                        colClasses = c("character", "character", "character"))

tabulation_table <- read.table("output/run4_rmoutliers2_asr_control/integrative_analyses/run9-figs/cohort_tabulations/tabulation_table.csv",
                               sep = ",", row.names = 1, header = TRUE)

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
            RACE %in% c("American Indian or Alaska Native",
                       "Native Hawaiian or Other Pacific Islander",
                       "Multiple Races") ~ "Other",
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

controls_list <- list(
    demographics = c("AGE_RAND", "SEX", "RACE", "ETHNIC"),
    risks = c("AGE_RAND", "SEX", "RACE", "ETHNIC", "DIABETES", "EGFR", "NO_DIALYSIS")
)

group_columns_for_testing <- c("rna_4cluster", "meth_3cluster", "composite_subtype", "composites_of_interest", "cs1_v_rest", "cs2_v_rest")
variables_to_test <- grep("BMI|SMOK|CKD|HYP|DIAB|DIA|SYS|LVEF|HEMOGLOB|HDLC|LDLC|TRIGLYC|TOTCHOL|EGFR|DEGISCH|CTMULT70|CTMULT50|CTNDV70|CTNDV50|IMGDEGIS|CUSTOM_IMGDEGIS_NONEMILD|IMGDEGIS_01_2_3", colnames(tabulation_table), value = TRUE)
censors <- c("C_PRIMARY", "C_CVDMI_P", "C_ACDMI_P", "C_MIPRIM", "C_CVD", "C_ACD")

biomarkers <- read.csv("data/biomarkers_20250906/20241025_biomarker_data_cleaned.csv", row.names = 1)
biomarker_variables <- colnames(biomarkers)[-1]

dat <- left_join(tabulation_table, biomarkers)

# ======================== Biomarkers ========================
dir.create(file.path(outdir, "biomarkers_stats"), showWarnings = FALSE)

biomarker_tabulations <- map(group_columns_for_testing, ~ {
    d <- drop_na(dat, all_of(c(.x, biomarker_variables)))
    t <- stats_table(d, .x, biomarker_variables)
    return(t)
})
names(biomarker_tabulations) <- group_columns_for_testing

map2(biomarker_tabulations, ~ write.csv(.x, file = file.path(outdir, "biomarkers_stats", paste0(.y, "_biomarker_tabulation.csv"))), .y = names(biomarker_tabulations))

# ======================== Heatmaps ========================
biomarker_heatmap_data <- dat %>%
    select(all_of(biomarker_variables), all_of(group_columns_for_testing[1:3])) %>%
    drop_na() %>%
    arrange(composite_subtype)

biomarker_heatmap_annotation <- HeatmapAnnotation(
    df = biomarker_heatmap_data %>% select(all_of(group_columns_for_testing[1:3])),
    col = color_mapping(biomarker_heatmap_data %>% select(all_of(group_columns_for_testing[1:3]))),
    annotation_name_side = "left"
)

biomarker_heatmap <- Heatmap(
    t(scale(biomarker_heatmap_data %>% select(all_of(biomarker_variables)))),
    col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
    name = "Biomarkers",
    top_annotation = biomarker_heatmap_annotation,
    column_split = biomarker_heatmap_data$composite_subtype,
    cluster_column_slices = FALSE,
    column_title = "Biomarkers by Subtype",
    cluster_rows = TRUE,
    show_row_dend = FALSE,
    show_column_names = FALSE,
    row_names_side = "left",
    width = unit(16, "cm"),
    height = unit(6, "cm"),
)

pdf(file = file.path(outdir, "biomarker_heatmap.pdf"), width = 12, height = 10)
draw(biomarker_heatmap, heatmap_legend_side = "right")
dev.off()

# Calculate average biomarker expression and ranks by subtype
rna_biomarkers_avg_expression <- dat %>%
    select(all_of(biomarker_variables), rna_4cluster) %>%
    drop_na() %>%
    group_by(rna_4cluster) %>%
    summarise(across(all_of(biomarker_variables), mean, na.rm = TRUE)) %>%
    column_to_rownames("rna_4cluster")
write.csv(rna_biomarkers_avg_expression, file = file.path(outdir, "rna_biomarkers_avg_expression.csv"))

rna_biomarkers_ranks <- rna_biomarkers_avg_expression %>%
    mutate(across(all_of(biomarker_variables), ~ rank(-.x, ties.method = "min")))
write.csv(rna_biomarkers_ranks, file = file.path(outdir, "rna_biomarkers_ranks.csv"))

meth_biomarkers_avg_expression <- dat %>%
    select(all_of(biomarker_variables), meth_3cluster) %>%
    drop_na() %>%
    group_by(meth_3cluster) %>%
    summarise(across(all_of(biomarker_variables), mean, na.rm = TRUE)) %>%
    column_to_rownames("meth_3cluster")
write.csv(meth_biomarkers_avg_expression, file = file.path(outdir, "meth_biomarkers_avg_expression.csv"))

meth_biomarkers_ranks <- meth_biomarkers_avg_expression %>%
    mutate(across(all_of(biomarker_variables), ~ rank(-.x, ties.method = "min")))
write.csv(meth_biomarkers_ranks, file = file.path(outdir, "meth_biomarkers_ranks.csv"))

heatmap_ranks_rna <- Heatmap(
    t(rna_biomarkers_ranks),
    col = colorRamp2(1:4, c("white", "grey70", "grey30", "black")),
    name = "Biomarker Ranks",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    width = unit(2, "cm"),
    height = unit(6, "cm"),
    border = TRUE,
    rect_gp = gpar(col = "black", lwd = 1),
)
pdf(file = file.path(outdir, "biomarker_ranks_heatmap_rna.pdf"), width = 12, height = 10)
draw(heatmap_ranks_rna, heatmap_legend_side = "right")
dev.off()

heatmap_ranks_meth <- Heatmap(
    t(meth_biomarkers_ranks),
    col = colorRamp2(1:3, c("white", "grey70", "black")),
    name = "Biomarker Ranks",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    width = unit(2, "cm"),
    height = unit(6, "cm"),
    border = TRUE,
    rect_gp = gpar(col = "black", lwd = 1),
)
pdf(file = file.path(outdir, "biomarker_ranks_heatmap_meth.pdf"), width = 12, height = 10)
draw(heatmap_ranks_meth, heatmap_legend_side = "right")
dev.off()

# ======================== Odds Ratios ========================
dir.create(file.path(outdir, "biomarkers_odds_ratios"), showWarnings = FALSE)

# Calculate odds ratios using one-vs-rest approach for each subtype level
calculate_odds_ratios_one_vs_rest <- function(data, group_var, biomarker_vars) {
    d <- drop_na(data, all_of(c(group_var, biomarker_vars)))

    unique_levels <- unique(d[[group_var]])
    one_vs_rest_results <- map(unique_levels, ~ {
        level <- .x
        d_ovr <- d %>%
            mutate(!!glue("{group_var}_ovr") := ifelse(!!sym(group_var) == level, 1, 0))

        ovr_results <- map(biomarker_vars, ~ {
            biomarker <- .x
            fmla <- as.formula(glue("{group_var}_ovr ~ {biomarker} + AGE_RAND + SEX + RACE + ETHNIC"))
            model <- glm(fmla, data = d_ovr, family = binomial)
            model_tidy <- tidy(model)
            model_tidy <- model_tidy %>%
                filter(term == biomarker) %>%
                mutate(
                    biomarker = biomarker,
                    odds_ratio = exp(estimate),
                    lower_ci = exp(estimate - 1.96 * std.error),
                    upper_ci = exp(estimate + 1.96 * std.error),
                    level = level
                ) %>%
                select(biomarker, level, odds_ratio, lower_ci, upper_ci, p.value)
            return(model_tidy)
        })
        ovr_results <- bind_rows(ovr_results)
        return(ovr_results)
    })
    one_vs_rest_results <- bind_rows(one_vs_rest_results)
    one_vs_rest_results$group <- group_var

    return(one_vs_rest_results)
}

odds_ratios <- map(group_columns_for_testing, ~ {
    calculate_odds_ratios_one_vs_rest(dat, .x, biomarker_variables)
})
odds_ratios <- bind_rows(odds_ratios)
write.csv(odds_ratios, file = file.path(outdir, "biomarkers_odds_ratios", "odds_ratios.csv"), row.names = FALSE)

create_odds_ratios_heatmap <- function(odds_ratios_data, group_name, output_dir) {
    odds_ratios_heatmap_data <- odds_ratios_data %>%
        mutate(
            odds_ratio = ifelse(p.value < 0.05, odds_ratio, NA),
        ) %>%
        filter(group == group_name) %>%
        select(biomarker, level, odds_ratio) %>%
        pivot_wider(names_from = level, values_from = odds_ratio) %>%
        column_to_rownames("biomarker")

    odds_ratios_heatmap <- Heatmap(
        odds_ratios_heatmap_data,
        col = colorRamp2(c(0.5, 1, 2), c("blue", "white", "red")),
        name = "Odds Ratio",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_order = sort(colnames(odds_ratios_heatmap_data)),
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        row_names_side = "left",
        column_names_side = "top",
        width = unit(1 * ncol(odds_ratios_heatmap_data), "cm"),
        height = unit(6, "cm"),
        border = TRUE,
        rect_gp = gpar(col = "black", lwd = 0.5),
    )

    pdf(file = file.path(output_dir, paste0(group_name, "_odds_ratios_heatmap.pdf")), width = 12, height = 10)
    draw(odds_ratios_heatmap, heatmap_legend_side = "right")
    dev.off()
}

create_odds_ratios_heatmap(odds_ratios, "rna_4cluster", file.path(outdir, "biomarkers_odds_ratios"))
create_odds_ratios_heatmap(odds_ratios, "meth_3cluster", file.path(outdir, "biomarkers_odds_ratios"))
create_odds_ratios_heatmap(odds_ratios, "composite_subtype", file.path(outdir, "biomarkers_odds_ratios"))

# ======================== Risk Modeling ========================
dir.create(file.path(outdir, "biomarkers_risk_modeling"), showWarnings = FALSE)

dat_adding <- read_csv("output/run4_rmoutliers2_asr_control/integrative_analyses/run9-figs/bioreptable_waddons.csv")
dat_adding <- dat_adding %>%
    select(
        PATNUM, starts_with("T_"),
    )
dat <- left_join(dat, dat_adding, by = "PATNUM")

# Convert biomarkers to tertiles for risk modeling
biomarker_risk_modeling_dat <- dat %>%
    mutate(
        across(
            all_of(biomarker_variables),
            ~ factor(ntile(.x, 3), levels = 1:3, labels = c("Low", "Medium", "High")),
            .names = "{col}_tertile"
        ),
    )

biomarker_risk_modeling_stats <- map(group_columns_for_testing, ~ {
    d <- drop_na(biomarker_risk_modeling_dat, all_of(c(.x, paste0(biomarker_variables, "_tertile"))))
    t <- stats_table(d, .x, paste0(biomarker_variables, "_tertile"))
    write.csv(t, file = file.path(outdir, "biomarkers_risk_modeling", paste0(.x, "_biomarker_risk_modeling_stats.csv")))
    return(t)
})
names(biomarker_risk_modeling_stats) <- group_columns_for_testing

# Sensitivity analysis: hazard ratios stratified by biomarker tertiles
biomarker_tertiles <- grep("_tertile", colnames(biomarker_risk_modeling_dat), value = TRUE)
biomarker_risk_modeling_results <- map(group_columns_for_testing[1:3], ~ {
    group_var <- .x
    filtered_res <- filtered_hazard_ratio_table(
        biomarker_risk_modeling_dat,
        group_var,
        biomarker_tertiles,
        censors = censors,
        controls = controls_list$demographics,
        censor_prefix = "C_",
        time_prefix = "T_"
    )
})
names(biomarker_risk_modeling_results) <- group_columns_for_testing[1:3]
biomarker_risk_modeling_results <- bind_rows(biomarker_risk_modeling_results, .id = "group")
write.csv(biomarker_risk_modeling_results, file = file.path(outdir, "biomarkers_risk_modeling", "biomarker_risk_modeling_results.csv"), row.names = FALSE)

filtered_heatmap_dat <- map(unique(biomarker_risk_modeling_results$group), ~ {
    group_name <- .x
    biomarker_risk_modeling_results %>%
        filter(group == group_name) %>%
        mutate(
            named_tertiles = glue("{x}_{y}")
        ) %>%
        select(
            term, named_tertiles, censor, hazard.ratio
        ) %>%
        pivot_wider(
            names_from = c(term, named_tertiles),
            values_from = c(hazard.ratio),
        )
})
names(filtered_heatmap_dat) <- unique(biomarker_risk_modeling_results$group)

# Compare subtypes to biomarker-based risk model using logistic regression
run_logistic_model <- function(data, group_var, outcome, controls = c("AGE_RAND", "SEX", "RACE", "ETHNIC")) {
    controls_str <- paste(controls, collapse = " + ")
    fmla <- as.formula(glue("{outcome} ~ {group_var} + {controls_str}"))
    print(fmla)
    print(glue("Data dimensions: {nrow(data)} rows, Outcome events: {sum(data[[outcome]], na.rm = TRUE)}"))

    model <- glm(fmla, data = data, family = binomial)
    return(model)
}

models <- map(censors, ~ {
    outcome <- .x
    controls <- controls_list$risks
    controls <- c(
        controls,
        "TROPT_UCR_clean_tr_tertile",
        "GDF_UCR_clean_tr_tertile",
        "NTPBNP_UCR_clean_tr_tertile",
        "CD40L_clean_tr_tertile"
    )

    d <- biomarker_risk_modeling_dat %>%
        drop_na(all_of(c(
            controls,
            outcome,
            "TROPT_UCR_clean_tr_tertile",
            "GDF_UCR_clean_tr_tertile",
            "NTPBNP_UCR_clean_tr_tertile",
            "CD40L_clean_tr_tertile",
            group_columns_for_testing[1:3]
        )))

    groups <- c(
        "rna_4cluster",
        "meth_3cluster",
        "composite_subtype",
        "TROPT_UCR_clean_tr_tertile + GDF_UCR_clean_tr_tertile + NTPBNP_UCR_clean_tr_tertile + CD40L_clean_tr_tertile"
    )

    models <- map(groups[1:3], ~ run_logistic_model(d, .x, outcome, controls = controls))
    models2 <- map(groups[4], ~ run_logistic_model(d, .x, outcome, controls = controls_list$risks))

    models <- c(models, models2)
    names(models) <- c(
        "RNA Subtypes",
        "Methylation Subtypes",
        "Composite Subtypes",
        "Biomarkers Model"
    )

    rocs <- map(models, ~ {
        preds <- predict(.x, type = "response")
        roc_obj <- roc(d[[outcome]], preds, ci = TRUE)
        return(roc_obj)
    })
    names(rocs) <- names(models)

    roc_labels <- map_chr(names(rocs), ~ {
        roc_obj <- rocs[[.x]]
        auc_val <- round(auc(roc_obj), 3)
        ci_lower <- round(ci(roc_obj)[1], 3)
        ci_upper <- round(ci(roc_obj)[3], 3)
        paste0(.x, " (AUC: ", auc_val, ", 95% CI: ", ci_lower, "-", ci_upper, ")")
    })
    names(roc_labels) <- names(rocs)

    g <- ggroc(rocs, legacy.axes = TRUE) +
        theme_bw() +
        theme(legend.position = "bottom") +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.7) +
        scale_color_discrete(name = "Model") +
        annotate("text",
            x = 1, y = 0, label = paste(roc_labels, collapse = "\n"),
            hjust = 1, vjust = 0, size = 3
        )
    plot_file <- file.path(outdir, "biomarkers_risk_modeling", paste0("roc_curve_", outcome, ".pdf"))
    ggsave(plot_file, plot = g, width = 8, height = 6)

    roc_comparisons <- combn(names(rocs), 2, function(x) {
        roc1 <- rocs[[x[1]]]
        roc2 <- rocs[[x[2]]]
        roc_test <- roc.test(roc1, roc2)
        tibble(
            model1 = x[1],
            model2 = x[2],
            p_value = roc_test$p.value,
            auc1 = auc(roc1)[1],
            auc2 = auc(roc2)[1],
            ci1 = paste0(round(ci(roc1)[1], 3), "-", round(ci(roc1)[3], 3)),
            ci2 = paste0(round(ci(roc2)[1], 3), "-", round(ci(roc2)[3], 3))
        )
    }, simplify = FALSE)
    roc_comparisons <- bind_rows(roc_comparisons)
    write.csv(roc_comparisons, file = file.path(outdir, "biomarkers_risk_modeling", paste0("roc_comparisons_", outcome, ".csv")), row.names = FALSE)

    return(list(models = models, rocs = rocs, roc_plot = g, roc_comparisons = roc_comparisons))
})
names(models) <- censors

# Create forest plot comparing composite subtypes to biomarker model
forest_plot_data <- bind_rows(
    map(models, ~ .x$roc_comparisons),
    .id = "censor"
)

forest_plot_data <- forest_plot_data %>%
    filter(
        model1 == "Composite Subtypes" & model2 == "Biomarkers Model"
    ) %>%
    mutate(
        ci.lower1 = as.numeric(str_extract(ci1, "\\d+\\.\\d+")),
        ci.upper1 = as.numeric(str_extract(ci1, "\\d+\\.\\d+$")),
        ci.lower2 = as.numeric(str_extract(ci2, "\\d+\\.\\d+\\d+")),
        ci.upper2 = as.numeric(str_extract(ci2, "\\d+\\.\\d+$"))
    ) %>%
    pivot_longer(
        cols = c(auc1, auc2, ci.lower1, ci.lower2, ci.upper1, ci.upper2),
        names_to = c(".value", "model"),
        names_pattern = "(.*)([12])$"
    ) %>%
    mutate(
        model = recode(model, `1` = "Composite Subtypes", `2` = "Biomarkers Model"),
        censor = factor(censor, levels = c("C_PRIMARY", "C_CVDMI_P", "C_ACDMI_P", "C_MIPRIM", "C_CVD", "C_ACD"), labels = c("MACE", "CVD or MI", "ACD or MI", "MI", "CVD", "ACD"))
    ) %>%
    select(
        censor, model, p_value,
        auc, ci.lower, ci.upper
    )

forest_plot <- ggplot(forest_plot_data, aes(x = auc, y = fct_rev(censor),
    xmin = ci.lower, xmax = ci.upper, color = model)) +
    geom_pointrange(position = position_dodge(width = 0.4)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey") +
    geom_text(aes(label = paste0("p=", round(p_value, 3))),
              x = 0.9, position = position_dodge(width = 0.4),
              color = "black", size = 3, hjust = 0) +
    scale_x_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, by = 0.1)) +
    labs(
        title = "AUC Comparison of Composite Subtypes vs Biomarkers",
        x = "AUC (95% CI)",
        y = "Censor",
        color = "Model"
    ) +
    theme_matt() +
    theme(
        legend.position = "bottom"
    )
ggsave(
    file = file.path(outdir, "biomarkers_risk_modeling", "forest_plot_auc_comparison.pdf"),
    plot = forest_plot, width = 6, height = 4
)
