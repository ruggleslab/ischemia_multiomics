###########################################################################
#
#              Multiomics Subtypes Additive Risk Models
#
###########################################################################
# Purpose: Evaluate whether composite subtypes add incremental predictive
#          value beyond traditional risk scores using ROC curve analysis
#
# Inputs:  - Cohort tabulation table with subtype assignments
#          - Clinical risk scores (SYNTAX, ASCVD, etc.)
#          - Biomarker, plaque, deconvolution data
#          - Outcome data
#
# Outputs: - Base vs composite subtype model comparisons
#          - ROC curves and AUC comparisons
#          - Forest plots of model performance
###########################################################################

experiment <- "muller_multiomic_subtypes_additive_models"
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

# ======================== Load and Process Data ========================
# TODO: Move these paths to config.R
colorguidefile <- "data/ischemia2021_colorguide.csv"
colorguide <- read.table(colorguidefile, sep = ",", header = TRUE, comment.char = "", colClasses = c("character", "character", "character"))

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

controls_list <- list(
    demographics = c("AGE_RAND", "SEX", "RACE", "ETHNIC"),
    risks = c("AGE_RAND", "SEX", "RACE", "ETHNIC", "DIABETES", "EGFR", "NO_DIALYSIS")
)

group_columns_for_testing <- c("rna_4cluster", "meth_3cluster", "composite_subtype", "composites_of_interest", "cs1_v_rest", "cs2_v_rest")
variables_to_test <- grep("BMI|SMOK|CKD|HYP|DIAB|DIA|SYS|LVEF|HEMOGLOB|HDLC|LDLC|TRIGLYC|TOTCHOL|EGFR|DEGISCH|CTMULT70|CTMULT50|CTNDV70|CTNDV50|IMGDEGIS|CUSTOM_IMGDEGIS_NONEMILD|IMGDEGIS_01_2_3", colnames(tabulation_table), value = TRUE)
censors <- c("C_PRIMARY", "C_CVDMI_P", "C_ACDMI_P", "C_MIPRIM", "C_CVD", "C_ACD")

biomarkers <- read.csv("data/biomarkers_20250906/20241025_biomarker_data_cleaned.csv", row.names = 1)
biomarker_variables <- colnames(biomarkers)[-1]

plaque_data <- read.csv("data/20250615-plaque.csv", row.names = 1)
plaque_variables <- colnames(plaque_data)[-c(1, ncol(plaque_data))]

deconvolution_data <- read.csv("data/20241107_whole_blood_cibersort_abs.csv", row.names = 1)
deconvolution_data$PATNUM <- rownames(deconvolution_data)
deconvolution_variables <- colnames(deconvolution_data)[-ncol(deconvolution_data)]

risk_scores <- read.csv("data/20210511_ISCHEMIA_riskscores_v1.csv", row.names = 1)
risk_scores$PATNUM <- rownames(risk_scores)
risk_score_variables <- colnames(risk_scores)[-ncol(risk_scores)]

meta <- read.csv("data/biorep_10_27.csv", row.names = 1)
meta <- meta %>% select(PATNUM, REGIONC, COUNTRY)
meta_variables <- c("REGIONC", "COUNTRY")

dat <- tabulation_table %>%
    left_join(biomarkers) %>%
    left_join(plaque_data) %>%
    left_join(deconvolution_data) %>%
    left_join(risk_scores) %>%
    left_join(meta)

# ======================== Risk Model Functions ========================
dir.create(file.path(outdir, "risk_factors"), showWarnings = FALSE)

run_logistic_model <- function(data, group_var, outcome, controls) {
    controls_str <- paste(controls, collapse = " + ")
    fmla <- as.formula(glue("{outcome} ~ {group_var} + {controls_str}"))
    print(fmla)
    print(glue("Data dimensions: {nrow(data)} rows, Outcome events: {sum(data[[outcome]], na.rm = TRUE)}"))

    model <- glm(fmla, data = data, family = binomial)
    return(model)
}

prepare_modeling_data <- function(data, outcome, controls, group_var) {
    required_cols <- c(controls, group_var, outcome)
    data %>%
        drop_na(all_of(required_cols))
}

# Fit base model (risk score alone) and composite model (risk score + subtypes)
fit_outcome_models <- function(data, outcome, risk_score, controls) {
    base_model <- run_logistic_model(data, risk_score, outcome, controls = controls)
    composite_model <- run_logistic_model(data, "composite_subtype", outcome, controls = c(controls, risk_score))

    models <- list(
        "Base Model" = base_model,
        "Composite Subtypes" = composite_model
    )

    return(models)
}

generate_roc_objects <- function(models, data, outcome) {
    rocs <- map(models, ~ {
        preds <- predict(.x, type = "response")
        roc_obj <- roc(data[[outcome]], preds, ci = TRUE)
        return(roc_obj)
    })
    names(rocs) <- names(models)
    return(rocs)
}

create_roc_plot <- function(rocs) {
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

    return(g)
}

perform_roc_comparisons <- function(rocs) {
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

    return(bind_rows(roc_comparisons))
}

save_outcome_results <- function(roc_plot, roc_comparisons, outcome, output_dir) {
    plot_file <- file.path(output_dir, paste0("roc_curve_", outcome, ".pdf"))
    ggsave(plot_file, plot = roc_plot, width = 8, height = 6)

    comp_file <- file.path(output_dir, paste0("roc_comparisons_", outcome, ".csv"))
    write.csv(roc_comparisons, file = comp_file, row.names = FALSE)
}

# Process all outcomes for a given risk score
process_outcomes <- function(dat, censors, risk_score, controls_list, outdir) {
    models <- map(censors, ~ {
        outcome <- .x

        d <- prepare_modeling_data(
            dat,
            outcome,
            controls_list$risks,
            "composite_subtype"
        )

        models <- fit_outcome_models(d, outcome, risk_score, controls_list$risks)

        rocs <- generate_roc_objects(models, d, outcome)

        roc_plot <- create_roc_plot(rocs)

        roc_comparisons <- perform_roc_comparisons(rocs)

        save_outcome_results(roc_plot, roc_comparisons, outcome, outdir)

        return(list(
            models = models,
            rocs = rocs,
            roc_plot = roc_plot,
            roc_comparisons = roc_comparisons
        ))
    })
    names(models) <- censors

    return(models)
}

create_auc_forest_plot <- function(models, outdir, censors) {
    forest_plot_data <- bind_rows(
        map(models, ~ .x$roc_comparisons),
        .id = "censor"
    )

    forest_plot_data <- forest_plot_data %>%
        filter(
            model2 == "Composite Subtypes" & model1 == "Base Model"
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
            model = recode(model, `1` = "Base Model", `2` = "Composite Subtypes"),
            censor = factor(censor,
                levels = c("C_PRIMARY", "C_CVDMI_P", "C_ACDMI_P", "C_MIPRIM", "C_CVD", "C_ACD"),
                labels = c("MACE", "CVD or MI", "ACD or MI", "MI", "CVD", "ACD")
            )
        ) %>%
        select(
            censor, model, p_value,
            auc, ci.lower, ci.upper
        )

    forest_plot <- ggplot(forest_plot_data, aes(
        x = auc, y = fct_rev(censor),
        xmin = ci.lower, xmax = ci.upper, color = model
    )) +
        geom_pointrange(position = position_dodge(width = 0.4)) +
        geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey") +
        geom_text(aes(label = paste0("p=", round(p_value, 3))),
            x = 0.85, position = position_dodge(width = 0.4),
            color = "black", size = 3, hjust = 0
        ) +
        scale_x_continuous(limits = c(0.5, 0.9), breaks = seq(0.5, 0.9, by = 0.1)) +
        labs(
            x = "AUC (95% CI)",
            y = NULL,
            color = "Model"
        ) +
        theme_matt() +
        theme(
            legend.position = "bottom"
        )

    ggsave(
        file = file.path(outdir, "forest_plot_auc_comparison.pdf"),
        plot = forest_plot, width = 6, height = 4
    )

    return(forest_plot)
}

# ======================== Run Analysis ========================
# Evaluate all risk scores to determine incremental predictive value
risk_score_variables <- c(risk_score_variables, "DEGRISCH")
models <- map(risk_score_variables, ~ {
    risk_score <- .x
    cat("Processing risk score:", risk_score, "\n")

    risk_score_dir <- file.path(outdir, risk_score)
    dir.create(risk_score_dir, showWarnings = FALSE, recursive = TRUE)

    models <- process_outcomes(dat, censors, risk_score, controls_list, risk_score_dir)

    forest_plot <- create_auc_forest_plot(models, risk_score_dir, censors)

    return(models)
})
names(models) <- risk_score_variables
saveRDS(models, file = file.path(outdir, "models.rds"))
