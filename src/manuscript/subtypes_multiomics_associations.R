###########################################################################
#
#           Multiomics Subtypes and Multi-omics Associations
#
###########################################################################
# Purpose: Analyze associations between RNA/methylation subtypes and
#          multiple data modalities (biomarkers, plaque, deconvolution,
#          risk scores, clinical regions)
#
# Inputs:  - Cohort tabulation table with subtype assignments
#          - Biomarker, plaque, deconvolution, risk score data
#          - Regional metadata
#
# Outputs: - Odds ratio analyses across all data types
#          - Regional distribution analyses
#          - Comprehensive multi-modal associations
###########################################################################

experiment <- "muller_multiomic_subtypes_and_omics"
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

# Scale deconvolution variables for consistent analysis
dat <- dat %>%
    mutate(
        across(deconvolution_variables, ~ scale(.x, center = TRUE, scale = TRUE)),
    )

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
        col = colorRamp2(c(0.1, 1, 2), c("blue", "white", "red")),
        name = "Odds Ratio",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_order = sort(colnames(odds_ratios_heatmap_data)),
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        row_names_side = "left",
        column_names_side = "top",
        width = unit(1 * ncol(odds_ratios_heatmap_data), "cm"),
        height = unit(8, "cm"),
        border = TRUE,
        rect_gp = gpar(col = "black", lwd = 0.5),
    )

    pdf(file = file.path(output_dir, paste0(group_name, "_odds_ratios_heatmap.pdf")), width = 12, height = 10)
    draw(odds_ratios_heatmap, heatmap_legend_side = "right")
    dev.off()
}

# Calculate odds ratios for each variable type
variable_groups <- list(
    biomarkers = biomarker_variables,
    plaque = plaque_variables,
    deconvolution = deconvolution_variables,
    risk_scores = risk_score_variables
)

for (var_type in names(variable_groups)) {
    var_list <- variable_groups[[var_type]]

    type_outdir <- file.path(outdir, paste0(var_type, "_odds_ratios"))
    dir.create(type_outdir, showWarnings = FALSE)

    odds_ratios <- map(group_columns_for_testing, ~ {
        calculate_odds_ratios_one_vs_rest(dat, .x, var_list)
    })
    odds_ratios <- bind_rows(odds_ratios)
    write.csv(odds_ratios, file = file.path(type_outdir, "odds_ratios.csv"), row.names = FALSE)

    create_odds_ratios_heatmap(odds_ratios, "rna_4cluster", type_outdir)
    create_odds_ratios_heatmap(odds_ratios, "meth_3cluster", type_outdir)
    create_odds_ratios_heatmap(odds_ratios, "composite_subtype", type_outdir)
}

# ======================== Regional Analysis ========================
dir.create(file.path(outdir, "regions"), showWarnings = FALSE)

cluster_types <- c("rna_4cluster", "meth_3cluster", "composite_subtype")

for (cluster_type in cluster_types) {
    dir.create(file.path(outdir, "regions", cluster_type), showWarnings = FALSE)

    dat_regions <- dat %>%
        select(PATNUM, REGIONC, !!sym(cluster_type)) %>%
        drop_na(REGIONC, !!sym(cluster_type))

    region_table <- with(dat_regions, table(REGIONC, get(cluster_type)))
    region_proportions <- prop.table(region_table, margin = 2)

    write.csv(region_table, file = file.path(outdir, "regions", cluster_type, "region_table.csv"))
    write.csv(region_proportions, file = file.path(outdir, "regions", cluster_type, "region_proportions.csv"))

    region_proportions_df <- as.data.frame(region_proportions)
    colnames(region_proportions_df) <- c("Region", "Subtype", "Proportion")
    region_proportions_df$Region <- factor(region_proportions_df$Region, levels = unique(region_proportions_df$Region))
    region_proportions_df <- region_proportions_df %>%
        mutate(
            N = round(Proportion * rowSums(region_table)),
            Proportion = paste0(round(Proportion * 100, 1), "%"),
            Label = paste0(N, " (", Proportion, ")")
        )
    region_proportions_df <- region_proportions_df %>%
        select(Region, Subtype, Label) %>%
        pivot_wider(names_from = Subtype, values_from = Label)
    write.csv(region_proportions_df, file = file.path(outdir, "regions", cluster_type, "region_proportions_labels.csv"), row.names = FALSE)

    chisq_result <- chisq.test(region_table)

    region_proportions_df <- as.data.frame(region_proportions)
    colnames(region_proportions_df) <- c("Region", "Subtype", "Proportion")
    region_proportions_df$Region <- factor(region_proportions_df$Region, levels = unique(region_proportions_df$Region))

    barplot <- ggplot(region_proportions_df, aes(x = Subtype, y = Proportion, fill = Region)) +
        geom_bar(stat = "identity", position = "fill") +
        labs(title = paste("Region Proportions for", cluster_type), x = "Region", y = "Proportion") +
        theme_minimal()
    ggsave(file.path(outdir, "regions", cluster_type, "region_proportions_barplot.pdf"),
           plot = barplot, width = length(unique(region_proportions_df$Subtype)) * 0.5, height = 2)

    write.csv(tidy(chisq_result), file = file.path(outdir, "regions", cluster_type, "chi_square_results.csv"))
}
