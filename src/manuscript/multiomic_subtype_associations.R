###########################################################################
#
#                       Multiomic Subtype Associations
#
###########################################################################
# Author: ISCHEMIA Study Team
# Date: Updated 2025-08-08
# Description: Multiomic Subtype Associations analysis for ISCHEMIA study
#

# Load configuration and utilities
source("config.R")
source("utils.R")

# Load required packages
# TODO: Update this list based on actual packages used in the script
# load_packages(c("package1", "package2"))

experiment <- "muller_multiomic_subtype_associations_reviewer_comments"
outdir <- file.path("output", experiment)
dir.create(outdir, showWarnings = FALSE)

# ======================== LIBRARIES ========================
library(tidyverse)
library(glue)
library(rmatt)
library(broom)
library(survival)
library(pROC)

set.seed(420)

# ======================== CODE ========================
dat <- read_csv("output/run4_rmoutliers2_asr_control/integrative_analyses/run9-figs/bioreptable_waddons.csv")

# set any values of 2 in the censors to 0 for the purposes of this analysis
dat <- dat %>%
    mutate_at(vars(starts_with("C_")), ~ ifelse(.x == 2, NA, .x)) %>%
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
        composite_subtype = interaction(
            gsub("meth_3cluster_", "MS", meth_3cluster),
            gsub("nmf_cluster_", "RS", rna_4cluster),
            drop = TRUE, sep = "_"
        ),
        composites_of_interest = factor(
            composite_subtype,
            levels = c("MS3_RS2", "MS3_RS4"),
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

        # Create composite subtypes
        composite_subtype = interaction(
            gsub("meth_3cluster_", "MS", meth_3cluster),
            gsub("nmf_cluster_", "RS", rna_4cluster),
            drop = TRUE, sep = "_"
        ),
        composites_of_interest = factor(
            composite_subtype,
            levels = c("MS3_RS2", "MS3_RS4"),
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

# --------------------------------- In group vs outgroup enrichment and summary ---------------------------------
group_columns_for_testing <- c("rna_4cluster", "meth_3cluster", "composite_subtype")
variables_to_test <- grep("BMI|SMOK|CKD|HYP|DIAB|DIA|SYS|LVEF|HEMOGLOB|HDLC|LDLC|TRIGLYC|TOTCHOL|EGFR|DEGISCH|CTMULT70|CTMULT50|CTNDV70|CTNDV50|IMGDEGIS|CUSTOM_IMGDEGIS_NONEMILD|IMGDEGIS_01_2_3", colnames(tabulation_table), value = TRUE)
censors <- c("C_PRIMARY", "C_CVDMI_P", "C_ACDMI_P", "C_MIPRIM", "C_CVD", "C_ACD")

# # We need to add adjusted (asr) odds ratios to the tabulation table
# groups_odds <- purrr::map(group_columns_for_testing, ~ {
#     values <- unique(na.omit(tabulation_table[[.x]]))
#     gs <- paste0(.x, "_", values)
#     tmp <- one_hot_encode_ovr(tabulation_table, .x, binary = FALSE)
#     purrr::map(gs, ~ {
#         g <- .x
#         purrr::map(variables_to_test, ~ {
#             v <- .x
#             # one hot encode the group column
#             tmp <- tmp %>%
#                 mutate(!!g := as.factor(!!sym(g)))
#             fmla <- as.formula(glue("{g} ~ {v} + AGE_RAND + SEX + RACE + ETHNIC"))
#             fit <- glm(fmla, data = tmp, family = binomial(link = "logit"))
#             res <- broom::tidy(fit)
#             res <- res %>%
#                 filter(grepl(v, term)) %>%
#                 mutate(
#                     group = g,
#                     variable = v,
#                     adjustments = paste0(controls_list$demographics, collapse = ", "),
#                     estimate = exp(estimate),
#                     p.value = p.value
#                 )
#             res
#         }) %>% bind_rows()
#     }) %>% bind_rows()
# })
# names(groups_odds) <- group_columns_for_testing
# saveRDS(groups_odds, file.path(outdir, "groups_odds.rds"))
# lapply(names(groups_odds), function(x) {
#     write_csv(groups_odds[[x]], file.path(outdir, paste0(x, "_odds_ratios.csv")))
# })

# groups_hazards <- purrr::map(group_columns_for_testing, ~ {
#     values <- unique(na.omit(dat[[.x]]))
#     gs <- paste0(.x, "_", values)
#     tmp <- one_hot_encode_ovr(dat, .x, binary = FALSE)
#     d <- purrr::map(gs, ~ {
#         g <- .x
#         res <- hazard_ratios_table(tmp, g, censors, controls_list$demographics, censor_prefix = "C_", time_prefix = "T_")
#         res <- res %>%
#             mutate(
#                 group = g,
#                 adjustments = "AGE, SEX, RACE, ETHNIC"
#             )
#         res
#     })
#     d <- bind_rows(d)
# })
# names(groups_hazards) <- group_columns_for_testing
# saveRDS(groups_hazards, file.path(outdir, "groups_hazards.rds"))
# lapply(names(groups_hazards), function(x) {
#     write_csv(groups_hazards[[x]], file.path(outdir, paste0(x, "_harzards_ratios.csv")))
# })

# # let's also make a logrank p value table
# groups_logranks <- purrr::map(group_columns_for_testing, ~ {
#     g <- .x
#     purrr::map(censors, ~ {
#         time <- gsub("C_", "T_", .x)
#         formula <- as.formula(glue("Surv({time}, {.x}) ~ {g} + AGE_RAND + SEX + RACE + ETHNIC"))
#         cox_model <- coxph(formula, data = dat)
#         summary_cox <- summary(cox_model)
#         res <- data.frame(
#             group = g,
#             censor = .x,
#             chisq = summary_cox$waldtest["test"],
#             p.value = summary_cox$waldtest["pvalue"]
#         )
#         res
#     }) %>% bind_rows()
# })
# names(groups_logranks) <- group_columns_for_testing
# saveRDS(groups_logranks, file.path(outdir, "groups_logranks.rds"))
# lapply(names(groups_logranks), function(x) {
#     write_csv(groups_logranks[[x]], file.path(outdir, paste0(x, "_logranks.csv")))
# })

# for (group_column_selected in group_columns_for_testing) {
#     COI_tabulation_subpath <- file.path(outdir, "cohort_tabulations/", group_column_selected, "/")
#     dir.create(COI_tabulation_subpath, showWarnings = FALSE, recursive = TRUE)
#     out1 <- ingroup_vs_outgroup_cohort_enrichment_tests(tabulation_table, group_column = group_column_selected, cohort_tabulations_outpath = COI_tabulation_subpath)
#     out1_ratiotable <- out1[[1]]
#     out1_citable <- out1[[2]]

#     # Add Ns to the col labels
#     N_annot_table <- table(tabulation_table[,group_column_selected])[colnames(out1_ratiotable)]
#     colnames_w_N <- paste0(names(N_annot_table), paste0("__N=", N_annot_table))

#     # heatmapcolorparam <- colorRamp2(breaks = c(-0.05, -0.000000001, 0.000000001, 0.05), c("white", "darkblue", "darkred", "white"))
#     # heatmapcolorparam <- colorRamp2(breaks = c(-5, -0.5, 0, 0.5, 5), c("darkblue", "blue", "white", "red", "darkred"))
#     heatmapcolorparam <- colorRamp2(breaks = c(0, 0.99999, 1, 1.00001, 5), c("#000099", "#ccccff", "white", "#ffb2b2", "#b20000"))
#     rowmetatable <- data.frame(category = unlist(lapply(strsplit(rownames(out1_ratiotable), split = "__"), function(x) x[1])), row.names = rownames(out1_ratiotable))[COI[COI %in% rownames(out1_ratiotable)],,drop=FALSE]
#     rowannotationlist <- annotationlist_builder(rowmetatable)
#     plottable <- out1_ratiotable[COI[COI %in% rownames(out1_ratiotable)],]
#     colnames(plottable) <- colnames_w_N
#     hm1 <- create_heatmap(
#         counttab = plottable,
#         # counttab = out1_ratiotable[COI[COI %in% rownames(out1_ratiotable)],],
#         subsetnum = FALSE, scale_data = FALSE,
#         rowmetatable = rowmetatable,
#         # rowmetatable = rowmetatable[COI[COI %in% rownames(out1_ratiotable)],,drop=FALSE],
#         rowannotationlist = rowannotationlist,
#         separate_legend = TRUE, heatmapcolorparam = heatmapcolorparam, addborders = TRUE)
#     pdf(paste0(COI_tabulation_subpath, "signed_filtered_sigfeature_cluster_hm.pdf"), 10, 10, useDingbats = FALSE)
#     draw(hm1[[1]])
#     junk <- dev.off()

#     write.table(out1_ratiotable, paste0(COI_tabulation_subpath, "signed_filtered_sigfeature_cluster_ratio_table.csv"), col.names = NA, row.names = TRUE, sep = ",")
#     write.table(out1_citable, paste0(COI_tabulation_subpath, "signed_filtered_sigfeature_cluster_ci_table.csv"), col.names = NA, row.names = TRUE, sep = ",")
# }

# ======================== Association Modelings ========================
dir.create(file.path(outdir, "association_modelings"), showWarnings = FALSE)

# Get some variables of interest
outcomes <- c("C_PRIMARY", "C_CVDMI_P", "C_ACDMI_P", "C_MIPRIM", "C_CVD", "C_ACD")
group_columns_for_testing <- c("rna_4cluster", "meth_3cluster", "composite_subtype")

# now let's make models for each of the outcomes
models <- purrr::map(outcomes, ~ {
    dir.create(file.path(outdir, "association_modelings", .x), showWarnings = FALSE, recursive = TRUE)
    outcome <- .x

    # Model 1: age, sex, race, ethnicity
    fmla_demos <- as.formula(glue("{outcome} ~ AGE_RAND + SEX + RACE + ETHNIC"))

    # Model 2: age, sex, race, ethnicity, diabetes, eGFR, LVEF
    fmla_risk <- as.formula(glue("{outcome} ~ AGE_RAND + SEX + RACE + ETHNIC + DIABETES + EGFR * NO_DIALYSIS"))

    # Create models for each group column (subtypes)
    purrr::map(group_columns_for_testing, ~ {
        group_col <- .x

        # make sure we have no missing data
        dat <- tabulation_table %>% drop_na(outcome, "composite_subtype", all_of(c("AGE_RAND", "SEX", "RACE", "ETHNIC", "DIABETES", "EGFR", "NO_DIALYSIS")))

        # Model 3: age, sex, race, ethnicity, diabetes, eGFR, LVEF, Subtype
        fmla_risk_group <- as.formula(glue("{outcome} ~ AGE_RAND + SEX + RACE + ETHNIC + DIABETES + EGFR * NO_DIALYSIS + {group_col}"))

        # Fit all three models
        model1 <- glm(fmla_demos, data = dat, family = binomial(link = "logit"))
        model2 <- glm(fmla_risk, data = dat, family = binomial(link = "logit"))
        model3 <- glm(fmla_risk_group, data = dat, family = binomial(link = "logit"))
 
        # Save models and comparisons
        model_list <- list(model1 = model1, model2 = model2, model3 = model3)
        saveRDS(model_list, file.path(outdir, "association_modelings", outcome, paste0(group_col, "_models.rds")))

        # Also compare models with ANOVA
        model_anova <- anova(model1, model2, model3, test = "Chisq")
        saveRDS(model_anova, file.path(outdir, "association_modelings", outcome, paste0(group_col, "_anova.rds")))

        # Calculate ROC curves for all models
        roc1 <- roc(dat[[outcome]], predict(model1, type = "response"), smooth = TRUE, ci = TRUE)
        roc2 <- roc(dat[[outcome]], predict(model2, type = "response"), smooth = TRUE, ci = TRUE)
        roc3 <- roc(dat[[outcome]], predict(model3, type = "response"), smooth = TRUE, ci = TRUE)

        # plot ROC curves
        # Create a data frame for plotting ROC curves
        roc_list <- list(
            "Model 1 (Demographics)" = roc1,
            "Model 2 (Risk Factors)" = roc2,
            "Model 3 (+ Subtypes)" = roc3
        )

        # Plot using ggroc
        g <- ggroc(roc_list, legacy.axes = TRUE) +
            theme_bw() +
            theme(legend.position = "bottom") +
            geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.7) +
            scale_color_manual(
                values = c("Model 1 (Demographics)" = "red", "Model 2 (Risk Factors)" = "blue", "Model 3 (+ Subtypes)" = "green"),
                name = "Models",
                labels = c(
                    paste0("Model 1 (Demographics) - AUC: ", round(auc(roc1), 3)),
                    paste0("Model 2 (Risk Factors) - AUC: ", round(auc(roc2), 3)),
                    paste0("Model 3 (+ Subtypes) - AUC: ", round(auc(roc3), 3))
                )
            )

        # Save the plot
        ggsave(file.path(outdir, "association_modelings", outcome, paste0(outcome, "_", group_col, "_roc.pdf")), g, width = 6, height = 6)

        # Save ROC objects
        roc_list <- list(roc1 = roc1, roc2 = roc2, roc3 = roc3)
        saveRDS(roc_list, file.path(outdir, "association_modelings", outcome, paste0(outcome, "_", group_col, "_roc.rds")))

        # get the aucs
        aucs <- c(
            outcome = outcome,
            group = group_col,
            model1 = auc(roc1),
            model1_ci_lower = roc1$ci[1],
            model1_ci_upper = roc1$ci[3],
            model2 = auc(roc2),
            model2_ci_lower = roc2$ci[1],
            model2_ci_upper = roc2$ci[3],
            model3 = auc(roc3),
            model3_ci_lower = roc3$ci[1],
            model3_ci_upper = roc3$ci[3]
        )
        out1 <- data.frame(
            outcome = outcome,
            group = group_col,
            model1_auc = aucs["model1"],
            model2_auc = aucs["model2"],
            model3_auc = aucs["model3"]
        )
        list(res = out1, aucs = aucs, roc_list = roc_list)
    }) %>% set_names(group_columns_for_testing)
}) %>% set_names(outcomes)
saveRDS(models, file.path(outdir, "association_models.rds"))

models_df_aucs <- bind_rows(map(models, map_dfr, "aucs"))
colnames(models_df_aucs) <- c("outcome", "group", "model1_auc", "model1_auc_ci_low", "model1_auc_ci_high", "model2_auc", "model2_auc_ci_low", "model2_auc_ci_high", "model3_auc", "model3_auc_ci_low", "model3_auc_ci_high")
write_csv(models_df_aucs, file.path(outdir, "association_models_aucs.csv"))

models_df_aucs_plotdata <- models_df_aucs %>%
    pivot_longer(cols = starts_with("model")) %>%
    mutate(
        model = gsub("_auc(_ci_low|_ci_high)?", "", name),
        metric = gsub("model[0-9]_", "", name),
        value = as.numeric(value)
    ) %>%
    # Remove unwanted groups
    filter(group != "cs1_v_rest", group != "cs2_v_rest") %>%
    # Set the outcome order
    mutate(
        outcome = factor(
            outcome,
            levels = c("C_ACD", "C_CVD", "C_MIPRIM", "C_ACDMI_P", "C_CVDMI_P", "C_PRIMARY"),
            labels = c("ACD", "CVD", "MI", "ACD / MI", "CVD / MI", "Primary")
        ),
        model = factor(
            model,
            levels = c("model3", "model2", "model1"),
            labels = c("Model 3 (+ Subtypes)", "Model 2 (Risk Factors)", "Model 1 (Demographics)")
        )
    ) %>%
    select(outcome, group, model, metric, value) %>%
    pivot_wider(names_from = metric, values_from = value)
models_df_aucs_plotdata

# Create a forest plot of AUCs
auc_forest_plot <- ggplot(models_df_aucs_plotdata, aes(x = auc, y = outcome, color = model, group = model)) +
    geom_point(size = 3, position = position_dodge(width = 0.6)) +
    geom_linerange(aes(xmin = auc_ci_low, xmax = auc_ci_high), position = position_dodge(width = 0.6)) +
    theme_bw() +
    theme(legend.position = "bottom") +
    facet_wrap(~group, scales = "free_y", nrow = 1) +
    labs(
        x = "AUC",
        y = "Outcome",
        color = "Model",
        title = "AUCs for Different Models Across Outcomes"
    )
ggsave(file.path(outdir, "association_models_auc_forest_plot.pdf"), auc_forest_plot, width = 12, height = 4)

# ======================== DeLong Testing ========================
# I'd like to do delong testing for the AUC curves now
# First, let's get the roc objects
roc_list <- map(models, map, "roc_list")

get_roc_tests <- function(l) {
    l <- lapply(l, function(x) {
        lapply(x, function(y) {
            roc.test(y$roc2, y$roc3, method = "bootstrap")
        })
    })
    l
}
model2_v_model3_roc_tests <- get_roc_tests(roc_list)
model2_v_model3_roc_tests_tidy <- bind_rows(map(model2_v_model3_roc_tests, map_dfr, tidy)) %>%
    mutate(
        outcome = rep(names(model2_v_model3_roc_tests), map_dbl(model2_v_model3_roc_tests, length)),
        group = rep(names(model2_v_model3_roc_tests[[1]]), length(model2_v_model3_roc_tests))
    )
write_csv(model2_v_model3_roc_tests_tidy, file.path(outdir, "model2_v_model3_roc_tests_tidy.csv"))

# ======================== Model Stats ========================
generate_p_values <- function(models) {
    # First map over the outcomes
    res1 <- map(models, function(outcome) {
        # Then map over the group columns
        res2 <- map(outcome, function(subtype) {
            roc_clinical <- subtype$roc_list$roc2
            roc_subtype <- subtype$roc_list$roc3
            test <- roc.test(roc_clinical, roc_subtype, method = "bootstrap", boot.n = 2000)
            res3 <- data.frame(
                roc_clinical = roc_clinical$auc,
                roc_subtype = roc_subtype$auc,
                p.value = test$p.value,
                method = "bootstrap"
            )
            res3
        })
        res2
    })
}

model_res <- generate_p_values(models)
saveRDS(model_res, file.path(outdir, "model_interaction_testing.rds"))

results_table <- purrr::map_dfr(model_res, ~ {
    purrr::map_dfr(.x, ~ {
        data.frame(
            group = deparse(substitute(.x)),
            p.value = .x$p.value
        )
    }, .id = "group")
}, .id = "outcome")
write_csv(results_table, file.path(outdir, "model_interaction_testing.csv"))

results_table %>% filter(p.value < 0.05)

# ======================== END ========================
writeLines(capture.output(sessionInfo()), file.path(outdir, "session.log"))

# ======================== OLD CODE ========================
# # Load in the metatable and bioreptable
# metatable_file <- "output/run4_rmoutliers2_asr_control/rna_processing/metatable_in.csv"
# metatable <- read.table(metatable_file, sep = ",", header = TRUE, row.names = 1)
# biorepfile <- "data/isch_BiorepData_3_23_2022_CUSTOM.csv"
# bioreptable <- read.table(biorepfile, header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)
# rownames(bioreptable) <- bioreptable[,"PATNUM"]

# rna_clustermembership_table_file <- "data/MULTIOMICS_SUBTYPE_LABELS/nmf_cluster_labels_with_subcluster_1_20220208.csv"
# rna_clustermembership_table <- read.table(rna_clustermembership_table_file, sep = ",", header = TRUE, row.names = 1)
# rna_clustermembership_table <- rna_clustermembership_table[!is.na(rna_clustermembership_table[,"rna_4cluster_w3AB"]),]
# rna_clustermembership_table <- rna_clustermembership_table[order(rna_clustermembership_table[,"rna_4cluster_w3AB"]),]

# meth_clustermembership_table_file <- "data/MULTIOMICS_SUBTYPE_LABELS/meth_cluster_membership_table_4_nodup_20220606.csv"
# meth_clustermembership_table <- read.table(meth_clustermembership_table_file, sep = ",", header = TRUE, row.names = 1)
# meth_clustermembership_table <- meth_clustermembership_table[order(meth_clustermembership_table[,"meth_4cluster"]),]

# combined_clustermembership_table <- merge(rna_clustermembership_table[,c("rna_4cluster_w3AB", "rna_4cluster"), drop = FALSE],
#                                           meth_clustermembership_table[,c("meth_4cluster", "meth_3cluster"),drop=FALSE],
#                                           by.x = "row.names", by.y = "row.names", all = TRUE)
# colnames(combined_clustermembership_table)[1] <- c("PATNUM")
# rownames(combined_clustermembership_table) <- combined_clustermembership_table[,"PATNUM"]

# ## Making conenience variables for all our samples, all rna, and all meth
# SOIall <- rownames(bioreptable)
# SOIrna <- unique(rownames(rna_clustermembership_table))
# SOImeth <- unique(rownames(meth_clustermembership_table))
# SOIboth <- intersect(rownames(rna_clustermembership_table), rownames(meth_clustermembership_table))
# SOIeither <- unique(c(rownames(rna_clustermembership_table), rownames(meth_clustermembership_table)))
