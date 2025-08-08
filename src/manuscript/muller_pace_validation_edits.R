###########################################################################
#
#                            muller_pace_validation_edits
#
###########################################################################
# Author: Matthew Muller
# Date: 2024-02-04
# Script Name: muller_pace_validation_edits
# Output directory:
experiment <- "muller_pace_validation_edits"
outdir <- file.path("output", paste0(experiment))
dir.create(outdir, showWarnings = F)

# ======================== LIBRARIES ========================#
library(tidyverse)
library(magrittr)
library(rmatt)
library(readxl)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(WGCNA)
library(ComplexHeatmap)

# LOAD FUNCTIONS
# space reserved for sourcing in functions
# source('https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R')
# source('https://raw.githubusercontent.com/mattmuller0/Rtools/main/stats_functions.R')

# ======================== PACE Data ========================#
# load in the clusters
clusters <- read.csv("output/run4_rmoutliers2_asr_control/validation/PACE/run4/singscore_to_cluster_classifier/50Up_singscoregenes/event_analysis/cluster_intable.csv", row.names = 1)
with(clusters, table(pace_predict_labels))

pace_wb_normcounttable <- read.table("data/validation_data/pace_grab/normcounttab.txt", header = TRUE, row.names = 1, sep = "\t")
paceONLY_wb_normcounttable <- pace_wb_normcounttable[, grepl("PACE", colnames(pace_wb_normcounttable))]
thrONLY_wb_normcounttable <- pace_wb_normcounttable[, grepl("THR", colnames(pace_wb_normcounttable))]

# get the metadata
pace_wb_metatable <- read.table("data/validation_data/pace_grab/metatable_in.csv", header = TRUE, row.names = 1, sep = ",")
pace_wb_full_metatable <- read.table("data/validation_data/pace_grab/pace_comb_allmetadata_20201014.csv", sep = ",", header = TRUE, row.names = 1)
# Add Age cat on to this for later calculations # 24-52 53-63 64-91
pace_wb_full_metatable[, "age_cat"] <- ifelse(
    pace_wb_full_metatable[, "current_age"] >= 71 &
        rownames(pace_wb_full_metatable) %in% colnames(paceONLY_wb_normcounttable), "71-97",
    ifelse(pace_wb_full_metatable[, "current_age"] >= 60 &
        rownames(pace_wb_full_metatable) %in% colnames(paceONLY_wb_normcounttable), "61-70",
    ifelse(pace_wb_full_metatable[, "current_age"] >= 25 &
        rownames(pace_wb_full_metatable) %in% colnames(paceONLY_wb_normcounttable), "25-60", NA)
    )
)

# We want to calc an ANY EVENT column: (should be MACLE2 effectively, but lets just double check) - its slightly different, so good to have
# pace_wb_full_metatable[,"censor_ANYEVENT"] <- ifelse(rowSums(pace_wb_full_metatable[,grepl("censor_", colnames(pace_wb_full_metatable)) & !grepl("_30", colnames(pace_wb_full_metatable))]) == 0, 0, 1)
pace_wb_full_metatable[, c("censor_ANYEVENT", "time_to_ANYEVENT")] <- t(apply(pace_wb_full_metatable[
    ,
    grepl("censor_|time_to", colnames(pace_wb_full_metatable)) & !grepl("_30", colnames(pace_wb_full_metatable))
], 1, function(x) {
    censor_vars <- x[grepl("censor", names(x))]
    timeto_vars <- x[grepl("time_to", names(x))]
    censor_out <- ifelse(sum(censor_vars) > 0, 1, 0)
    if (censor_out == 1) {
        timeto_out <- min(timeto_vars, na.rm = TRUE)
    } else {
        timeto_out <- max(timeto_vars, na.rm = TRUE)
    }
    c(censor_ANYEVENT = censor_out, time_to_ANYEVENT = timeto_out)
}))
pace_wb_full_metatable <- pace_wb_full_metatable %>%
    mutate(
        censor_ANYAMPU = ifelse(censor_ampu == 1 | censor_minor_ampu == 1, 1, 0),
        censor_ANYREOP = ifelse(censor_reop == 1 | censor_minor_reop == 1, 1, 0),
        time_to_ANYAMPU = pmin(time_to_ampu, time_to_minor_ampu, na.rm = TRUE),
        time_to_ANYREOP = pmin(time_to_reop, time_to_minor_reop, na.rm = TRUE)
    )


## These were the SOI from the PACE Paper
# SOI <- unique(c(rownames(pace_wb_metatable[pace_wb_metatable[,"Cohort"] %in% "PACE",,drop=FALSE]),
#                 rownames(na.omit(pace_wb_metatable[,"comp_THRGrOr_PACEg4_v_THRg4",drop=FALSE]))))

# So the final run that was used was run 4, therefore we will base our analysis on this -- that is the clusters we loaded in at the start
pace_wb_metatable <- pace_wb_metatable[rownames(clusters), ]
pace_wb_full_metatable <- pace_wb_full_metatable[rownames(clusters), ]

# add the eGFR and other metabolic data into the metadata
pace_metabolics <- readxl::read_xlsx("data/validation_data/pace_grab/CKDLTA_vMaster.xlsx")
pace_wb_full_metatable[, c("eGFR", "creatinine")] <- pace_metabolics[match(rownames(pace_wb_full_metatable), pace_metabolics$subject_name), c("eGFR", "creat_b")]


# save the data
metadata_with_clusters <- cbind(pace_wb_metatable, clusters)
metadata_full_with_clusters <- cbind(pace_wb_full_metatable, clusters)
write.csv(metadata_with_clusters, file.path(outdir, "metadata_with_clusters.csv"))

with(pace_wb_full_metatable, table(censor_ANYEVENT, censor_MACLE2, useNA = "always"))
pace_wb_full_metatable %>%
    filter(censor_MACE == 1) %>%
    select(starts_with("censor_"))

# --------------------------------- Read in our metatable and bioreptable ---------------------------------
isch_genestocolorstab_file <- "output/run4_rmoutliers2_asr_control/WGCNA/WGCNA_power14_size30/wgcna_genestocolors.csv"
isch_genestocolorstab <- read.table(isch_genestocolorstab_file, ",", header = TRUE, row.names = 1)
eigengene_order <- paste0(
    "ME",
    c(
        "salmon", "magenta", "purple", "midnightblue", "turquoise",
        "orange", "pink", "red",
        "blue", "tan", "brown",
        "greenyellow", "black", "cyan", "green"
    )
)

SOI <- unique(c(rownames(pace_wb_metatable[pace_wb_metatable[, "Cohort"] %in% "PACE", , drop = FALSE]), rownames(na.omit(pace_wb_metatable[, "comp_THRGrOr_PACEg4_v_THRg4", drop = FALSE]))))

# Let's get the eigengenes for the PACE data
sharedgenes <- intersect(rownames(isch_genestocolorstab), rownames(pace_wb_normcounttable))
pace_eigengenetable <- orderMEs(moduleEigengenes(t(pace_wb_normcounttable[sharedgenes, SOI]), isch_genestocolorstab[sharedgenes, "moduleColors"])$eigengenes)
pace_eigengenetable[, "MEorange"] <- NA

pace_eigengenetable_ordered <- pace_eigengenetable[, eigengene_order]
write.csv(pace_eigengenetable_ordered, file.path(outdir, "pace_eigengenetable_ordered.csv"))
colnames(pace_eigengenetable_ordered) <- gsub("ME", "", colnames(pace_eigengenetable_ordered))

# Heatmap of the eigengenes
top_anno <- HeatmapAnnotation(
    df = metadata_with_clusters %>% select(pace_predict_labels),
    col = list(
        pace_predict_labels = c("nmf_cluster_1" = "blue", "nmf_cluster_2" = "red", "nmf_cluster_3" = "green", "nmf_cluster_4" = "purple")
    ),
    annotation_label = "",
    show_legend = FALSE
)
heatmap <- Heatmap(
    t(pace_eigengenetable_ordered),
    name = "Eigengene",
    top_annotation = top_anno,
    show_row_names = TRUE,
    row_names_side = "left",
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    cluster_column_slices = FALSE,
    column_split = metadata_with_clusters$pace_predict_labels,
    column_title = c("RS1", "RS2", "RS3", "RS4"),
    width = unit(14, "cm"),
    height = unit(10, "cm")
)
pdf(file.path(outdir, "pace_eigengene_heatmap.pdf"), width = 10, height = 10)
draw(heatmap)
dev.off()

# summarize the expression of the eigengenes by cluster
dat <- cbind(pace_eigengenetable_ordered, metadata_with_clusters)
eigengene_avg_expr_stats <- dat %>%
    select(all_of(gsub("ME", "", eigengene_order)), pace_predict_labels) %>%
    gather(key = "eigengene", value = "value", -pace_predict_labels) %>%
    group_by(pace_predict_labels, eigengene) %>%
    summarize(
        mean = median(value, na.rm = TRUE),
        sd = sd(value, na.rm = TRUE),
        n = sum(!is.na(value))
    )
eigengene_order_avg <- eigengene_avg_expr_stats %>%
    select(pace_predict_labels, eigengene, mean) %>%
    pivot_wider(names_from = pace_predict_labels, values_from = mean) %>%
    column_to_rownames("eigengene")

eigengene_order_ranked <- apply(eigengene_order_avg, 1, function(x) factor(rank(x, ties.method = "min")))
eigengene_order_ranked <- t(eigengene_order_ranked)[gsub("ME", "", eigengene_order), ]

colors = structure(c("#E5E5E5", "#A1A0A0", "#5E5E5E", "#181818"), names = c("1", "2", "3", "4"))
order_hm <- Heatmap(
    eigengene_order_ranked,
    name = "Rank",
    border = TRUE,
    rect_gp = gpar(col = "black"),
    col = colors,
    show_row_names = TRUE,
    row_names_side = "left",
    show_column_names = TRUE,
    column_labels = c("RS1", "RS2", "RS3", "RS4"),
    column_names_side = "top",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    width = unit(4, "cm"),
    height = unit(10, "cm")
)
# Let's do them both here
pdf(file.path(outdir, "pace_eigengene_heatmap_order.pdf"), width = 10, height = 10)
draw(heatmap + order_hm)
dev.off()

# ======================== Clustering Stats Table ========================#
dir.create(file.path(outdir, "clustering_stats"), showWarnings = F)
# let's make a clustering stats table
colnames(metadata_full_with_clusters)
colnames(metadata_with_clusters)
voi <- c(
    # demographics
    "age", "sex", "race", "ethnicity",
    # risk factors
    "bmi", "smoking", "diabetes", "hypertension", "CAD", "bmi",
    # labs
    "ldl_b_n", "hdl_b_n", "hba1c_b_n", "tchol_b_n", "tg_b_n",
    "eGFR", "creatinine"
)
clust_stats <- stats_table(metadata_with_clusters, "pace_predict_labels", voi)
clust_full_stats <- stats_table(metadata_full_with_clusters, "pace_predict_labels", voi)
write.csv(clust_stats, file.path(outdir, "clustering_stats", "pace_clust_stats.csv"))
write.csv(clust_full_stats, file.path(outdir, "clustering_stats", "pace_clust_full_stats.csv"))
censors <- c("censor_death", "censor_ampu", "censor_MACE", "censor_MALE2", "censor_ANYEVENT")
clust_events <- stats_table(metadata_full_with_clusters, "pace_predict_labels", censors)
write.csv(clust_events, file.path(outdir, "clustering_stats", "pace_clust_events.csv"))

# make a barplot of the clusters for the events
censors <- grep("censor_", colnames(metadata_full_with_clusters), value = TRUE)
tmpDat <- metadata_full_with_clusters[, c("pace_predict_labels", censors)]
tmpDat %<>% mutate_at(censors, ~ factor(., levels = c(0, 1), labels = c("Censored", "Event")))

clust_events <- stats_table(tmpDat, "pace_predict_labels", censors)
write.csv(clust_events, file.path(outdir, "clustering_stats", "pace_clust_events.csv"))

tmpDat_long <- gather(tmpDat, key = "event", value = "censor", -pace_predict_labels) %>%
    mutate(
        pace_predict_labels = factor(pace_predict_labels, labels = c("RS1", "RS2", "RS3", "RS4")),
        event = factor(
            case_when(
                event == "censor_death" ~ "Death",
                event == "censor_ampu" ~ "Amputation",
                event == "censor_MACE" ~ "MACE",
                event == "censor_MALE2" ~ "MALE",
                event == "censor_MACLE2" ~ "MACLE",
                event == "censor_ANYEVENT" ~ "Any Event",
                TRUE ~ NA_character_
            ),
            levels = c("Death", "Amputation", "MACE", "MALE", "MACLE", "Any Event")
        )
    )
p <- ggplot(tmpDat_long, aes(x = pace_predict_labels, fill = censor)) +
    geom_bar(position = "fill") +
    facet_wrap(~event, ncol = 5) +
    theme_bw() +
    theme(
        legend.position = "none",
        panel.grid.major = element_blank()
    ) +
    labs(x = NULL, y = "Proportion", title = NULL) +
    scale_fill_manual(values = c("Censored" = "white", "Event" = "red"))
ggsave(file.path(outdir, "clustering_stats", "pace_clust_events.pdf"), p, width = 12, height = 3)
p


# ======================== END ========================#
writeLines(capture.output(sessionInfo()), file.path(outdir, "session.log"))
save.image(file.path(outdir, "image.RData"))
