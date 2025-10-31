###########################################################################
#
#                       PACE Validation
#
###########################################################################
# Author: ISCHEMIA Study Team
# Date: Updated 2025-08-08
# Description: PACE Validation analysis for ISCHEMIA study
#

# Load configuration and utilities
source("config.R")
source("utils.R")

# Load required packages
# TODO: Update this list based on actual packages used in the script
# load_packages(c("package1", "package2"))

source("../code/overlap_finder_function.R")
source("../code/summarize_table_function.R")
source("code/ischemia2021_exploratory_analysis_functions.R") 
source("../code/cohort_matching_functions.R")
source("../code/rnaseq_scripts/deseq_functions.R")

library(survival) #nolint
library(survminer) #nolint

# Outfolder
outfilepathmaster <- "output/run4_rmoutliers2_asr_control/validation/PACE/run4_mimic_MM_edits2/"
dir.create(outfilepathmaster, recursive = TRUE, showWarnings = FALSE)

# save.image(file = paste0(outfilepathmaster, "/ischemia_exploratory.RData"))
# load(file = paste0(outfilepathmaster, "/ischemia_exploratory.RData"))

# --------------------------------- Read in the ischemia genemodule information ---------------------------------
isch_genestocolorstab_file <- "output/run3_rmoutliers2/WGCNA/WGCNA_power14_size30/wgcna_genestocolors.csv"
isch_genestocolorstab_file <- "data/WGCNA_power14_size30/wgcna_genestocolors.csv"
isch_genestocolorstab <- read.table(isch_genestocolorstab_file, ",", header = TRUE, row.names = 1)
eigengene_order <- paste0("ME",
                     c("magenta", "salmon", "purple", "midnightblue", "turquoise", 
                       "yellow", "pink", "red", 
                       "blue", "tan", "brown",
                       "greenyellow", "black", "cyan", "green"))
colorguidefile <- "data/ischemia2021_colorguide.csv"
colorguide <- read.table(colorguidefile, sep = ",", header = TRUE, comment.char = "", colClasses = c("character", "character", "character"))



# --------------------------------- Create PACE Eigengene and Singscore tables ---------------------------------
## PACE WB counttable and metatable
# pace_wb_normcounttable <- read.table("../../projects/newman-pace-2020/output_hg19_newcontrols/run1_nooutlier1/rna_processing/normcounttab.txt", header = TRUE, row.names = 1, sep = "\t")
pace_wb_normcounttable <- read.table("data/validation_data/pace_grab/normcounttab.txt", header = TRUE, row.names = 1, sep = "\t")
paceONLY_wb_normcounttable <- pace_wb_normcounttable[,grepl("PACE", colnames(pace_wb_normcounttable))]
thrONLY_wb_normcounttable <- pace_wb_normcounttable[,grepl("THR", colnames(pace_wb_normcounttable))]

### USE THIS FOR 
# pace_wb_metatable <- read.table("../../projects/newman-pace-2020/output_hg19_newcontrols/run1_nooutlier1/rna_processing/metatable_in.csv", header = TRUE, row.names = 1, sep = ",")
# pace_wb_full_metatable <- read.table("../../projects/newman-pace-2020/data/pace_comb_allmetadata_20201014.csv",
#                                      sep = ",", header = TRUE, row.names = 1)
pace_wb_metatable <- read.table("data/validation_data/pace_grab/metatable_in.csv", header = TRUE, row.names = 1, sep = ",")
pace_wb_full_metatable <- read.table("data/validation_data/pace_grab/pace_comb_allmetadata_20201014.csv",
                                     sep = ",", header = TRUE, row.names = 1)
# Add Age cat on to this for later calculations # 24-52 53-63 64-91
pace_wb_full_metatable[,"age_cat"] <- ifelse(pace_wb_full_metatable[,"current_age"] >= 71 & 
                                          rownames(pace_wb_full_metatable) %in% colnames(paceONLY_wb_normcounttable), "71-97",
                                      ifelse(pace_wb_full_metatable[,"current_age"] >= 60 & 
                                          rownames(pace_wb_full_metatable) %in% colnames(paceONLY_wb_normcounttable), "61-70",
                                      ifelse(pace_wb_full_metatable[,"current_age"] >= 25 & 
                                          rownames(pace_wb_full_metatable) %in% colnames(paceONLY_wb_normcounttable), "25-60", NA)))

# We want to calc an ANY EVENT column: (should be MACLE2 effectively, but lets just double check) - its slightly different, so good to have
# pace_wb_full_metatable[,"censor_ANYEVENT"] <- ifelse(rowSums(pace_wb_full_metatable[,grepl("censor_", colnames(pace_wb_full_metatable)) & !grepl("_30", colnames(pace_wb_full_metatable))]) == 0, 0, 1)
pace_wb_full_metatable[,c("censor_ANYEVENT", "time_to_ANYEVENT")] <- t(apply(pace_wb_full_metatable[,
    grepl("censor_|time_to", colnames(pace_wb_full_metatable)) & !grepl("_30", colnames(pace_wb_full_metatable))], 1, function(x) {
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
pace_wb_full_metatable <- pace_wb_full_metatable %>% ## ADDED BY MATTHEW MULLER 2024-03-05
    mutate(
        censor_MALE2Death = ifelse(censor_MALE2 == 1 | censor_death == 1, 1, 0),
        time_to_MALE2Death = pmin(time_to_MALE2, time_to_death, na.rm = TRUE),
        censor_ANYAMPU = ifelse(censor_ampu == 1 | censor_minor_ampu == 1, 1, 0),
        time_to_ANYAMPU = pmin(time_to_ampu, time_to_minor_ampu, na.rm = TRUE),
        censor_ANYREOP = ifelse(censor_reop == 1 | censor_minor_reop == 1, 1, 0),
        time_to_ANYREOP = pmin(time_to_reop, time_to_minor_reop, na.rm = TRUE)
    )

## These were the SOI from the PACE Paper
SOI <- unique(c(rownames(pace_wb_metatable[pace_wb_metatable[,"Cohort"] %in% "PACE",,drop=FALSE]),
                rownames(na.omit(pace_wb_metatable[,"comp_THRGrOr_PACEg4_v_THRg4",drop=FALSE]))))

## Need quick heatmaps for this:
colmetatable <- pace_wb_metatable[SOI, "Cohort", drop = FALSE]
colmetatable[,"Cohort"] <- ifelse(colmetatable[,"Cohort"] %in% "PACE", "PACE", "Control")
colannotationlist = annotationlist_builder(colmetatable)
heatmapcolorparam = colorRamp2(breaks = c(-4, 0, 4), colors = c("blue", "white", "red"))
hmout1 <- create_heatmap(counttab = pace_wb_normcounttable[,SOI],
                         colmetatable = colmetatable, colannotationlist = colannotationlist,
                         colclusterparam = TRUE, rowclusterparam = TRUE,
                         columnsplitparam = colmetatable,
                         heatmapcolorparam = heatmapcolorparam, separate_legend = TRUE)
pdf(paste0(outfilepathmaster, "pace_wb_normcounttable_hm.pdf"), 10, 10, useDingbats = FALSE)
draw(hmout1[[1]])
junk <- dev.off()

hmout2 <- create_heatmap(counttab = paceONLY_wb_normcounttable[,SOI[SOI %in% colnames(paceONLY_wb_normcounttable)]],
                         colmetatable = NULL, colannotationlist = NULL,
                         colclusterparam = TRUE, rowclusterparam = TRUE,
                         heatmapcolorparam = heatmapcolorparam, separate_legend = TRUE)
pdf(paste0(outfilepathmaster, "paceONLY_wb_normcounttable.pdf"), 10, 10, useDingbats = FALSE)
draw(hmout2[[1]])
junk <- dev.off()


# (A) Looking at modules:
sharedgenes <- intersect(rownames(isch_genestocolorstab), rownames(pace_wb_normcounttable))
pace_eigengenetable <- orderMEs(moduleEigengenes(t(pace_wb_normcounttable[sharedgenes,SOI]),
                                                 isch_genestocolorstab[sharedgenes,"moduleColors"])$eigengenes)
# (B) Looking at singscore scores
# Looking at singscore values: ## TOSH WORKING HERE - ADD THIS TO ABOVE TABLE AND RUN IT ALL AT ONCE!
## NEED TOO ADD THE DIFF SINGSCORE GENESETS HERE! - read in all combos of rna_subtypes and numbers here:
singscore_table_list <- singscore_genelists_list <- list()
singscore_path <- paste0("output/run4_rmoutliers2_asr_control/integrative_analyses/run7/signature_scoring/")
singscore_list <- c("rna_subtype1", "rna_subtype2", "rna_subtype3", "rna_subtype3A", "rna_subtype3B", "rna_subtype4", "rna_subtype3Bv3A")
scoring_genesetsize = c("all", 1000, 500, 100, 50)
singscore_grid <- expand.grid(singscore_list, scoring_genesetsize)
for (singscore_param_num in seq_len(nrow(singscore_grid))) {
    singscore_result <- singscore_grid[singscore_param_num,"Var1"]
    genesetsize <- singscore_grid[singscore_param_num,"Var2"]
    outlabel <- paste0(singscore_result, "__", genesetsize)

    singscore_basepath <- paste0(singscore_path, singscore_result, "/singscore/")
    singscore_subpath <- list.files(singscore_basepath)[grepl(paste0(genesetsize, "UP"), list.files(singscore_basepath)) | grepl(paste0(genesetsize, "DOWN"), list.files(singscore_basepath))]
    upfilename <- list.files(paste0(singscore_basepath, singscore_subpath, "/"))[grepl("GOIup", list.files(paste0(singscore_basepath, singscore_subpath, "/")))]
    downfilename <- list.files(paste0(singscore_basepath, singscore_subpath, "/"))[grepl("GOIdown", list.files(paste0(singscore_basepath, singscore_subpath, "/")))]

    singscore_genelists_list[[outlabel]][["Upset"]] <- read.table(paste0(singscore_basepath, singscore_subpath, "/", upfilename), sep = ",", header = TRUE)[,1]
    singscore_genelists_list[[outlabel]][["Downset"]] <- read.table(paste0(singscore_basepath, singscore_subpath, "/", downfilename), sep = ",", header = TRUE)[,1]
}

## For every score - read in the sets, and then apply singscore to get a value for every sample in our PACE cohort - repeat for all cluster scores
# pace_ss_table_list <- list()
# for (singscore_genelist_num in seq_along(singscore_genelists_list)) {
#     ss_upset <- singscore_genelists_list[[singscore_genelist_num]][["Upset"]]
#     ss_downset <- singscore_genelists_list[[singscore_genelist_num]][["Downset"]]
#     ss_label <- names(singscore_genelists_list)[singscore_genelist_num]
#     # singscore_table <- singscore_sample_scoring(GOIup=ss_upset, GOIdown=ss_downset, counttable=pace_wb_normcounttable)
#     singscore_table <- singscore_sample_scoring(GOIup=ss_upset, GOIdown=ss_downset, counttable=pace_wb_normcounttable)
#     pace_ss_table_list[[singscore_genelist_num]] <- cbind(Row.names = rownames(singscore_table), singscore_table)
#     names(pace_ss_table_list)[singscore_genelist_num] <- ss_label
# }
# pace_combined_singscore_table <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Row.names", all = TRUE, sort = FALSE),
#                                         lapply(pace_ss_table_list, function(x) x[,c("Row.names", "TotalScore", "UpScore", "DownScore")]))
# dimnames(pace_combined_singscore_table) <- list(pace_combined_singscore_table[,"Row.names"], c("Row.names", names(pace_ss_table_list)))
# pace_combined_singscore_table <- pace_combined_singscore_table[,!grepl("Row.names", colnames(pace_combined_singscore_table))]
pace_ss_table_list <- list()
for (singscore_genelist_num in seq_along(singscore_genelists_list)) {
    ss_upset <- singscore_genelists_list[[singscore_genelist_num]][["Upset"]]
    ss_downset <- singscore_genelists_list[[singscore_genelist_num]][["Downset"]]
    ss_label <- names(singscore_genelists_list)[singscore_genelist_num]
    singscore_table <- singscore_sample_scoring(GOIup=ss_upset, GOIdown=ss_downset, counttable=pace_wb_normcounttable)
    # singscore_table <- singscore_sample_scoring(GOIup=ss_upset, GOIdown=ss_downset, counttable=paceONLY_wb_normcounttable)
    pace_ss_table_list[[singscore_genelist_num]] <- cbind(Row.names = rownames(singscore_table), singscore_table)
    names(pace_ss_table_list)[singscore_genelist_num] <- ss_label
}
pace_combined_singscore_table <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Row.names", all = TRUE, sort = FALSE), 
                                        lapply(pace_ss_table_list, function(x) x[,c("Row.names", "TotalScore", "UpScore", "DownScore")]))
# dimnames(pace_combined_singscore_table) <- list(pace_combined_singscore_table[,"Row.names"], c("Row.names", names(pace_ss_table_list)))
dimnames(pace_combined_singscore_table) <- list(pace_combined_singscore_table[,"Row.names"], 
                                                c("Row.names", apply(expand.grid(paste0(c("TotalScore", "UpScore", "DownScore"), "__"), names(singscore_genelists_list)), 1, paste0, collapse = "")))
pace_combined_singscore_table <- pace_combined_singscore_table[,!grepl("Row.names", colnames(pace_combined_singscore_table))]
# write.table(pace_combined_singscore_table, paste0(subanalysis_outfilepath, "pace_", genesetsize, scoring_type_label, "_singscore_table.csv"),
#             sep = ",", col.names = NA, row.names = TRUE)


# --------------------------------- Viz our eigengene and singscore tables ---------------------------------
## So for every sample - do a simple bar plot for every sample of singscore vaues



# --------------------------------- eigen_v_annot ---------------------------------
dir.create(paste0(outfilepathmaster, "eigen_v_annot/"), showWarnings = FALSE, recursive = TRUE)

# Also create a singscore table for our various singscores? ## TOSH WORKING HERE
boxplot_intable_temp <- merge(pace_eigengenetable[SOI,eigengene_order], 
                              pace_combined_singscore_table[SOI,
        grepl(paste(c("rna_subtype1__", "rna_subtype2__", "rna_subtype3__", "rna_subtype4__"), collapse = "|"), colnames(pace_combined_singscore_table))],
    by = "row.names", all.x = TRUE)
boxplot_intable <- data.frame(t(data.frame(boxplot_intable_temp[,!grepl("Row.names", colnames(boxplot_intable_temp))], row.names = SOI)))
eigen_and_singscore_order <- colnames(boxplot_intable_temp[,2:ncol(boxplot_intable_temp)])

# Create filtered event table
# Function for taking events with continuous time to event, and creating new columns filtered by day cutoffs - censoring appropriately
toe_filter_inlist <- list(
    death = pace_wb_full_metatable[,c("time_to_death", "censor_death")], MI = pace_wb_full_metatable[,c("time_to_MI", "censor_MI")],
    stroke = pace_wb_full_metatable[,c("time_to_stroke", "censor_stroke")], ampu = pace_wb_full_metatable[,c("time_to_ampu", "censor_ampu")],
    minor_ampu = pace_wb_full_metatable[,c("time_to_minor_ampu", "censor_minor_ampu")], DeathMI = pace_wb_full_metatable[,c("time_to_DeathMI", "censor_DeathMI")],
    MALEDeath = pace_wb_full_metatable[,c("time_to_MALEDeath", "censor_MALEDeath")], MACEampu = pace_wb_full_metatable[,c("time_to_MACEampu", "censor_MACEampu")],
    ampuDeath = pace_wb_full_metatable[,c("time_to_ampuDeath", "censor_ampuDeath")], MACE = pace_wb_full_metatable[,c("time_to_MACE", "censor_MACE")],
    MALE = pace_wb_full_metatable[,c("time_to_MALE", "censor_MALE")], MALE2 = pace_wb_full_metatable[,c("time_to_MALE2", "censor_MALE2")],
    MACLE = pace_wb_full_metatable[,c("time_to_MACLE", "censor_MACLE")], MACLE2 = pace_wb_full_metatable[,c("time_to_MACLE2", "censor_MACLE2")],
    ANYEVENT = pace_wb_full_metatable[,c("time_to_ANYEVENT", "censor_ANYEVENT")], MALE2Death = pace_wb_full_metatable[,c("time_to_MALE2Death", "censor_MALE2Death")],
    ANYAMPU = pace_wb_full_metatable[,c("time_to_ANYAMPU", "censor_ANYAMPU")], ANYREOP = pace_wb_full_metatable[,c("time_to_ANYREOP", "censor_ANYREOP")]
)

time_filters <- c(182.625 * c(1,2,3,4,6))
toe_table <- pace_wb_full_metatable[,grepl("censor_|time_to", colnames(pace_wb_full_metatable)) & !grepl("_30", colnames(pace_wb_full_metatable))]
toe_table_outlist <- list()
for (toe_table in toe_filter_inlist) {
    colnames(toe_table) <- gsub("censor_", "C_", colnames(toe_table))
    colnames(toe_table) <- gsub("time_to_", "T_", colnames(toe_table))
    filtered_toe_table <- filter_time_to_event_table(time_to_event_table=toe_table, time_filters=time_filters)
    toe_label <- gsub("C_", "", colnames(toe_table)[1])
    toe_table_outlist[[toe_label]] <- filtered_toe_table
}
EOI_table <- do.call(cbind.data.frame, toe_table_outlist)
colnames(EOI_table) <- unlist(lapply(toe_table_outlist, colnames))


# Create the annotation table that we need
## Combine our custom filtered event table and metatable
boxplot_annotationtable <- merge(EOI_table[,grepl("filtered", colnames(EOI_table))], pace_wb_metatable[!grepl("censor_|time_to_", colnames(pace_wb_metatable))], by = "row.names", all.y = TRUE)
rownames(boxplot_annotationtable) <- boxplot_annotationtable[,"Row.names"]

## First annotation table - is Cohort or any other metadata annotations
boxplot_annotationtable_meta <- boxplot_annotationtable[SOI, c("Cohort", "sex", "age_cat"), drop = FALSE]
boxplot_annotationtable_meta[,"Cohort"] <- factor(ifelse(boxplot_annotationtable_meta[,"Cohort"] %in% "PACE", "PACE", "THR"), levels = c("PACE", "THR"))
# OHE non binary variables
if (sum(apply(boxplot_annotationtable_meta, 2, function(x) length(unique(x)) > 2) > 0)) {
    newcols <- data.frame(make_comparison_columns(
        boxplot_annotationtable_meta[,apply(boxplot_annotationtable_meta, 2, function(x) length(unique(x)) > 2), drop = FALSE]))
    newcols[] <- lapply(newcols, function(x) {
        factor(x, levels = c(unique(x)[!grepl("not_", unique(x))], unique(x)[grepl("not_", unique(x))]))
    })
    boxplot_annotationtable_meta <- cbind(boxplot_annotationtable_meta[,!apply(boxplot_annotationtable_meta, 2, function(x) length(unique(x)) > 2),drop = FALSE], newcols[rownames(boxplot_annotationtable_meta),])
}

## Second is our custom filtered event comp columns
boxplot_annotationtable_filtevent <- boxplot_annotationtable[SOI, c(colnames(boxplot_annotationtable)[grepl("filtered", colnames(boxplot_annotationtable))])]
w1 <- which(boxplot_annotationtable_filtevent==1,arr.ind=TRUE)
w2 <- which(boxplot_annotationtable_filtevent==2,arr.ind=TRUE)
boxplot_annotationtable_filtevent[w1] <- names(boxplot_annotationtable_filtevent)[w1[,"col"]]
boxplot_annotationtable_filtevent[w2] <- NA
boxplot_annotationtable_filtevent[] <- lapply(boxplot_annotationtable_filtevent, function(x) {
    factor(x, levels = c(unique(x)[unique(x) != "0"], unique(x)[unique(x) == "0"]))
})

## Third is custom filtered events, but the no events are all the same of absolutely NO event in that time frame
outlist <- list()
for (time_filter_sel in time_filters) {
    grabtab <- boxplot_annotationtable_filtevent[,grepl(time_filter_sel, colnames(boxplot_annotationtable_filtevent))]
    # So the key here is to turn the 0s for each column into the 0s for the ANYEVENT column
    grabtab[] <- apply(grabtab, 2, function(x) ifelse(!x %in% c(0, NA), x, ifelse(grabtab[,grepl("_ANYEVENT", colnames(grabtab))] == 0, 0, NA)))
    outlist[[paste0("time__", time_filter_sel)]] <- grabtab
}
boxplot_annotationtable_filteventVnoneevent <- do.call(cbind, unname(outlist))
boxplot_annotationtable_filteventVnoneevent[] <- lapply(boxplot_annotationtable_filteventVnoneevent, function(x) {
    treat = unique(x)[unique(x) != "0"]
    control = unique(x)[unique(x) == "0"]
    factor(x, levels = c(treat, control))
})

## Fourth - I can do my matched PACE THR analyses as well...
## 4a - match within PACE vs noevent
MATCH_boxplot_annotationtable_filtevent <- data.frame(apply(boxplot_annotationtable_filtevent, 2, function(x) {

    treatlabel <- as.character(unique(x)[!unique(x) %in% c("0", NA)])
    ## Match out PACE and THR samples
    matchtable <- na.omit(pace_wb_metatable[c(names(x[x %in% treatlabel]), names(x[x %in% "0"])), c("age", "sex")])
    matchtable[,"matchcohort"] <- ifelse(rownames(matchtable) %in% names(x[x %in% treatlabel]), "PACE_EVENT", "PACE_NOEVENT")
    ## First function creates a "matchtable" that will give you all of your possible matches
    set.seed(12345) ## NEED THIS TO PRESERVE MATCHING BETWEEN RUNS... OOPS ## BEST
    out1 <- create_match_table(intable = matchtable, groupcat = "matchcohort", treatvar = "PACE_EVENT", controlvar = "PACE_NOEVENT",
                               controlcat = c("age", "sex"), discretevar = "sex", discrete_DF = 0, contvar = list(age = c(var = "age", range = 5)))
    ## This function helps you select the matches from your table # rownames - treat, colnames - control
    out2 <- create_matches_v2(matchtab = out1, control_to_treat_ratio = 1, seedparam = 12345)$final_fullratiomatch_table
    out4 <- data.frame(factor(c(rep(treatlabel, length(rownames(out2))), rep(0, length(colnames(out2))), 
                                rep(NA, (length(x) - length(rownames(out2)) - length(colnames(out2))))), levels = c(treatlabel, 0)), 
                       row.names = c(rownames(out2), colnames(out2), names(x)[!names(x) %in% c(rownames(out2), colnames(out2))]))[names(x),]
    out4
}))
dimnames(MATCH_boxplot_annotationtable_filtevent) <- list(rownames(boxplot_annotationtable_filtevent), paste0("MATCHED_", colnames(boxplot_annotationtable_filtevent)))
MATCH_boxplot_annotationtable_filtevent[] <- lapply(MATCH_boxplot_annotationtable_filtevent, function(x) {
    factor(x, levels = c(unique(x)[unique(x) != "0"], unique(x)[unique(x) == "0"]))
})

## 4b - match within PACE vs NONEevent
MATCH_boxplot_annotationtable_filteventVnoneevent <- data.frame(apply(boxplot_annotationtable_filteventVnoneevent, 2, function(x) {
    treatlabel <- as.character(unique(x)[!unique(x) %in% c("0", NA)])
    print(treatlabel)
    
    ## Match out PACE and THR samples
    matchtable <- na.omit(pace_wb_metatable[c(names(x[x %in% treatlabel]), names(x[x %in% "0"])), c("age", "sex")])
    matchtable[,"matchcohort"] <- ifelse(rownames(matchtable) %in% names(x[x %in% treatlabel]), "PACE_EVENT", "PACE_NOEVENT")
    ## First function creates a "matchtable" that will give you all of your possible matches
    set.seed(12345) ## NEED THIS TO PRESERVE MATCHING BETWEEN RUNS... OOPS ## BEST
    out1 <- create_match_table(intable = matchtable, groupcat = "matchcohort", treatvar = "PACE_EVENT", controlvar = "PACE_NOEVENT",
                               controlcat = c("age", "sex"), discretevar = "sex", discrete_DF = 0, contvar = list(age = c(var = "age", range = 5)))
    ## This function helps you select the matches from your table # rownames - treat, colnames - control
    out2 <- create_matches_v2(matchtab = out1, control_to_treat_ratio = 1, seedparam = 12345)$final_fullratiomatch_table
    out4 <- data.frame(factor(c(rep(treatlabel, length(rownames(out2))), rep(0, length(colnames(out2))), 
                                rep(NA, (length(x) - length(rownames(out2)) - length(colnames(out2))))), levels = c(treatlabel, 0)), 
                       row.names = c(rownames(out2), colnames(out2), names(x)[!names(x) %in% c(rownames(out2), colnames(out2))]))[names(x),]
    out4
}))
dimnames(MATCH_boxplot_annotationtable_filteventVnoneevent) <- list(rownames(boxplot_annotationtable_filteventVnoneevent), paste0("MATCHED_", colnames(boxplot_annotationtable_filteventVnoneevent)))
MATCH_boxplot_annotationtable_filteventVnoneevent[] <- lapply(MATCH_boxplot_annotationtable_filteventVnoneevent, function(x) {
    factor(x, levels = c(unique(x)[unique(x) != "0"], unique(x)[unique(x) == "0"]))
})

## 4c - match between PACE and THR
MATCH_boxplot_annotationtable_filteventVthr <- data.frame(apply(boxplot_annotationtable_filtevent, 2, function(x) {
    treatlabel <- as.character(unique(x)[!unique(x) %in% c("0", NA)])
    
    ## Match out PACE and THR samples
    matchtable <- na.omit(pace_wb_metatable[c(names(x[x %in% treatlabel]), names(x[grepl("THR", names(x))])), c("age", "sex")])
    matchtable[,"matchcohort"] <- ifelse(rownames(matchtable) %in% names(x[x %in% treatlabel]), "PACE_EVENT", "THR")
    ## First function creates a "matchtable" that will give you all of your possible matches
    set.seed(12345) ## NEED THIS TO PRESERVE MATCHING BETWEEN RUNS... OOPS ## BEST
    out1 <- create_match_table(intable = matchtable, groupcat = "matchcohort", treatvar = "PACE_EVENT", controlvar = "THR",
                               controlcat = c("age", "sex"), discretevar = "sex", discrete_DF = 0, contvar = list(age = c(var = "age", range = 5)))
    ## This function helps you select the matches from your table # rownames - treat, colnames - control
    out2 <- create_matches_v2(matchtab = out1, control_to_treat_ratio = 1, seedparam = 12345)$final_fullratiomatch_table
    out4 <- data.frame(factor(c(rep(treatlabel, length(rownames(out2))), rep(0, length(colnames(out2))), 
                                rep(NA, (length(x) - length(rownames(out2)) - length(colnames(out2))))), levels = c(treatlabel, 0)), 
                       row.names = c(rownames(out2), colnames(out2), names(x)[!names(x) %in% c(rownames(out2), colnames(out2))]))[names(x),]
    out4
}))
dimnames(MATCH_boxplot_annotationtable_filteventVthr) <- list(rownames(boxplot_annotationtable_filtevent), paste0("MATCHED_", colnames(boxplot_annotationtable_filtevent)))
MATCH_boxplot_annotationtable_filteventVthr[] <- lapply(MATCH_boxplot_annotationtable_filteventVthr, function(x) {
    factor(x, levels = c(unique(x)[unique(x) != "0"], unique(x)[unique(x) == "0"]))
})


# I want to do this twice - once with the event vs no event, and then again with event vs NONE event
## Three times? And do one with cohort, one with event vs no event, and then once more with event vs NONE event (then I can add additional annotation tables if I want them)
bp_annotationtable_list <- list(bp_annot_1 = list(label = "eigen_v_meta", annottable = boxplot_annotationtable_meta),
                                bp_annot_2 = list(label = "eigen_v_filtevent", annottable = boxplot_annotationtable_filtevent),
                                bp_annot_3 = list(label = "eigen_v_filteventVnoneevent", annottable = boxplot_annotationtable_filteventVnoneevent),
                                
                                bp_annot_4 = list(label = "eigen_v_MATCHfiltevent", annottable = MATCH_boxplot_annotationtable_filtevent),
                                bp_annot_5 = list(label = "eigen_v_MATCHfilteventVnoneevent", annottable = MATCH_boxplot_annotationtable_filteventVnoneevent),
                                bp_annot_6 = list(label = "eigen_v_MATCHfilteventVthr", annottable = MATCH_boxplot_annotationtable_filteventVthr)
                                )

for (bp_annotation in bp_annotationtable_list) {
    boxplot_annotationtable <- bp_annotation[["annottable"]]
    boxplot_label <- bp_annotation[["label"]]
    
    outfilepath <- paste0(outfilepathmaster, "eigen_v_annot/", boxplot_label, "/")
    dir.create(outfilepath, showWarnings = FALSE, recursive = TRUE)
    write.table(boxplot_annotationtable, paste0(outfilepathmaster, "eigen_v_annot/", boxplot_label, "/", boxplot_label, "_boxplot_annotationtable.csv"), sep = ",", col.names = NA, row.names = TRUE)
    
    ## Failsafe here:
    boxplot_annotationtable <- boxplot_annotationtable[,!colSums(is.na(boxplot_annotationtable)) == nrow(boxplot_annotationtable)]
    
    boxplots_from_counttable_by_annotation_out <- boxplots_from_counttable_by_annotation(counttable = boxplot_intable,
                                                                                         boxplot_annotationtable = boxplot_annotationtable, outfilepath = outfilepath, 
                                                                                         calculate_corrvalues = FALSE, bp_color_ref = NULL)
    bp_stat_table <- data.frame(boxplots_from_counttable_by_annotation_out[[1]])
    bp_stat_table[,"signed_wilcox"] <- as.numeric(bp_stat_table[,"wilcox_pval"]) * apply(bp_stat_table[,c("Feature1_mean", "Feature2_mean")], 1, function(x) ifelse(as.numeric(x[1]) > as.numeric(x[2]), 1, -1))
    write.table(bp_stat_table, paste0(outfilepathmaster, "eigen_v_annot/", boxplot_label, "/", boxplot_label, "_stattable.csv"),
                sep = ",", col.names = TRUE, row.names = FALSE)
    
    # Now I want a table to make a heatmap viz out of
    wilcox_pval_table <- reshape2::dcast(bp_stat_table, Cohort ~ Feature, value.var = "signed_wilcox")
    rownames(wilcox_pval_table) <- wilcox_pval_table[,"Cohort"]
    wilcox_pstat_table <- -log10(abs(wilcox_pval_table[colnames(boxplot_annotationtable),eigen_and_singscore_order]))
    wilcox_sign_table <- sign(wilcox_pval_table[colnames(boxplot_annotationtable),eigen_and_singscore_order])
    wilcox_plottable <- data.frame(t(wilcox_pstat_table * 
                                     wilcox_sign_table))[eigen_and_singscore_order, colnames(boxplot_annotationtable), drop = FALSE]
    
    write.table(wilcox_plottable[eigen_and_singscore_order,], 
                paste0(outfilepathmaster, "eigen_v_annot/", boxplot_label, "/", boxplot_label, "_plottable.csv"),
                sep = ",", col.names = NA, row.names = TRUE)
    
    # Plot out a heatmap - and add a reference for the eigengenes as being "risk" ones of "protective" ones
    rowmetatable <- eigenegene_to_cluster_table <- data.frame(eigengene = rownames(wilcox_plottable), 
        rna_4cluster = ifelse(rownames(wilcox_plottable) %in% c("MEmagenta", "MEsalmon", "MEpurple", "MEmidnightblue", "MEturquoise") |
                        grepl("rna_subtype1", rownames(wilcox_plottable)), "nmf_cluster_1",
                 ifelse(rownames(wilcox_plottable) %in% c("MEyellow", "MEpink", "MEred") |
                        grepl("rna_subtype2", rownames(wilcox_plottable)), "nmf_cluster_2", 
                 ifelse(rownames(wilcox_plottable) %in% c("MEblue", "MEtan") |
                        grepl("rna_subtype3", rownames(wilcox_plottable)), "nmf_cluster_3", 
                 ifelse(rownames(wilcox_plottable) %in% c("MEbrown", "MEgreenyellow", "MEblack", "MEcyan", "MEgreen") |
                        grepl("rna_subtype4", rownames(wilcox_plottable)), "nmf_cluster_4", NA)))), 
                                                              row.names = rownames(wilcox_plottable))
    
    rowmetatable <- rowmetatable[order(rowmetatable[,2]),,drop=FALSE]
    
    # genecustomcolorlist <- list(eigengene = gsub("ME", "", eigengene_order))
    # names(genecustomcolorlist[["eigengene"]]) <- eigengene_order
    
    genecustomcolorlist <- list(eigengene = ifelse(grepl("ME", rownames(rowmetatable)), gsub("ME", "", rownames(rowmetatable)), "white"))
    names(genecustomcolorlist[["eigengene"]]) <- rownames(rowmetatable)
    
    rowannotationlist <- c(genecustomcolorlist, custom_annotation_list_from_colorguide("rna_4cluster", colorguide = colorguide))
    
    heatmapcolorparam <- colorRamp2(breaks = c(-3, log10(0.05), 0, -log10(0.05), 3), c("blue", "white", "white", "white", "red"))
    hmout <- create_heatmap(counttab = wilcox_plottable[rownames(rowmetatable),,drop=FALSE], subsetnum = FALSE, scale_data = FALSE,
                            colmetatable = NULL, colannotationlist = NULL,
                            rowmetatable = rowmetatable, rowannotationlist = rowannotationlist,
                            addborders = TRUE, heatmapcolorparam = heatmapcolorparam
    )
    pdf(paste0(outfilepathmaster, "eigen_v_annot/", boxplot_label, "/", boxplot_label, "_heatmap.pdf"), 
        useDingbats = FALSE, width = 15, height = 10)
    draw(hmout[[1]])
    junk <- dev.off()
}




# --------------------------------- eigen_quant_km_analysis ---------------------------------
# (2) For those that we have time-to-event data for - do a time to event with quantiles of eigengene
dir.create(paste0(outfilepathmaster, "eigen_quant_km_analysis/"), showWarnings = FALSE, recursive = TRUE)

## Suppose I could also do singscore quants too huh? ## TOSH WORKING HERE
# boxplot_intable_temp <- merge(pace_eigengenetable[SOI,eigengene_order], pace_combined_singscore_table[SOI,], by = "row.names", all.x = TRUE)
# boxplot_intable <- data.frame(t(data.frame(boxplot_intable_temp[,!grepl("Row.names", colnames(boxplot_intable_temp))], row.names = SOI)))



param_table <- expand.grid(eigengene = eigengene_order, quantiles = c(2,3,4,5))
cluster_intable_list <- list()
for (param_num in seq_len(nrow(param_table))) {
    eigen_sel <- param_table[param_num, "eigengene"]
    quantile_sel <- param_table[param_num, "quantiles"]
    quantile_out <- continuous_to_named_quantile(pace_eigengenetable[grepl("PACE", rownames(pace_eigengenetable)),eigen_sel,drop=FALSE], 
                                                 number_of_quantiles = quantile_sel)[[1]]
    colnames(quantile_out) <- paste0(eigen_sel, "_Quant", quantile_sel)
    cluster_intable_list[[param_num]] <- quantile_out
}
cluster_intable <- do.call(cbind, cluster_intable_list)

pace_wb_full_metatable
event_intable <- pace_wb_full_metatable[rownames(cluster_intable), grepl("censor_|time_to_", colnames(pace_wb_full_metatable))]
event_intable <- event_intable[, !grepl("_30", colnames(event_intable))]
# event_intable <- pace_wb_metatable[grepl("PACE", rownames(pace_wb_metatable)), grepl("censor_|time_to_", colnames(pace_wb_metatable))]
event_intable <- data.frame(cbind(PATNUM = rownames(event_intable), event_intable))
colnames(event_intable) <- gsub("censor_", "C_", colnames(event_intable))
colnames(event_intable) <- gsub("time_to_", "T_", colnames(event_intable))

plottable_outlist <- list()
for (group_column_selected in colnames(cluster_intable)) {
    
    selected_cluster_table <- cbind(PATNUM = rownames(cluster_intable[order(cluster_intable[,group_column_selected]), group_column_selected,drop=FALSE]),
                                    cluster_intable[order(cluster_intable[,group_column_selected]), group_column_selected,drop=FALSE])
    
    # COIref_tablist <- list(
    #     # IMGDEGIS = bioreptable_waddons[,"IMGDEGIS", drop=FALSE],
    #     IMGDEGIS_01_3 = bioreptable_waddons[,"IMGDEGIS_01_3", drop=FALSE],
    #     # CTNDV50 = bioreptable_waddons[,"CTNDV50", drop=FALSE],
    #     # CTNDV50_13_clean = bioreptable_waddons[,"CTNDV50_13_clean", drop=FALSE],
    #     CTNDV70_13_clean = bioreptable_waddons[,"CTNDV70_13_clean", drop=FALSE],
    #     CKD = bioreptable_waddons[,"CKD",drop=FALSE],
    #     highlowage = highlowage
    # )
    
    # plotCOI <- c(names(COIref_tablist),
    #              paste0(group_column_selected, "_",
    #                     unique(na.omit(selected_cluster_table[,group_column_selected]))[order(unique(na.omit(selected_cluster_table[,group_column_selected])))])
    # )
    COIref_tablist <- NULL
    
    plotCOI <- paste0(group_column_selected, "_", unique(na.omit(selected_cluster_table[,group_column_selected]))[
        order(unique(na.omit(selected_cluster_table[,group_column_selected])))])
    
    survivalplot_outfilepath <- paste0(outfilepathmaster, "eigen_quant_km_analysis/", group_column_selected, "/")
    dir.create(survivalplot_outfilepath, showWarnings = FALSE, recursive = TRUE)
    summary_heatmap_out <- cluster_vs_event_kmanalysis_summary_heatmap(selected_cluster_table = selected_cluster_table,
                                                                   COIref_tablist = COIref_tablist,
                                                                   survivalplot_outfilepath = survivalplot_outfilepath, 
                                                                   bioreptable_waddons = event_intable,
                                                                   eventpvalcutoff = 0.05, return_output_tables = TRUE,
                                                                   plotCOI)
    pdf(paste0(outfilepathmaster, "eigen_quant_km_analysis/", group_column_selected, "/", group_column_selected, "_event_summary_heatmap.pdf"),
        12, 10, useDingbats = FALSE)
    draw(summary_heatmap_out[[1]])
    junk <- dev.off()
    
    plottable_outlist[[group_column_selected]] <- plottable_out <- summary_heatmap_out[[2]]
    
    
    # I also want to do the KM analysis with all of the groups together for each group
    # if (group_column_selected == "rna_3A_v_3B_only") { next }
    
    # survivalplot_outfilepath <- paste0(outfilepathmaster, "cluster_v_event_survival_NONBINARY/", group_column_selected, "/")
    # dir.create(survivalplot_outfilepath, showWarnings = FALSE, recursive = TRUE)
    # # COIref_tablist <- list(
    # #     CUSTOM_IMGDEGIS_NONEMILD = bioreptable_waddons[,"CUSTOM_IMGDEGIS_NONEMILD", drop=FALSE],
    # #     CTNDV70 = bioreptable_waddons[,"CTNDV70", drop=FALSE],
    # #     CKD = bioreptable_waddons[,"CKD",drop=FALSE],
    # #     age_cat = bioreptable_waddons[,"CAGE_RND",drop=FALSE]
    # # )
    # # plotCOI <- c(names(COIref_tablist), group_column_selected)
    # 
    # coxph_control_table <- factorize_metatable(bioreptable_waddons[,c("CAGE_RND", "SEX", "RACE")])
    # 
    # summary_heatmap <- cluster_NONBINARY_vs_event_kmanalysis_summary_heatmap(selected_cluster_table = selected_cluster_table,
    #                                                                          COIref_tablist = COIref_tablist,
    #                                                                          survivalplot_outfilepath = survivalplot_outfilepath, 
    #                                                                          bioreptable_waddons = bioreptable_waddons,
    #                                                                          eventpvalcutoff = 0.05,
    #                                                                          return_pairwise_coxph_values = TRUE, coxph_control_table = NULL,
    #                                                                          plotCOI)
    # pdf(paste0(outfilepathmaster, "cluster_v_event_survival_NONBINARY/", group_column_selected, "/", group_column_selected, "_event_summary_heatmap.pdf"),
    #     12, 10, useDingbats = FALSE)
    # draw(summary_heatmap)
    # junk <- dev.off()
    
}

## TOSH - NEED TO SUMMARIZE THIS BETTER

library(stringr)
library(stringi)



plottable_full <- do.call(cbind, plottable_outlist)
write.table(plottable_full, paste0(outfilepathmaster, "eigen_quant_km_analysis/", "eigen_quant_km_summary_table.csv"),
            sep = ",", col.names = NA, row.names = TRUE)
columnname_info <- data.frame(stri_split_fixed(gsub(".*\\.|_pval", "", colnames(plottable_full)), "_", 2, simplify = TRUE), 
                              row.names = colnames(plottable_full))
# colnames(columnname_info) <- c("eigengene", "quantile")
# columnname_info[,c(3,4)] <- stri_split_fixed(gsub(".*\\.|_pval", "", columnname_info[,2]), "_", 2, simplify = TRUE)

colmetatable_temp <- data.frame(cbind(eigengene = rep(eigengene_order, each = 14), quantile = unique(columnname_info[,2])))
colmetatable_temp[,c(3,4)] <- stri_split_fixed(gsub(".*\\.|_pval", "", colmetatable_temp[,2]), "_", 2, simplify = TRUE)
dimnames(colmetatable_temp) <- list(apply(colmetatable_temp, 1, function(x) paste0(x[1], "_", x[3], ".", x[1], "_", x[2], "_pval")), 
                                    c("eigengene", "quantile_label", "quantile_category", "HighLow"))

colmetatable <- colmetatable_temp[,c(1,3)]
colmetatable[,1] <- factor(colmetatable[,1], levels = eigengene_order)
colmetatable[,"HighLow"] <- ifelse(grepl("_Low_", rownames(colmetatable)), "Low", ifelse(grepl("_High_", rownames(colmetatable)), "High", NA))
colmetatable[,"rna_4cluster"] <- 
    ifelse(colmetatable[,"eigengene"] %in% c("MEmagenta", "MEsalmon", "MEpurple", "MEmidnightblue", "MEturquoise"), "nmf_cluster_1",
    ifelse(colmetatable[,"eigengene"] %in% c("MEyellow", "MEpink", "MEred"), "nmf_cluster_2", 
    ifelse(colmetatable[,"eigengene"] %in% c("MEblue", "MEtan"), "nmf_cluster_3", 
    ifelse(colmetatable[,"eigengene"] %in% c("MEbrown", "MEgreenyellow", "MEblack", "MEcyan", "MEgreen"), "nmf_cluster_4", NA))))
colmetatable <- colmetatable[,c("eigengene", "rna_4cluster", "quantile_category", "HighLow")]
write.table(plottable_full, paste0(outfilepathmaster, "eigen_quant_km_analysis/", "eigen_quant_km_summary_colmetatable.csv"),
            sep = ",", col.names = NA, row.names = TRUE)

genecustomcolorlist <-  list(eigengene = gsub("ME", "", eigengene_order))
names(genecustomcolorlist[["eigengene"]]) <- eigengene_order
# rowannotationlist <- c(genecustomcolorlist, custom_annotation_list_from_colorguide("rna_4cluster", colorguide = colorguide))
genecustomcolorlist <-  c(genecustomcolorlist, custom_annotation_list_from_colorguide("rna_4cluster", colorguide = colorguide))
# names(genecustomcolorlist[["eigengene"]]) <- eigengene_order
colannotationlist <- annotationlist_builder(colmetatable, customcolorlist = genecustomcolorlist)
# Reorder the counttable based off off of the new metatable
plottable_full <- plottable_full[,rownames(colmetatable)]

#       quantile = colmetatable[match(colmetatable[,2], unique(colmetatable[,2])),2])

# heatmapcolorparam <- colorRamp2(breaks = c(-1, log10(0.05) - 0.000001, log10(0.05), -0.000001, 0.000001, -log10(0.05), log10(0.05) + 0.000001, 1), 
#                                 colors = c("white", "white", "lightblue", "blue", "red", "lightred", "white", "white"))
eventpvalcutoff <- 0.05
heatmapcolorparam <- colorRamp2(c(-1, -eventpvalcutoff - 0.0001, -eventpvalcutoff, -0.0001, -1e-100, 0,
             1e-100, 0.0001, eventpvalcutoff, eventpvalcutoff + 0.0001, 1),
           c("grey", "grey", "#bcb2fd", "#11007d", "#11007d", "white",
             "#7d0000", "#7d0000", "#ffb2b2", "grey", "grey"))
hmout <- create_heatmap(counttab = plottable_full, colmetatable = colmetatable, colannotationlist = colannotationlist,
                        scale_data = FALSE, heatmapcolorparam = heatmapcolorparam, addborders = TRUE)
pdf(paste0(outfilepathmaster, "eigen_quant_km_analysis/", "eigen_quant_km_summary_heatmap.pdf"), 15, 10, useDingbats = FALSE)
draw(hmout[[1]])
junk <- dev.off()

# Ok but really - I want to know if the top and bottom quantiles are sig, so select for only top and bottom
# And then sort by top and bottom by eigengene, I think thats a cleaner picture
plottable_sel <- plottable_full[,grepl("_Low_|_High_", colnames(plottable_full))]
colmetatable_sel <- colmetatable[colnames(plottable_sel),]
colmetatable_sel <- colmetatable_sel[order(colmetatable_sel[,"HighLow"], colmetatable_sel[,"quantile_category"], colmetatable_sel[,"eigengene"]),]
plottable_sel<- plottable_sel[,rownames(colmetatable_sel)]
colannotationlist <- annotationlist_builder(colmetatable_sel, customcolorlist = genecustomcolorlist)
hmout_sel <- create_heatmap(counttab = plottable_sel, colmetatable = colmetatable_sel, colannotationlist = colannotationlist,
                        scale_data = FALSE, heatmapcolorparam = heatmapcolorparam, addborders = TRUE)
pdf(paste0(outfilepathmaster, "eigen_quant_km_analysis/", "eigen_quant_km_summary_heatmap_sel.pdf"), 15, 10, useDingbats = FALSE)
draw(hmout_sel[[1]])
junk <- dev.off()

# Select our hits and only plot those:
whenpos_prot_h <- rownames(colmetatable_sel[colmetatable_sel[,"rna_4cluster"] %in% c("nmf_cluster_1") & 
                                                  colmetatable_sel[,"HighLow"] %in% c("High"),])
whenpos_risk_h <- rownames(colmetatable_sel[colmetatable_sel[,"rna_4cluster"] %in% c("nmf_cluster_2", "nmf_cluster_4") & 
                                                  colmetatable_sel[,"HighLow"] %in% c("High"),])

whenpos_risk_l <- rownames(colmetatable_sel[colmetatable_sel[,"rna_4cluster"] %in% c("nmf_cluster_1") & 
                                                  colmetatable_sel[,"HighLow"] %in% c("Low"),])
whenpos_prot_l <- rownames(colmetatable_sel[colmetatable_sel[,"rna_4cluster"] %in% c("nmf_cluster_2", "nmf_cluster_4") & 
                                                  colmetatable_sel[,"HighLow"] %in% c("Low"),])


hittable <- plottable_sel
# hittable[,whenpos_prot_h] <- hittable[,whenpos_prot_h] < 0.05 & sign(hittable[,whenpos_prot_h]) == 1
# hittable[,whenpos_risk_h] <- hittable[,whenpos_risk_h] < 0.05 & sign(hittable[,whenpos_risk_h]) == 1
# hittable[,whenpos_risk_l] <- hittable[,whenpos_risk_l] > -0.05 & sign(hittable[,whenpos_risk_l]) == -1
# hittable[,whenpos_prot_l] <- hittable[,whenpos_prot_l] > -0.05 & sign(hittable[,whenpos_prot_l]) == -1
hittable[,whenpos_prot_h] <- hittable[,whenpos_prot_h] > -0.05 & sign(hittable[,whenpos_prot_h]) == -1
hittable[,whenpos_risk_h] <- hittable[,whenpos_risk_h] < 0.05 & sign(hittable[,whenpos_risk_h]) == 1
hittable[,whenpos_risk_l] <- hittable[,whenpos_risk_l] < 0.05 & sign(hittable[,whenpos_risk_l]) == 1
hittable[,whenpos_prot_l] <- hittable[,whenpos_prot_l] > -0.05 & sign(hittable[,whenpos_prot_l]) == -1
hittable[,!colnames(hittable) %in% c(whenpos_prot_h, whenpos_risk_h, whenpos_risk_l, whenpos_prot_l)] <- FALSE

plottable_sel_hits <- plottable_sel
plottable_sel_hits[!hittable] <- 1

hmout_sel_hits <- create_heatmap(counttab = plottable_sel_hits, colmetatable = colmetatable_sel, colannotationlist = colannotationlist,
                            scale_data = FALSE, heatmapcolorparam = heatmapcolorparam, addborders = TRUE)
pdf(paste0(outfilepathmaster, "eigen_quant_km_analysis/", "eigen_quant_km_summary_heatmap_sel_hits.pdf"), 15, 10, useDingbats = FALSE)
draw(hmout_sel_hits[[1]])
junk <- dev.off()

# Helper table to actually figure this out
## I need to remove some events - because I think they arent right


# genefinder_table <- merge(reshape2::melt(cbind(event = rownames(plottable_full), plottable_full)),
#                           colmetatable_temp, by.x = "variable", by.y = "row.names")
# genefinder_table <- merge(genefinder_table, eigenegene_to_cluster_table, by = "eigengene", all.x = TRUE)
# genefinder_table[,"clustertype"] <- ifelse(genefinder_table[,"rna_4cluster"] %in% "nmf_cluster_1", "protect", 
#                                     ifelse(genefinder_table[,"rna_4cluster"] %in% c("nmf_cluster_2", "nmf_cluster_4"), "risk", NA))
# genefinder_table <- genefinder_table[!genefinder_table[,3] %in% c("MACLE2", "noevent_CUST"),]
# 
# # Ok - the noevent_CUST variable is actually the opposite direction - so PURELY for the hit finder - Im going to flip the sign on all of the noevent_CUST vals
# # genefinder_table[genefinder_table[,"event"] == "noevent_CUST", "value"] <- -genefinder_table[genefinder_table[,"event"] == "noevent_CUST","value"]
# 
# genefinder_table[,"cluster_hit"] <- NA
# genefinder_table[(genefinder_table[,"rna_4cluster"] %in% "nmf_cluster_1" & abs(genefinder_table[,"value"]) < 0.05 &
#     ((grepl("High", genefinder_table[,"quantile"]) & sign(genefinder_table[,"value"]) == -1) | 
#      (grepl("Low", genefinder_table[,"quantile"]) & sign(genefinder_table[,"value"]) == 1))), "cluster_hit"] <- "RS1_hit"
# genefinder_table[(genefinder_table[,"rna_4cluster"] %in% "nmf_cluster_2" & abs(genefinder_table[,"value"]) < 0.05 &
#                       ((grepl("High", genefinder_table[,"quantile"]) & sign(genefinder_table[,"value"]) == 1) | 
#                            (grepl("Low", genefinder_table[,"quantile"]) & sign(genefinder_table[,"value"]) == -1))), "cluster_hit"] <- "RS2_hit"
# genefinder_table[(genefinder_table[,"rna_4cluster"] %in% "nmf_cluster_4" & abs(genefinder_table[,"value"]) < 0.05 &
#                       ((grepl("High", genefinder_table[,"quantile"]) & sign(genefinder_table[,"value"]) == 1) | 
#                            (grepl("Low", genefinder_table[,"quantile"]) & sign(genefinder_table[,"value"]) == -1))), "cluster_hit"] <- "RS4_hit"
# # Lets simplify by looking for pval < 0.05 first
# genefinder_table_filt1 <- na.omit(genefinder_table)
# 
# 
# xx1 <- merge(eigenegene_to_cluster_table, data.frame(table(genefinder_table_filt1[,1])), by.x = "row.names", by.y = "Var1", all = TRUE)
# xx1[order(xx1[,3]),]

# ## Outfile
# # outfilepathsurvival = paste0("../../projects/newman-pace-2020/output_hg19_newcontrols/run1_nooutlier1/survival_analysis/run1/tertilesplit/")
# # outfilepathsurvival = paste0("../../projects/newman-pace-2020/output_hg19_newcontrols/run1_nooutlier1/survival_analysis/run1/binarysplit/")
# # outfilepathsurvival = paste0("../../projects/newman-pace-2020/output_hg19_newcontrols/run1_nooutlier1/survival_analysis/run1/05sdtertile/")
# 
# # outfilepathsurvival = paste0("../../projects/newman-pace-2020/output_hg19_newcontrols/run1_nooutlier1/survival_analysis/run1/topbotquartile/")
# # outfilepathsurvival = paste0("../../projects/newman-pace-2020/output_hg19_newcontrols/run1_nooutlier1/survival_analysis/run1/topbottertile/")
# 
# outfilepathsurvival = paste0("../../projects/newman-pace-2020/output_hg19_newcontrols/run1_nooutlier1/survival_analysis/run2/topbotquartile/")
# # outfilepathsurvival = paste0("../../projects/newman-pace-2020/output_hg19_newcontrols/run1_nooutlier1/survival_analysis/run3_testpacegroup/topbotquartile/")
# dir.create(outfilepathsurvival, recursive = TRUE, showWarnings = FALSE)
# 
# ## Ok - so we want to iterate over every GOI, and every event, and see what comes out
# # survival_intable <- merge(inmetatable[grepl("PACE", inmetatable[,"PATNUM"]),c("PATNUM", EOIcols)], 
# #                           GOIcat_counttable, by = "PATNUM")
# survival_intable <- merge(inmetatable[grepl("PACE", inmetatable[,"PATNUM"]),c("PATNUM", EOIcols)],
#                           GOIcat_counttable, by = "PATNUM")
# survival_intable <- survival_intable[survival_intable[,1] %in% SOI,]
# 
# fullpvaloutlist <- fullcumhazoutlist <- list()
# for (genenum in seq_len(nrow(GOItab))) {
#     ## Select gene
#     GOIsel <- rownames(GOItab)[genenum]
#     
#     survpvallist <- survHRlist <- list()
#     for (eventnum in seq_len(length(EOI))) {
#         ## Select event
#         eventsel <- EOI[eventnum]
#         
#         ## Create the subsurvivaltab
#         survivaldata <- data.frame(na.omit(survival_intable[,c(paste0(c("time_to_", "censor_"), eventsel),GOIsel)]))
#         rownames(survivaldata) <- survival_intable[,"PATNUM"]
#         
#         ## Run the analysis
#         survivalanalysisout <- create_survival_plot(survivaldata = survivaldata, timebreakparam = NULL)
#         outsurvtable <- survivalanalysisout$outsurvtable
#         outsurvplot <- survivalanalysisout$outsurvplot
#         outsurvpvalue <- survivalanalysisout$outsurvpvalue
#         
#         ## Save out the pvalue of the log-rank test
#         survpvallist[[eventnum]] <- outsurvpvalue[,"pval"]
#         names(survpvallist)[eventnum] <- eventsel
#         
#         # Ok, I think we need to return the HR as well for coxph to get some kind of effect size for this analysis
#         outcoxphobject <- survivalanalysisout$outcoxphobject
#         hrvalue <- summary(outcoxphobject)[["conf.int"]][,c("exp(coef)", "lower .95", "upper .95")]
#         
#         ## Save out the coxph HR of the coxph test
#         survHRlist[[eventnum]] <- hrvalue
#         names(survHRlist)[eventnum] <- eventsel
#         
#         ## Write out the plot and table
#         outsubdir <- paste0(outfilepathsurvival, GOIsel, "/", eventsel, "/")
#         dir.create(outsubdir, showWarnings = FALSE, recursive = TRUE)
#         
#         write.table(outsurvtable, paste0(outsubdir, GOIsel, "_", eventsel, "_survival_stat_table.csv"),
#                     sep = ",", col.names = TRUE, row.names = FALSE)
#         pdf(paste0(outsubdir, GOIsel, "_", eventsel, "_survival_plot.pdf"), width = 8, height = 5)
#         print(outsurvplot)
#         junk <- dev.off()
#         
#     }
#     survpvaltab <- do.call(rbind, survpvallist)
#     colnames(survpvaltab) <- paste0(GOIsel, "_pval")
#     survhazardtab <- do.call(rbind, survHRlist)
#     colnames(survhazardtab) <- paste0(GOIsel, "_", colnames(survhazardtab))
#     
#     fullpvaloutlist[[genenum]] <- survpvaltab
#     names(fullpvaloutlist[genenum]) <- paste0(GOIsel, "_", eventsel)
#     fullcumhazoutlist[[genenum]] <- survhazardtab
#     names(fullcumhazoutlist)[genenum] <- paste0(GOIsel, "_", eventsel)
#     
# }
# fullpvalouttable <- t(do.call(cbind, fullpvaloutlist))
# write.table(fullpvalouttable, paste0(outfilepathsurvival, "full_survival_pval_table.csv"), sep = ",", col.names = NA, row.names = TRUE)
# fullcumhazouttable <- t(do.call(cbind, fullcumhazoutlist))
# write.table(fullcumhazouttable, paste0(outfilepathsurvival, "full_survival_cumhaz_table.csv"), sep = ",", col.names = NA, row.names = TRUE)
# 
# ## To create a ranking metric for the importance of a gene in the prediction of events - can we combine -log10pval of the deseq with the pval of the KM analysis...
# # log10fullstatouttable <- data.frame(-log10(fullstatouttable))
# # log10deseqpadjtable <- data.frame(-log10(indeseqtable[rownames(fullstatouttable), "padj",drop=FALSE]))
# # # log10deseqpadjtable <- do.call(cbind, rep(-log10(indeseqtable[rownames(fullstatouttable), "padj",drop=FALSE]), 6))
# # tt1 <- log10fullstatouttable * t(log10deseqpadjtable) * t(indeseqtable[rownames(fullstatouttable), "log2FoldChange",drop=FALSE])
# # 
# # 
# 
# 
# # fullstatouttablefile <- "../../projects/newman-pace-2020/output_hg19_newcontrols/run1_nooutlier1/survival_analysis/run1/topbotquartile/full_survival_stat_table.csv"
# fullstatouttablefile <- "../../projects/newman-pace-2020/output_hg19_newcontrols/run1_nooutlier1/survival_analysis/run2/topbotquartile/full_survival_stat_table.csv"
# # fullstatouttablefile <- "../../projects/newman-pace-2020/output_hg19_newcontrols/run1_nooutlier1/survival_analysis/run3_testpacegroup/topbotquartile/full_survival_stat_table.csv"
# fullstatouttable <- read.table(fullstatouttablefile, sep = ",", header = TRUE, row.names = 1)
# 
# ## Turn the event outcome into a heatmap
# # eventpvalcutoff <- 0.05
# eventpvalcutoff <- 0.1
# 
# # GOI heatmap with events
# ## Lets create a heatmap with the GOI, and then use a side annotation for if its indicative of events
# # eventpvaltabfile <- "../../projects/newman-pace-2020/output_hg19_newcontrols/run1_nooutlier1/survival_analysis/topbotquartile/full_survival_stat_table.csv"
# # outeventpvaltab <- read.table(eventpvaltabfile, sep = ",", header = TRUE, row.names = 1)
# outeventpvaltab <- fullstatouttable
# 
# ## Create the initial maptab
# GOI <- rownames(GOItab)
# outeventpvaltab <- outeventpvaltab[GOI,]
# maptab <- data.frame(outeventpvaltab[rowSums(outeventpvaltab[,EOIvec] < eventpvalcutoff) > 0, EOIvec])
# 
# ## Ok, so what if we do multihyp correction
# # p.adjust(unlist(maptab), method = "fdr") ## Sanity check - it works
# # “holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”
# correctedmat <- matrix(p.adjust(as.vector(as.matrix(maptab)), method='fdr'),ncol=ncol(maptab))
# rownames(correctedmat) <- rownames(maptab)
# colnames(correctedmat) <- colnames(maptab)
# maptab <- correctedmat
# 
# # Create a heatmap of just the event pvalues
# EOIvec <- colnames(outeventpvaltab)
# ## Add in geneordering
# # geneorder <- indeseqtable[rownames(outeventpvaltab)[rowSums(outeventpvaltab[,EOIvec] < eventpvalcutoff) > 0],]
# # geneorder <- rownames(geneorder[order(geneorder[,"padj"], decreasing = FALSE),]) ## order by padj
# # geneorder <- rownames(geneorder[order(geneorder[,"log2FoldChange"], decreasing = TRUE),]) ## order by log2fc
# geneorder <- rownames(correctedmat[order(indeseqtable[rownames(correctedmat),"log2FoldChange"], decreasing = TRUE),]) ## order by log2fc
# ## Order by avg pval across all events
# # geneorder <- names(rowMeans(correctedmat[rownames(geneorder),])[order(rowMeans(correctedmat[rownames(geneorder),]), decreasing = FALSE)])
# # geneorder <- names(rowMeans(correctedmat[order(rowMeans(correctedmat), decreasing = FALSE),]))
# 
# # maptab <- data.frame(outeventpvaltab[rowSums(outeventpvaltab[,EOIvec] < eventpvalcutoff) > 0, EOIvec])
# maptab <- maptab[geneorder,][rowSums(maptab[geneorder,] < eventpvalcutoff) > 0,]
# 
# 
# heatmapcolorparam = colorRamp2(c(1, eventpvalcutoff + 0.0001, eventpvalcutoff, 0), c("white", "white", "#c9e0dc", "#00705E"))
# heatmapcolorLEGEND = colorRamp2(c(eventpvalcutoff, 0), c("#c9e0dc", "#00705E"))
# # heatmapcolorparam = colorRamp2(c(1, 0.1 + 0.0001, 0.1, 0), c("white", "white", "#c9e0dc", "#00705E"))
# # heatmapcolorLEGEND = colorRamp2(c(0.1, 0), c("#c9e0dc", "#00705E"))
# 
# rowmetatable <- indeseqtable[rownames(maptab),c("padj", "log2FoldChange")]
# rowannotationlist <- annotationlist_builder(rowmetatable, customcolorlist = list(
#     # log2FoldChange = colorRamp2(c(ceiling(max(rowmetatable[,"log2FoldChange"])), 0, floor(min(rowmetatable[,"log2FoldChange"]))), brewer.pal(3,"RdBu")),
#     # log2FoldChange = colorRamp2(c(2.5, 0, -1.5), brewer.pal(3,"RdBu")), ##CUSTOM
#     log2FoldChange = colorRamp2(c(2, 0, -1), c("#A30000", "White", "#0000A7")), ##CUSTOM 
#     padj = colorRamp2(seq(0.05,0,length = 5), brewer.pal(5,"Reds"))
# ))
# ## Define parameters for each of the labels on the annotation bars
# temp1 <- vector("list", length(rowannotationlist))
# names(temp1) = names(rowannotationlist)
# annotlegendlist = lapply(temp1, function(x)
#     x[[1]] = list(title_gp=gpar(fontsize=8, fontface="bold"), labels_gp=gpar(fontsize=8)))
# haside = rowAnnotation(df = rowmetatable,
#                        col = rowannotationlist,
#                        na_col = "white",
#                        gp = gpar(fontsize = 0.5),
#                        show_annotation_name=TRUE,
#                        annotation_name_gp = gpar(fontsize = 8, fontface="bold"),
#                        annotation_name_side = "top",
#                        simple_anno_size = unit(min(60/length(rowannotationlist), 5),"mm"),
#                        show_legend = TRUE,
#                        annotation_legend_param = annotlegendlist)
# 
# 
# ht1 = Heatmap(as.matrix(maptab), 
#               col = heatmapcolorparam,    ## Define the color scale for the heatmap
#               row_title = "Genes",                                       ## Name the rows
#               column_title = "Events",                                  ## Name the columns
#               border = TRUE,
#               na_col = "white",
#               rect_gp = gpar(col = "black", lwd = 0.5),
#               
#               cluster_columns = FALSE,                         ## Cluster the columns or leave as is
#               cluster_rows = FALSE,                            ## Cluster the rows or leave as is
#               clustering_method_rows = "ward.D2",
#               clustering_method_columns = "ward.D2",
#               #"ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", centroid"
#               
#               show_column_names = TRUE,                                  ## Show the Column Names
#               column_names_gp = gpar(fontsize = 6),                      ## Change the size of the column names
#               show_row_names = TRUE,                                    ## Show the row names
#               row_names_side = "left",                                   ## Place the row names on the Left side of the heatmap
#               row_names_gp = gpar(fontsize=6),
#               
#               show_row_dend = TRUE, #nrow(maptab) <=500,                                     ## Show the dendrogram on the rows
#               show_column_dend = TRUE,                                   ## Show the dendrogram on the columns
#               
#               heatmap_legend_param = list(
#                   col_fun = heatmapcolorLEGEND,
#                   at = seq(0, eventpvalcutoff, 0.01),
#                   # at = seq(0, 0.1, 0.01),
#                   # title = ifelse(samplesample==FALSE, "Zscore", "Spearman\nCorrelation"),
#                   title = "zscore",
#                   legend_height = unit(2.5, "cm"),
#                   title_gp = gpar(fontsize = 8, fontface = "bold")),
#               height = unit(min((nrow(maptab)/2), 12),"cm"),
#               width = unit(min(ncol(maptab), 18),"cm")
#               
# )
# draw(ht1+haside)
# 
# # hmoutfile <- paste0(outfilepathsurvival, "event_heatmap_with_annot.pdf")
# # hmoutfile <- paste0(outfilepathsurvival, "CORRECTED_event_heatmap_log2fcsort_with_annot_005.pdf")
# hmoutfile <- paste0(outfilepathsurvival, "CORRECTED_event_heatmap_log2fcsort_with_annot_01.pdf")
# # hmoutfile <- paste0(outfilepathsurvival, "CORRECTED_event_heatmap_log2fcsort_with_annot.pdf")
# pdf(hmoutfile, 11, 9)
# draw(ht1+haside)
# junk <- dev.off()
# 
# 
# 
# 
# ## Lets output a cleaned version of the cumhazouttable, and a heatmap with the CI and text within it
# # tt1 <- data.frame(fullcumhazouttable)
# genestoagg <- unique(gsub("_.*", "", rownames(fullcumhazouttable)))
# ## I can definitely aggregate this, but idk how, so I am gonna use a forloop instead.
# HRoutlist <- textoutlist <- list()
# for (genenum in seq_len(length(genestoagg))) {
#     ## Grab the rows we need
#     genesel <- genestoagg[genenum]
#     selrows <- fullcumhazouttable[grepl(genesel, rownames(fullcumhazouttable)),]
#     
#     ## Output to a new table
#     # textoutlist[[genenum]] <- apply(selrows, 2, function(x){paste0(round(x[1], 2), " (",round(x[2], 2), " - ", round(x[3], 2), ")")})
#     ## Lets try just the CIs
#     textoutlist[[genenum]] <- apply(selrows, 2, function(x){paste0("(",round(x[2], 2), " - ", round(x[3], 2), ")")})
#     names(textoutlist)[genenum] <- genesel
#     
#     ## Write out just the HR to a list
#     HRoutlist[[genenum]] <- selrows[grepl("_exp", rownames(selrows)),]
#     names(HRoutlist)[genenum] <- genesel
# }
# HRtable <- do.call(rbind, HRoutlist)
# HRtexttable <- do.call(rbind, textoutlist)
# 
# heatmapcolorparam <- colorRamp2(c(0, 1, 6), c( "#0000A7", "White", "#A30000"))
# hrmap = Heatmap(HRtable[rownames(maptab),], 
#                 border = TRUE,
#                 rect_gp = gpar(col = "black", lwd = 1),
#                 col = heatmapcolorparam,    ## Define the color scale for the heatmap
#                 row_title = "Gene",                                       ## Name the rows
#                 column_title = "Event",                                  ## Name the columns
#                 
#                 cluster_columns = FALSE,                         ## Cluster the columns or leave as is
#                 cluster_rows = FALSE,                            ## Cluster the rows or leave as is
#                 
#                 show_column_names = TRUE  ,                               ## Show the Column Names
#                 column_names_gp = gpar(fontsize = 6),                      ## Change the size of the column names
#                 show_row_names = TRUE,                                   ## Show the row names
#                 row_names_side = "left",                                   ## Place the row names on the Left side of the heatmap
#                 row_names_gp = gpar(fontsize=6),
#                 
#                 heatmap_legend_param = list(
#                     title = "Hazard Ratio",
#                     legend_height = unit(2.5, "cm"),
#                     title_gp = gpar(fontsize = 8, fontface = "bold")),
#                 
#                 height = unit(20,"cm"),
#                 width = unit(10,"cm"),
#                 
#                 ## Adds in the cor value from the maptab
#                 # cell_fun = function(j, i, x, y, width, height, fill) {
#                 #     grid.text(paste0(sprintf("%.2f", t(moduleTraitCor)[,colnames(maptab),drop=FALSE][i, j]), "\n", "(",
#                 #                      sprintf("%.2f", t(moduleTraitPvalue)[,colnames(maptab),drop=FALSE][i, j]), ")"),
#                 #               x, y, gp = gpar(fontsize = 8))
#                 # }
#                 cell_fun = function(j, i, x, y, width, height, fill) {
#                     # grid.text(paste0(sprintf("%.2f", moduleTraitCor[i, j]), "\n", "(",
#                     #                  sprintf("%.2f", moduleTraitPvalue[i, j]), ")"),
#                     grid.text(HRtexttable[rownames(maptab),][i, j],
#                               x, y, gp = gpar(fontsize = 6))
#                 }
#                 
# )
# hmoutfile2 <- paste0(outfilepathsurvival, "HR_heatmap.pdf")
# # hmoutfile <- paste0(outfilepathsurvival, "CORRECTED_event_heatmap_log2fcsort_with_annot.pdf")
# pdf(hmoutfile2, 11, 9)
# draw(hrmap)
# junk <- dev.off()



# --------------------------------- rnaclusters -> diff exp genes -> singscore of genes -> new cluster labels in other pops ---------------------------------
# Ok brand new idea
## (1) With your 4 rna clusters - use diffexp to get genes that are diff exp for that cluster vs not that cluster
## (2) with those genes, create gene sets
## (3) Then run singscore with those genesets to get single clusterscores for each sample
## (4) Then use that info to train and test a model that should have very high prediction of the original cluster based off of those singscore values
## (5) Then we can use those genes to calcuate new singscore values for each clustertype of our NEW population (PACE) and then see if those clusters associate with risk

# Will need this for the ischemia dataset training
bioreptable_waddons <- read.table("output/run4_rmoutliers2_asr_control/integrative_analyses/run7/bioreptable_waddons.csv", sep = ",", header = TRUE)
rownames(bioreptable_waddons) <- bioreptable_waddons[,"PATNUM"]
cluster_membership_table <- na.omit(bioreptable_waddons[,"rna_4cluster",drop=FALSE])
eventpvalcutoff <- 0.1


## OK - I THINK I NEED TO WRAP FROM HERE ON DOWN WITH A GENESETSIZE FOR LOOP
scoring_genesetsize = c("all", 1000, 500, 100, 50)
scoring_type = list("TotalScore", "UpScore", "DownScore", c("TotalScore", "UpScore", "DownScore"))
param_grid <- expand.grid(scoring_genesetsize, scoring_type)
# for (genesetsize in scoring_genesetsize) {
for (param_row in seq_len(nrow(param_grid))) {
    # Select parameters
    genesetsize = as.character(param_grid[param_row,1])
    scoring_type = param_grid[param_row,2][[1]]
    scoring_type_label = ifelse(scoring_type == "TotalScore", "Tot",
                         ifelse(scoring_type == "UpScore", "Up", 
                         ifelse(scoring_type == "DownScore", "Down", NA)))
    if (length(scoring_type_label) == 3) {scoring_type_label = "All"}
    
    subanalysis_outfilepath <- paste0(outfilepathmaster, "singscore_to_cluster_classifier/", genesetsize, scoring_type_label, "_singscoregenes/")
    dir.create(subanalysis_outfilepath, showWarnings = FALSE, recursive = TRUE)
    ## (1-3) Diffexp genes, singscore, values
    singscore_list <- c(
        # "rna_subtype1", "rna_subtype2", "rna_subtype3", "rna_subtype4"
        "rna_subtype1", "rna_subtype2", "rna_subtype3", "rna_subtype4", "rna_subtype3A", "rna_subtype3B", "rna_subtype3Bv3A"
        # "meth_subtype1", "meth_subtype2", "meth_subtype3"
    )
    singscore_table_list <- singscore_genelists_list <- list()
    singscore_path <- paste0("output/run4_rmoutliers2_asr_control/integrative_analyses/run7/signature_scoring/")
    for (singscore_result in singscore_list) {
        singscore_basepath <- paste0(singscore_path, singscore_result, "/singscore/")
        singscore_subpath <- list.files(singscore_basepath)[grepl(paste0(genesetsize, "UP"), list.files(singscore_basepath)) | grepl(paste0(genesetsize, "DOWN"), list.files(singscore_basepath))]
        upfilename <- list.files(paste0(singscore_basepath, singscore_subpath, "/"))[grepl("GOIup", list.files(paste0(singscore_basepath, singscore_subpath, "/")))]
        downfilename <- list.files(paste0(singscore_basepath, singscore_subpath, "/"))[grepl("GOIdown", list.files(paste0(singscore_basepath, singscore_subpath, "/")))]
        
        intable_singscore <- read.table(paste0(singscore_basepath, singscore_subpath, "/", singscore_result, "_singscore_outtable.csv"),
                                        sep = ",", header = TRUE)
        singscore_table_list[[singscore_result]] <- intable_singscore
        
        singscore_genelists_list[[singscore_result]][["Upset"]] <- read.table(paste0(singscore_basepath, singscore_subpath, "/", upfilename),
                                                                              sep = ",", header = TRUE)[,1]
        singscore_genelists_list[[singscore_result]][["Downset"]] <- read.table(paste0(singscore_basepath, singscore_subpath, "/", downfilename),
                                                                                sep = ",", header = TRUE)[,1]
    }
    combined_singscore_table <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Row.names", all = TRUE, sort = FALSE), 
                                       lapply(singscore_table_list, function(x) x[,c("Row.names", scoring_type)]))
    dimnames(combined_singscore_table) <- list(combined_singscore_table[,"Row.names"], 
        c("Row.names", apply(expand.grid(paste0(scoring_type, "__"), names(singscore_table_list)), 1, paste0, collapse = "")))
    combined_singscore_table <- combined_singscore_table[,!grepl("Row.names", colnames(combined_singscore_table))]
    write.table(combined_singscore_table, paste0(subanalysis_outfilepath, "isch_", genesetsize, scoring_type_label, "_singscore_table.csv"),
                sep = ",", col.names = NA, row.names = TRUE)
    
    ## Now using this, train our model using the 4 subtypes ( we will do 3A,B after this - I think its a better way to do it)
    featuretable <- combined_singscore_table[,
        grepl(paste(c("rna_subtype1$", "rna_subtype2$", "rna_subtype3$", "rna_subtype4$"), collapse = "|"), colnames(combined_singscore_table))]
    outcometable <- na.omit(bioreptable_waddons[,"rna_4cluster",drop=FALSE])[rownames(combined_singscore_table),,drop=FALSE]
    outcomelevels <- c("nmf_cluster_1", "nmf_cluster_2", "nmf_cluster_3", "nmf_cluster_4")
    cohort_multi_classification_ml_analysis_out <- cohort_multi_classification_ml_analysis(featuretable, outcometable, outcomelevels, seedparam=11111,
                                                                                           presplitdata = NULL, train_partition_percent = 0.7,
                                                                                           subsampleparam = NULL, 
                                                                                           # models_to_run = c("xgbTree", "svmRadial", "rf", "glmnet", "nnet", "kknn")) 
                                                                                           models_to_run = "glmnet")
    saveRDS(cohort_multi_classification_ml_analysis_out, paste0(subanalysis_outfilepath, "cohort_multi_classification_ml_analysis_out.RDS"))
    write.table(rbind(cohort_multi_classification_ml_analysis_out[[3]], cohort_multi_classification_ml_analysis_out[[5]]), 
                paste0(subanalysis_outfilepath, "isch_", genesetsize, scoring_type_label, "_modeling_results.csv"),
                sep = ",", col.names = NA, row.names = TRUE)
    pdf(paste0(subanalysis_outfilepath, "isch_", genesetsize, scoring_type_label, "_modeling_testrocplot.pdf"))
    print(cohort_multi_classification_ml_analysis_out[[10]])
    junk <- dev.off()
    # cohort_multi_classification_ml_analysis_out <- readRDS(paste0(subanalysis_outfilepath, "cohort_multi_classification_ml_analysis_out.RDS"))
    # glmnet is the best!
    ischemia_RS_subtype_modfit <- cohort_multi_classification_ml_analysis_out[["modellist"]][["glmnet_modfit"]]
    
    
    # Write out a heatmap here of the featuretable with the outcome as an annotation
    colannotationlist <- custom_annotation_list_from_colorguide(COI = "rna_4cluster", colorguide)
    isch_model_out <- create_heatmap(counttab = t(featuretable), subsetnum = FALSE, scale_data = TRUE,
                                     colmetatable = outcometable, colannotationlist = colannotationlist,
                                     rowmetatable = NULL, rowannotationlist = NULL,
                                     colclusterparam = TRUE, rowclusterparam = FALSE, separate_legend = FALSE,
                                     columnsplitparam = outcometable
                                     )
    pdf(paste0(subanalysis_outfilepath, "isch_", genesetsize, scoring_type_label, "_singscore_hm.pdf"), 15, 10, useDingbats = FALSE)
    draw(isch_model_out[[1]])
    junk <- dev.off()
    
    
    
    ## Now with this - we also want to then FURTHER split our 3s into 3A and 3B using a predictive algorithm:
    ## Now using this, train our model using the 4 subtypes ( we will do 3A,B after this - I think its a better way to do it)
    subsplitsamples <- rownames(bioreptable_waddons)[bioreptable_waddons[,"rna_4cluster_w3AB"] %in% c("nmf_cluster_3A", "nmf_cluster_3B")]
    featuretable <- combined_singscore_table[subsplitsamples,]
    outcometable <- na.omit(bioreptable_waddons[subsplitsamples,"rna_4cluster_w3AB",drop=FALSE])
    outcomelevels <- c("nmf_cluster_3A", "nmf_cluster_3B")
    RS3AB_subcohort_multi_classification_ml_analysis_out <- cohort_multi_classification_ml_analysis(featuretable, outcometable, outcomelevels, seedparam=11111,
                                                                                           presplitdata = NULL, train_partition_percent = 0.6,
                                                                                           subsampleparam = NULL, 
                                                                                           models_to_run = "glmnet") 
    saveRDS(RS3AB_subcohort_multi_classification_ml_analysis_out, paste0(subanalysis_outfilepath, "RS3AB_subcohort_multi_classification_ml_analysis_out.RDS"))
    # RS3AB_subcohort_multi_classification_ml_analysis_out <- readRDS(paste0(subanalysis_outfilepath, "RS3AB_subcohort_multi_classification_ml_analysis_out.RDS"))
    # glmnet is the best! (again!!)
    ischemia_RS3AB_subtype_modfit <- RS3AB_subcohort_multi_classification_ml_analysis_out[["modellist"]][["glmnet_modfit"]]
    
    
    
    # Then return the model we want - and use that model to predict the classess of our validation cohort
    ## First we have to generate the singscore values for each sample in our PACE cohort
    ## For every score - read in the sets, and then apply singscore to get a value for every sample in our PACE cohort - repeat for all cluster scores
    pace_ss_table_list <- list()
    for (singscore_genelist_num in seq_along(singscore_genelists_list)) {
        ss_upset <- singscore_genelists_list[[singscore_genelist_num]][["Upset"]]
        ss_downset <- singscore_genelists_list[[singscore_genelist_num]][["Downset"]]
        ss_label <- names(singscore_genelists_list)[singscore_genelist_num]
        # singscore_table <- singscore_sample_scoring(GOIup=ss_upset, GOIdown=ss_downset, counttable=pace_wb_normcounttable)
        singscore_table <- singscore_sample_scoring(GOIup=ss_upset, GOIdown=ss_downset, counttable=paceONLY_wb_normcounttable)
        pace_ss_table_list[[singscore_genelist_num]] <- cbind(Row.names = rownames(singscore_table), singscore_table)
        names(pace_ss_table_list)[singscore_genelist_num] <- ss_label
    }
    pace_combined_singscore_table <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Row.names", all = TRUE, sort = FALSE), 
                                            lapply(pace_ss_table_list, function(x) x[,c("Row.names", scoring_type)]))
    # dimnames(pace_combined_singscore_table) <- list(pace_combined_singscore_table[,"Row.names"], c("Row.names", names(pace_ss_table_list)))
    dimnames(pace_combined_singscore_table) <- list(pace_combined_singscore_table[,"Row.names"], 
        c("Row.names", apply(expand.grid(paste0(scoring_type, "__"), names(singscore_table_list)), 1, paste0, collapse = "")))
    pace_combined_singscore_table <- pace_combined_singscore_table[,!grepl("Row.names", colnames(pace_combined_singscore_table))]
    write.table(pace_combined_singscore_table, paste0(subanalysis_outfilepath, "pace_", genesetsize, scoring_type_label, "_singscore_table.csv"),
                sep = ",", col.names = NA, row.names = TRUE)
    
    # Now we have our new input data - we grab our most successful model from the ischemia data, and apply it here to give these samples an RNA subtype
    # ischemia_RS_subtype_modfit
    # grepl(paste(c("rna_subtype1$", "rna_subtype2$", "rna_subtype3$", "rna_subtype4$"), collapse = "|"), colnames(combined_singscore_table))
    pace_inmodeling_singscore_table <- pace_combined_singscore_table[,grepl(paste(c("rna_subtype1$", "rna_subtype2$", "rna_subtype3$", "rna_subtype4$"), collapse = "|"), colnames(pace_combined_singscore_table))]
    pace_predict_labels <- predict(ischemia_RS_subtype_modfit, pace_inmodeling_singscore_table)
    pace_predict_values <- predict(ischemia_RS_subtype_modfit, pace_inmodeling_singscore_table, type = "prob")
    pace_predicted_values_and_labels <- cbind(pace_predict_values, pace_predict_labels)
    pace_subsplit_samples <- rownames(pace_predicted_values_and_labels[pace_predicted_values_and_labels[,"pace_predict_labels"] %in% "nmf_cluster_3",])
    write.table(pace_predicted_values_and_labels, paste0(subanalysis_outfilepath, "pace_predicted_values_and_labels.csv"),
                sep = ",", col.names = NA, row.names = TRUE)
    
    # Write out a heatmap here of the featuretable with the outcome as an annotation
    pace_colmetatable <- data.frame(rna_4cluster = pace_predicted_values_and_labels[,"pace_predict_labels"], row.names = rownames(pace_predicted_values_and_labels))
    colannotationlist <- custom_annotation_list_from_colorguide(COI = "rna_4cluster", colorguide)
    pace_model_out <- create_heatmap(counttab = t(pace_inmodeling_singscore_table), subsetnum = FALSE, scale_data = TRUE,
                                     colmetatable = pace_colmetatable, colannotationlist = colannotationlist,
                                     rowmetatable = NULL, rowannotationlist = NULL,
                                     colclusterparam = TRUE, rowclusterparam = FALSE, separate_legend = FALSE,
                                     columnsplitparam = pace_colmetatable
    )
    pdf(paste0(subanalysis_outfilepath, "pace_", genesetsize, scoring_type_label, "_singscore_hm.pdf"), 15, 10, useDingbats = FALSE)
    draw(pace_model_out[[1]])
    junk <- dev.off()
    
    
    # But also do a subprediction for our 3A and 3Bs
    subsplit_pace_predict_labels <- predict(ischemia_RS3AB_subtype_modfit, pace_combined_singscore_table[pace_subsplit_samples,])
    subsplit_pace_predict_values <- predict(ischemia_RS3AB_subtype_modfit, pace_combined_singscore_table[pace_subsplit_samples,], type = "prob")
    rownames(subsplit_pace_predict_values) <- pace_subsplit_samples
    subsplit_pace_predicted_values_and_labels <- cbind(subsplit_pace_predict_values, subsplit_pace_predict_labels)
    
    
    
    # Now with this - do they have differences in events?
    ## I think I want 3A and 3B in there
    
    source("code/ischemia2021_exploratory_analysis_functions.R")
    dir.create(paste0(subanalysis_outfilepath, "event_analysis/"), showWarnings = FALSE, recursive = TRUE)
    
    cluster_intable <- cbind(pace_predicted_values_and_labels[,"pace_predict_labels",drop=FALSE], pace_w_subsplit = as.character(pace_predicted_values_and_labels[,"pace_predict_labels"]))
    cluster_intable[rownames(cluster_intable[,"pace_w_subsplit",drop=FALSE]) %in% rownames(subsplit_pace_predicted_values_and_labels),"pace_w_subsplit"] <- as.character(subsplit_pace_predicted_values_and_labels[,"subsplit_pace_predict_labels"])
    cluster_intable[,"pace_w_subsplit"] <- factor(cluster_intable[,"pace_w_subsplit"], levels = c("nmf_cluster_1", "nmf_cluster_2", "nmf_cluster_3A", "nmf_cluster_3B", "nmf_cluster_4"))
    write.table(cluster_intable, paste0(subanalysis_outfilepath, "event_analysis/", "cluster_intable.csv"), sep = ",", col.names = NA, row.names = TRUE)
    
    event_intable <- pace_wb_full_metatable[rownames(cluster_intable), grepl("censor_|time_to_", colnames(pace_wb_full_metatable))]
    event_intable <- event_intable[, !grepl("_30", colnames(event_intable))]
    # event_intable <- pace_wb_metatable[grepl("PACE", rownames(pace_wb_metatable)), grepl("censor_|time_to_", colnames(pace_wb_metatable))]
    event_intable <- data.frame(cbind(PATNUM = rownames(event_intable), event_intable))
    colnames(event_intable) <- gsub("censor_", "C_", colnames(event_intable))
    colnames(event_intable) <- gsub("time_to_", "T_", colnames(event_intable))
    write.table(event_intable, paste0(subanalysis_outfilepath, "event_analysis/", "event_intable.csv"), sep = ",", col.names = NA, row.names = TRUE)
    
    plottable_outlist <- list()
    for (group_column_selected in colnames(cluster_intable)) {
        
        selected_cluster_table <- cbind(PATNUM = rownames(cluster_intable[order(cluster_intable[,group_column_selected]), group_column_selected,drop=FALSE]),
                                        cluster_intable[order(cluster_intable[,group_column_selected]), group_column_selected,drop=FALSE])
        
        # COIref_tablist <- list(
        #     # IMGDEGIS = bioreptable_waddons[,"IMGDEGIS", drop=FALSE],
        #     IMGDEGIS_01_3 = bioreptable_waddons[,"IMGDEGIS_01_3", drop=FALSE],
        #     # CTNDV50 = bioreptable_waddons[,"CTNDV50", drop=FALSE],
        #     # CTNDV50_13_clean = bioreptable_waddons[,"CTNDV50_13_clean", drop=FALSE],
        #     CTNDV70_13_clean = bioreptable_waddons[,"CTNDV70_13_clean", drop=FALSE],
        #     CKD = bioreptable_waddons[,"CKD",drop=FALSE],
        #     highlowage = highlowage
        # )
        
        # plotCOI <- c(names(COIref_tablist),
        #              paste0(group_column_selected, "_",
        #                     unique(na.omit(selected_cluster_table[,group_column_selected]))[order(unique(na.omit(selected_cluster_table[,group_column_selected])))])
        # )
        # COIref_tablist <- NULL
        
        plotCOI <- paste0(group_column_selected, "_", unique(na.omit(selected_cluster_table[,group_column_selected]))[
            order(unique(na.omit(selected_cluster_table[,group_column_selected])))])
        
        survivalplot_outfilepath <- paste0(subanalysis_outfilepath, "event_analysis/", group_column_selected, "/")
        dir.create(survivalplot_outfilepath, showWarnings = FALSE, recursive = TRUE)
        summary_heatmap_out <- cluster_vs_event_kmanalysis_summary_heatmap(selected_cluster_table = selected_cluster_table,
                                                                           COIref_tablist = NULL,
                                                                           survivalplot_outfilepath = survivalplot_outfilepath, 
                                                                           bioreptable_waddons = event_intable,
                                                                           eventpvalcutoff = eventpvalcutoff, return_output_tables = TRUE,
                                                                           plotCOI)
        pdf(paste0(survivalplot_outfilepath, group_column_selected, "_event_summary_heatmap.pdf"),
            12, 10, useDingbats = FALSE)
        draw(summary_heatmap_out[[1]])
        junk <- dev.off()
        
        
        ## Repeat with NO FILTER
        heatmapcolorparam = colorRamp2(c(-1, -1 - 0.0001, -1, -0.0001, -1e-100, 0, 1e-100, 0.0001, 1, 1 + 0.0001, 1),
                                       c("grey", "grey", "#bcb2fd", "#11007d", "#11007d", "white", "#7d0000", "#7d0000", "#ffb2b2", "grey", "grey"))
        ht1_nofilt = Heatmap(as.matrix(summary_heatmap_out[[2]]),
                      col = heatmapcolorparam, row_title = "Events", column_title = "Cohorts",
                      border = TRUE, na_col = "white", rect_gp = gpar(col = "black", lwd = 0.5),
                      cluster_columns = FALSE, cluster_rows = FALSE,clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
                      show_column_names = TRUE, column_names_gp = gpar(fontsize = 6), show_row_names = TRUE, 
                      row_names_side = "left", row_names_gp = gpar(fontsize=6),show_row_dend = TRUE, show_column_dend = TRUE,
                      heatmap_legend_param = list(title = "KM pval", at = c(-eventpvalcutoff, 0, eventpvalcutoff)),
                      height = unit(min((nrow(summary_heatmap_out[[2]])/2), 12),"cm"), width = unit(min(ncol(summary_heatmap_out[[2]]), 18),"cm")
        )
        pdf(paste0(survivalplot_outfilepath, group_column_selected, "_event_nopvalfilt_summary_heatmap.pdf"),12, 10, useDingbats = FALSE)
        draw(ht1_nofilt)
        junk <- dev.off()
        
        
        plottable_outlist[[group_column_selected]] <- plottable_out <- summary_heatmap_out[[2]]
        
        
        # NONBINARY ANALYSIS
        survivalplot_outfilepath <- paste0(subanalysis_outfilepath, "event_analysis/", group_column_selected, "_NONBINARY/")
        dir.create(survivalplot_outfilepath, showWarnings = FALSE, recursive = TRUE)
        # COIref_tablist <- list(
        #     CUSTOM_IMGDEGIS_NONEMILD = bioreptable_waddons[,"CUSTOM_IMGDEGIS_NONEMILD", drop=FALSE],
        #     CTNDV70 = bioreptable_waddons[,"CTNDV70", drop=FALSE],
        #     CKD = bioreptable_waddons[,"CKD",drop=FALSE],
        #     age_cat = bioreptable_waddons[,"CAGE_RND",drop=FALSE]
        # )
        plotCOI <- c(names(COIref_tablist), group_column_selected)
        coxph_control_table <- pace_wb_full_metatable[rownames(selected_cluster_table), c("sex1", "race1", "age_cat")]
        #
        summary_heatmap_NONBINARY <- cluster_NONBINARY_vs_event_kmanalysis_summary_heatmap(selected_cluster_table = selected_cluster_table,
                                                                                 COIref_tablist = NULL,
                                                                                 survivalplot_outfilepath = survivalplot_outfilepath,
                                                                                 bioreptable_waddons = event_intable,
                                                                                 eventpvalcutoff = eventpvalcutoff, return_output_tables = TRUE, 
                                                                                 return_pairwise_coxph_values = TRUE, 
                                                                                 coxph_control_table = coxph_control_table,
                                                                                 EOI = NULL, plotCOI)
        pdf(paste0(survivalplot_outfilepath, group_column_selected, "_event_summary_heatmap.pdf"),
            12, 10, useDingbats = FALSE)
        draw(summary_heatmap_NONBINARY[[1]])
        junk <- dev.off()
        
        ## Repeat with NO FILTER
        heatmapcolorparam = colorRamp2(c(-1, -1 - 0.0001, -1, -0.0001, -1e-100, 0, 1e-100, 0.0001, 1, 1 + 0.0001, 1),
                                       c("grey", "grey", "#bcb2fd", "#11007d", "#11007d", "white", "#7d0000", "#7d0000", "#ffb2b2", "grey", "grey"))
        ht2_nofilt = Heatmap(as.matrix(summary_heatmap_NONBINARY[[2]]),
                             col = heatmapcolorparam, row_title = "Events", column_title = "Cohorts",
                             border = TRUE, na_col = "white", rect_gp = gpar(col = "black", lwd = 0.5),
                             cluster_columns = FALSE, cluster_rows = FALSE,clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
                             show_column_names = TRUE, column_names_gp = gpar(fontsize = 6), show_row_names = TRUE, 
                             row_names_side = "left", row_names_gp = gpar(fontsize=6),show_row_dend = TRUE, show_column_dend = TRUE,
                             heatmap_legend_param = list(title = "KM pval", at = c(-eventpvalcutoff, 0, eventpvalcutoff)),
                             height = unit(min((nrow(summary_heatmap_out[[2]])/2), 12),"cm"), width = unit(min(ncol(summary_heatmap_out[[2]]), 18),"cm")
        )
        pdf(paste0(survivalplot_outfilepath, group_column_selected, "_event_nopvalfilt_summary_heatmap.pdf"),12, 10, useDingbats = FALSE)
        draw(ht2_nofilt)
        junk <- dev.off()
        
    }
    
}



# --------------------------------- PACE WB - matching diffexp for DeathMI ---------------------------------
# dir.create(paste0(outfilepathmaster, "diffexp_matched_analysis/"), showWarnings = FALSE, recursive = TRUE)
# ## PACE WB counttable and metatable
# pace_wb_normcounttable <- read.table("../../projects/newman-pace-2020/output_hg19_newcontrols/run1_nooutlier1/rna_processing/normcounttab.txt", header = TRUE, row.names = 1, sep = "\t")
# paceONLY_wb_normcounttable <- pace_wb_normcounttable[,grepl("PACE", colnames(pace_wb_normcounttable))]
# thrONLY_wb_normcounttable <- pace_wb_normcounttable[,grepl("THR", colnames(pace_wb_normcounttable))]
# pace_wb_rawcounttable <- read.table("../../projects/newman-pace-2020/output_hg19_newcontrols/run1_nooutlier1/rna_processing/filtrawcounttab.txt", header = TRUE, row.names = 1, sep = "\t")
# 
# ### USE THIS FOR 
# pace_wb_metatable <- read.table("../../projects/newman-pace-2020/output_hg19_newcontrols/run1_nooutlier1/rna_processing/metatable_in.csv", header = TRUE, row.names = 1, sep = ",")
# pace_wb_full_metatable <- read.table("../../projects/newman-pace-2020/data/pace_comb_allmetadata_20201014.csv",
#                                      sep = ",", header = TRUE, row.names = 1)
# # Add Age cat on to this for later calculations
# pace_wb_full_metatable[,"age_cat"] <- ifelse(pace_wb_full_metatable[,"current_age"] >= 71 & 
#                                                  rownames(pace_wb_full_metatable) %in% colnames(paceONLY_wb_normcounttable), "71-97",
#                                              ifelse(pace_wb_full_metatable[,"current_age"] >= 60 & 
#                                                         rownames(pace_wb_full_metatable) %in% colnames(paceONLY_wb_normcounttable), "61-70",
#                                                     ifelse(pace_wb_full_metatable[,"current_age"] >= 25 & 
#                                                                rownames(pace_wb_full_metatable) %in% colnames(paceONLY_wb_normcounttable), "25-60", 
#                                                            NA)))
# # 24-52 53-63 64-91
# SOI <- unique(c(rownames(pace_wb_metatable[pace_wb_metatable[,"Cohort"] %in% "PACE",,drop=FALSE]), rownames(na.omit(pace_wb_metatable[,"comp_THRGrOr_PACEg4_v_THRg4",drop=FALSE]))))
# 
# # Ok - we are doing a simple Diffexp here - so we need different day counts here for event and then go from there
# EOI <- c("DeathMI", "death", "MI")
# daycutoff_list <- c(182.625 * c(1,2,3,4,6,10))
# params_grid <- expand.grid(EOI, daycutoff_list)
# statoutlist <- list()
# for (param_num in seq_len(nrow(params_grid))) {
#     EOI <- params_grid[param_num,"Var1"]
#     daycutoff <- params_grid[param_num,"Var2"]
#     
#     subanalysis_outfolder <- paste0(outfilepathmaster, "diffexp_matched_analysis/", EOI, "_matched_analysis/", EOI, round(daycutoff), "/")
#     dir.create(subanalysis_outfolder, showWarnings = FALSE, recursive = TRUE)
#     # First - grab the PACE samples that we need
#     pace_samples_grab <- rownames(pace_wb_metatable[pace_wb_metatable[,paste0("censor_", EOI)] %in% 1 & 
#                                                     pace_wb_metatable[,paste0("time_to_", EOI)] < daycutoff,])
#     pace_samples_grab <- pace_samples_grab[pace_samples_grab %in% SOI]
#     
#     # thr_samples_grab <- rownames(pace_wb_metatable[pace_wb_metatable[,"Cohort"] %in% "THR_GreatControl",])
#     thr_samples_grab <- rownames(pace_wb_metatable[!pace_wb_metatable[,"Cohort"] %in% "PACE",])
#     thr_samples_grab <- thr_samples_grab[thr_samples_grab %in% SOI]
#     
#     ## Match out PACE and THR samples
#     matchtable <- na.omit(pace_wb_metatable[c(pace_samples_grab, thr_samples_grab),c("Cohort", "age", "sex")])
#     # matchtable <- na.omit(pace_wb_metatable[SOI,c("Cohort", "age", "sex")])
#     matchtable[,"Cohort"] <- ifelse(matchtable[,"Cohort"] %in% "PACE", "PACE", "THR")
#     ## First function creates a "matchtable" that will give you all of your possible matches
#     set.seed(12345) ## NEED THIS TO PRESERVE MATCHING BETWEEN RUNS... OOPS ## BEST
#     out1 <- create_match_table(intable = matchtable, groupcat = "Cohort", treatvar = "PACE", controlvar = "THR",
#                                controlcat = c("age", "sex"), discretevar = "sex", discrete_DF = 0, contvar = list(age = c(var = "age", range = 5)))
#     ## This function helps you select the matches from your table
#     out2 <- create_matches(matchtab = out1, treatcontrolratiomax = 1)
#     
#     ## Gotta at least have 3 matches to pass - or else we go next:
#     if(nrow(out2) < 3) {next}
#     matched_table <- pace_wb_metatable[c(rownames(out2), colnames(out2)), c("Cohort", "age", "age_cat", "sex", "race",
#                                         "ethnicity", "smoking", "bmi", "diabetes", "hypertension", "hyperlipidemia", "CAD")]
#     
#     # First do summary of table, w group of interest
#     testdata_wgroup_summarystats <- summarize_table(intable = cbind(SampleID = rownames(matched_table), matched_table),
#                                                     groupvar = "Cohort", calc_stats = TRUE, calc_ci = FALSE,
#                                                     outfile = paste0(subanalysis_outfolder, "matcheddata_wgroup_summarystats.csv"))
#     testdata_wgroup_cleaned_summarystats <- clean_summarize_table_output(sumstattable_input = testdata_wgroup_summarystats,
#                                                 addpercents = "vertical", contsummary = c("mean", "sd"), roundpvaldigits = 3)
#     write.table(testdata_wgroup_cleaned_summarystats, paste0(subanalysis_outfolder, "matcheddata_wgroup_cleaned_summarystats.csv"), 
#                 sep = ",", col.names = TRUE, row.names = FALSE)
# 
#     
#     # Now - (1) diffexp here, (2) check gene module differences, (3) check subtype scores amongst this group
#     comp_label <- paste0("comp_", EOI, "__", EOI, round(daycutoff), "_v_NO", EOI, round(daycutoff))
#     compcols <- data.frame(c(rep(1, nrow(out2)), rep(0, ncol(out2))), row.names = c(rownames(out2), colnames(out2)))
#     colnames(compcols) <- comp_label
#     deseq_table_out <- DESeq_analysis(compcols = compcols, controlcols = NULL, rawcounttab = pace_wb_rawcounttable)[[1]]
#     write.table(deseq_table_out, paste0(subanalysis_outfolder, comp_label, "_deseq_table.csv"), 
#                 sep = ",", col.names = NA, row.names = TRUE)
#     
#     # (A) Looking at modules:
#     sharedgenes <- intersect(rownames(isch_genestocolorstab), rownames(pace_wb_normcounttable))
#     pace_eigengenetable <- orderMEs(moduleEigengenes(t(pace_wb_normcounttable[sharedgenes,SOI]),
#                                                      isch_genestocolorstab[sharedgenes,"moduleColors"])$eigengenes)
#     
#     
#     # (B) Looking at singscore scores
#     # Looking at singscore values: ## TOSH WORKING HERE - ADD THIS TO ABOVE TABLE AND RUN IT ALL AT ONCE!
#     ## NEED TOO ADD THE DIFF SINGSCORE GENESETS HERE! - read in all combos of rna_subtypes and numbers here:
#     singscore_table_list <- singscore_genelists_list <- list()
#     singscore_path <- paste0("output/run4_rmoutliers2_asr_control/integrative_analyses/run7/signature_scoring/")
#     singscore_list <- c("rna_subtype1", "rna_subtype2", "rna_subtype3", "rna_subtype3A", "rna_subtype3B", "rna_subtype4", "rna_subtype3Bv3A")
#     scoring_genesetsize = c("all", 1000, 500, 100, 50)
#     singscore_grid <- expand.grid(singscore_list, scoring_genesetsize)
#     for (singscore_param_num in seq_len(nrow(singscore_grid))) {
#         singscore_result <- singscore_grid[singscore_param_num,"Var1"]
#         genesetsize <- singscore_grid[singscore_param_num,"Var2"]
#         outlabel <- paste0(singscore_result, "__", genesetsize)
# 
#         singscore_basepath <- paste0(singscore_path, singscore_result, "/singscore/")
#         singscore_subpath <- list.files(singscore_basepath)[grepl(paste0(genesetsize, "UP"), list.files(singscore_basepath)) | grepl(paste0(genesetsize, "DOWN"), list.files(singscore_basepath))]
#         upfilename <- list.files(paste0(singscore_basepath, singscore_subpath, "/"))[grepl("GOIup", list.files(paste0(singscore_basepath, singscore_subpath, "/")))]
#         downfilename <- list.files(paste0(singscore_basepath, singscore_subpath, "/"))[grepl("GOIdown", list.files(paste0(singscore_basepath, singscore_subpath, "/")))]
# 
#         singscore_genelists_list[[outlabel]][["Upset"]] <- read.table(paste0(singscore_basepath, singscore_subpath, "/", upfilename), sep = ",", header = FALSE)[,1]
#         singscore_genelists_list[[outlabel]][["Downset"]] <- read.table(paste0(singscore_basepath, singscore_subpath, "/", downfilename), sep = ",", header = FALSE)[,1]
#     }
#     
#     ## For every score - read in the sets, and then apply singscore to get a value for every sample in our PACE cohort - repeat for all cluster scores
#     pace_ss_table_list <- list()
#     for (singscore_genelist_num in seq_along(singscore_genelists_list)) {
#         ss_upset <- singscore_genelists_list[[singscore_genelist_num]][["Upset"]]
#         ss_downset <- singscore_genelists_list[[singscore_genelist_num]][["Downset"]]
#         ss_label <- names(singscore_genelists_list)[singscore_genelist_num]
#         # singscore_table <- singscore_sample_scoring(GOIup=ss_upset, GOIdown=ss_downset, counttable=pace_wb_normcounttable)
#         singscore_table <- singscore_sample_scoring(GOIup=ss_upset, GOIdown=ss_downset, counttable=pace_wb_normcounttable)
#         pace_ss_table_list[[singscore_genelist_num]] <- cbind(Row.names = rownames(singscore_table), singscore_table)
#         names(pace_ss_table_list)[singscore_genelist_num] <- ss_label
#     }
#     pace_combined_singscore_table <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Row.names", all = TRUE, sort = FALSE), 
#                                             lapply(pace_ss_table_list, function(x) x[,c("Row.names", "TotalScore")]))
#     dimnames(pace_combined_singscore_table) <- list(pace_combined_singscore_table[,"Row.names"], c("Row.names", names(pace_ss_table_list)))
#     pace_combined_singscore_table <- pace_combined_singscore_table[,!grepl("Row.names", colnames(pace_combined_singscore_table))]
#     # write.table(pace_combined_singscore_table, paste0(outfilepathmaster, "singscore_to_cluster_classifier/", "pace_combined_singscore_table.csv"), sep = ",", col.names = NA, row.names = TRUE)
#     
#     ## Plot out the boxplots for each analysis
#     boxplot_counttable <- merge(pace_eigengenetable, pace_combined_singscore_table, by = "row.names")
#     rownames(boxplot_counttable) <- boxplot_counttable[,"Row.names"]
#     boxplot_counttable <- boxplot_counttable[SOI, !grepl("Row.names", colnames(boxplot_counttable))]
#     boxplot_annotationtable <- compcols
#     boxplot_annotationtable[,1] <- ifelse(compcols[,1] == 1, paste0(EOI, round(daycutoff)), paste0("NO", EOI, round(daycutoff)))
#     
#     bpanalysisout <- boxplots_from_counttable_by_annotation(counttable = t(boxplot_counttable), boxplot_annotationtable = boxplot_annotationtable,
#                                            outfilepath = paste0(subanalysis_outfolder, "eigengene_boxplots/"), calculate_corrvalues = FALSE)
#     write.table(bpanalysisout[["wilcox_summarytable"]], paste0(subanalysis_outfolder, "eigengene_boxplots/", "wilcox_summarytable.csv"), 
#                 sep = ",", col.names = TRUE, row.names = FALSE)
#     
#     statoutlist[[param_num]] <- bpanalysisout[["wilcox_summarytable"]]
# 
# }
# fullstat_outtable <- data.frame(do.call(rbind, statoutlist))
# write.table(fullstat_outtable, paste0(outfilepathmaster, "diffexp_matched_analysis/", "fullstat_outtable.csv"),
#                                       sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
# fullstat_outtable[,"pstat"] <- -log10(as.numeric(fullstat_outtable[,"wilcox_pval"])) * 
#     ifelse(fullstat_outtable[,"Feature1_mean"] > fullstat_outtable[,"Feature2_mean"], 1, -1)
# # Viz the table as a heatmap
# hmplottab <- reshape2::dcast(data = fullstat_outtable, Cohort ~ Feature, value.var = "pstat")
# rownames(hmplottab) <- hmplottab[,"Cohort"]
# hmplottab <- hmplottab[,c(eigengene_order, names(singscore_genelists_list))]
# heatmapcolorparam <- colorRamp2(breaks = c(-2, log10(0.05), 0, -log10(0.05), 2), 
#                                      c("#000099", "#ccccff", "white", "#ffb2b2", "#b20000"))
# heatmapcolorparam_filt <- colorRamp2(breaks = c(-2, log10(0.05), log10(0.05)+0.0000001, 0, -log10(0.05)-0.0000001, -log10(0.05), 2), 
#                                 c("#000099", "#ccccff", "grey", "grey", "grey", "#ffb2b2", "#b20000"))
# 
# hmout <- create_heatmap(counttab = hmplottab, subsetnum = FALSE, scale_data = FALSE, colmetatable = NULL,colannotationlist = NULL, rowmetatable = NULL, rowannotationlist = NULL, colclusterparam = FALSE, rowclusterparam = FALSE, separate_legend = FALSE, heatmapcolorparam = heatmapcolorparam, addborders = TRUE, columnsplitparam = NULL, rowsplitparam = NULL)
# pdf(paste0(outfilepathmaster, "diffexp_matched_analysis/", "fullstat_outtable_hm.pdf"), width = 12, height = 8, useDingbats = FALSE)
# draw(hmout[[1]])
# junk <- dev.off()
# 
# 
# hmout_filt <- create_heatmap(counttab = hmplottab, subsetnum = FALSE, scale_data = FALSE, colmetatable = NULL,colannotationlist = NULL, rowmetatable = NULL, rowannotationlist = NULL, colclusterparam = FALSE, rowclusterparam = FALSE, separate_legend = FALSE, heatmapcolorparam = heatmapcolorparam_filt, addborders = TRUE, columnsplitparam = NULL, rowsplitparam = NULL)
# pdf(paste0(outfilepathmaster, "diffexp_matched_analysis/", "fullstat_outtable_hm_filt.pdf"), width = 12, height = 8, useDingbats = FALSE)
# draw(hmout_filt[[1]])
# junk <- dev.off()

# --------------------------------- END ---------------------------------



