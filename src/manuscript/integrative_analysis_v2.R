###########################################################################
#
#                       Integrative2
#
###########################################################################
# Author: ISCHEMIA Study Team
# Date: Updated 2025-08-08
# Description: Integrative2 analysis for ISCHEMIA study
#

# Load configuration and utilities
source("config.R")
source("utils.R")

# Load required packages
# TODO: Update this list based on actual packages used in the script
# load_packages(c("package1", "package2"))

packagelist = c("Hmisc", "tools", "psych")
junk <- lapply(packagelist, function(xxx) suppressMessages(
  require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))

# source("# EXTERNAL_FUNCTION: mgc_plotting_functions.R")
# source("# EXTERNAL_FUNCTION: rnaseq_scripts/geneset_analysis_functions.R")
# source("# EXTERNAL_FUNCTION: WGCNA_functions.R")
# source("# EXTERNAL_FUNCTION: process_metadata_functions.R")
# source("# EXTERNAL_FUNCTION: mgc_survival_functions.R")
# source("# EXTERNAL_FUNCTION: summarize_table_function.R")
# source("# EXTERNAL_FUNCTION: overlap_finder_function.R")


## Outpath
# outfilepathintegration = create_output_dir("analysis_output")
outfilepathintegration = create_output_dir("analysis_output")
dir.create(outfilepathintegration, recursive = TRUE, showWarnings = FALSE)

## Infiles:
metatable_file <- file.path(DATA_DIR, "metatable_in.csv")
metatable <- read.table(metatable_file, sep = ",", header = TRUE, row.names = 1)
biorepfile <- "# PATH_UPDATED: docs/isch_BiorepData_10_26_2021.csv"
bioreptable <- read.table(biorepfile, header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)
rownames(bioreptable) <- bioreptable[,"PATNUM"]
## Selecting here the metaCOI we want to grab in later analyses.
# metaCOI <- c(colnames(metatable)[colnames(metatable) %in% colnames(bioreptable)])
metaCOI <- c(colnames(metatable)[colnames(metatable) %in% colnames(bioreptable)], "SMOKSTAT", "RVBPDIA", "RVBPSYS")

# rna_clustermembership_table_file <- "# PATH_UPDATED: output/run3_rmoutliers2/NMF/rnaseq_only_run2/nmf_clustermembership_table.csv"
rna_clustermembership_table_file <- "# PATH_UPDATED: data/MULTIOMICS_SUBTYPE_LABELS/nmf_cluster_labels_with_subcluster_1_20220208.csv"
rna_clustermembership_table <- read.table(rna_clustermembership_table_file, sep = ",", header = TRUE, row.names = 1)
rna_clustermembership_table <- rna_clustermembership_table[!is.na(rna_clustermembership_table[,"combined"]),]
# rna_clustermembership_table[,"combined"] <- apply(rna_clustermembership_table, 1, function(x) paste(na.omit(x)))
# write.table(rna_clustermembership_table, paste0(outfilepathintegration, "RNA_cluster_membership_table.csv"), sep = ",", row.names = TRUE, col.names = NA)

# meth_clustermembership_table_file <- "# PATH_UPDATED: data/K_means_cluster_group.csv"
meth_clustermembership_table_file <- "# PATH_UPDATED: data/MULTIOMICS_SUBTYPE_LABELS/meth_cluster_membership_table_2_20220210.csv"
meth_clustermembership_table <- read.table(meth_clustermembership_table_file, sep = ",", header = TRUE, row.names = 1)

## Have to add sex back to this - because those are separate for male and female............. so I guess make them F1 and M1..?
# meth_clustermembership_table_wsex <- unique(merge(bioreptable[,c("PATNUM", "SEX")], meth_clustermembership_table[,c("Sample_Patnum", "combat")],
#                                            by.x = "PATNUM", by.y = "Sample_Patnum"))
# meth_clustermembership_table_wsex[,"cluster_wsex"] <- apply(meth_clustermembership_table_wsex[,c("combat", "SEX")], 1, function(x) 
#   paste0(x[1], "_", ifelse(x[2] == "Male", "M", "F")))
# rownames(meth_clustermembership_table_wsex) <- meth_clustermembership_table_wsex[,"PATNUM"]

combined_clustermembership_table <- merge(rna_clustermembership_table[,"combined", drop = FALSE], 
                                          meth_clustermembership_table[,c("kpp_cluster"),drop=FALSE], 
                                          by.x = "row.names", by.y = "row.names", all = TRUE)
colnames(combined_clustermembership_table) <- c("PATNUM", "rnacluster", "methcluster")
rownames(combined_clustermembership_table) <- combined_clustermembership_table[,"PATNUM"]
table(combined_clustermembership_table[,c(2,3)])

## Making conenience variables for all our samples, all rna, and all meth
SOIall <- unique(c(rownames(rna_clustermembership_table), rownames(meth_clustermembership_table)))
SOIrna <- unique(rownames(rna_clustermembership_table))
SOImeth <- unique(rownames(meth_clustermembership_table))
SOIboth <- intersect(rownames(rna_clustermembership_table), rownames(meth_clustermembership_table))

# Counttable:
normcounttable_file <- file.path(DATA_DIR, "normcounttab.txt")
normcounttable <- read.table(normcounttable_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
# WGCNA
eigengenecounttable_file <- "# PATH_UPDATED: output/run3_rmoutliers2/WGCNA/WGCNA_power14_size30/wgcna_eigengenes.csv"
eigengenecounttable <- read.table(eigengenecounttable_file, sep = ",", header = TRUE, row.names = 1)
genestocolorstab_file <- "# PATH_UPDATED: output/run3_rmoutliers2/WGCNA/WGCNA_power14_size30/wgcna_genestocolors.csv"
genestocolorstab <- read.table(genestocolorstab_file, ",", header = TRUE, row.names = 1)
eigengenes <- paste0("ME",
                     c("magenta", "salmon", "purple", "midnightblue", "turquoise", 
                       "yellow", "pink", "red", 
                       "blue", "tan", "brown",
                       "greenyellow", "black", "cyan", "green"))








## First viz - overview of biorep (number of rnaseq, meth, and then after processing, how many we have left)
overlapinlist <- list(
  rnaseq = rownames(rna_clustermembership_table),
  meth = rownames(meth_clustermembership_table)
)
overlapout <- overlap_finder(overlapinlist)
overlaptable <- overlapout$overlaptable
vennplot <- overlapout$vennplot
overlapgrouptab <- overlapout$overlapgrouptab
colnames(overlapgrouptab)[1] <- "Row.names"
overlapsummary <- overlapout$overlapsummary

dir.create(paste0(outfilepathintegration, "omics_sample_comparison/"), showWarnings = FALSE, recursive = TRUE)
write.table(overlaptable, paste0(outfilepathintegration, "omics_sample_comparison/", "testoverlap_overlaptable.csv"), 
            sep = ",", col.names = TRUE, row.names = FALSE)
write.table(overlapgrouptab, paste0(outfilepathintegration, "omics_sample_comparison/", "testoverlap_overlapgrouptable.csv"), 
            sep = ",", col.names = TRUE, row.names = FALSE)
write.table(overlapsummary, paste0(outfilepathintegration, "omics_sample_comparison/", "testoverlap_overlapsummary.csv"), 
            sep = ",", col.names = TRUE, row.names = TRUE)
pdf(paste0(outfilepathintegration, "omics_sample_comparison/", "testoverlap_overlapvenn.pdf"))
grid.draw(vennplot)
junk <- dev.off()




## Second viz - the overlap of the ischemia/CT variables - I really think that is a useful viz
## I would like to get an idea on classifiers - so simple heatmap of the labels for imaging, anatomy, and duke
align_metatable_reflist <- list(
  "IMGDEGIS" = list("0" = "None", "1" = "Mild", "2" = "Moderate", "3" = "Severe", "NA" = "Uninterpretable"),
  "CTNDV50" = list("0" = "Non-evaluable", "1" = "Mild", "2" = "Moderate", "3" = "Severe"),
  "DUKESCORE" = list("7" = "Left Main >=50%", 
                     "6" = "3 Vessels with Severe (>=70%) Plaque or 2 Severe (>=70%) Plaque Including Proximal LAD", 
                     "5" = "3 Vessels with at least Moderate (>=50%) Plaque or 2 Severe (>=70%) Plaque or Severe (>=70%) Proximal LAD Plaque", 
                     "4" = "2 Vessels with at least Moderate (>=50%) Plaque or 1 Severe (>=70%) Plaque", 
                     "3" = "1 Vessel with at least Moderate (>=50%) Plaque"))
labeltable_numeric <- align_metatable(intable = bioreptable[bioreptable[,"PATNUM"] %in% SOIall,
                                                          c("IMGDEGIS", "CTNDV50", "DUKESCORE"),drop=FALSE], align_metatable_reflist)
labeltable_numeric <- labeltable_numeric[do.call(order, labeltable_numeric),]

# Grab the original labels into a table, and sort by the sorted numeric table
labeltable <- bioreptable[rownames(labeltable_numeric), c("IMGDEGIS", "CTNDV50", "DUKESCORE")]

# ## Add on an annotation if its RNAseq and meth?
# labeltable <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "PATNUM", all = TRUE, sort = FALSE), list(
#   cbind.data.frame(PATNUM = rownames(labeltable), labeltable),
#   cbind.data.frame(PATNUM = combined_clustermembership_table[!is.na(combined_clustermembership_table[,"rnacluster"]),"PATNUM"], 
#                    rnasample = "rnasample"),
#   cbind.data.frame(PATNUM = combined_clustermembership_table[!is.na(combined_clustermembership_table[,"methcluster"]),"PATNUM"], 
#                    methsample = "methsample")
# ))
# rownames(labeltable) <- labeltable[,"PATNUM"]
# labeltable <- labeltable[rownames(labeltable_numeric),!grepl("PATNUM", colnames(labeltable))]

## Make the annotation (the main plot we want)
IMGDEGIS_color_f <- colorRampPalette(colors = c("#e4b2b2", "#A70000"))
CTNDV50_color_f <- colorRampPalette(colors = c("#cf99d9", "#8700A1"))
DUKESCORE_color_f <- colorRampPalette(colors = c("#99c7a5", "#005115"))
annotationlist1 = annotationlist_builder(labeltable, 
                                         customcolorlist = list(IMGDEGIS = IMGDEGIS_color_f(4),
                                                                CTNDV50 = CTNDV50_color_f(4),
                                                                DUKESCORE = DUKESCORE_color_f(5)
                                                                # rnasample = c("rnasample" = "black"),
                                                                # methsample = c("methsample" = "black")
                                         ))
names(annotationlist1[["IMGDEGIS"]]) <- c("None", "Mild", "Moderate", "Severe")
names(annotationlist1[["CTNDV50"]]) <- c("Non-evaluable", "1", "2", "3")
names(annotationlist1[["DUKESCORE"]]) <- c("1 Vessel with at least Moderate (>=50%) Plaque", 
                                           "2 Vessels with at least Moderate (>=50%) Plaque or 1 Severe (>=70%) Plaque", 
                                           "3 Vessels with at least Moderate (>=50%) Plaque or 2 Severe (>=70%) Plaque or Severe (>=70%) Proximal LAD Plaque", 
                                           "3 Vessels with Severe (>=70%) Plaque or 2 Severe (>=70%) Plaque Including Proximal LAD",
                                           "Left Main >=50%")
## Dummy table to work with function.
dummytab <- labeltable
dummytab[] <- 1
labelplot <- create_heatmap(counttab = dummytab, subsetnum = FALSE, scale_data = FALSE, colmetatable = NULL, colannotationlist = NULL, 
                            rowmetatable = labeltable, rowannotationlist = annotationlist1, 
                            colclusterparam = FALSE, rowclusterparam = FALSE, separate_legend = FALSE, heatmapcolorparam = NULL)
pdf(paste0(outfilepathintegration, "ischemia_test_labels_hm.pdf"), 11.5, 8, useDingbats = FALSE)
draw(labelplot[[1]])
junk <- dev.off()

# Do summary stats by datatype:
dir.create(paste0(outfilepathintegration, "cohort_tabulations/"), showWarnings = FALSE, recursive = TRUE)
## Create tabulation table, with characterize of certain cols, and selection of certain samples
tabulation_table <- merge(bioreptable[,c("PATNUM", metaCOI)], combined_clustermembership_table, by = "PATNUM")
rownames(tabulation_table) <- tabulation_table[,"PATNUM"]
colstocharacterize <- c(colnames(tabulation_table)[grepl("C_", colnames(tabulation_table))])
tabulation_table[,colstocharacterize] <- apply(tabulation_table[,colstocharacterize], 2, as.character)

summarystatfile <- paste0(outfilepathintegration, "cohort_tabulations/", "COI_tabulation_SOImeth.csv")
sumstattab_seq <- summarize_table(intable = tabulation_table[SOImeth,], groupvar = NULL, outfile = summarystatfile, calc_stats = FALSE)
summarystatfile <- paste0(outfilepathintegration, "cohort_tabulations/", "COI_tabulation_SOIrna.csv")
sumstattab_seq <- summarize_table(intable = tabulation_table[SOIrna,], groupvar = NULL, outfile = summarystatfile, calc_stats = FALSE)
summarystatfile <- paste0(outfilepathintegration, "cohort_tabulations/", "COI_tabulation_SOIall.csv")
sumstattab_seq <- summarize_table(intable = tabulation_table[SOIall,], groupvar = NULL, outfile = summarystatfile, calc_stats = FALSE)
summarystatfile <- paste0(outfilepathintegration, "cohort_tabulations/", "COI_tabulation_SOIboth.csv")
sumstattab_seq <- summarize_table(intable = tabulation_table[SOIboth,], groupvar = NULL, outfile = summarystatfile, calc_stats = FALSE)

# Now also add on the respsective clusters and rerun with tabulations by cluster
summarystatfile <- paste0(outfilepathintegration, "cohort_tabulations/", "COI_tabulation_SOImeth_methcluster.csv")
sumstattab_seq <- summarize_table(intable = tabulation_table[SOImeth,], groupvar = "methcluster", outfile = summarystatfile, calc_stats = TRUE)
summarystatfile <- paste0(outfilepathintegration, "cohort_tabulations/", "COI_tabulation_SOIrna_rnacluster.csv")
sumstattab_seq <- summarize_table(intable = tabulation_table[SOIrna,], groupvar = "rnacluster", outfile = summarystatfile, calc_stats = TRUE)




## Use this function to do an A vs not A comparison for each cluster
summary_table_pval_outlist <- summary_table_sign_outlist <- list()
for (nmfnum in seq_len(length(sort(unique(na.omit(tabulation_table[,"rnacluster"])))))) {
  nmfsel <- sort(unique(na.omit(tabulation_table[,"rnacluster"])))[nmfnum]
  tabulation_table_sel <- tabulation_table[SOIrna,]
  tabulation_table_sel[,"rnacluster"] <- ifelse(tabulation_table_sel[,"rnacluster"] == nmfsel, nmfsel, paste0("NOT_", nmfsel))
  
  summarystatfile <- paste0(outfilepathintegration, "cohort_tabulations/", "COI_tabulation_SOIrna_rnacluster_", nmfsel, ".csv")
  sumstattab_seq <- summarize_table(intable = tabulation_table_sel[SOIrna,], groupvar = "rnacluster", 
                                    outfile = summarystatfile, calc_stats = TRUE)
  
  collapsed_summary_table <- cbind(category = rep(names(sumstattab_seq), lapply(sumstattab_seq, nrow)), do.call(rbind.fill, sumstattab_seq))
  collapsed_summary_table[collapsed_summary_table[,"rnacluster_cat"] == "Min.","rnacluster_cat"] <- "NUMERIC"
  collapsed_summary_table[,"combined_catval"] <- apply(collapsed_summary_table[,c("category", "rnacluster_cat")], 1, function(x) paste(x, collapse = "__"))
  ## First one - is a pval table
  summary_table_pval <- collapsed_summary_table[!collapsed_summary_table[,"statval"] %in% c(NA, ""),c("category", "rnacluster_cat", "statval", "combined_catval")]
  summary_table_pval_outlist[[nmfnum]] <- cbind(summary_table_pval[,c("combined_catval", "statval")], cluster = nmfsel)
  names(summary_table_pval_outlist)[nmfnum] <- nmfsel
  
  ## For each category - 
  signoutlist <- list()
  for (combined_catval_num in seq_len(length(summary_table_pval[,"combined_catval"]))) {
    combined_catval_sel <- summary_table_pval[,"combined_catval"][combined_catval_num]
    if (grepl("__NUMERIC", combined_catval_sel)) {
      catsel <- gsub("NUMERIC", "Mean", combined_catval_sel)
      datasel <- collapsed_summary_table[collapsed_summary_table[,"combined_catval"] == catsel,]
      signout <- ifelse(datasel[,nmfsel] > datasel[,paste0("NOT_", nmfsel)], 1, -1)
    } else {
      datasel <- collapsed_summary_table[collapsed_summary_table[,"combined_catval"] == combined_catval_sel,]
      outratios <- as.numeric(unlist(datasel[,c(nmfsel, paste0("NOT_", nmfsel))])) / 
        as.numeric(unlist(collapsed_summary_table[collapsed_summary_table[,"combined_catval"] == "PATNUM__total",
                                                  c(nmfsel, paste0("NOT_", nmfsel))]))
      signout <- ifelse(outratios[1] > outratios[2], 1, -1)
    }
    signoutlist[[combined_catval_num]] <- data.frame(combined_catval = combined_catval_sel, cluster = nmfsel, statsign = as.numeric(signout))
    names(signoutlist)[combined_catval_num] <- combined_catval_sel
  }
  signouttable <- do.call(rbind, signoutlist)
  summary_table_sign_outlist[[nmfnum]] <- signouttable
  names(summary_table_sign_outlist)[nmfnum] <- nmfsel
  # summary_table_value <- collapsed_summary_table[!summary_table_pval[,"statval"] %in% c(NA, ""),c("category", "rnacluster_cat", "statval")]
  # summary_table_pval[summary_table_pval[,"rnacluster_cat"] == "Min.","rnacluster_cat"] <- "NUMERIC"
}
# Summary sign out table
summary_table_sign_temp <- do.call(rbind, summary_table_sign_outlist)
summary_table_sign_temp2 <- dcast(summary_table_sign_temp, combined_catval ~ cluster, value.var = "statsign")
summary_table_sign_FULL <- apply(summary_table_sign_temp2[,2:ncol(summary_table_sign_temp2)], 2, function(x) as.numeric(as.character(x)))
rownames(summary_table_sign_FULL) <- summary_table_sign_temp2[,"combined_catval"]
summary_table_sign_FULL <- summary_table_sign_FULL[summary_table_pval_outlist[[1]][,1], ]

summary_table_pval_temp <- do.call(rbind, summary_table_pval_outlist)
summary_table_pval_temp2 <- dcast(summary_table_pval_temp, combined_catval ~ cluster, value.var = "statval")
summary_table_pval_FULL <- apply(summary_table_pval_temp2[,2:ncol(summary_table_pval_temp2)], 2, function(x) as.numeric(as.character(x)))
rownames(summary_table_pval_FULL) <- summary_table_pval_temp2[,"combined_catval"]
summary_table_pval_FULL <- summary_table_pval_FULL[summary_table_pval_outlist[[1]][,1], ]
summary_table_pval_FULL_signed <- summary_table_pval_FULL * summary_table_sign_FULL

summary_table_pval_FULL_filtered <- summary_table_pval_FULL
summary_table_pval_FULL_filtered[summary_table_pval_FULL > 0.05] <- NA
summary_table_pval_FULL_filtered_signed <- summary_table_pval_FULL_filtered * summary_table_sign_FULL

# COI <- rownames(summary_table_pval_FULL_filtered_signed)
COI <- c("SEX__Female", "AGE_ENRL__NUMERIC",  "ETHNIC__Hispanic or Latino", 
         "RACE__American Indian or Alaska Native",  "RACE__Asian",  "RACE__Black or African American",  "RACE__White", "RACE__Multiple Races", 
         "BMI__NUMERIC", "HEMOGLOB__NUMERIC",  "HDLC__NUMERIC",  "LDLC__NUMERIC",  "TRIGLYC__NUMERIC",  "TOTCHOL__NUMERIC",
         "RVBPDIA__NUMERIC", "RVBPSYS__NUMERIC", 
         "MESTATIN__Yes",  "MEAPASP__Yes", 
         "SMOKSTAT__Current Smoker", "SMOKSTAT__Former Smoker", "SMOKSTAT__Never Smoked",
         "CKD__Yes",  "HYPTENSE__Yes",  "DIABETES__Yes", 
         "IMGDEGIS__Mild",  "IMGDEGIS__Moderate",  "IMGDEGIS__None",  "IMGDEGIS__Severe", 
         "CTNDV50__1",  "CTNDV50__2",  "CTNDV50__3",  "CTMULT50__Yes", 
         "DUKESCORE__1 Vessel with at least Moderate (>=50%) Plaque", 
         "DUKESCORE__2 Vessels with at least Moderate (>=50%) Plaque or 1 Severe (>=70%) Plaque", 
         "DUKESCORE__3 Vessels with at least Moderate (>=50%) Plaque or 2 Severe (>=70%) Plaque or Severe (>=70%) Proximal LAD Plaque", 
         "DUKESCORE__3 Vessels with Severe (>=70%) Plaque or 2 Severe (>=70%) Plaque Including Proximal LAD",  "DUKESCORE__Left Main >=50%", 
         "C_PRIMARY__0",  "C_PRIMARY__1", 
         "C_CVDMI_P__0",  "C_CVDMI_P__1", 
         
         # "methcluster__1_F",  "methcluster__1_M",  "methcluster__2_F",  "methcluster__2_M", 
         # "methcluster__3_F",  "methcluster__3_M",  "methcluster__4_F",  "methcluster__4_M"
         "methcluster__I", "methcluster__II", "methcluster__III", "methcluster__IV"  
         )

heatmapcolorparam <- colorRamp2(breaks = c(-0.05, -0.000000001, 0.000000001, 0.05), c("white", "darkblue", "darkred", "white"))
rowmetatable <- data.frame(category = unlist(lapply(strsplit(rownames(summary_table_pval_FULL_filtered), split = "__"), function(x) x[1])),
                                             row.names = rownames(summary_table_pval_FULL_filtered))
rowannotationlist <- annotationlist_builder(rowmetatable) 
hm1 <- create_heatmap(counttab = summary_table_pval_FULL_filtered_signed[COI,], subsetnum = FALSE, scale_data = FALSE,
               rowmetatable = rowmetatable[COI,,drop=FALSE], rowannotationlist = rowannotationlist,
               separate_legend = TRUE, heatmapcolorparam = heatmapcolorparam, addborders = TRUE)
pdf(paste0(outfilepathintegration, "cohort_tabulations/", "signed_filtered_sigfeature_rnacluster_hm.pdf"), 10, 10, useDingbats = FALSE)
draw(hm1[[1]])
junk <- dev.off()

heatmapcolorparam <- colorRamp2(breaks = c(-1, -0.000000001, 0.000000001, 1), c("white", "darkblue", "darkred", "white"))
hm1 <- create_heatmap(counttab = summary_table_pval_FULL_signed[COI,], subsetnum = FALSE, scale_data = FALSE,
                      rowmetatable = rowmetatable[COI,,drop=FALSE], rowannotationlist = rowannotationlist,
                      separate_legend = TRUE, heatmapcolorparam = heatmapcolorparam, addborders = TRUE)
pdf(paste0(outfilepathintegration, "cohort_tabulations/", "signed_sigfeature_rnacluster_hm.pdf"), 10, 10, useDingbats = FALSE)
draw(hm1[[1]])
junk <- dev.off()







## REPEAT FOR METH CLUSTER
clustersubgroups <- c("I", "II", "III", "IV")
summary_table_pval_outlist <- summary_table_sign_outlist <- list()
for (clusternum in seq_len(length(clustersubgroups))) {
  nmfsel <- clustersubgroups[clusternum]
  tabulation_table_sel <- tabulation_table[tabulation_table[,"PATNUM"] %in% SOImeth,]
  tabulation_table_sel[,"methcluster"] <- ifelse(tabulation_table_sel[,"methcluster"] == nmfsel, nmfsel, paste0("NOT_", nmfsel))
  
  summarystatfile <- paste0(outfilepathintegration, "cohort_tabulations/", "COI_tabulation_SOImethM_methcluster_", nmfsel, ".csv")
  sumstattab_seq <- summarize_table(intable = tabulation_table_sel[SOImeth,], groupvar = "methcluster", 
                                    outfile = summarystatfile, calc_stats = TRUE)
  
  collapsed_summary_table <- cbind(category = rep(names(sumstattab_seq), lapply(sumstattab_seq, nrow)), do.call(rbind.fill, sumstattab_seq))
  collapsed_summary_table[collapsed_summary_table[,"methcluster_cat"] == "Min.","methcluster_cat"] <- "NUMERIC"
  collapsed_summary_table[,"combined_catval"] <- apply(collapsed_summary_table[,c("category", "methcluster_cat")], 1, function(x) paste(x, collapse = "__"))
  ## First one - is a pval table
  summary_table_pval <- collapsed_summary_table[!collapsed_summary_table[,"statval"] %in% c(NA, ""),c("category", "methcluster_cat", "statval", "combined_catval")]
  summary_table_pval_outlist[[clusternum]] <- cbind(summary_table_pval[,c("combined_catval", "statval")], cluster = nmfsel)
  names(summary_table_pval_outlist)[clusternum] <- nmfsel
  
  ## For each category - 
  signoutlist <- list()
  for (combined_catval_num in seq_len(length(summary_table_pval[,"combined_catval"]))) {
    combined_catval_sel <- summary_table_pval[,"combined_catval"][combined_catval_num]
    if (grepl("__NUMERIC", combined_catval_sel)) {
      catsel <- gsub("NUMERIC", "Mean", combined_catval_sel)
      datasel <- collapsed_summary_table[collapsed_summary_table[,"combined_catval"] == catsel,]
      signout <- ifelse(datasel[,nmfsel] > datasel[,paste0("NOT_", nmfsel)], 1, -1)
    } else {
      datasel <- collapsed_summary_table[collapsed_summary_table[,"combined_catval"] == combined_catval_sel,]
      outratios <- as.numeric(unlist(datasel[,c(nmfsel, paste0("NOT_", nmfsel))])) / 
        as.numeric(unlist(collapsed_summary_table[collapsed_summary_table[,"combined_catval"] == "PATNUM__total",
                                                  c(nmfsel, paste0("NOT_", nmfsel))]))
      signout <- ifelse(outratios[1] > outratios[2], 1, -1)
    }
    signoutlist[[combined_catval_num]] <- data.frame(combined_catval = combined_catval_sel, cluster = nmfsel, statsign = as.numeric(signout))
    names(signoutlist)[combined_catval_num] <- combined_catval_sel
  }
  signouttable <- do.call(rbind, signoutlist)
  summary_table_sign_outlist[[clusternum]] <- signouttable
  names(summary_table_sign_outlist)[clusternum] <- nmfsel
  # summary_table_value <- collapsed_summary_table[!summary_table_pval[,"statval"] %in% c(NA, ""),c("category", "rnacluster_cat", "statval")]
  # summary_table_pval[summary_table_pval[,"rnacluster_cat"] == "Min.","rnacluster_cat"] <- "NUMERIC"
}
# Summary sign out table
summary_table_sign_temp <- do.call(rbind, summary_table_sign_outlist)
summary_table_sign_temp2 <- dcast(summary_table_sign_temp, combined_catval ~ cluster, value.var = "statsign")
summary_table_sign_FULL <- apply(summary_table_sign_temp2[,2:ncol(summary_table_sign_temp2)], 2, function(x) as.numeric(as.character(x)))
rownames(summary_table_sign_FULL) <- summary_table_sign_temp2[,"combined_catval"]
summary_table_sign_FULL <- summary_table_sign_FULL[summary_table_pval_outlist[[1]][,1], ]

summary_table_pval_temp <- do.call(rbind, summary_table_pval_outlist)
summary_table_pval_temp2 <- dcast(summary_table_pval_temp, combined_catval ~ cluster, value.var = "statval")
summary_table_pval_FULL <- apply(summary_table_pval_temp2[,2:ncol(summary_table_pval_temp2)], 2, function(x) as.numeric(as.character(x)))
rownames(summary_table_pval_FULL) <- summary_table_pval_temp2[,"combined_catval"]
summary_table_pval_FULL <- summary_table_pval_FULL[summary_table_pval_outlist[[1]][,1], ]
summary_table_pval_FULL_signed <- summary_table_pval_FULL * summary_table_sign_FULL

summary_table_pval_FULL_filtered <- summary_table_pval_FULL
summary_table_pval_FULL_filtered[summary_table_pval_FULL > 0.05] <- NA
summary_table_pval_FULL_filtered_signed <- summary_table_pval_FULL_filtered * summary_table_sign_FULL

# COI <- rownames(summary_table_pval_FULL_filtered_signed)
COI <- c("AGE_ENRL__NUMERIC",  "ETHNIC__Hispanic or Latino", 
         "RACE__American Indian or Alaska Native",  "RACE__Asian",  "RACE__Black or African American",  "RACE__White", "RACE__Multiple Races", 
         "BMI__NUMERIC", "HEMOGLOB__NUMERIC",  "HDLC__NUMERIC",  "LDLC__NUMERIC",  "TRIGLYC__NUMERIC",  "TOTCHOL__NUMERIC",
         "RVBPDIA__NUMERIC", "RVBPSYS__NUMERIC", 
         "MESTATIN__Yes",  "MEAPASP__Yes", 
         "SMOKSTAT__Current Smoker", "SMOKSTAT__Former Smoker", "SMOKSTAT__Never Smoked",
         "CKD__Yes",  "HYPTENSE__Yes",  "DIABETES__Yes", 
         "IMGDEGIS__Mild",  "IMGDEGIS__Moderate",  "IMGDEGIS__None",  "IMGDEGIS__Severe", 
         "CTNDV50__1",  "CTNDV50__2",  "CTNDV50__3",  "CTMULT50__Yes", 
         "DUKESCORE__1 Vessel with at least Moderate (>=50%) Plaque", 
         "DUKESCORE__2 Vessels with at least Moderate (>=50%) Plaque or 1 Severe (>=70%) Plaque", 
         "DUKESCORE__3 Vessels with at least Moderate (>=50%) Plaque or 2 Severe (>=70%) Plaque or Severe (>=70%) Proximal LAD Plaque", 
         "DUKESCORE__3 Vessels with Severe (>=70%) Plaque or 2 Severe (>=70%) Plaque Including Proximal LAD",  "DUKESCORE__Left Main >=50%", 
         "C_PRIMARY__0",  "C_PRIMARY__1", 
         "C_CVDMI_P__0",  "C_CVDMI_P__1",
         "rnacluster__nmf_cluster_1", "rnacluster__nmf_cluster_2", "rnacluster__nmf_cluster_3A","rnacluster__nmf_cluster_3A", "rnacluster__nmf_cluster_4"  
)

heatmapcolorparam <- colorRamp2(breaks = c(-0.05, -0.000000001, 0.000000001, 0.05), c("white", "darkblue", "darkred", "white"))
rowmetatable <- data.frame(category = unlist(lapply(strsplit(rownames(summary_table_pval_FULL_filtered), split = "__"), function(x) x[1])),
                           row.names = rownames(summary_table_pval_FULL_filtered))
rowannotationlist <- annotationlist_builder(rowmetatable) 
hm1 <- create_heatmap(counttab = summary_table_pval_FULL_filtered_signed[COI,], subsetnum = FALSE, scale_data = FALSE,
                      rowmetatable = rowmetatable[COI,,drop=FALSE], rowannotationlist = rowannotationlist,
                      separate_legend = TRUE, heatmapcolorparam = heatmapcolorparam, addborders = TRUE)
pdf(paste0(outfilepathintegration, "cohort_tabulations/", "signed_filtered_sigfeature_methcluster_hm.pdf"), 10, 10, useDingbats = FALSE)
draw(hm1[[1]])
junk <- dev.off()

heatmapcolorparam <- colorRamp2(breaks = c(-1, -0.000000001, 0.000000001, 1), c("white", "darkblue", "darkred", "white"))
hm1 <- create_heatmap(counttab = summary_table_pval_FULL_signed[COI,], subsetnum = FALSE, scale_data = FALSE,
                      rowmetatable = rowmetatable[COI,,drop=FALSE], rowannotationlist = rowannotationlist,
                      separate_legend = TRUE, heatmapcolorparam = heatmapcolorparam, addborders = TRUE)
pdf(paste0(outfilepathintegration, "cohort_tabulations/", "signed_sigfeature_methclustersummary_table_pval_outlist_hm.pdf"), 10, 10, useDingbats = FALSE)
draw(hm1[[1]])
junk <- dev.off()










# Figure 5 - NMF cluster viz with WGCNA Modules:
## Read in eigengene table
GOI <- rownames(genestocolorstab[genestocolorstab[,"moduleColors"] != "grey",]) ## Only plotting on genes that are NOT grey

# Now lets make a summary - where we have a split heatmap - with modules labeled and split on rows, and then sample labeled and split on columns
hmplottab <- normcounttable[GOI,]
colsplitparam <- factor(combined_clustermembership_table[SOIrna, "rnacluster"],
                        levels = c("nmf_cluster_1", "nmf_cluster_2", "nmf_cluster_3A", "nmf_cluster_3B", "nmf_cluster_4"))
rowsplitparam <- factor(genestocolorstab[GOI,"moduleColors"], levels = c("magenta", "salmon", "purple", "midnightblue", "turquoise", 
                                                                         "yellow", "pink", "red", 
                                                                         "blue", "tan", "brown",
                                                                         "greenyellow", "black", "cyan", "green"))

sampleannottable <- merge(metatable[SOIrna,c("SEX", "AGE_CAT", "IMGDEGIS", "CTNDV50", "DUKESCORE")], 
                          combined_clustermembership_table[SOIrna, "rnacluster", drop=FALSE], by = "row.names")
rownames(sampleannottable) <- sampleannottable[,"Row.names"]
sampleannottable <- sampleannottable[SOIrna,!grepl("Row.names", colnames(sampleannottable))]

## Make the annotation (the main plot we want)
IMGDEGIS_color_f <- colorRampPalette(colors = c("#e4b2b2", "#A70000"))
CTNDV50_color_f <- colorRampPalette(colors = c("#cf99d9", "#8700A1"))
DUKESCORE_color_f <- colorRampPalette(colors = c("#99c7a5", "#005115"))
AGECAT_color_f <- colorRampPalette(colors = c("#e0d1b2", "#946300"))
sampleannotationlist = annotationlist_builder(sampleannottable, 
                                         customcolorlist = list(SEX = c("Male" = "brown", "Female" = "green"),
                                                                AGE_CAT = AGECAT_color_f(3),
                                                                IMGDEGIS = IMGDEGIS_color_f(4),
                                                                CTNDV50 = CTNDV50_color_f(4),
                                                                DUKESCORE = DUKESCORE_color_f(5)
                                                                ))
names(sampleannotationlist[["AGE_CAT"]]) <- c("33_63", "63_70", "70_93")
names(sampleannotationlist[["IMGDEGIS"]]) <- c("None", "Mild", "Moderate", "Severe")
names(sampleannotationlist[["CTNDV50"]]) <- c("Non-evaluable", "1", "2", "3")
names(sampleannotationlist[["DUKESCORE"]]) <- c("1 Vessel with at least Moderate (>=50%) Plaque", 
                                           "2 Vessels with at least Moderate (>=50%) Plaque or 1 Severe (>=70%) Plaque", 
                                           "3 Vessels with at least Moderate (>=50%) Plaque or 2 Severe (>=70%) Plaque or Severe (>=70%) Proximal LAD Plaque", 
                                           "3 Vessels with Severe (>=70%) Plaque or 2 Severe (>=70%) Plaque Including Proximal LAD",
                                           "Left Main >=50%")


geneannottable <- genestocolorstab[GOI,"moduleColors", drop=FALSE]
genecustomcolorlist <- list(moduleColors = unique(genestocolorstab[,"moduleColors"]))
names(genecustomcolorlist[["moduleColors"]]) <- unique(genestocolorstab[,"moduleColors"])
geneannotationlist <- annotationlist_builder(geneannottable, customcolorlist = genecustomcolorlist)

## Top annotation
temp1 <- vector("list", length(sampleannotationlist))
names(temp1) = names(sampleannotationlist)
annotlegendlist = lapply(temp1, function(x) x[[1]] =
                           list(title_gp=gpar(fontsize=5, fontface="bold"), labels_gp=gpar(fontsize=4)))
## Define a param that will go through each annotation - and keep a legend if its continuous or has less than 10 discrete terms, otherwise hide the legend
showlegendparam = unname(unlist(lapply(sampleannotationlist, function(x) {
  numterms = tryCatch(length(na.omit(unique(x))), error=function(e) NULL)
  is.null(numterms) || numterms <= 10})))
hatop = HeatmapAnnotation(df = sampleannottable,
                          col = sampleannotationlist,
                          na_col = "white",
                          show_annotation_name = TRUE,
                          annotation_name_gp = gpar(fontsize = 8, fontface="bold"),
                          annotation_name_side = "left",
                          simple_anno_size = unit(min(60/length(sampleannotationlist), 5),"mm"),
                          show_legend = TRUE,
                          annotation_legend_param = annotlegendlist)


## Side annotation
## Define parameters for each of the labels on the annotation bars
temp1 <- vector("list", length(geneannotationlist))
names(temp1) = names(geneannotationlist)
annotlegendlist = lapply(temp1, function(x)
  x[[1]] = list(title_gp=gpar(fontsize=8, fontface="bold"), labels_gp=gpar(fontsize=8)))

## Define a param that will go through each annotation - and keep a legend if its continuous or has less than 10 discrete terms, otherwise hide the legend
showlegendparam = unname(unlist(lapply(geneannotationlist, function(x) {
  numterms = tryCatch(length(na.omit(unique(x))), error=function(e) NULL)
  is.null(numterms) || numterms <= 10})))

## Look for any empty annotations - fill them with white, and later, make sure to hide their legend
emptyannots = names(sapply(geneannotationlist, length)[sapply(geneannotationlist, length)==0])
if (length(emptyannots) > 0){
  for (i in 1:length(emptyannots)) {
    temp1 = "white"
    names(temp1) = emptyannots[i]
    geneannotationlist[[emptyannots[i]]] = temp1
  }
  showlegendparam[which(names(geneannotationlist) %in% emptyannots)] = FALSE
}
## Add param that will bolden the side annotation bars if it is <100, and omit the grid lines if more
if (nrow(geneannottable) < 100) {
  sideannotation_linebold_param <- gpar(fontsize = 0.5)} else {sideannotation_linebold_param <- NULL}
haside = rowAnnotation(df = geneannottable,
                       col = geneannotationlist,
                       na_col = "white",
                       # gp = gpar(fontsize = 0.01),
                       gp = sideannotation_linebold_param,
                       show_annotation_name=TRUE,
                       annotation_name_gp = gpar(fontsize = 8, fontface="bold"),
                       annotation_name_side = "top",
                       simple_anno_size = unit(min(60/length(geneannotationlist), 5),"mm"),
                       show_legend = TRUE,
                       annotation_legend_param = annotlegendlist)


hmplottab = as.matrix(t(apply(hmplottab, 1, function(x) zscore(x))))
heatmapcolorparam <- colorRamp2(breaks = c(4, 0, -4), c("darkred", "white", "darkblue"))
summary_outhm <- Heatmap(as.matrix(hmplottab),
                         col = heatmapcolorparam,
                         row_title = "Genes",column_title = "Samples",
                         row_split = rowsplitparam, column_split = colsplitparam,
                         top_annotation = hatop,
                         
                         cluster_columns = TRUE, cluster_rows = TRUE,
                         cluster_row_slices = FALSE, cluster_column_slices = FALSE,
                         
                         show_row_dend = FALSE, show_column_dend = FALSE,
                         
                         # rect_gp = gpar(col = "black", lwd = 0.5),
                         
                         show_column_names = FALSE, column_names_gp = gpar(fontsize = 6), 
                         show_row_names = FALSE,
                         
                         heatmap_legend_param = list(
                           title = "Zscore",
                           title_gp = gpar(fontsize = 8, fontface = "bold")),
                         height = unit(min((nrow(hmplottab)/2), 12),"cm"),
                         width = unit(min(ncol(hmplottab), 18),"cm")
                         
)

pdf(paste0(outfilepathintegration, "nmf_eigen_meta_summary_heatmap.pdf"), useDingbats = FALSE, width = 20, height = 15)
draw(summary_outhm + haside)
junk <- dev.off()










# ## For each eigengenecluster - we want to see if the clusters have different values
dir.create(paste0(outfilepathintegration, "wgcna_nmf_integration/"), showWarnings = FALSE, recursive = TRUE)
# 
# ## The cleanest way to run this is over a list of cohort TABLES, and then work with each of those
COItablist <- list(
    nmfcluster_allsamples = rna_clustermembership_table[,"combined",drop=FALSE]
#   nmfcluster_coresamples = nmf_clustercores_table[,"combined",drop=FALSE],
#   hemoglobin = addonmetatable[,"HEMOGLOB", drop=FALSE],
#   platelet = addonmetatable[,"PLATELET", drop=FALSE],
#   wbc = addonmetatable[,"WBC", drop=FALSE],
#   IMGDEGIS = addonmetatable[,"IMGDEGIS", drop=FALSE],
#   IMGDEGIS_01_2_3 = addonmetatable[,"IMGDEGIS_01_2_3", drop=FALSE],
#   IMGDEGIS_01_23 = addonmetatable[,"IMGDEGIS_01_23", drop=FALSE],
#   CLUSTER1_IMGDEGIS_01_2_3 = addonmetatable[,"CLUSTER1_IMGDEGIS_01_2_3",drop=FALSE],
#   CLUSTER2_IMGDEGIS_01_2_3 = addonmetatable[,"CLUSTER2_IMGDEGIS_01_2_3",drop=FALSE],
#   CLUSTER3_IMGDEGIS_01_2_3 = addonmetatable[,"CLUSTER3_IMGDEGIS_01_2_3",drop=FALSE],
#   CLUSTER4_IMGDEGIS_01_2_3 = addonmetatable[,"CLUSTER4_IMGDEGIS_01_2_3",drop=FALSE],
#   CTNDV50 = addonmetatable[,"CTNDV50",drop=FALSE]
)
# 
## Run a for loop over all cohort tables

ALL_wilcox_statoutlist <- ALL_spearman_statoutlist <- ALL_bptablist <- list()
for (COItabnum in seq_len(length(COItablist))) {
  ## For each cohort, grab the table we need
  COItab_sel <- COItablist[[COItabnum]]
  COItab_label <- names(COItablist)[COItabnum]

  COIbpoutpath <- paste0(outfilepathintegration, "wgcna_nmf_integration/boxplots_", COItab_label, "/")
  dir.create(COIbpoutpath, showWarnings = FALSE, recursive = TRUE)

  wilcox_statoutlist <- spearman_statoutlist <- bptablist <- eigenranklist <- list()
  for (eigengenenum in seq_len(length(eigengenes))) {
    ## Select the category to analyze
    eigen_sel <- eigengenes[eigengenenum]
    wgcna_sel <- eigengenecounttable[,eigen_sel,drop=FALSE]
    # nmf_sel <- nmf_clustermembership_table[,"combined",drop=FALSE]
    bptab <- na.omit(merge(wgcna_sel, COItab_sel, by = "row.names")[,c(3,1,2)])
    bptab[,2] <- c("cohort")

    ## IF THE FIRST COL IS CONT (like for labs) - then split into quartiles
    rawspearmanstatout <- NULL
    if (is.numeric(bptab[,1])) {
      ## Also sneak in an actual corr with the raw values here jsut to see
      scattertab <- bptab[,c(1,3)]
      scattertab[,3] <- gsub("ME", "", colnames(scattertab)[2])

      scatterout <- scatter_plotter(indata = scattertab[,c(1,2)], colorvar = scattertab[,3,drop=FALSE],
                                    labsparam = list(title = paste0(eigen_sel, " values for samples by ", COItab_label), x = COItab_label, y = eigen_sel),
                                    plotstats = TRUE)
      pdf(paste0(COIbpoutpath, eigen_sel, "_rawval_scatter.pdf"), useDingbats = FALSE)
      print(scatterout)
      junk <- dev.off()

      ## And capture this spearman value
      rawspearmanstatval <- cor.test(scattertab[,1], scattertab[,2], method = "spearman", exact = FALSE)
      rawspearmanstatout <- c(paste0(COItab_label, "_raw"), eigen_sel, rawspearmanstatval$p.value, rawspearmanstatval$estimate)

      ## Create the bptab
      cattab <- cut2(t(bptab[,1]),g=4)
      levels(cattab)[match(levels(cattab)[1],levels(cattab))] <- paste0(COItab_label, "_Q1min")
      levels(cattab)[match(levels(cattab)[2],levels(cattab))] <- paste0(COItab_label, "_Q2")
      levels(cattab)[match(levels(cattab)[3],levels(cattab))] <- paste0(COItab_label, "_Q3")
      levels(cattab)[match(levels(cattab)[4],levels(cattab))] <- paste0(COItab_label, "_Q4max")
      bptab[,1] <- as.character(cattab)
    }

    ## Pre calc the avg per group so we can sort in that order:
    avgeigenvalue <- aggregate(bptab[,3], by = list(bptab[,1]), mean)
    avgeigenvalue <- avgeigenvalue[order(avgeigenvalue[,2]),]
    avgeigenvalue[,"Eigen_Score_Rank"] <- seq(1:nrow(avgeigenvalue))
    colnames(avgeigenvalue) <- c("Group", "Eigengene", "Eigen_Score_Rank")
    eigenranklist[[eigengenenum]] <- avgeigenvalue

    ## OUTPUT THIS
    bpout <- boxplot_plotter(boxplottable = bptab, xsplit = "category",
                             labsparam = list(title = paste0(eigen_sel, " values for samples by ", COItab_label), x = COItab_label, y = eigen_sel,
                                              catorder = "cohort", featorder = avgeigenvalue[,1]),
                             plotstats = "intra", testtypeparam = "wilcox.test"
    )
    pdf(paste0(COIbpoutpath, eigen_sel, "_boxplot.pdf"), useDingbats = FALSE)
    print(bpout)
    junk <- dev.off()
    ## Save out the bptab
    bptablist[[eigengenenum]] <- bptab


    scattertab <- bptab[,c(1,3)]
    scattertab[,3] <- avgeigenvalue[,3][match(scattertab[,1], avgeigenvalue[,1])]
    scattertab[,4] <- gsub("ME", "", colnames(scattertab)[2])

    scatterlabelparam <- as.vector(avgeigenvalue[,1])
    names(scatterlabelparam) <- avgeigenvalue[,3]
    scatterout <- scatter_plotter(indata = scattertab[,c(3,2)], colorvar = scattertab[,4,drop=FALSE],
                                  labsparam = list(title = paste0(eigen_sel, " values for samples by ", COItab_label), x = COItab_label, y = eigen_sel),
                                  plotstats = TRUE, addjitterparam = 0.1)
    scatterout <- scatterout + scale_x_continuous(breaks = avgeigenvalue[,3], labels = avgeigenvalue[,1])
    # scatterout <- scatterout + geom_point(alpha = 0.7, shape = 16, position = position_jitter(width = 0.1))
    pdf(paste0(COIbpoutpath, eigen_sel, "_scatter.pdf"), useDingbats = FALSE)
    print(scatterout)
    junk <- dev.off()

    # Grab the stats into a table
    # Create the combinatorial table for all combos of the groups
    combtab = combn(as.character(unique(bptab[!is.na(bptab[,1]),1])), 2, simplify=F)
    # Apply a functiona cross each combo using the bptab and pulling out each group
    wilcoxstatsout <- lapply(combtab, function(x){
      group1 <- bptab[bptab[,1] %in% x[1], 3]
      group2 <- bptab[bptab[,1] %in% x[2], 3]
      out1 <- c(COItab_label, x[1], x[2], eigen_sel, wilcox.test(group1, group2)$p.value)
      out2 <- cohens_d(group1, group2)
      c(out1, out2)
    })
    wilcox_statoutlist[[eigengenenum]] <- do.call(rbind, wilcoxstatsout)

    # Repeat for the spearman correlation


    spearmanstatout <- cor.test(scattertab[,3], scattertab[,2], method = "spearman", exact = FALSE)
    if (!is.null(rawspearmanstatout)) { ## Attach the raw spearman corr out as well if applicable
      spearman_statoutlist[[eigengenenum]] <- rbind(rawspearmanstatout, c(paste0(COItab_label, "_quartile"), eigen_sel, spearmanstatout$p.value, spearmanstatout$estimate))
    } else {
      spearman_statoutlist[[eigengenenum]] <- c(COItab_label, eigen_sel, spearmanstatout$p.value, spearmanstatout$estimate)
    }

    ## Name all of our outlists
    names(bptablist)[eigengenenum] <- names(wilcox_statoutlist)[eigengenenum] <-
      names(spearman_statoutlist)[eigengenenum] <- names(eigenranklist)[eigengenenum] <- eigen_sel

  }
  # Cat our lists
  wilcox_summarytable <- do.call(rbind, wilcox_statoutlist)
  colnames(wilcox_summarytable) <- c("Cohort", "Feature1", "Feature2", "Eigengene", "wilcox_pval", "cohens_d")
  spearman_summarytable <- do.call(rbind, spearman_statoutlist)
  colnames(spearman_summarytable) <- c("Cohort", "Eigengene", "spearman_pval", "spearman_rval")
  eigenrank_summarytable <- do.call(rbind, eigenranklist)
  eigenrank_summarytable[,"eigengene"] <- gsub("\\..*", "", rownames(eigenrank_summarytable))
  eigenrank_summarytable <- eigenrank_summarytable[order(eigenrank_summarytable[,"Group"], eigenrank_summarytable[,"Eigen_Score_Rank"], decreasing = TRUE),]

  # Write out the results
  write.table(wilcox_summarytable, paste0(COIbpoutpath, COItab_label,"_wilcox_summary_table.csv"),
              sep = ",", col.names = TRUE, row.names = FALSE)
  write.table(spearman_summarytable, paste0(COIbpoutpath, COItab_label, "_spearman_summary_table.csv"),
              sep = ",", col.names = TRUE, row.names = FALSE)
  write.table(eigenrank_summarytable, paste0(COIbpoutpath, COItab_label, "_eigenrank_summary_table.csv"),
              sep = ",", col.names = TRUE, row.names = FALSE)

  # Specifically for the eigenrank_summarytable - I want to output a quick HM
  eigenrank_hmplottab <-dcast(eigenrank_summarytable, eigengene ~ Group, value.var = "Eigen_Score_Rank")
  rownames(eigenrank_hmplottab) <- eigenrank_hmplottab[,"eigengene"]
  eigenrank_hmplottab <- eigenrank_hmplottab[do.call(order, eigenrank_hmplottab[,2:ncol(eigenrank_hmplottab)]),2:ncol(eigenrank_hmplottab)]
  eigenrank_hmplot <- create_heatmap(counttab = eigenrank_hmplottab, subsetnum = FALSE, scale_data = FALSE,
                                     rowclusterparam = FALSE, colclusterparam = FALSE)
  pdf(paste0(COIbpoutpath, COItab_label, "_eigenrank_summary_heatmap.pdf"), useDingbats = FALSE)
  draw(eigenrank_hmplot[[1]])
  junk <- dev.off()
  write.table(eigenrank_summarytable, paste0(COIbpoutpath, COItab_label, "_eigenrank_summary_table.csv"), sep = ",", row.names = TRUE, col.names = NA)

  #     ## Save out the table
  #     bptablist[[eigengenenum]] <- bptab

  ALL_wilcox_statoutlist[[COItabnum]] <- wilcox_summarytable
  names(ALL_wilcox_statoutlist[COItabnum]) <- COItab_label
  ALL_spearman_statoutlist[[COItabnum]] <- spearman_summarytable
  names(ALL_spearman_statoutlist[COItabnum]) <- COItab_label

}
ALL_wilcox_sumarytable <- do.call(rbind, ALL_wilcox_statoutlist)
ALL_spearman_sumarytable <- do.call(rbind, ALL_spearman_statoutlist)


## Lets get a viz of this to see what it looks like... (maybe with the WGCNA custom func? Could be useful)
# Need to unmelt tables to fo this
temptab <- data.frame(ALL_wilcox_sumarytable)
temptab[,"combined_label"] <- apply(temptab[,c("Cohort", "Feature1","Feature2")], 1, function(x) paste(x, collapse = "_"))

wilcoxpval_plottab <- dcast(temptab[,c("combined_label", "Eigengene", "wilcox_pval")], Eigengene ~ combined_label, value.var = "wilcox_pval")
rownames(wilcoxpval_plottab) <- wilcoxpval_plottab[,"Eigengene"]
wilcoxpval_plottab <- wilcoxpval_plottab[,!grepl("Eigengene", colnames(wilcoxpval_plottab))]
keepeigennames <- rownames(wilcoxpval_plottab)
wilcoxpval_plottab <- t(apply(wilcoxpval_plottab, 2, function(x) as.numeric(as.character(x))))
colnames(wilcoxpval_plottab) <- keepeigennames

wilcoxdval_plottab <- dcast(temptab[,c("combined_label", "Eigengene", "cohens_d")], Eigengene ~ combined_label, value.var = "cohens_d")
rownames(wilcoxdval_plottab) <- wilcoxdval_plottab[,"Eigengene"]
wilcoxdval_plottab <- wilcoxdval_plottab[,!grepl("Eigengene", colnames(wilcoxdval_plottab))]
wilcoxdval_plottab <- t(apply(wilcoxdval_plottab, 2, function(x) as.numeric(as.character(x))))
colnames(wilcoxdval_plottab) <- keepeigennames

wilcoxhmout <- WGCNA_custom_heatmap(moduleTraitPvalue = wilcoxpval_plottab, moduleTraitCor = wilcoxdval_plottab)
pdf(paste0(outfilepathintegration, "wgcna_nmf_integration/", "wilcox_values.pdf"), 10, 25, useDingbats = FALSE)
print(wilcoxhmout)
junk <- dev.off()




## Lets get a viz of this to see what it looks like... (maybe with the WGCNA custom func? Could be useful)
# Need to unmelt tables to fo this
temptab <- data.frame(ALL_spearman_sumarytable)

spearmanpval_plottab <- dcast(temptab[,c("Cohort", "Eigengene", "spearman_pval")], Eigengene ~ Cohort, value.var = "spearman_pval")
rownames(spearmanpval_plottab) <- spearmanpval_plottab[,"Eigengene"]
spearmanpval_plottab <- spearmanpval_plottab[,!grepl("Eigengene", colnames(spearmanpval_plottab)),drop=FALSE]
keepeigennames <- rownames(spearmanpval_plottab)
spearmanpval_plottab <- t(apply(spearmanpval_plottab, 2, function(x) as.numeric(as.character(x))))
colnames(spearmanpval_plottab) <- keepeigennames

spearmanrval_plottab <- dcast(temptab[,c("Cohort", "Eigengene", "spearman_rval")], Eigengene ~ Cohort, value.var = "spearman_rval")
rownames(spearmanrval_plottab) <- spearmanrval_plottab[,"Eigengene"]
spearmanrval_plottab <- spearmanrval_plottab[,!grepl("Eigengene", colnames(spearmanrval_plottab)),drop=FALSE]
spearmanrval_plottab <- t(apply(spearmanrval_plottab, 2, function(x) as.numeric(as.character(x))))
colnames(spearmanrval_plottab) <- keepeigennames

spearmanhmout <- WGCNA_custom_heatmap(moduleTraitPvalue = spearmanpval_plottab, moduleTraitCor = spearmanrval_plottab)
pdf(paste0(outfilepathintegration, "wgcna_nmf_integration/", "spearman_values.pdf"), 10, 20, useDingbats = FALSE)
print(spearmanhmout)
junk <- dev.off()








## Heatmap of eigengene values - with annotation of nmf cluster (and maybe annotation for nmf cluster driver for modules??)
rnacluster_eigen_ranktable_temp <- read.table("# PATH_UPDATED: output/run3_rmoutliers2/integrative_analyses/run3/wgcna_nmf_integration/boxplots_nmfcluster_allsamples/nmfcluster_allsamples_eigenrank_summary_table.csv", sep = ",", header = TRUE, row.names = 1)
# Rank table
rnacluster_eigen_ranktable_temp2 <- dcast(rnacluster_eigen_ranktable_temp, Group ~ eigengene, value.var = "Eigen_Score_Rank")
rnacluster_eigen_ranktable <- data.frame(t(rnacluster_eigen_ranktable_temp2[,2:ncol(rnacluster_eigen_ranktable_temp2)]))
colnames(rnacluster_eigen_ranktable) <- rnacluster_eigen_ranktable_temp2[,"Group"]
# Avg eigengene table
rnacluster_eigen_avgtable_temp2 <- dcast(rnacluster_eigen_ranktable_temp, Group ~ eigengene, value.var = "Eigengene")
rnacluster_eigen_avgtable <- data.frame(t(rnacluster_eigen_avgtable_temp2[,2:ncol(rnacluster_eigen_avgtable_temp2)]))
colnames(rnacluster_eigen_avgtable) <- rnacluster_eigen_avgtable_temp2[,"Group"]


## If I want to only do drivers:
clusterOrder <- c("nmf_cluster_1", "nmf_cluster_2", "nmf_cluster_3A", "nmf_cluster_3B", "nmf_cluster_4")
moduleOrder <- c("turquoise", "salmon", "magenta", "midnightblue", "purple",
                 "red", "yellow", "pink",
                 "blue", "tan", "brown",
                 "green", "black", "greenyellow", "cyan")
plottab <- rnacluster_eigen_avgtable[paste0("ME", moduleOrder),clusterOrder] ## Need to remove the grey eigengene - not important
rowmeta <- cbind(rnacluster_eigen_ranktable[rownames(plottab),],
                 data.frame(moduleColor = gsub("ME", "", rownames(plottab)), row.names = rownames(plottab)))
genecustomcolorlist <- list(moduleColor = gsub("ME", "", rownames(plottab)))
names(genecustomcolorlist[["moduleColor"]]) <- gsub("ME", "", rownames(plottab))
rowannotationlist <- annotationlist_builder(rowmeta, customcolorlist = genecustomcolorlist)
heatmapcolorparam <- colorRamp2(
  c(-0.05, 0, 0.05), c("blue", "white", "red"))
outhm1 <- create_heatmap(counttab = plottab, scale_data = FALSE, separate_legend = TRUE,
                         colmetatable = NULL, colannotationlist = NULL,
                         rowmetatable = rowmeta, rowannotationlist = rowannotationlist,
                         colclusterparam = FALSE, rowclusterparam = FALSE, heatmapcolorparam = heatmapcolorparam, addborders = TRUE)
outheatmapfile1 = paste0(outfilepathintegration, "nmf_eigenavg_heatmap.pdf")
pdf(outheatmapfile1, useDingbats = FALSE, height = 7, width = 13)
draw(outhm1$heatmap)
junk <- dev.off()


























# ### Maybe do a survival analysis for EVERY EVENT with our clusters??? Maybe... Worth a shot?
# # Do a survival analysis for each cohort over every event....
# # And I guess we use the max # of samples per cohort... thats probably the most fair.
# outfilepathsurvival = paste0(outfilepathintegration, "nmf_survival_analysis/run3_nmfcluster_IMGDEGIS_nmfclustervrest_wlabs/")
# dir.create(outfilepathsurvival, recursive = TRUE, showWarnings = FALSE)
# 
# EOI <- gsub("^C_", "", colnames(bioreptable)[grepl("^C_", colnames(bioreptable))])
# EOIcols <- apply(expand.grid(c("T_", "C_"), EOI, stringsAsFactors = FALSE), 1, paste, collapse = "")
# 
# ## I am also going to compare each cluster against the rest of samples as well - that may also be interesting (and make the HR more understandable)
# expand_rnacomps <- make_comparison_columns(rna_clustermembership_table[,"combined",drop=FALSE])
# expand_methcomps <- make_comparison_columns(meth_clustermembership_table_wsex[,"cluster_wsex",drop=FALSE])
# 
# ## Also want to output the labs data to see if its predictive for anything (hopefully not....)
# CBCmetatable_quartiled <- bioreptable[,c("HEMOGLOB", "PLATELET", "WBC")]
# for (metavarnum in seq_len(ncol(CBCmetatable_quartiled))) {
#   colsel <- CBCmetatable_quartiled[,metavarnum,drop=FALSE]
#   collable <- colnames(colsel)
#   cattab <- cut2(t(colsel[,1]),g=4)
#   levels(cattab)[match(levels(cattab)[1],levels(cattab))] <- paste0(collable, "_Q1min")
#   levels(cattab)[match(levels(cattab)[2],levels(cattab))] <- paste0(collable, "_Q2")
#   levels(cattab)[match(levels(cattab)[3],levels(cattab))] <- paste0(collable, "_Q3")
#   levels(cattab)[match(levels(cattab)[4],levels(cattab))] <- paste0(collable, "_Q4max")
#   CBCmetatable_quartiled[,metavarnum] <- as.character(cattab)
# }
# expand_CBCmetatable_quartiled <- make_comparison_columns(na.omit(CBCmetatable_quartiled))
# 
# IMGDEGIS_addontable <- bioreptable[,"IMGDEGIS", drop=FALSE]
# IMGDEGIS_addontable[,"IMGDEGIS_01_2_3"] <- ifelse(IMGDEGIS_addontable[,"IMGDEGIS"] %in% c("None", "Mild"), "NoneMild", 
#                                            ifelse(IMGDEGIS_addontable[,"IMGDEGIS"] %in% "Moderate", "Moderate",
#                                            ifelse(IMGDEGIS_addontable[,"IMGDEGIS"] %in% "Severe", "Severe",
#                                            NA)))
# IMGDEGIS_addontable[,"IMGDEGIS_01_3"] <- ifelse(IMGDEGIS_addontable[,"IMGDEGIS"] %in% c("None", "Mild"), "NoneMild",
#                                          ifelse(IMGDEGIS_addontable[,"IMGDEGIS"] %in% "Severe", "Severe",
#                                          NA))
# 
# CTNDV50_addontable <- bioreptable[,"CTNDV50", drop=FALSE]
# CTNDV50_addontable[,"CTNDV50_123_clean"] <- ifelse(CTNDV50_addontable[,"CTNDV50"] %in% "1", "1", 
#                                   ifelse(CTNDV50_addontable[,"CTNDV50"] %in% "2", "2",
#                                   ifelse(CTNDV50_addontable[,"CTNDV50"] %in% "3", "3",
#                                   NA)))
# CTNDV50_addontable[,"CTNDV50_13_clean"] <- ifelse(CTNDV50_addontable[,"CTNDV50"] %in% "1", "1",
#                                   ifelse(CTNDV50_addontable[,"CTNDV50"] %in% "3", "3",
#                                   NA))
# 
# 
# combined_clustermembership_table[,"rna_AND_methcluster"] <- apply(combined_clustermembership_table[,c("rnacluster", "methcluster")], 1, function(x) paste0(x, collapse = "_AND_"))
# 
# COItablist <- list(
# 
#   IMGDEGIS = IMGDEGIS_addontable[,"IMGDEGIS", drop=FALSE],
#   IMGDEGIS_01_2_3 = IMGDEGIS_addontable[,"IMGDEGIS_01_2_3", drop=FALSE],
#   IMGDEGIS_01_3 = IMGDEGIS_addontable[,"IMGDEGIS_01_3", drop=FALSE],
#   
#   CTNDV50 = CTNDV50_addontable[,"CTNDV50", drop=FALSE],
#   CTNDV50_123_clean = CTNDV50_addontable[,"CTNDV50_123_clean", drop=FALSE],
#   CTNDV50_13_clean = CTNDV50_addontable[,"CTNDV50_13_clean", drop=FALSE],
#   
#   nmfcluster_allsamples = rna_clustermembership_table[,"combined",drop=FALSE],
#   nmfcluster1vrest = data.frame(expand_rnacomps[,"combined_nmf_cluster_1"]),
#   nmfcluster2vrest = data.frame(expand_rnacomps[,"combined_nmf_cluster_2"]),
#   nmfcluster3vrest = data.frame(expand_rnacomps[,"combined_nmf_cluster_3"]),
#   nmfcluster4vrest = data.frame(expand_rnacomps[,"combined_nmf_cluster_4"]),
#   
#   methcluster_allsamples = meth_clustermembership_table_wsex[,"cluster_wsex", drop = FALSE],
#   methcluster1Mvrest = data.frame(expand_methcomps[,"cluster_wsex_1_M"]),
#   methcluster1Fvrest = data.frame(expand_methcomps[,"cluster_wsex_1_F"]),
#   methcluster2Mvrest = data.frame(expand_methcomps[,"cluster_wsex_2_M"]),
#   methcluster2Fvrest = data.frame(expand_methcomps[,"cluster_wsex_2_F"]),
#   methcluster3Mvrest = data.frame(expand_methcomps[,"cluster_wsex_3_M"]),
#   methcluster3Fvrest = data.frame(expand_methcomps[,"cluster_wsex_3_F"]),
#   methcluster4Mvrest = data.frame(expand_methcomps[,"cluster_wsex_4_M"]),
#   methcluster4Fvrest = data.frame(expand_methcomps[,"cluster_wsex_4_F"]),
#   
#   # hemoglobin = addonmetatable[,"HEMOGLOB", drop=FALSE],
#   # platelet = addonmetatable[,"PLATELET", drop=FALSE],
#   # wbc = addonmetatable[,"WBC", drop=FALSE],
#   # CLUSTER1_IMGDEGIS_01_2_3 = addonmetatable[,"CLUSTER1_IMGDEGIS_01_2_3",drop=FALSE],
#   # CLUSTER2_IMGDEGIS_01_2_3 = addonmetatable[,"CLUSTER2_IMGDEGIS_01_2_3",drop=FALSE],
#   # CLUSTER3_IMGDEGIS_01_2_3 = addonmetatable[,"CLUSTER3_IMGDEGIS_01_2_3",drop=FALSE],
#   # CLUSTER4_IMGDEGIS_01_2_3 = addonmetatable[,"CLUSTER4_IMGDEGIS_01_2_3",drop=FALSE],
#   # CTNDV50 = addonmetatable[,"CTNDV50",drop=FALSE]
#   
#   HEMOGLOBQ1vrest = data.frame(expand_CBCmetatable_quartiled[,"HEMOGLOB_HEMOGLOB_Q1min"]),
#   HEMOGLOBQ2vrest = data.frame(expand_CBCmetatable_quartiled[,"HEMOGLOB_HEMOGLOB_Q2"]),
#   HEMOGLOBQ3vrest = data.frame(expand_CBCmetatable_quartiled[,"HEMOGLOB_HEMOGLOB_Q3"]),
#   HEMOGLOBQ4vrest = data.frame(expand_CBCmetatable_quartiled[,"HEMOGLOB_HEMOGLOB_Q4max"]),
#   PLATELETQ1vrest = data.frame(expand_CBCmetatable_quartiled[,"PLATELET_PLATELET_Q1min"]),
#   PLATELETQ2vrest = data.frame(expand_CBCmetatable_quartiled[,"PLATELET_PLATELET_Q2"]),
#   PLATELETQ3vrest = data.frame(expand_CBCmetatable_quartiled[,"PLATELET_PLATELET_Q3"]),
#   PLATELETQ4vrest = data.frame(expand_CBCmetatable_quartiled[,"PLATELET_PLATELET_Q4max"]),
#   WBCQ1vrest = data.frame(expand_CBCmetatable_quartiled[,"WBC_WBC_Q1min"]),
#   WBCQ2vrest = data.frame(expand_CBCmetatable_quartiled[,"WBC_WBC_Q2"]),
#   WBCQ3vrest = data.frame(expand_CBCmetatable_quartiled[,"WBC_WBC_Q3"]),
#   WBCQ4vrest = data.frame(expand_CBCmetatable_quartiled[,"WBC_WBC_Q4max"])
# )
# 
# ## Ok - so we want to iterate over every GOI, and every event, and see what comes out
# fullpvaloutlist <- fullcumhazoutlist <- list()
# for (COInum in seq_len(length(COItablist))) {
#   ## Select COI
#   COIsel <- COItablist[[COInum]]
#   COIlabel <- names(COItablist)[COInum]
#   
#   survival_intable <- na.omit(merge(bioreptable[,c("PATNUM", EOIcols)], cbind(PATNUM = rownames(COIsel), COIsel), by = "PATNUM"))
#   colnames(survival_intable)[ncol(survival_intable)] <- COIlabel
#   
#   survpvallist <- survHRlist <- list()
#   for (eventnum in seq_len(length(EOI))) {
#     ## Select event
#     eventsel <- EOI[eventnum]
#     
#     ## Create the subsurvivaltab
#     survivaldata <- data.frame(na.omit(survival_intable[,c(paste0(c("T_", "C_"), eventsel), COIlabel)]), row.names = survival_intable[,"PATNUM"])
#     survivaldata <- survivaldata[!survivaldata[,3] %in% c(NA, "") &
#                                  !survivaldata[,2] %in% c(NA, "NotAvailable") & !survivaldata[,1] %in% c(NA, "NotAvailable"),]
#     survivaldata[,1] <- as.numeric(survivaldata[,1])
#     survivaldata[,2] <- as.numeric(survivaldata[,2])
#     
#     ## Ok, for the competing events - we have to clean those out cause it screws with the analysis - just remove I guess?
#     survivaldata <- survivaldata[!survivaldata[,grepl("^C_", colnames(survivaldata))] %in% 2,]
#     
#     # ## Apply a max time of 5 years - it cleans up the data a little - NOPE, dont need it
#     # maxtimelimit <- 5*365
#     # survivaldata <- survivaldata[survivaldata[,grepl("^T_", colnames(survivaldata))] < maxtimelimit,]
#     
#     ## Need to properly factorize our cohort data
#     if (COIlabel %in% c("nmfcluster_allsamples")) {survivaldata[,3] <- factor(survivaldata[,3], levels = c(
#       "nmf_cluster_1", "nmf_cluster_2", "nmf_cluster_3", "nmf_cluster_4"))}
#     if (COIlabel %in% c("methcluster_allsamples")) {survivaldata[,3] <- factor(survivaldata[,3], levels = c(
#       "1_F", "1_M", "2_F", "2_M", "3_F", "3_M", "4_F", "4_M"))}
#     if (COIlabel %in% c("IMGDEGIS")) {survivaldata[,3] <- factor(survivaldata[,3], levels = c("None", "Mild", "Moderate", "Severe"))}
#     if (COIlabel %in% c("IMGDEGIS_01_2_3")) {survivaldata[,3] <- factor(survivaldata[,3], levels = c("NoneMild", "Moderate", "Severe"))}
#     if (COIlabel %in% c("IMGDEGIS_01_3")) {survivaldata[,3] <- factor(survivaldata[,3], levels = c("NoneMild", "Severe"))}
#     if (COIlabel %in% c("CTNDV50")) {survivaldata[,3] <- factor(survivaldata[,3], levels = c("Non-evaluable", "1", "2", "3"))}
#     if (COIlabel %in% c("CTNDV50_123_clean")) {survivaldata[,3] <- factor(survivaldata[,3], levels = c("1", "2", "3"))}
#     if (COIlabel %in% c("CTNDV50_13_clean")) {survivaldata[,3] <- factor(survivaldata[,3], levels = c("1", "3"))}
#     if (COIlabel %in% c("nmfcluster1vrest", "nmfcluster2vrest", "nmfcluster3vrest", "nmfcluster4vrest")) {
#       cohortlabels <- unique(survivaldata[,3])
#       survivaldata[,3] <- factor(survivaldata[,3], levels = c(as.character(cohortlabels[grepl("not_", cohortlabels)]), 
#                                                               as.character(cohortlabels[!grepl("not_", cohortlabels)])))
#     }
#     if (COIlabel %in% c("methcluster1Mvrest", "methcluster1Fvrest", "methcluster2Mvrest", "methcluster2Fvrest",
#                         "methcluster3Mvrest", "methcluster3Fvrest","methcluster4Mvrest", "methcluster4Fvrest")) {
#       cohortlabels <- unique(survivaldata[,3])
#       survivaldata[,3] <- factor(survivaldata[,3], levels = c(as.character(cohortlabels[grepl("not_", cohortlabels)]), 
#                                                               as.character(cohortlabels[!grepl("not_", cohortlabels)])))
#     }
#     # if (COIlabel %in% c("HEMOGLOBQ1vrest", "HEMOGLOBQ2vrest", "HEMOGLOBQ3vrest", "HEMOGLOBQ4vrest",
#     #                     "PLATELETQ1vrest", "PLATELETQ2vrest", "PLATELETQ3vrest", "PLATELETQ4vrest",
#     #                      WBCQ1vrest", "WBCQ2vrest", "WBCQ3vrest", "WBCQ4vrest")) {
#     #     cohortlabels <- unique(survivaldata[,3])
#     #     survivaldata[,3] <- factor(survivaldata[,3], levels = c(as.character(cohortlabels[grepl("not_", cohortlabels)]), 
#     #                                                             as.character(cohortlabels[!grepl("not_", cohortlabels)])))
#     # }
#     
#     ## Run the analysis
#     survivalanalysisout <- create_survival_plot(survivaldata = survivaldata, timebreakparam = NULL, ylimitparam = c(0.25,1))
#     # survivalanalysisout <- create_survival_plot(survivaldata = survivaldata, timebreakparam = NULL)
#     outsurvtable <- survivalanalysisout$outsurvtable
#     outsurvplot <- survivalanalysisout$outsurvplot
#     outsurvpvalue <- survivalanalysisout$outsurvpvalue
#     
#     ## Save out the pvalue of the log-rank test
#     survpvallist[[eventnum]] <- outsurvpvalue[,"pval"]
#     names(survpvallist)[eventnum] <- eventsel
#     
#     # Ok, I think we need to return the HR as well for coxph to get some kind of effect size for this analysis
#     outcoxphobject <- survivalanalysisout$outcoxphobject
#     ## Put a hack in here to get the orientation correct. If theres only one comp, then coerce to a 1x3 dataframe:
#     coxtempout <- summary(outcoxphobject)[["conf.int"]][,c("exp(coef)", "lower .95", "upper .95")]
#     if(is.null(dim(coxtempout))){
#       coxtempout <- t(data.frame(coxtempout))
#       rownames(coxtempout) <- levels(survivaldata[,3])[2] ## I can do this because I know the first level is always the ref and I set this earlier above ^^
#     } 
#     hrvalue <- cbind(coxtempout, Event = eventsel, Cohort = COIlabel)
#     
#     ## Save out the coxph HR of the coxph test
#     survHRlist[[eventnum]] <- hrvalue
#     names(survHRlist)[eventnum] <- eventsel
#     
#     ## Write out the plot and table
#     outsubdir <- paste0(outfilepathsurvival, COIlabel, "/", eventsel, "/")
#     dir.create(outsubdir, showWarnings = FALSE, recursive = TRUE)
#     
#     write.table(outsurvtable, paste0(outsubdir, COIlabel, "_", eventsel, "_survival_stat_table.csv"),
#                 sep = ",", col.names = TRUE, row.names = FALSE)
#     pdf(paste0(outsubdir, COIlabel, "_", eventsel, "_survival_plot.pdf"), width = 8, height = 5, onefile=FALSE)
#     print(outsurvplot)
#     junk <- dev.off()
#     
#   }
#   survpvaltab <- do.call(rbind, survpvallist)
#   colnames(survpvaltab) <- paste0(COIlabel, "_pval")
#   survhazardtab <- do.call(rbind, survHRlist)
#   colnames(survhazardtab) <- paste0(COIlabel, "_", colnames(survhazardtab))
#   
#   fullpvaloutlist[[COInum]] <- survpvaltab
#   names(fullpvaloutlist[COInum]) <- paste0(COIlabel)
#   fullcumhazoutlist[[COInum]] <- survhazardtab
#   names(fullcumhazoutlist)[COInum] <- paste0(COIlabel)
#   
# }
# fullpvalouttable <- do.call(cbind, fullpvaloutlist)
# write.table(fullpvalouttable, paste0(outfilepathsurvival, "full_survival_pval_table.csv"), sep = ",", col.names = NA, row.names = TRUE)
# fullcumhazouttable <- data.frame(do.call(rbind, fullcumhazoutlist))
# colnames(fullcumhazouttable) <- c("HR", "lower_CI", "upper_CI", "Event", "Cohort")
# fullcumhazouttable[,c("HR", "lower_CI", "upper_CI")] <- apply(fullcumhazouttable[,c("HR", "lower_CI", "upper_CI")], 2, function(x) 
#   as.numeric(as.character(x)))
# write.table(fullcumhazouttable, paste0(outfilepathsurvival, "full_survival_cumhaz_table.csv"), sep = ",", col.names = NA, row.names = TRUE)
# 
# 
# ## Summary heatmap of the pvals
# # fullstatouttablefile <- "# PATH_UPDATED: output/run1c_rmoutliers2/integrative_analyses/nmf_survival_analysis/full_survival_pval_table.csv"
# # fullstatouttable <- read.table(fullstatouttablefile, sep = ",", header = TRUE, row.names = 1)
# fullstatouttable <- fullpvalouttable
# 
# ## Turn the event outcome into a heatmap
# eventpvalcutoff <- 0.05
# # eventpvalcutoff <- 0.99
# plotCOI <- c("IMGDEGIS_01_3", "CTNDV50_13_clean",
#              "nmfcluster1vrest", "nmfcluster2vrest", "nmfcluster3vrest", "nmfcluster4vrest",
#              "methcluster1Mvrest", "methcluster1Fvrest", "methcluster2Mvrest", "methcluster2Fvrest",
#              "methcluster3Mvrest", "methcluster3Fvrest", "methcluster4Mvrest", "methcluster4Fvrest")
# 
# # GOI heatmap with events
# ## Lets create a heatmap with the GOI, and then use a side annotation for if its indicative of events
# # eventpvaltabfile <- "/Users/tosh/Desktop/Ruggles_Lab/projects/newman-pace-2020/output_hg19_newcontrols/run1_nooutlier1/survival_analysis/topbotquartile/full_survival_stat_table.csv"
# # outeventpvaltab <- read.table(eventpvaltabfile, sep = ",", header = TRUE, row.names = 1)
# 
# ## TOSH - selecting just nmf for now, the cell stuff is important, but not worth it atm
# outeventpvaltab <- fullstatouttable[,paste0(plotCOI, "_pval")]
# 
# ## Create the initial maptab
# # maptab <- data.frame(outeventpvaltab[rowSums(outeventpvaltab[,EOIvec] < eventpvalcutoff) > 0, EOIvec])
# maptab <- data.frame(outeventpvaltab)
# 
# ## MULTIPLE HYPOTHESIS CORRECTION - OMITTING FOR NOW BECAUSE WE HAVE WAY TOO MANY BAD HYPOTHESES AND ITS CLOGGING IT UP
# # p.adjust(unlist(maptab), method = "fdr") ## Sanity check - it works
# # “holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”
# 
# # correctedmat <- matrix(p.adjust(as.vector(as.matrix(maptab)), method='fdr'),ncol=ncol(maptab))
# # rownames(correctedmat) <- rownames(maptab)
# # colnames(correctedmat) <- colnames(maptab)
# # maptab <- correctedmat
# 
# # Create a heatmap of just the event pvalues
# EOIvec <- colnames(outeventpvaltab)
# ## Add in geneordering
# # geneorder <- indeseqtable[rownames(outeventpvaltab)[rowSums(outeventpvaltab[,EOIvec] < eventpvalcutoff) > 0],]
# # geneorder <- rownames(geneorder[order(geneorder[,"padj"], decreasing = FALSE),]) ## order by padj
# # geneorder <- rownames(geneorder[order(geneorder[,"log2FoldChange"], decreasing = TRUE),]) ## order by log2fc
# # geneorder <- rownames(correctedmat[order(indeseqtable[rownames(correctedmat),"log2FoldChange"], decreasing = TRUE),]) ## order by log2fc
# ## Order by avg pval across all events
# # geneorder <- names(rowMeans(correctedmat[rownames(geneorder),])[order(rowMeans(correctedmat[rownames(geneorder),]), decreasing = FALSE)])
# # geneorder <- names(rowMeans(correctedmat[order(rowMeans(correctedmat), decreasing = FALSE),]))
# geneorder <- rownames(outeventpvaltab[order(outeventpvaltab[,1], decreasing = FALSE),])
# 
# # maptab <- data.frame(outeventpvaltab[rowSums(outeventpvaltab[,EOIvec] < eventpvalcutoff) > 0, EOIvec])
# 
# 
# # maptab <- maptab[geneorder,][rowSums(maptab[geneorder,] < eventpvalcutoff) > 0,]
# 
# ## Add in coloring here based on directionality:
# hazardmaptabtemp <- dcast(data.frame(fullcumhazouttable[fullcumhazouttable[,"Cohort"] %in% plotCOI, c("Event", "Cohort", "HR")]),
#                           Event ~ Cohort, value.var = "HR")
# rownames(hazardmaptabtemp) <- hazardmaptabtemp[,1]
# # hazardmaptab <- hazardmaptabtemp[rownames(maptab),!grepl("Event", colnames(hazardmaptabtemp))]
# hazardmaptabtemp <- hazardmaptabtemp[rownames(maptab),gsub("_pval", "", colnames(maptab))]
# # hazardmaptab <- hazardmaptabtemp[rownames(maptab),gsub("_pval", "", colnames(maptab))]
# hazardmaptab <- apply(hazardmaptabtemp, 2, as.numeric)
# rownames(hazardmaptab) <- hazardmaptabtemp[,1]
# hazardsigntab <- ifelse(hazardmaptab > 1, 1, -1)
# 
# maptab <- maptab * hazardsigntab
# 
# heatmapcolorparam = colorRamp2(c(-1, -eventpvalcutoff - 0.0001, -eventpvalcutoff, -0.0001, 0, 0.0001, eventpvalcutoff, eventpvalcutoff + 0.0001, 1), 
#                                c("grey", "grey", "#bcb2fd", "#11007d", "white", "#7d0000", "#fd9999", "grey", "grey"))
# heatmapcolorLEGEND = colorRamp2(c(-eventpvalcutoff, 0, eventpvalcutoff), c("blue", "white", "red"))
# 
# 
# 
# # heatmapcolorparam = colorRamp2(c(1, eventpvalcutoff + 0.0001, eventpvalcutoff, 0), c("white", "white", "#c9e0dc", "#00705E"))
# # heatmapcolorLEGEND = colorRamp2(c(eventpvalcutoff, 0), c("#c9e0dc", "#00705E"))
# 
# ## LEAVE IN ALL EVENTS:
# # heatmapcolorparam = colorRamp2(c(1, 0), c("white", "#00705E"))
# # heatmapcolorLEGEND = colorRamp2(c(1, 0), c("#white", "#00705E"))
# 
# # heatmapcolorparam = colorRamp2(c(1, 0.1 + 0.0001, 0.1, 0), c("white", "white", "#c9e0dc", "#00705E"))
# # heatmapcolorLEGEND = colorRamp2(c(0.1, 0), c("#c9e0dc", "#00705E"))
# 
# ht1 = Heatmap(as.matrix(maptab), 
#               col = heatmapcolorparam,    ## Define the color scale for the heatmap
#               row_title = "Events", column_title = "Cohorts",
#               border = TRUE, na_col = "white", rect_gp = gpar(col = "black", lwd = 0.5),
#               
#               cluster_columns = FALSE, cluster_rows = FALSE,clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
#               #"ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", centroid"
#               
#               show_column_names = TRUE, column_names_gp = gpar(fontsize = 6), 
#               show_row_names = TRUE, row_names_side = "left", row_names_gp = gpar(fontsize=6),
#               show_row_dend = TRUE, show_column_dend = TRUE,
#               
#               heatmap_legend_param = list(
#                 col_fun = heatmapcolorLEGEND,
#                 title = "pvalue", legend_height = unit(2.5, "cm"), title_gp = gpar(fontsize = 8, fontface = "bold")),
#               height = unit(min((nrow(maptab)/2), 12),"cm"), width = unit(min(ncol(maptab), 18),"cm")
#               
# )
# draw(ht1)
# 
# 
# hmoutfile <- paste0(outfilepathsurvival, "event_heatmap_005.pdf")
# # hmoutfile <- paste0(outfilepathsurvival, "event_heatmap_010.pdf")
# # hmoutfile <- paste0(outfilepathsurvival, "event_heatmap_099.pdf")
# # hmoutfile <- paste0(outfilepathsurvival, "event_heatmap_all_005.pdf")
# # hmoutfile <- paste0(outfilepathsurvival, "event_heatmap_all_005_wcolor.pdf")
# pdf(hmoutfile, 11, 9)
# draw(ht1)
# junk <- dev.off()














# ## Heatmap of eigengene values - with annotation of nmf cluster (and maybe annotation for nmf cluster driver for modules??)
# # wgcna_eigengenes_table ## table we need
# rnacluster_avgeigen_table <- aggregate(merge(rna_clustermembership_table[,"combined",drop=FALSE], 
#                                              eigengenecounttable, by = "row.names")[,colnames(eigengenecounttable)], 
#                                        by = list(rnacluster_eigen_ranktable[,"combined"]), mean)
# rownames(rnacluster_avgeigen_table) <- rnacluster_avgeigen_table[,"Group.1"]
# rnacluster_avgeigen_table <- t(rnacluster_avgeigen_table[,!colnames(rnacluster_avgeigen_table) %in% "Group.1"])
# rnacluster_eigenrank_table <- t(apply(tt1[,colnames(eigengenecounttable)], 2, function(x) rank(-x)))
# colnames(rnacluster_eigenrank_table) <- colnames(rnacluster_avgeigen_table)
# # rnacluster_eigenrank_table <- rnacluster_eigenrank_table[do.call(order, data.frame(rnacluster_eigenrank_table)),]
# 
# # geneannottable <- genestocolorstab[GOI,"moduleColors", drop=FALSE]
# 
# # geneannotationlist <- annotationlist_builder(geneannottable, customcolorlist = genecustomcolorlist)
# 
# 
# 
# ## If I want to only do drivers:
# clusterOrder = c("nmf_cluster_1", "nmf_cluster_2","nmf_cluster_3", "nmf_cluster_4")
# moduleOrder <- c("turquoise", "salmon", "magenta", "midnightblue", "purple", 
#                  "red", "yellow", "pink", 
#                  "blue", "tan", "brown",  
#                  "green", "black", "greenyellow", "cyan")
# plottab <- rnacluster_avgeigen_table[paste0("ME", moduleOrder),clusterOrder] ## Need to remove the grey eigengene - not important
# rowmeta <- cbind(rnacluster_eigenrank_table[rownames(plottab),], 
#                  data.frame(moduleColor = gsub("ME", "", rownames(plottab)), row.names = rownames(plottab)))
# genecustomcolorlist <- list(moduleColor = gsub("ME", "", rownames(plottab)))
# names(genecustomcolorlist[["moduleColor"]]) <- gsub("ME", "", rownames(plottab))
# rowannotationlist <- annotationlist_builder(rowmeta, customcolorlist = genecustomcolorlist)
# heatmapcolorparam <- colorRamp2(
#   c(-0.05, 0, 0.05), c("blue", "white", "red"))
# outhm1 <- create_heatmap(counttab = plottab, scale_data = FALSE, separate_legend = TRUE,
#                          colmetatable = NULL, colannotationlist = NULL,
#                          rowmetatable = rowmeta, rowannotationlist = rowannotationlist, 
#                          colclusterparam = FALSE, rowclusterparam = FALSE, heatmapcolorparam = heatmapcolorparam, addborders = TRUE)
# outheatmapfile1 = paste0(outfilepathintegration, "nmf_eigenavg_heatmap.pdf")
# pdf(outheatmapfile1, useDingbats = FALSE, height = 7, width = 13)
# draw(outhm1$heatmap)
# junk <- dev.off()







# ## For each eigengenecluster - we want to see if the clusters have different values
# dir.create(paste0(outfilepathintegration, "wgcna_nmf_integration/"), showWarnings = FALSE, recursive = TRUE)
# 
# ## The cleanest way to run this is over a list of cohort TABLES, and then work with each of those
# COItablist <- list(
#   nmfcluster_allsamples = nmf_clustermembership_table[,"combined",drop=FALSE],
#   nmfcluster_coresamples = nmf_clustercores_table[,"combined",drop=FALSE],
#   hemoglobin = addonmetatable[,"HEMOGLOB", drop=FALSE],
#   platelet = addonmetatable[,"PLATELET", drop=FALSE],
#   wbc = addonmetatable[,"WBC", drop=FALSE],
#   IMGDEGIS = addonmetatable[,"IMGDEGIS", drop=FALSE],
#   IMGDEGIS_01_2_3 = addonmetatable[,"IMGDEGIS_01_2_3", drop=FALSE],
#   IMGDEGIS_01_23 = addonmetatable[,"IMGDEGIS_01_23", drop=FALSE],
#   CLUSTER1_IMGDEGIS_01_2_3 = addonmetatable[,"CLUSTER1_IMGDEGIS_01_2_3",drop=FALSE],
#   CLUSTER2_IMGDEGIS_01_2_3 = addonmetatable[,"CLUSTER2_IMGDEGIS_01_2_3",drop=FALSE],
#   CLUSTER3_IMGDEGIS_01_2_3 = addonmetatable[,"CLUSTER3_IMGDEGIS_01_2_3",drop=FALSE],
#   CLUSTER4_IMGDEGIS_01_2_3 = addonmetatable[,"CLUSTER4_IMGDEGIS_01_2_3",drop=FALSE],
#   CTNDV50 = addonmetatable[,"CTNDV50",drop=FALSE]
# )
# 
# ## Run a for loop over all cohort tables
# 
# ALL_wilcox_statoutlist <- ALL_spearman_statoutlist <- ALL_bptablist <- list()
# for (COItabnum in seq_len(length(COItablist))) {
#   ## For each cohort, grab the table we need
#   COItab_sel <- COItablist[[COItabnum]]
#   COItab_label <- names(COItablist)[COItabnum]
#   
#   COIbpoutpath <- paste0(outfilepathintegration, "wgcna_nmf_integration/boxplots_", COItab_label, "/")
#   dir.create(COIbpoutpath, showWarnings = FALSE, recursive = TRUE)
#   
#   wilcox_statoutlist <- spearman_statoutlist <- bptablist <- eigenranklist <- list()
#   for (eigengenenum in seq_len(length(eigengenes))) {
#     ## Select the category to analyze
#     eigen_sel <- eigengenes[eigengenenum]
#     wgcna_sel <- wgcna_eigengenes_table[,eigen_sel,drop=FALSE]
#     # nmf_sel <- nmf_clustermembership_table[,"combined",drop=FALSE]
#     bptab <- na.omit(merge(wgcna_sel, COItab_sel, by = "row.names")[,c(3,1,2)])
#     bptab[,2] <- c("cohort")
#     
#     ## IF THE FIRST COL IS CONT (like for labs) - then split into quartiles
#     rawspearmanstatout <- NULL
#     if (is.numeric(bptab[,1])) {
#       ## Also sneak in an actual corr with the raw values here jsut to see
#       scattertab <- bptab[,c(1,3)]
#       scattertab[,3] <- gsub("ME", "", colnames(scattertab)[2])
#       
#       scatterout <- scatter_plotter(indata = scattertab[,c(1,2)], colorvar = scattertab[,3,drop=FALSE], 
#                                     labsparam = list(title = paste0(eigen_sel, " values for samples by ", COItab_label), x = COItab_label, y = eigen_sel), 
#                                     plotstats = TRUE)
#       pdf(paste0(COIbpoutpath, eigen_sel, "_rawval_scatter.pdf"), useDingbats = FALSE)
#       print(scatterout)
#       junk <- dev.off()
#       
#       ## And capture this spearman value
#       rawspearmanstatval <- cor.test(scattertab[,1], scattertab[,2], method = "spearman", exact = FALSE)
#       rawspearmanstatout <- c(paste0(COItab_label, "_raw"), eigen_sel, rawspearmanstatval$p.value, rawspearmanstatval$estimate)
#       
#       ## Create the bptab
#       cattab <- cut2(t(bptab[,1]),g=4)
#       levels(cattab)[match(levels(cattab)[1],levels(cattab))] <- paste0(COItab_label, "_Q1min")
#       levels(cattab)[match(levels(cattab)[2],levels(cattab))] <- paste0(COItab_label, "_Q2")
#       levels(cattab)[match(levels(cattab)[3],levels(cattab))] <- paste0(COItab_label, "_Q3")
#       levels(cattab)[match(levels(cattab)[4],levels(cattab))] <- paste0(COItab_label, "_Q4max")
#       bptab[,1] <- as.character(cattab)
#     }
#     
#     ## Pre calc the avg per group so we can sort in that order:
#     avgeigenvalue <- aggregate(bptab[,3], by = list(bptab[,1]), mean)
#     avgeigenvalue <- avgeigenvalue[order(avgeigenvalue[,2]),]
#     avgeigenvalue[,"Eigen_Score_Rank"] <- seq(1:nrow(avgeigenvalue))
#     colnames(avgeigenvalue) <- c("Group", "Eigengene", "Eigen_Score_Rank")
#     eigenranklist[[eigengenenum]] <- avgeigenvalue
#     
#     ## OUTPUT THIS
#     
#     
#     bpout <- boxplot_plotter(boxplottable = bptab, xsplit = "category", 
#                              labsparam = list(title = paste0(eigen_sel, " values for samples by ", COItab_label), x = COItab_label, y = eigen_sel, 
#                                               catorder = "cohort", featorder = avgeigenvalue[,1]), 
#                              plotstats = "intra", testtypeparam = "wilcox.test"
#     )
#     pdf(paste0(COIbpoutpath, eigen_sel, "_boxplot.pdf"), useDingbats = FALSE)
#     print(bpout)
#     junk <- dev.off()
#     ## Save out the bptab
#     bptablist[[eigengenenum]] <- bptab
#     
#     
#     scattertab <- bptab[,c(1,3)]
#     scattertab[,3] <- avgeigenvalue[,3][match(scattertab[,1], avgeigenvalue[,1])]
#     scattertab[,4] <- gsub("ME", "", colnames(scattertab)[2])
#     
#     scatterlabelparam <- as.vector(avgeigenvalue[,1])
#     names(scatterlabelparam) <- avgeigenvalue[,3]
#     scatterout <- scatter_plotter(indata = scattertab[,c(3,2)], colorvar = scattertab[,4,drop=FALSE], 
#                                   labsparam = list(title = paste0(eigen_sel, " values for samples by ", COItab_label), x = COItab_label, y = eigen_sel), 
#                                   plotstats = TRUE, addjitterparam = 0.1)
#     scatterout <- scatterout + scale_x_continuous(breaks = avgeigenvalue[,3], labels = avgeigenvalue[,1])
#     # scatterout <- scatterout + geom_point(alpha = 0.7, shape = 16, position = position_jitter(width = 0.1))
#     pdf(paste0(COIbpoutpath, eigen_sel, "_scatter.pdf"), useDingbats = FALSE)
#     print(scatterout)
#     junk <- dev.off()
#     
#     # Grab the stats into a table
#     # Create the combinatorial table for all combos of the groups
#     combtab = combn(as.character(unique(bptab[!is.na(bptab[,1]),1])), 2, simplify=F)
#     # Apply a functiona cross each combo using the bptab and pulling out each group
#     wilcoxstatsout <- lapply(combtab, function(x){
#       group1 <- bptab[bptab[,1] %in% x[1], 3]
#       group2 <- bptab[bptab[,1] %in% x[2], 3]
#       out1 <- c(COItab_label, x[1], x[2], eigen_sel, wilcox.test(group1, group2)$p.value)
#       out2 <- cohens_d(group1, group2)
#       c(out1, out2)
#     })
#     wilcox_statoutlist[[eigengenenum]] <- do.call(rbind, wilcoxstatsout)
#     
#     # Repeat for the spearman correlation
#     
#     
#     spearmanstatout <- cor.test(scattertab[,3], scattertab[,2], method = "spearman", exact = FALSE)
#     if (!is.null(rawspearmanstatout)) { ## Attach the raw spearman corr out as well if applicable
#       spearman_statoutlist[[eigengenenum]] <- rbind(rawspearmanstatout, c(paste0(COItab_label, "_quartile"), eigen_sel, spearmanstatout$p.value, spearmanstatout$estimate))
#     } else {
#       spearman_statoutlist[[eigengenenum]] <- c(COItab_label, eigen_sel, spearmanstatout$p.value, spearmanstatout$estimate)
#     }
#     
#     ## Name all of our outlists
#     names(bptablist)[eigengenenum] <- names(wilcox_statoutlist)[eigengenenum] <- 
#       names(spearman_statoutlist)[eigengenenum] <- names(eigenranklist)[eigengenenum] <- eigen_sel
#     
#   }
#   # Cat our lists
#   wilcox_summarytable <- do.call(rbind, wilcox_statoutlist)
#   colnames(wilcox_summarytable) <- c("Cohort", "Feature1", "Feature2", "Eigengene", "wilcox_pval", "cohens_d")
#   spearman_summarytable <- do.call(rbind, spearman_statoutlist)
#   colnames(spearman_summarytable) <- c("Cohort", "Eigengene", "spearman_pval", "spearman_rval")
#   eigenrank_summarytable <- do.call(rbind, eigenranklist)
#   eigenrank_summarytable[,"eigengene"] <- gsub("\\..*", "", rownames(eigenrank_summarytable))
#   eigenrank_summarytable <- eigenrank_summarytable[order(eigenrank_summarytable[,"Group"], eigenrank_summarytable[,"Eigen_Score_Rank"], decreasing = TRUE),]
#   
#   # Write out the results
#   write.table(wilcox_summarytable, paste0(COIbpoutpath, COItab_label,"_wilcox_summary_table.csv"),
#               sep = ",", col.names = TRUE, row.names = FALSE)
#   write.table(spearman_summarytable, paste0(COIbpoutpath, COItab_label, "_spearman_summary_table.csv"),
#               sep = ",", col.names = TRUE, row.names = FALSE)
#   write.table(eigenrank_summarytable, paste0(COIbpoutpath, COItab_label, "_eigenrank_summary_table.csv"),
#               sep = ",", col.names = TRUE, row.names = FALSE)
#   
#   # Specifically for the eigenrank_summarytable - I want to output a quick HM
#   eigenrank_hmplottab <-dcast(eigenrank_summarytable, eigengene ~ Group, value.var = "Eigen_Score_Rank")
#   rownames(eigenrank_hmplottab) <- eigenrank_hmplottab[,"eigengene"]
#   eigenrank_hmplottab <- eigenrank_hmplottab[do.call(order, eigenrank_hmplottab[,2:ncol(eigenrank_hmplottab)]),2:ncol(eigenrank_hmplottab)]
#   eigenrank_hmplot <- create_heatmap(counttab = eigenrank_hmplottab, subsetnum = FALSE, scale_data = FALSE, 
#                                      rowclusterparam = FALSE, colclusterparam = FALSE)
#   pdf(paste0(COIbpoutpath, COItab_label, "_eigenrank_summary_heatmap.pdf"), useDingbats = FALSE)
#   draw(eigenrank_hmplot[[1]])
#   junk <- dev.off()
#   write.table(eigenrank_summarytable, paste0(COIbpoutpath, COItab_label, "_eigenrank_summary_table.csv"), sep = ",", row.names = TRUE, col.names = NA)
#   
#   #     ## Save out the table
#   #     bptablist[[eigengenenum]] <- bptab
#   
#   ALL_wilcox_statoutlist[[COItabnum]] <- wilcox_summarytable
#   names(ALL_wilcox_statoutlist[COItabnum]) <- COItab_label
#   ALL_spearman_statoutlist[[COItabnum]] <- spearman_summarytable
#   names(ALL_spearman_statoutlist[COItabnum]) <- COItab_label
#   
# }
# ALL_wilcox_sumarytable <- do.call(rbind, ALL_wilcox_statoutlist)
# ALL_spearman_sumarytable <- do.call(rbind, ALL_spearman_statoutlist)
# 
# 
# ## Lets get a viz of this to see what it looks like... (maybe with the WGCNA custom func? Could be useful)
# # Need to unmelt tables to fo this
# temptab <- data.frame(ALL_wilcox_sumarytable)
# temptab[,"combined_label"] <- apply(temptab[,c("Cohort", "Feature1","Feature2")], 1, function(x) paste(x, collapse = "_"))
# 
# wilcoxpval_plottab <- dcast(temptab[,c("combined_label", "Eigengene", "wilcox_pval")], Eigengene ~ combined_label, value.var = "wilcox_pval")
# rownames(wilcoxpval_plottab) <- wilcoxpval_plottab[,"Eigengene"]
# wilcoxpval_plottab <- wilcoxpval_plottab[,!grepl("Eigengene", colnames(wilcoxpval_plottab))]
# keepeigennames <- rownames(wilcoxpval_plottab)
# wilcoxpval_plottab <- t(apply(wilcoxpval_plottab, 2, function(x) as.numeric(as.character(x))))
# colnames(wilcoxpval_plottab) <- keepeigennames
# 
# wilcoxdval_plottab <- dcast(temptab[,c("combined_label", "Eigengene", "cohens_d")], Eigengene ~ combined_label, value.var = "cohens_d")
# rownames(wilcoxdval_plottab) <- wilcoxdval_plottab[,"Eigengene"]
# wilcoxdval_plottab <- wilcoxdval_plottab[,!grepl("Eigengene", colnames(wilcoxdval_plottab))]
# wilcoxdval_plottab <- t(apply(wilcoxdval_plottab, 2, function(x) as.numeric(as.character(x))))
# colnames(wilcoxdval_plottab) <- keepeigennames
# 
# wilcoxhmout <- WGCNA_custom_heatmap(moduleTraitPvalue = wilcoxpval_plottab, moduleTraitCor = wilcoxdval_plottab)
# pdf(paste0(outfilepathintegration, "wgcna_nmf_integration/", "wilcox_values.pdf"), 10, 25, useDingbats = FALSE)
# print(wilcoxhmout)
# junk <- dev.off()
# 
# 
# 
# 
# ## Lets get a viz of this to see what it looks like... (maybe with the WGCNA custom func? Could be useful)
# # Need to unmelt tables to fo this
# temptab <- data.frame(ALL_spearman_sumarytable)
# 
# spearmanpval_plottab <- dcast(temptab[,c("Cohort", "Eigengene", "spearman_pval")], Eigengene ~ Cohort, value.var = "spearman_pval")
# rownames(spearmanpval_plottab) <- spearmanpval_plottab[,"Eigengene"]
# spearmanpval_plottab <- spearmanpval_plottab[,!grepl("Eigengene", colnames(spearmanpval_plottab))]
# keepeigennames <- rownames(spearmanpval_plottab)
# spearmanpval_plottab <- t(apply(spearmanpval_plottab, 2, function(x) as.numeric(as.character(x))))
# colnames(spearmanpval_plottab) <- keepeigennames
# 
# spearmanrval_plottab <- dcast(temptab[,c("Cohort", "Eigengene", "spearman_rval")], Eigengene ~ Cohort, value.var = "spearman_rval")
# rownames(spearmanrval_plottab) <- spearmanrval_plottab[,"Eigengene"]
# spearmanrval_plottab <- spearmanrval_plottab[,!grepl("Eigengene", colnames(spearmanrval_plottab))]
# spearmanrval_plottab <- t(apply(spearmanrval_plottab, 2, function(x) as.numeric(as.character(x))))
# colnames(spearmanrval_plottab) <- keepeigennames
# 
# spearmanhmout <- WGCNA_custom_heatmap(moduleTraitPvalue = spearmanpval_plottab, moduleTraitCor = spearmanrval_plottab)
# pdf(paste0(outfilepathintegration, "wgcna_nmf_integration/", "spearman_values.pdf"), 10, 20, useDingbats = FALSE)
# print(spearmanhmout)
# junk <- dev.off()
# 
# 




