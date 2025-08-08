################################
# Name: MacIntosh Cornwell
# Email: mcornwell1957@gmail.com
################################
## ISCHEMIA Figures

## Load in Libraries
packagelist = c("Hmisc", "tools")
junk <- lapply(packagelist, function(xxx) suppressMessages(
  require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))

source("/Users/tosh/Desktop/Ruggles_Lab/code/mgc_plotting_functions.R")
source("/Users/tosh/Desktop/Ruggles_Lab/code/rnaseq_scripts/geneset_analysis_functions.R")
source("/Users/tosh/Desktop/Ruggles_Lab/code/WGCNA_functions.R")
source("/Users/tosh/Desktop/Ruggles_Lab/code/process_metadata_functions.R")
source("/Users/tosh/Desktop/Ruggles_Lab/code/mgc_survival_functions.R")
source("/Users/tosh/Desktop/Ruggles_Lab/code/summarize_table_function.R")
source("/Users/tosh/Desktop/Ruggles_Lab/code/overlap_finder_function.R")

## Outpath
# outfilepathintegration = paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run1c_rmoutliers2/integrative_analyses/")
# outfilepathfigures = paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/figures/raw/run3_20210729/")
outfilepathnmf = paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/figures/raw/run3_20210729/")
dir.create(outfilepathfigures, recursive = TRUE, showWarnings = FALSE)

## Infiles
inmetafile <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run2b_rmoutliers2_controlage/rna_processing/metatable_in.csv"
inmetatable <- read.table(inmetafile, sep = ",", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
# SOI <- rownames(na.omit(inmetatable[,"comp_ischemia__Sev_v_MildNone",drop=FALSE]))
SOI <- rownames(inmetatable)

## Grab the matched samples
incountfile <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run2b_rmoutliers2_controlage/rna_processing/normcounttab.txt"
normcounttable <- read.table(incountfile, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

## What about comparing it to the biomarker data...
inbiomarkertable_file <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/biomarker_data_cleaned_20210325.csv"
# inbiomarkertable_file <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/biomarker_data_cleaned_20210211.csv"
inbiomarkertable <- read.table(inbiomarkertable_file, sep = ",", header = TRUE)
biomarkerlabels <- colnames(inbiomarkertable)[grepl("_clean", colnames(inbiomarkertable))]

inbiomarkertable <- inbiomarkertable[rowSums(is.na(inbiomarkertable)) != (ncol(inbiomarkertable)-1), 
                                     grepl("_clean|PATNUM", colnames(inbiomarkertable))]

# Metabalomics table
inmetabalomics_file <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/ischemia_metabalomics_simplified_withID_cleaned_duplicates.csv"
# inmetabalomics_file <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/ischemia_metabalomics_simplified_withID_table.csv"
inmetabalomics_table <- read.table(inmetabalomics_file, sep = ",", header = TRUE, check.names = FALSE)
## Ok, 230 is kind of arbitrary, but that seems to be the magic cut off number, and then the rest we can assess as a case by case basis
inmetabalomics_table <- inmetabalomics_table[rowSums(is.na(inmetabalomics_table))<230, ]
## Change all remaining "TAG" values to NA
inmetabalomics_table[inmetabalomics_table == "TAG"] <- NA

# Methylation table
# inmethylation_file <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/ischemia_methylation_dummydata.csv"
# inmethylationtable <- read.table(inmethylation_file, sep = ",", header = TRUE, row.names = 1, check.names = FALSE)

#inmethyl_file <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/ischemia_methylation_beta_top100k.RDS"
# load(file = inmethyl_file) # Loads methyltable_convertedIDs
# inmethyl_file <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/ischemia_methylation_DC_beta_top100k_nosnp_noxy.RDS"
# load(inmethyl_file)
# ## We are going to subsample this down to top 8000, that seems to be a magic number... but I guess we could do more?
# methylvariances <- apply(topmethyltable, 1, var)
# methylvariances <- methylvariances[order(methylvariances, decreasing = TRUE)]
# numbermethylsites <- 8000
# topmethyltable <- methyltable_convertedIDs[names(methylvariances)[1:numbermethylsites],]

inmethylfile <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/top100k_vars_probes_rm_dupe_XY.rds"
# inmethylfile <- "/gpfs/data/ischemialab/projects/nmf_multiomics_clustering/data/top100k_vars_probes_rm_dupe_XY.rds"
topmethyltable <- readRDS(file = inmethylfile)
#load(file = inmethyl_file) # Loads methyltable_convertedIDs
## We are going to subsample this down to top 8000, that seems to be a magic number... but I guess we could do more?
#methylvariances <- apply(methyltable_convertedIDs, 1, var)
methylvariances <- apply(topmethyltable, 1, var)
methylvariances <- methylvariances[order(methylvariances, decreasing = TRUE)]
numbermethylsites <- 8000
topmethyltable <- topmethyltable[names(methylvariances)[1:numbermethylsites],]




# Full biorep table
inbioreptablefile <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/biorep_10_27.csv"
inbioreptable <- read.table(inbioreptablefile, sep = ",", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("NA", "", NA))

## Ok, to start off with, lets get an idea of what I have to work with from each data set - and what the overlap of the datatypes are...
overlap_inlist = list(
                      RNAseq_PATNUMS = colnames(normcounttable),
                      Methyl_PATNUMS = colnames(topmethyltable),
                      Biomarker_PATNUMS = inbiomarkertable[,"PATNUM"],
                      Metabalomics_PATNUMS = inmetabalomics_table[,"PATNUM"]
                     )
overlap_analysisout <- overlap_finder(overlap_inlist)
overlaptable <- overlap_analysisout$overlaptable
vennplot <- overlap_analysisout$vennplot
overlapgrouptab <- overlap_analysisout$overlapgrouptab
overlapsummary <- overlap_analysisout$overlapsummary

dir.create(paste0(outfilepathnmf, "sample_datatype_overlap_analysis/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(outfilepathnmf, "sample_datatype_overlap_analysis/venndiagram.pdf"))
grid.draw(vennplot)
junk <- dev.off()
write.table(overlaptable, paste0(outfilepathnmf, "sample_datatype_overlap_analysis/overlaptable.csv"), sep = ",", row.names = FALSE, col.names = TRUE)
write.table(overlapgrouptab, paste0(outfilepathnmf, "sample_datatype_overlap_analysis/overlapgrouptab.csv"), sep = ",", row.names = FALSE, col.names = TRUE)
write.table(overlapsummary, paste0(outfilepathnmf, "sample_datatype_overlap_analysis/overlapsummary.csv"), sep = ",", row.names = TRUE, col.names = TRUE)



## I think an upset plot here would be much more attractive
upset_indata <- overlap_inlist
upsetout <- create_upset_plot(upset_indata, comboparam = "all", transposeparam = TRUE)
upsetobject <- upsetout$upsetobject
upsetplot <- upsetout$upsetplot

pdf(paste0(outfilepathnmf, "sample_datatype_overlap_analysis/upset_plot.pdf"))
draw(upsetplot)
junk <- dev.off()




## Thats just the sample overlap - lets get a sense for how many values are missing for each feature
dir.create(paste0(outfilepathnmf, "missing_sample_feature_analysis/"), showWarnings = FALSE, recursive = TRUE)

rough_combined_table1 <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "PATNUM", all = TRUE, sort = FALSE),
                                list(
                                  t(rbind.data.frame(PATNUM = colnames(normcounttable), normcounttable)),
                                  t(rbind.data.frame(PATNUM = colnames(topmethyltable), topmethyltable))
                                  # inbiomarkertable,
                                  # inmetabalomics_table[,!colnames(inmetabalomics_table) %in% "Inventory.Code"]
                                ))
rough_combined_table2 <- data.frame(t(rough_combined_table1))
colnames(rough_combined_table2) <- rough_combined_table2["PATNUM",]
rough_combined_table2 <- rough_combined_table2[!rownames(rough_combined_table2) %in% "PATNUM",]

rough_combined_table2[,"datatype"] <- c(rep("RNAseq", nrow(normcounttable)),
                                        rep("Methylation", nrow(topmethyltable))
                                        # rep("biomarker", ncol(inbiomarkertable[,!colnames(inbiomarkertable) %in% "PATNUM"])),
                                        # rep("metabolomics", ncol(inmetabalomics_table[,!colnames(inmetabalomics_table) %in% c("PATNUM", "Inventory.Code")]))
                                        )

## Ok, all samples with RNA are going to be filled RNA by definition - so lets not make more work
# full_combinedtable_RNAsel <- rough_combined_table2[,colSums(!is.na(rough_combined_table2[rough_combined_table2[,"datatype"] %in% "RNAseq",])) == nrow(normcounttable)]
# # Now with these, lets select for just the nonRNA rows
# noRNA_combinedtable_RNAsel <- full_combinedtable_RNAsel[!full_combinedtable_RNAsel[,"datatype"] %in% "RNAseq",]
# tt1 <- rowSums(is.na(noRNA_combinedtable_RNAsel))


# Lets just grab the samples we need from the previous overlap I did
# combinedtable_overlapsamples <- rough_combined_table2[,c(colnames(rough_combined_table2)[colnames(rough_combined_table2) %in% overlaptable[,ncol(overlaptable)]], "datatype")]
# write.table(combinedtable_overlapsamples, 
#             paste0(outfilepathnmf, "missing_sample_feature_analysis/", "combinedtable_nofiltering.csv"),
#             sep = ",", col.names = NA, row.names = TRUE)
# 
# # Ok Lets get a sense for (1) how many features each sample are missing
# samples_missing_features <- data.frame(NumberFeaturesMissing = colSums(is.na(combinedtable_overlapsamples)))
# histout1 <- plot_histogram(data = samples_missing_features, binparam = max(samples_missing_features),
#                labsparam = list(title = "Number of missing features for each Sample", 
#                                 x = "Number of Missing Features", y = "Number of samples with given quantity of missing features"))
# pdf(paste0(outfilepathnmf, "missing_sample_feature_analysis/", "samples_with_missing_features_hist.pdf"))
# print(histout1)
# junk <- dev.off()
# write.table(samples_missing_features[order(samples_missing_features[,1]),,drop=FALSE], 
#             paste0(outfilepathnmf, "missing_sample_feature_analysis/", "samples_with_missing_features_table.csv"),
#             sep = ",", col.names = NA, row.names = TRUE)
# 
# # and (2) how many samples each feature are missing
# features_missing_samples <- cbind.data.frame(datatype = combinedtable_overlapsamples[,"datatype"], missingvalues = rowSums(is.na(combinedtable_overlapsamples)))
# # RNA isnt helpful, so lets remove those for now
# features_missing_samples_noRNA <- features_missing_samples[!features_missing_samples[,"datatype"] %in% "RNAseq",]
# histout2_biomarker <- plot_histogram(data = features_missing_samples_noRNA[features_missing_samples_noRNA[,"datatype"] %in% "biomarker",2,drop=FALSE],
#                            binparam = max(samples_missing_features),
#                            labsparam = list(title = "Number of missing Samples for each Feature", 
#                                             x = "Number of Missing Samples", y = "Number of Features with given quantity of missing Samples"))
# pdf(paste0(outfilepathnmf, "missing_sample_feature_analysis/", "BIOMARKER_features_with_missing_samples_hist.pdf"))
# print(histout2_biomarker)
# junk <- dev.off()
# 
# histout2_metab <- plot_histogram(data = features_missing_samples_noRNA[features_missing_samples_noRNA[,"datatype"] %in% "metabolomics",2,drop=FALSE],
#                             binparam = max(samples_missing_features),
#                             labsparam = list(title = "Number of missing Samples for each Feature", 
#                                              x = "Number of Missing Samples", y = "Number of Features with given quantity of missing Samples"))
# pdf(paste0(outfilepathnmf, "missing_sample_feature_analysis/", "METAB_features_with_missing_samples_hist.pdf"))
# print(histout2_metab)
# junk <- dev.off()
# 
# write.table(features_missing_samples_noRNA[order(features_missing_samples_noRNA[,2]),,drop=FALSE], 
#             paste0(outfilepathnmf, "missing_sample_feature_analysis/", "features_with_missing_samples_table.csv"),
#             sep = ",", col.names = NA, row.names = TRUE)
# 
# ## Lets get a cross-sectional overview with a HM of missing features (excluding RNA because again, pointless)
# hmplottab <- combinedtable_overlapsamples[!combinedtable_overlapsamples[,"datatype"] %in% c("RNAseq", "Methylation"), 
#                                           !colnames(combinedtable_overlapsamples) %in% "datatype"]
# hmplottab[] <- ifelse(is.na(hmplottab), 1, 0)
# # hmplottab <- hmplottab[rowSums(hmplottab) !=0 , colSums(hmplottab) !=0 ]
# 
# heatmapcolorparam <- colorRamp2(breaks = c(0,1), c("white", "darkgreen"))
# missing_data_heatmap <- function(hmplottab, rowsplitparam = rowsplitparam) {
#   summary_outhm <- Heatmap(as.matrix(hmplottab),
#                            col = heatmapcolorparam,    ## Define the color scale for the heatmap
#                            row_title = paste0(dim(hmplottab), collapse=" x "),                                       ## Name the rows
#                            column_title = "Samples",                                  ## Name the columns
#                            
#                            row_split = rowsplitparam[,1],
#                            
#                            cluster_columns = TRUE,                         ## Cluster the columns or leave as is
#                            cluster_rows = TRUE,                            ## Cluster the rows or leave as is
#                            
#                            border = TRUE,
#                            na_col = "white",
#                            # rect_gp = gpar(col = "black", lwd = 0.5),
#                            
#                            show_column_names = TRUE,                                  ## Show the Column Names
#                            column_names_gp = gpar(fontsize = 4),                      ## Change the size of the column names
#                            show_row_names = nrow(hmplottab) < 200,                                    ## Show the row names
#                            row_names_side = "left",                                   ## Place the row names on the Left side of the heatmap
#                            row_names_gp = gpar(fontsize=6),
#                            
#                            heatmap_legend_param = list(
#                              title = "Missing Data",
#                              title_gp = gpar(fontsize = 8, fontface = "bold")),
#                            height = unit(min((nrow(hmplottab)/2), 12),"cm"),
#                            width = unit(min(ncol(hmplottab), 18),"cm")
#                            
#   )
#   return(summary_outhm)
# }
# 
# rowsplitparam <- combinedtable_overlapsamples[!combinedtable_overlapsamples[,"datatype"] %in% c("RNAseq", "Methylation"), 
#                                               colnames(combinedtable_overlapsamples) %in% "datatype",drop=FALSE]
# summary_outhm <- missing_data_heatmap(hmplottab, rowsplitparam = rowsplitparam)
# pdf(paste0(outfilepathnmf, "missing_sample_feature_analysis/", "missing_sample_feature_summary_ALLSAMPLES_heatmap.pdf"), 10, 10, useDingbats = FALSE)
# draw(summary_outhm)
# junk <- dev.off()
# 
# 
# rowsplitparam <- rowsplitparam[rownames(hmplottab[rowSums(hmplottab) !=0 , colSums(hmplottab) !=0 ]),,drop=FALSE]
# summary_outhm <- missing_data_heatmap(hmplottab[rowSums(hmplottab) !=0 , colSums(hmplottab) !=0 ], rowsplitparam = rowsplitparam)
# pdf(paste0(outfilepathnmf, "missing_sample_feature_analysis/", "missing_sample_feature_summary_ONLYMISSING_heatmap.pdf"), 10, 10)
# draw(summary_outhm)
# junk <- dev.off()
# 
# 
# 
# ## Now, let's (manually) remove some features and samples to maximize out population
# features_to_remove <- c("IL6_clean",
#                         "bOHbutyrate", "Gly", "Pyruvate", "Acetate", "Gln", "Creatinine",
#                         "XXL-VLDL-PL %", "XXL-VLDL-C %", "XXL-VLDL-CE %", "XXL-VLDL-FC %", "XXL-VLDL-TG %")
# samples_to_remove <- c("011006-516", "064001-502")
# 
# combinedtable_overlapsamples_manualfilter <- combinedtable_overlapsamples[!rownames(combinedtable_overlapsamples) %in% features_to_remove,
#                                                                           !colnames(combinedtable_overlapsamples) %in% samples_to_remove]
# 
# hmplottab <- combinedtable_overlapsamples_manualfilter[!combinedtable_overlapsamples_manualfilter[,"datatype"] %in% "RNAseq", 
#                                           !colnames(combinedtable_overlapsamples_manualfilter) %in% "datatype"]
# hmplottab[] <- ifelse(is.na(hmplottab), 1, 0)
# # hmplottab <- hmplottab[rowSums(hmplottab) !=0 , colSums(hmplottab) !=0 ]
# 
# rowsplitparam <- rowsplitparam[rownames(hmplottab[rowSums(hmplottab) !=0 , colSums(hmplottab) !=0 ]),,drop=FALSE]
# summary_outhm <- missing_data_heatmap(hmplottab[rowSums(hmplottab) !=0 , colSums(hmplottab) !=0 ], rowsplitparam = rowsplitparam)
# pdf(paste0(outfilepathnmf, "missing_sample_feature_analysis/", "missing_sample_feature_summary_onlymissing_manualremoval_heatmap.pdf"), 10, 10)
# draw(summary_outhm)
# junk <- dev.off()
# 
# samples_missing_features_manualfilter <- data.frame(NumberFeaturesMissing = colSums(is.na(combinedtable_overlapsamples_manualfilter)))
# features_missing_samples_manualfilter <- cbind.data.frame(datatype = combinedtable_overlapsamples_manualfilter[,"datatype"],
#                                              missingvalues = rowSums(is.na(combinedtable_overlapsamples_manualfilter)))
# 
# 
# write.table(features_missing_samples_manualfilter[order(features_missing_samples_manualfilter[,2]),,drop=FALSE], 
#             paste0(outfilepathnmf, "missing_sample_feature_analysis/", "features_with_missing_samples_manualremoval_table.csv"),
#             sep = ",", col.names = NA, row.names = TRUE)
# 
# write.table(samples_missing_features_manualfilter[order(samples_missing_features_manualfilter[,1]),,drop=FALSE], 
#             paste0(outfilepathnmf, "missing_sample_feature_analysis/", "samples_with_missing_features_manualremoval_table.csv"),
#             sep = ",", col.names = NA, row.names = TRUE)
# 
# 
# ## Imputing Values (trying median and KNN)
# dir.create(paste0(outfilepathnmf, "missing_sample_feature_analysis/imputation_analysis/"), showWarnings = FALSE, recursive = TRUE)
# imputation_intable <- combinedtable_overlapsamples_manualfilter[,!colnames(combinedtable_overlapsamples_manualfilter) %in% "datatype"]
# imputation_intable <- t(apply(imputation_intable, 2, as.numeric))
# colnames(imputation_intable) <- rownames(combinedtable_overlapsamples_manualfilter)
# 
# knn_imputation_out <- mgc_imputation(imputation_intable = imputation_intable, features_to_impute = NULL,
#                                      imputation_method = "knn", 
#                                      outplotpath = paste0(outfilepathnmf, "missing_sample_feature_analysis/imputation_analysis/knn_imputed/"))
# knn_imputed_data <- knn_imputation_out$imputed_outtable
# write.table(knn_imputed_data, paste0(outfilepathnmf, "missing_sample_feature_analysis/imputation_analysis/knn_imputed/", "knn.csv"), sep = ",", col.names = NA, row.names = TRUE)
# 
# median_imputation_out <- mgc_imputation(imputation_intable = imputation_intable, features_to_impute = NULL,
#                                      imputation_method = "median",
#                                      outplotpath = paste0(outfilepathnmf, "missing_sample_feature_analysis/imputation_analysis/median_imputed/"))
# median_imputed_data <- median_imputation_out$imputed_outtable
# write.table(median_imputed_data, paste0(outfilepathnmf, "missing_sample_feature_analysis/imputation_analysis/median_imputed/", "median_data.csv"), sep = ",", col.names = NA, row.names = TRUE)

## Now that we have this, lets do a PCA of just our biomarker/metabolomic data BEFORE imputation and AFTER imputation, to make sure there are no drastic global shifts
# pcaplottab_noimputation <- imputation_intable[,rownames(features_missing_samples_manualfilter[features_missing_samples_manualfilter[,2] > 0,])]
# pcaplotout_noimputation <- pca_plotter(pcadata = pcaplottab_noimputation, colorvar = NULL,
#                           scalecolor=FALSE, separatelegend = FALSE,
#                           labelpoints = FALSE, labsparam = list(title = "PCA before Imputation", x = "PC1", y = "PC2"))



## Need LOTS of formatting for Multiomic NMF (to even the playing field)
# To enable integrative multi-omics clustering, we required all data types (and converted if necessary) to represent ratios to either a common reference measured in each TMT plex (proteome, phosphoproteome, acetylproteome) or 
# DONE - an in-silico common reference calculated as the median abundance across all samples (mRNA, see “RNA quantification”). 
# DONE - All data tables were then concatenated and only features quantified in all tumors were used for subsequent analysis. 
# DONE - Features with the lowest standard deviation (bottom 5th percentile) across all samples were deemed uninformative and were removed from the dataset. 
# DONE - Each row in the data matrix was further scaled and standardized such that all features from different data types were represented as z-scores.

## 1 - create combined table
# We need all tables to be features on x, and samples on y, but the easiest is to merge on the x-axis, and then flip and adjust
# rough_combined_table1 <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "PATNUM", all = TRUE, sort = FALSE),
#                                 list(t(rbind.data.frame(PATNUM = colnames(normcounttable), normcounttable)),
#                                      nmfcluster_biomarker_intab,
#                                      inmetabalomics_table))
# rough_combined_table2 <- t(rough_combined_table1)
# colnames(rough_combined_table2) <- rough_combined_table2["PATNUM",]
# rough_combined_table2 <- rough_combined_table2[!rownames(rough_combined_table2) %in% "PATNUM",]
# 
# rough_combined_table3 <- rough_combined_table2[,colSums(is.na(rough_combined_table2)) == 0]

# rough_combined_table3 <- t(knn_imputed_data)
meth_and_rna_samples <- overlapgrouptab[grepl("RNAseq_PATNUMS", overlapgrouptab[,2]) & grepl("Methyl_PATNUMS", overlapgrouptab[,2]),1]
rough_combined_table3 <- rough_combined_table2[,meth_and_rna_samples]

## Ok, so with this we only have (306)... samples that have ALL info... That is substantially less than before............ but let run with it for now, and no NAs!


## 2 - Turn everything into a ratio of medians for each feature AND scale
# using forloop for now because theres a bug somewhere...
normalizationoutlist <- sdoutlist <- list()
for(rownum in seq_len(nrow(rough_combined_table3))) {
  rowsel <- as.numeric(as.character(rough_combined_table3[rownum,]))
  rowlabel <- rownames(rough_combined_table3)[rownum]
  # First median center all data
  outrow <- rowsel/median(rowsel)
  
  ## Here - we grab the SD to later filter out our least informative features
  sdoutlist[[rownum]] <- sd(outrow)
  names(sdoutlist)[rownum] <- rowlabel
  
  # Then scale to zscore
  outrow2 <- scale(outrow)
  
  if (sum(is.na(outrow2)) > 0){break}
  
  normalizationoutlist[[rownum]] <- outrow2
  names(normalizationoutlist)[rownum] <- rowlabel
}
rough_combined_table4 <- t(do.call(cbind.data.frame, normalizationoutlist))
colnames(rough_combined_table4) <- colnames(rough_combined_table3)

# 3 - SD reference table
percentilecutoffparam <- 0.05
percentilecutoffvalue <- quantile(unlist(sdoutlist), percentilecutoffparam)
percentile_cutoff_features <- names(unlist(sdoutlist)[unlist(sdoutlist) < percentilecutoffvalue])
# remove bottom 5 percetile of SD features to clean up our data
rough_combined_table5 <- rough_combined_table4[!rownames(rough_combined_table4) %in% percentile_cutoff_features,]
write.table(rough_combined_table5, paste0(outfilepathnmf, "nmfintable_notsplit.csv"), sep = ",", col.names = NA, row.names = TRUE)


# DONE - Since NMF requires a non-negative input matrix, the data matrix of z-scores was further converted into a non-negative matrix as follows:
# 1) Create one data matrix with all negative numbers zeroed.
# 2) Create another data matrix with all positive numbers zeroed and the signs of all negative numbers removed.
# 3) Concatenate both matrices resulting in a data matrix twice as large as the original, but with positive values only and zeros and hence appropriate for NMF.
posvaluetable <- rough_combined_table5
posvaluetable[posvaluetable < 0] <- 0

negvaluetable <- rough_combined_table5
negvaluetable[negvaluetable > 0] <- 0
negvaluetable <- -negvaluetable

rough_combined_table6 <- rbind(posvaluetable, negvaluetable)


# Set this as our NMF table
nmfintable <- rough_combined_table6
write.table(nmfintable, paste0(outfilepathnmf, "nmfintable_posneg.csv"), sep = ",", col.names = NA, row.names = TRUE)

## Step 1 - clean data as needed
## NOTHING YET
# DONE - all done above

## Step 2 test rank
dir.create(paste0(outfilepathnmf, "ranktestanalysis/"), showWarnings = FALSE, recursive = TRUE)

nmf_ranktest_out <- nmf_ranktest(incounttable = nmfintable, ranksparam = 2:8, nrunparam = 50,
                                 seedparam = 123456, plotconsensusmapparam = paste0(outfilepathnmf, "ranktestanalysis/"))
ranktest_object <- nmf_ranktest_out$ranktestobject
rankmetricplot <- nmf_ranktest_out$rankmetricplot
ranktest_measures <- nmf_ranktest_out$rankmetrictable
suggested_rank <- nmf_ranktest_out$suggestedrank

save(object = ranktest_object, file = paste0(outfilepathnmf, "ranktestanalysis/", "ranktest_object.rds"))
pdf(paste0(outfilepathnmf, "ranktestanalysis/", "rankmetricplot.pdf"), useDingbats = FALSE)
print(rankmetricplot)
junk <- dev.off()
write.table(ranktest_measures, paste0(outfilepathnmf, "ranktestanalysis/", "ranktest_measures.csv"), sep = ",", col.names = TRUE, row.names = FALSE)

## Step 3 run nmf
dir.create(paste0(outfilepathnmf, "rank_4_nmf/"), showWarnings = FALSE, recursive = TRUE)
nmfout <- run_nmf(incounttable = nmfintable, ranksparam = 4, methodparam = NULL, nrun = 200)
nmfout_object <- nmfout$nmfobject
nmf_w <- nmfout$nmf_w_mat
nmf_h <- nmfout$nmf_h_mat
nmf_h_membership_scores <- nmfout$nmf_h_membership_scores
feature_scores <- nmfout$nmf_featurescores
driving_features <- nmfout$nmf_drivingfeatures
driving_features_named <- lapply(driving_features, function(x) names(feature_scores[x]))
driving_features_table <- do.call(cbind.fill, driving_features_named)
colnames(driving_features_table) <- paste0("nmf_cluster_", seq(1,ncol(driving_features_table)))

## Define the cluster each sample belongs to
clustermembershipmat <- t(apply(nmf_h_membership_scores, 1, function(x) x == max(x)))
colnames(clustermembershipmat) <- paste0("nmf_cluster_", seq(1,ncol(clustermembershipmat)))
# HAX: https://stackoverflow.com/questions/26703363/replace-value-with-the-name-of-its-respective-column
w <- which(clustermembershipmat == TRUE,arr.ind=TRUE)
clustermembershipmat[w] <- colnames(clustermembershipmat)[w[,"col"]]
clustermembershipmat[clustermembershipmat == "FALSE"] <- NA

## Define our clustercores
clustercores <- clustermembershipmat
clustercores[nmf_h_membership_scores <= 0.5] <- NA

saveRDS(object = nmfout_object, file = paste0(outfilepathnmf, "rank_4_nmf/", "nmfout_object.rds"))
write.table(nmf_w, paste0(outfilepathnmf, "rank_4_nmf/", "nmf_W.csv"), sep = ",", col.names = FALSE, row.names = TRUE)
write.table(nmf_h, paste0(outfilepathnmf, "rank_4_nmf/", "nmf_H.csv"), sep = ",", col.names = TRUE, row.names = FALSE)
write.table(nmf_h_membership_scores, paste0(outfilepathnmf, "rank_4_nmf/", "nmf_H_membership_scores.csv"), sep = ",", col.names = FALSE, row.names = TRUE)
write.table(data.frame(feature_scores), paste0(outfilepathnmf, "rank_4_nmf/", "nmf_feature_scores.csv"), sep = ",", col.names = FALSE, row.names = TRUE)
write.table(driving_features_table, paste0(outfilepathnmf, "rank_4_nmf/", "nmf_driving_features_table.csv"), sep = ",", col.names = TRUE, row.names = FALSE)
write.table(clustermembershipmat, paste0(outfilepathnmf, "rank_4_nmf/", "nmf_clustermembership_table.csv"), sep = ",", col.names = NA, row.names = TRUE)
write.table(clustercores, paste0(outfilepathnmf, "rank_4_nmf/", "nmf_clustercores_table.csv"), sep = ",", col.names = NA, row.names = TRUE)






## Write all of this out ## TOSH - THIS WILL ALL HAVE TO BE ADJUSTED WHEN YOU RERUN THE FUNCTIONALIZED FORM OF THE CODE