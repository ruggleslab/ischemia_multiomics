################################
# Name: MacIntosh Cornwell
# Email: mcornwell1957@gmail.com
################################
## ISCHEMIA Subtype Analysis

packagelist = c("psych")
junk <- lapply(packagelist, function(xxx) suppressMessages(
    require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))
source("/Users/tosh/Desktop/Ruggles_Lab/code/overlap_finder_function.R")
source("/Users/tosh/Desktop/Ruggles_Lab/code/summarize_table_function.R")
source("/Users/tosh/Desktop/Ruggles_Lab/code/rnaseq_scripts/deseq_functions.R")
source("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/code/ischemia2021_exploratory_analysis_functions.R")


# Outfolder
outfilepathmaster <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/integrative_analyses/run9-figs/"
dir.create(outfilepathmaster, recursive = TRUE, showWarnings = FALSE)

# save.image(file = paste0(outfilepathmaster, "/ischemia_exploratory.RData"))
# load(file = paste0(outfilepathmaster, "/ischemia_exploratory.RData"))


# --------------------------------- Utility tables for colors and statistical cutoffs ---------------------------------
# Colorguide
colorguidefile <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/ischemia2021_colorguide.csv"
colorguide <- read.table(colorguidefile, sep = ",", header = TRUE, comment.char = "", colClasses = c("character", "character", "character"))

# Differential analysis cutoffs
## Ugly - but allows for easily adjusting cutoffs if we want
diff_feat_analysis_cutoffs <- suppressWarnings(read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/integrative_analyses/manual_annotation_tables/diff_feat_analysis_cutoffs.csv",sep = ",", header = TRUE, row.names = 1))
comp_ischemia__Sev_v_MildNone__PVALTEST <- diff_feat_analysis_cutoffs["comp_ischemia__Sev_v_MildNone", "pval_test"]
comp_ischemia__Sev_v_MildNone__PVALCUTOFF <- diff_feat_analysis_cutoffs["comp_ischemia__Sev_v_MildNone", "pval_cutoff"]
comp_ischemia__Sev_v_MildNone__FCCUTOFF <- diff_feat_analysis_cutoffs["comp_ischemia__Sev_v_MildNone", "fc_cutoff"]
comp_anatomy70__3v_v_1v__PVALTEST <- diff_feat_analysis_cutoffs["comp_anatomy70__3v_v_1v", "pval_test"]
comp_anatomy70__3v_v_1v__PVALCUTOFF <- diff_feat_analysis_cutoffs["comp_anatomy70__3v_v_1v", "pval_cutoff"]
comp_anatomy70__3v_v_1v__FCCUTOFF <- diff_feat_analysis_cutoffs["comp_anatomy70__3v_v_1v", "fc_cutoff"]
DMP_ischemia__Sev_v_MildNone__PVALTEST <- diff_feat_analysis_cutoffs["DMP_ischemia__Sev_v_MildNone", "pval_test"]
DMP_ischemia__Sev_v_MildNone__PVALCUTOFF <- diff_feat_analysis_cutoffs["DMP_ischemia__Sev_v_MildNone", "pval_cutoff"]
DMP_ischemia__Sev_v_MildNone__FCCUTOFF <- diff_feat_analysis_cutoffs["DMP_ischemia__Sev_v_MildNone", "fc_cutoff"]
DMP_anatomy70__3v_v_1v__PVALTEST <- diff_feat_analysis_cutoffs["DMP_anatomy70__3v_v_1v", "pval_test"]
DMP_anatomy70__3v_v_1v__PVALCUTOFF <- diff_feat_analysis_cutoffs["DMP_anatomy70__3v_v_1v", "pval_cutoff"]
DMP_anatomy70__3v_v_1v__FCCUTOFF <- diff_feat_analysis_cutoffs["DMP_anatomy70__3v_v_1v", "fc_cutoff"]

# --------------------------------- Read in our metatable and bioreptable ---------------------------------
metatable_file <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/rna_processing/metatable_in.csv"
metatable <- read.table(metatable_file, sep = ",", header = TRUE, row.names = 1)
biorepfile <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/isch_BiorepData_3_23_2022_CUSTOM.csv"
bioreptable <- read.table(biorepfile, header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)
rownames(bioreptable) <- bioreptable[,"PATNUM"]

rna_clustermembership_table_file <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/MULTIOMICS_SUBTYPE_LABELS/nmf_cluster_labels_with_subcluster_1_20220208.csv"
rna_clustermembership_table <- read.table(rna_clustermembership_table_file, sep = ",", header = TRUE, row.names = 1)
rna_clustermembership_table <- rna_clustermembership_table[!is.na(rna_clustermembership_table[,"rna_4cluster_w3AB"]),]
rna_clustermembership_table <- rna_clustermembership_table[order(rna_clustermembership_table[,"rna_4cluster_w3AB"]),]

meth_clustermembership_table_file <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/MULTIOMICS_SUBTYPE_LABELS/meth_cluster_membership_table_4_nodup_20220606.csv"
meth_clustermembership_table <- read.table(meth_clustermembership_table_file, sep = ",", header = TRUE, row.names = 1)
meth_clustermembership_table <- meth_clustermembership_table[order(meth_clustermembership_table[,"meth_4cluster"]),]

combined_clustermembership_table <- merge(rna_clustermembership_table[,c("rna_4cluster_w3AB", "rna_4cluster"), drop = FALSE], 
                                          meth_clustermembership_table[,c("meth_4cluster", "meth_3cluster"),drop=FALSE], 
                                          by.x = "row.names", by.y = "row.names", all = TRUE)
colnames(combined_clustermembership_table)[1] <- c("PATNUM")
rownames(combined_clustermembership_table) <- combined_clustermembership_table[,"PATNUM"]

## Making conenience variables for all our samples, all rna, and all meth
SOIall <- rownames(bioreptable)
SOIrna <- unique(rownames(rna_clustermembership_table))
SOImeth <- unique(rownames(meth_clustermembership_table))
SOIboth <- intersect(rownames(rna_clustermembership_table), rownames(meth_clustermembership_table))
SOIeither <- unique(c(rownames(rna_clustermembership_table), rownames(meth_clustermembership_table)))

# Normalized counts
normcounttable_file <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/rna_processing/normcounttab.txt"
normcounttable <- read.table(normcounttable_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

# biomarkertable attached to bioreptable:
biomarkerfile <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/docs/ischemia_biomarkers/Marker_Meta_Joint_Dropdup_release.V1.2.csv"
biomarkertable <- read.table(biomarkerfile, sep = ",", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE,
                             na.strings = c("", NA, "NA"))
biomarkertable_sel <- biomarkertable[,c("PATNUM", colnames(biomarkertable)[grepl("_clean$", colnames(biomarkertable))])]

# WGCNA
eigengenecounttable_file <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run3_rmoutliers2/WGCNA/WGCNA_power14_size30/wgcna_eigengenes.csv"
eigengenecounttable <- read.table(eigengenecounttable_file, sep = ",", header = TRUE, row.names = 1)
genestocolorstab_file <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run3_rmoutliers2/WGCNA/WGCNA_power14_size30/wgcna_genestocolors.csv"
genestocolorstab <- read.table(genestocolorstab_file, ",", header = TRUE, row.names = 1)
eigengenes <- paste0("ME",
                     c("magenta", "salmon", "purple", "midnightblue", "turquoise", 
                       "yellow", "pink", "red", 
                       "blue", "tan", "brown",
                       "greenyellow", "black", "cyan", "green"))



# --------------------------------- meth WGCNA read in and analysis ---------------------------------
dir.create(paste0(outfilepathmaster, "meth_WGCNA/"), showWarnings = FALSE, recursive = TRUE)

# Read in meth count table and metatable
methmetatable_file <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/final_meth_data_v1_20230207/kpp_metadata.rds"
methmetatable <- readRDS(methmetatable_file)
methcounttable_file <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/final_meth_data_v1_20230207/limma_noCovariate.rds"
methcounttable <- readRDS(methcounttable_file)
colnames(methcounttable) <- methmetatable[match(methmetatable[,"Sample_Name"], colnames(methcounttable)), "Sample_Patnum"]
methcounttable <- methcounttable[,-578] # Need to remove the duplicate of this sample: "048003-005", columns 204 and 578

# Read in format and output the meth annotation table
meth_wgcna_annotation_goterm_table_list <- readRDS("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/methylation_data_20221020/WGCNA_pathway_annotation_hypergeo.rds")
meth_wgcna_genestocolors <- read.table(paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/integrative_analyses/run7/", "meth_WGCNA/", "meth_wgcna_genestocolors.csv"), sep = ",", header = TRUE, row.names = 1)

# Add the adjusted p value and add the gene ratio value
meth_wgcna_annotation_goterm_table_list <- lapply(meth_wgcna_annotation_goterm_table_list, function(x){
    x[,"padj"] <- x[,"padj"]
    x[,"GeneRatio"] <- apply(x, 1, function(x2) gsub(" ", "", paste0(x2["numDEInCat"], "/",  x2["numInCat"])))
    x
})
# reformat to match other (just to be consistent)
# module, size, pval, p.adjust, GeneRatio, ont, ID
genes_to_color_reftable <- merge(unique(meth_wgcna_genestocolors), data.frame(table(meth_wgcna_genestocolors[,2])),
                                 by.x = "moduleColors", by.y = "Var1")
colnames(genes_to_color_reftable) <- c("moduleColors", "moduleLabels", "size")
# 1+2 - add color back and module size
meth_wgcna_annotation_table <- merge(bind_rows(meth_wgcna_annotation_goterm_table_list, .id = "groups"),
                                     genes_to_color_reftable, by.x = "groups", by.y = "moduleLabels", all.x = TRUE)
# 5 - reorder columns
meth_wgcna_annotation_table <- meth_wgcna_annotation_table[,c("moduleColors", "size", "over_represented_pvalue", "padj", "GeneRatio", "category", "category")]
colnames(meth_wgcna_annotation_table) <- c("module", "size", "pvalue", "p.adjust", "GeneRatio", "ont", "ID")
rownames(meth_wgcna_annotation_table) <- apply(meth_wgcna_annotation_table[,c("module", "ID")], 1, function(x) paste(x, collapse = "."))

write.table(meth_wgcna_annotation_table, paste0(outfilepathmaster, "meth_WGCNA/", "GO_meth_wgcna_module_annotations.csv"),
            sep = ",", col.names = NA, row.names = TRUE)

write.table(meth_wgcna_annotation_table[meth_wgcna_annotation_table[,"p.adjust"] < 0.05,], 
            paste0(outfilepathmaster, "meth_WGCNA/", "GO_meth_wgcna_module_annotations_padjval05filt.csv"),
            sep = ",", col.names = NA, row.names = TRUE)



# 119 entries long - DMP tables
# they are fisher tests for another gene annotation, but idk what they are... so im gonna ignore for now...
# meth_wgcna_annotation_fisher_table_file <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/methylation_data_20220608/WGCNA/WGCNA_pathway_annotation_fisher.rds"
# meth_wgcna_annotation_fisher_table <- readRDS(meth_wgcna_annotation_fisher_table_file)


# # Ok well only the top (15?) have any annotation, so thats actuall easier to deal with...
# # ok this looks more like an actual GO term annoation... for each module
# meth_wgcna_annotation_goterm_table_file <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/methylation_data_20220608/WGCNA/WGCNA_pathway_annotation_Goterm.rds"
# meth_wgcna_annotation_goterm_table_list <- readRDS(meth_wgcna_annotation_goterm_table_file)# also 119 long
# 
# 
# 
# # List where each entry is a module and the values in it are the probes in that module
# meth_wgcna_ht_molecular_order_no0_file <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/methylation_data_20220608/WGCNA/WGCNA_ht_molecular_order_without0.rds"
# meth_wgcna_ht_molecular_order_no0 <- readRDS(meth_wgcna_ht_molecular_order_no0_file)
# 
# 
# 
# 
# 
# probe_to_module_list <- meth_wgcna_ht_molecular_order_no0
# temp_meth_genestocolors <- merge(data.frame(methprobe = rownames(methcounttable)),
#                                  cbind.data.frame(methprobe = unlist(probe_to_module_list), 
#                                                   moduleLabels = rep(seq_len(length(probe_to_module_list)), unlist(lapply(probe_to_module_list, length)))),
#                                  by = "methprobe", all = TRUE)
# temp_meth_genestocolors[is.na(temp_meth_genestocolors[,"moduleLabels"]),"moduleLabels"] <- 0
# meth_wgcna_genestocolors <- data.frame(moduleLabels = as.numeric(temp_meth_genestocolors[,"moduleLabels"]),
#                                        row.names = temp_meth_genestocolors[,"methprobe"])
# meth_wgcna_genestocolors[,"moduleColors"] <- labels2colors(meth_wgcna_genestocolors[,"moduleLabels"])
# 
# write.table(meth_wgcna_genestocolors, paste0(outfilepathmaster, "meth_WGCNA/", "meth_wgcna_genestocolors.csv"), sep = ",", col.names = NA, row.names = TRUE)
# # meth_wgcna_genestocolors <- read.table(paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/integrative_analyses/run7/", "meth_WGCNA/", "meth_wgcna_genestocolors.csv"), sep = ",", header = TRUE, row.names = 1)
# 
# # 2 - use that to get an eigengene table (733x119) - NOT RUN FOR TIME AND SPACE REASONS
# # meth_eigengenecounttable <- orderMEs(moduleEigengenes(t(methcounttable),
# #                                             meth_wgcna_genestocolors[rownames(methcounttable), "moduleColors"])$eigengenes)
# # Error: vector memory exhausted (limit reached?)
# # Lets try with the NON 0 genes
# # filtered_meth_wgcna_genestocolors <- meth_wgcna_genestocolors[!meth_wgcna_genestocolors[,"moduleColors"] %in% "grey",]
# # meth_eigengenecounttable <- orderMEs(moduleEigengenes(t(methcounttable[rownames(filtered_meth_wgcna_genestocolors),]),
# # filtered_meth_wgcna_genestocolors[rownames(filtered_meth_wgcna_genestocolors), "moduleColors"])$eigengenes)
# # write.table(meth_eigengenecounttable, paste0(outfilepathmaster, "meth_WGCNA/", "meth_eigengene_nogrey_table.csv"), sep = ",", col.names = NA, row.names = TRUE)
# meth_eigengenecounttable <- read.table(paste0(outfilepathmaster, "meth_WGCNA/", "meth_eigengene_nogrey_table.csv"),
#                                        sep = ",", header = TRUE, row.names = 1)
# 
# # 3 - grab only the modules with any kind of annotation (which is limited): 1,2,3,4,5,6,7,8,9,10,12,13,14,22,23
# # meth_wgcna_annotation_goterm_table
# 
# # Add the adjusted p value and add the gene ratio value
# meth_wgcna_annotation_goterm_table_list <- lapply(meth_wgcna_annotation_goterm_table_list, function(x){
#     x[,"over_represented_p.adjust"] <- p.adjust(x[,"over_represented_pvalue"])
#     x[,"GeneRatio"] <- apply(x, 1, function(x2) gsub(" ", "", paste0(x2["numDEInCat"], "/",  x2["numInCat"])))
#     x
# })
# # reformat to match other (just to be consistent)
# # module, size, pval, p.adjust, GeneRatio, ont, ID
# genes_to_color_reftable <- merge(unique(meth_wgcna_genestocolors), data.frame(table(meth_wgcna_genestocolors[,2])),
#                                  by.x = "moduleColors", by.y = "Var1")
# colnames(genes_to_color_reftable) <- c("moduleColors", "moduleLabels", "size")
# # 1+2 - add color back and module size
# meth_wgcna_annotation_table <- merge(bind_rows(meth_wgcna_annotation_goterm_table_list, .id = "groups"),
#                                      genes_to_color_reftable, by.x = "groups", by.y = "moduleLabels", all.x = TRUE)
# # 5 - reorder columns
# meth_wgcna_annotation_table <- meth_wgcna_annotation_table[,c("moduleColors", "size", "over_represented_pvalue", "over_represented_p.adjust", "GeneRatio", "ontology", "term")]
# colnames(meth_wgcna_annotation_table) <- c("module", "size", "pvalue", "p.adjust", "GeneRatio", "ont", "ID")
# rownames(meth_wgcna_annotation_table) <- apply(meth_wgcna_annotation_table[,c("module", "ID")], 1, function(x) paste(x, collapse = "."))
# 
# write.table(meth_wgcna_annotation_table, paste0(outfilepathmaster, "meth_WGCNA/", "GO_meth_wgcna_module_annotations.csv"),
#             sep = ",", col.names = NA, row.names = TRUE)

# 3b - maybe we should try our modules for a eigengene avg heatmap instead?
# 3c - if i want an actual feature heatmap - i need to select a subset, probably start with the most varied? maybe the ones form the DMP?





# meth_wgcna_annotation_goterm_table_file <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/WGCNA_pathway_annotation_hypergeo.rds"
# meth_wgcna_annotation_goterm_table_list <- readRDS(meth_wgcna_annotation_goterm_table_file)# also 119 long
# 
# # Add the adjusted p value and add the gene ratio value
# meth_wgcna_annotation_goterm_table_list <- lapply(meth_wgcna_annotation_goterm_table_list, function(x){
#     x[,"over_represented_p.adjust"] <- p.adjust(x[,"over_represented_pvalue"])
#     x[,"GeneRatio"] <- apply(x, 1, function(x2) gsub(" ", "", paste0(x2["numDEInCat"], "/",  x2["numInCat"])))
#     x
# })
# # reformat to match other (just to be consistent)
# # module, size, pval, p.adjust, GeneRatio, ont, ID
# genes_to_color_reftable <- merge(unique(meth_wgcna_genestocolors), data.frame(table(meth_wgcna_genestocolors[,2])),
#                                  by.x = "moduleColors", by.y = "Var1")
# colnames(genes_to_color_reftable) <- c("moduleColors", "moduleLabels", "size")
# # 1+2 - add color back and module size
# meth_wgcna_annotation_table <- merge(bind_rows(meth_wgcna_annotation_goterm_table_list, .id = "groups"),
#                                      genes_to_color_reftable, by.x = "groups", by.y = "moduleLabels", all.x = TRUE)
# # 5 - reorder columns
# meth_wgcna_annotation_table <- meth_wgcna_annotation_table[,c("moduleColors", "size", "over_represented_pvalue", "over_represented_p.adjust", "GeneRatio", "ontology", "term")]
# colnames(meth_wgcna_annotation_table) <- c("module", "size", "pvalue", "p.adjust", "GeneRatio", "ont", "ID")
# rownames(meth_wgcna_annotation_table) <- apply(meth_wgcna_annotation_table[,c("module", "ID")], 1, function(x) paste(x, collapse = "."))
# 
# write.table(meth_wgcna_annotation_table, paste0(outfilepathmaster, "meth_WGCNA/", "GO_meth_wgcna_module_annotations.csv"),
#             sep = ",", col.names = NA, row.names = TRUE)



# --------------------------------- WES CHIP Analysis read in and incorporation ---------------------------------
# # Read in WES info
# CHIP_sampleID_conversion_table <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/farheen/chip_output_20230406/Sample_ID_with_PATNUM_ischemia_CUSTOM.csv",
#                                              sep = ",", header = TRUE, na.strings = c("", "NA", NA, " "), check.names = FALSE, row.names = 1)
# ## This sample is duplicated and causing issues - im gonna remove for now:
# # 8002310725	11/16/22	036007-001	8002310725
# # S-191208-05937	6/3/22	036007-001	8002301553
# CHIP_sampleID_conversion_table <- CHIP_sampleID_conversion_table[!rownames(CHIP_sampleID_conversion_table) %in% "8002310725",]
# 
# # CHIP mutation table
# CHIP_mutation_table <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/farheen/chip_output_20230406/CHIP_mutation_calls.csv",
#                                              sep = ",", header = TRUE, na.strings = c("", "NA", NA, " "), check.names = FALSE, row.names = 1)
# CHIP_mutation_table_PASS <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/farheen/chip_output_20230406/CHIP_mutation_calls_PASS.csv",
#                                              sep = ",", header = TRUE, na.strings = c("", "NA", NA, " "), check.names = FALSE, row.names = 1)
# CHIP_mutation_table_combined <- merge(CHIP_mutation_table[,c("freq_of_CHIP_mutations", "PATNUM")],
#                                       setNames(CHIP_mutation_table_PASS[,c("freq_of_CHIP_mutations", "PATNUM")], 
#                                                c("freq_of_CHIP_mutations_PASS", "PATNUM")),
#                                       by = "PATNUM", all = TRUE)
# CHIP_mutation_table_combined[is.na(CHIP_mutation_table_combined[,"freq_of_CHIP_mutations_PASS"]),"freq_of_CHIP_mutations_PASS"] <- 0
# 
# # GAM variant tables
# GAM_all_variant_table <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/farheen/chip_output_20230406/GAM_all_variants.csv",
#                                              sep = ",", header = TRUE, na.strings = c("", "NA", NA, " "), check.names = FALSE, row.names = 1)
# GAM_all_variant_table <- GAM_all_variant_table[!rownames(GAM_all_variant_table) %in% "8002310725",]
# rownames(GAM_all_variant_table) <- CHIP_sampleID_conversion_table[
#     rownames(GAM_all_variant_table)[rownames(GAM_all_variant_table) %in% rownames(CHIP_sampleID_conversion_table)], "PATNUM"]
# GAM_all_variant_table_PASS <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/farheen/chip_output_20230406/GAM_PASS_variants.csv",
#                                              sep = ",", header = TRUE, na.strings = c("", "NA", NA, " "), check.names = FALSE, row.names = 1)
# GAM_all_variant_table_PASS <- GAM_all_variant_table_PASS[!rownames(GAM_all_variant_table_PASS) %in% "8002310725",]
# rownames(GAM_all_variant_table_PASS) <- CHIP_sampleID_conversion_table[
#     rownames(GAM_all_variant_table_PASS)[rownames(GAM_all_variant_table_PASS) %in% rownames(CHIP_sampleID_conversion_table)], "PATNUM"]
# GAM_all_variant_table_combined <- merge(GAM_all_variant_table,
#                                       setNames(GAM_all_variant_table_PASS, paste0(colnames(GAM_all_variant_table_PASS), "_PASS")),
#                                       by = "row.names", all = TRUE)
# rownames(GAM_all_variant_table_combined) <- GAM_all_variant_table_combined[,"Row.names"]
# GAM_all_variant_table_combined <- GAM_all_variant_table_combined[,sort(c(colnames(GAM_all_variant_table), paste0(colnames(GAM_all_variant_table_PASS), "_PASS")))]
# 
# # WES full table
# WES_fulltable <- merge(CHIP_mutation_table_combined, GAM_all_variant_table_combined, by.x = "PATNUM", by.y = "row.names")
# rownames(WES_fulltable) <- WES_fulltable[,"PATNUM"]


# CHIP mutation tables
# Theres a duplicate error with these - so I am removing both of these for now: "036007-001", "36007-001"
total_chip_mutant_table <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/farheen/chip_output_20230413/correct_no_of_mutations_all.csv",
                                      sep = ",", header = TRUE)
total_chip_mutant_table_PASS <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/farheen/chip_output_20230413/correct_no_of_mutations_PASS.csv",
                                           sep = ",", header = TRUE)
total_chip_mutant_table_combined <- merge(total_chip_mutant_table[,c("Freq", "PATNUM")],
                                          setNames(total_chip_mutant_table_PASS[,c("Freq", "PATNUM")], c("Freq_PASS", "PATNUM")),
                                          by = "PATNUM", all = TRUE)
total_chip_mutant_table_combined[is.na(total_chip_mutant_table_combined[,"Freq_PASS"]),"Freq_PASS"] <- 0


chip_variant_table <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/farheen/chip_output_20230413/correct_variants_all.csv",
                                      sep = ",", header = TRUE)
chip_variant_full_table <- dcast(chip_variant_table, PATNUM ~ Gene.refGene, fun.aggregate = length)
chip_variant_table_PASS <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/farheen/chip_output_20230413/correct_variants_PASS.csv",
                                           sep = ",", header = TRUE)
chip_variant_full_table_PASS <- dcast(chip_variant_table_PASS, PATNUM ~ Gene.refGene, fun.aggregate = length)
colnames(chip_variant_full_table_PASS) <- c("PATNUM", paste0(colnames(chip_variant_full_table_PASS)[!colnames(chip_variant_full_table_PASS) %in% "PATNUM"], "_PASS"))

chip_variant_full_table_combined <- merge(chip_variant_full_table, chip_variant_full_table_PASS, by = "PATNUM", all = TRUE)
chip_variant_full_table_combined[is.na(chip_variant_full_table_combined)] <- 0
rownames(chip_variant_full_table_combined) <- chip_variant_full_table_combined[,"PATNUM"]
chip_variant_full_table_combined <- chip_variant_full_table_combined[,sort(colnames(chip_variant_full_table_combined))]






# --------------------------------- Adding to and expanding the bioreptable --------------------------------- 
## Dummy heatmap to compare the disease labels we have (CTNDV50, IMGDEGIS, DUKESCORE) in our cohort
bioreptable_waddons <- addon_and_edit_biorep_table(bioreptable = bioreptable,
                                                   addclustertypes = TRUE,
                                                   biomarkertable_sel, SOIall, rna_clustermembership_table, meth_clustermembership_table)
write.table(bioreptable_waddons, paste0(outfilepathmaster, "bioreptable_waddons.csv"), sep = ",", col.names = NA, row.names = TRUE)

# --------------------------------- Omics Overlap ---------------------------------
## Overlap of our omics types
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

dir.create(paste0(outfilepathmaster, "omics_sample_comparison/"), showWarnings = FALSE, recursive = TRUE)
write.table(overlaptable, paste0(outfilepathmaster, "omics_sample_comparison/", "testoverlap_overlaptable.csv"), 
            sep = ",", col.names = TRUE, row.names = FALSE)
write.table(overlapgrouptab, paste0(outfilepathmaster, "omics_sample_comparison/", "testoverlap_overlapgrouptable.csv"), 
            sep = ",", col.names = TRUE, row.names = FALSE)
write.table(overlapsummary, paste0(outfilepathmaster, "omics_sample_comparison/", "testoverlap_overlapsummary.csv"), 
            sep = ",", col.names = TRUE, row.names = TRUE)
pdf(paste0(outfilepathmaster, "omics_sample_comparison/", "testoverlap_overlapvenn.pdf"))
grid.draw(vennplot)
junk <- dev.off()

# --------------------------------- Disease label plot --------------------------------- 
## Dummy heatmap to compare the disease labels we have (CTNDV50, IMGDEGIS, DUKESCORE) in our cohort
labelplot <- ischemia_disease_metrics_plot(bioreptable_waddons, colorguide, plotmetrics = c("IMGDEGIS_01_2_3", "CTNDV70_combo0noneval"))
dir.create(paste0(outfilepathmaster, "cohort_tabulations/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(outfilepathmaster, "cohort_tabulations/", "ischemia_test_labels_hm.pdf"), 11.5, 8, useDingbats = FALSE)
draw(labelplot[[1]])
junk <- dev.off()

# Comparing association 
quicktest1 <- bioreptable_waddons[SOIeither,c("PATNUM", "IMGDEGIS_01_2_3", "CTNDV70_combo0noneval")]
# quicktest1[is.na(quicktest1)] <- "NotAv"
# quicktest1 <- na.omit(quicktest1)
quicktest1 <- quicktest1[quicktest1[,"IMGDEGIS_01_2_3"] %in% c("NoneMild", "Moderate", "Severe") & quicktest1[,"CTNDV70_combo0noneval"] %in% c("1", "2", "3"),]
chisq.test(quicktest1[,"IMGDEGIS_01_2_3"], quicktest1[,"CTNDV70_combo0noneval"], simulate.p.value = TRUE)
sumstattab_seq <- summarize_table(intable = quicktest1[SOIeither,c("PATNUM", "IMGDEGIS_01_2_3", "CTNDV70_combo0noneval")], groupvar = "IMGDEGIS_01_2_3", 
                                  outfile = paste0(outfilepathmaster, "cohort_tabulations/", "COI_tabulation_",  "CUSTOM_IMGDEGIS_NONEMILD_vs_CTNDV70", ".csv"), calc_stats = TRUE)

# --------------------------------- Cohort Tabulations --------------------------------- 
## Tabulated out stats of interest per group comparison - RNAseq subtype, Methsubtype, etc.
dir.create(paste0(outfilepathmaster, "cohort_tabulations/"), showWarnings = FALSE, recursive = TRUE)
metaCOI <- c("CUSTOM_maxfollow", 
             "AGE_RAND", "CAGE_RND", "SEX", "RACE", "CUSTOM_RACE_AGGR", "ETHNIC",
             "BMI", "CUSTOM_BMI_obese", "HYPTENSE", "DIABETES", "DIABTRT", "SMOKSTAT", "CUSTOM_SMOKSTAT_currentsmoker", 
             "FAMHXPHD", "PRIORMI", "PRIPCIYN", "HXCABG", "CUSTOM_PCIorCABG", "XPRIFALYN", "PRICEREB", "PRIPVDYN",
             "LVEFCONT", "CUSTOM_LVEFCONT_LVEFund55", "EFLT45", "EGFR", "CUSTOM_EGFR_EGFRund60", "DIALYSIS",
             "HDLC", "TRIGLYC", "TOTCHOL", "RVBPDIA",  "HEMOGLOB",
             "LDLC", "CUSTOM_LDLC_LDLCund70", "RVBPSYS", "CUSTOM_RVBPSYS_RVBPSYSund140", "HEMOGA1C",
             "MESTATIN", "CUSTOM_RVBPSYS_RVBPSYSund140", "CUSTOM_LDLCund70andMESTATIN", "MEAPAC",
             "DEGRISCH", "CUSTOM_DEGRISCH_NONEMILD", "IMGDEGIS", "CUSTOM_IMGDEGIS_NONEMILD", "IMGDEGIS_01_2_3",
             "CTNDV70", "CTMULT70", "CTNDV70_combo0noneval",
             "CTNDV50", "CTMULT50",
             "DUKESCORE",
             "TREATMNT", "CKD", "CUSTOM_CKD_TRT",
             "SAQ7_AF", "SAQ7_PL", "SAQ7_QL", "SAQ7_SS",
             "MEAPASP", "MEAPCLOP", "HIGHSTAT", "MELLEZET", "MEAHACE", "MEAHARB", "CUSTOM_ACEIandARB", "INSULDBT",
             "C_PRIMARY", "C_CVDMI_P", "C_ACD", "C_CVD", "C_ACDMI_P", "C_MIPRIM", "C_RCA", "C_HOSP_UA", "C_HOSP_HF",
             colnames(biomarkertable)[grepl("_clean$", colnames(biomarkertable))],
             "rna_4cluster_w3AB", "rna_4cluster", "meth_4cluster", "meth_3cluster"
)


## Create tabulation table, with characterize of certain cols, and selection of certain samples
# tabulation_table <- merge(bioreptable_waddons[,c("PATNUM", metaCOI)], combined_clustermembership_table, by = "PATNUM", all = TRUE)
# rownames(tabulation_table) <- tabulation_table[,"PATNUM"]
tabulation_table <- bioreptable_waddons[, c("PATNUM", metaCOI)]
write.table(tabulation_table, paste0(outfilepathmaster, "cohort_tabulations/", "tabulation_table.csv"), sep =",", row.names = TRUE, col.names = NA)

tab_param_list <- list(list(SOImeth = SOImeth, groupvar = NULL), list(SOIrna = SOIrna, groupvar = NULL),
                       list(SOIall = SOIall, groupvar = NULL), list(SOIboth = SOIboth, groupvar = NULL),
                       list(SOIeither = SOIeither, groupvar = NULL),
                       list(SOIrna = SOIrna, groupvar = "rna_4cluster_w3AB"), list(SOIrna = SOIrna, groupvar = "rna_4cluster"), 
                       list(SOImeth = SOImeth, groupvar = "meth_4cluster"), list(SOImeth = SOImeth, groupvar = "meth_3cluster"),
                       list(SOIall_CKD = SOIall, groupvar = "CKD"))
for (param_sel in tab_param_list) {
    filenameparam <- if(is.null(param_sel[[2]])) NULL else {paste0("__", param_sel[[2]], "_stats")}
    calc_stats_param <- ifelse(is.null(param_sel[[2]]), FALSE, TRUE)
    summarystatfile <- paste0(outfilepathmaster, "cohort_tabulations/", "COI_tabulation_",  names(param_sel)[1], filenameparam, ".csv")
    sumstattab_seq <- summarize_table(intable = tabulation_table[param_sel[[1]],], groupvar = param_sel[[2]], 
                                      outfile = summarystatfile, calc_stats = calc_stats_param)
    clean_sumstattab_seq <- clean_summarize_table_output(sumstattable_input = sumstattab_seq, addpercents = "vertical", 
                                                         contsummary = c("median", "iqr"), roundpvaldigits = 3)
    cleansummarystatfile <- paste0(outfilepathmaster, "cohort_tabulations/", "COI_tabulation_",  names(param_sel)[1], filenameparam, "_clean.csv")
    write.table(clean_sumstattab_seq, cleansummarystatfile, sep = ",", col.names = TRUE, row.names = FALSE)
    
    ## Adding in here a specific WES section
    dir.create(paste0(outfilepathmaster, "cohort_tabulations/", "WES_tabulation/"), showWarnings = FALSE, recursive = TRUE)
    WES_tabulation_table <- merge(chip_variant_full_table_combined,
                                  bioreptable_waddons[, c("PATNUM", "rna_4cluster", "meth_3cluster", "meth_4cluster", "rna_4cluster_w3AB", "CKD")],
                                  by = "PATNUM", all = TRUE)
    rownames(WES_tabulation_table) <- WES_tabulation_table[,"PATNUM"]
    
    summarystatfile <- paste0(outfilepathmaster, "cohort_tabulations/", "WES_tabulation/", "WES_tabulation_",  names(param_sel)[1], filenameparam, ".csv")
    sumstattab_seq <- summarize_table(intable = WES_tabulation_table[param_sel[[1]],], groupvar = param_sel[[2]], 
                                      outfile = summarystatfile, calc_stats = calc_stats_param)
    clean_sumstattab_seq <- clean_summarize_table_output(sumstattable_input = sumstattab_seq, addpercents = "vertical", 
                                                         contsummary = c("median", "iqr"), roundpvaldigits = 3)
    cleansummarystatfile <- paste0(outfilepathmaster, "cohort_tabulations/", "WES_tabulation/", "WES_tabulation_",  names(param_sel)[1], filenameparam, "_clean.csv")
    write.table(clean_sumstattab_seq, cleansummarystatfile, sep = ",", col.names = TRUE, row.names = FALSE)

}

## Quick tab for comparing IMGDEGIS custom and CTNDV70 custom
apply(bioreptable_waddons[SOIeither,c("IMGDEGIS_01_2_3", "CTNDV70_combo0noneval")], 2, table)
apply(bioreptable_waddons[SOIeither,c("IMGDEGIS_01_2_3", "CTNDV70_combo0noneval")], 2, function(x) {
    rbind(table(x, useNA = "ifany"), table(x, useNA = "ifany")/length(x))
})


# --------------------------------- Rectangle plots for descriptor characteristics --------------------------------- 
## Draw a "rectangle plot" for each of our descriptor characteristics
dir.create(paste0(outfilepathmaster, "cohort_tabulations/", "rectplots/"), showWarnings = FALSE, recursive = TRUE)
COI <- c("SEX", "CAGE_RND", "DIABETES", "IMGDEGIS", "CTNDV70", "DUKESCORE", "IMGDEGIS_01_2_3", "CTNDV70_combo0noneval")
COI_colref_list <- create_color_reference(COI=COI, colorguide=colorguide)
# Create a rect plot for each COI and each omics type
param_intable <- expand.grid(COI = COI, omicstype = c("RNA", "METH"))
param_inlist <- split(param_intable, seq(nrow(param_intable)))

for (selected_parameters in param_inlist) {
    COI_sel <- as.character(selected_parameters[,"COI"])
    omicstype_sel <- as.character(selected_parameters[,"omicstype"])
    if (omicstype_sel == "RNA") {sample_sel <- SOIrna}
    if (omicstype_sel == "METH") {sample_sel <- SOImeth}
    bioreptable_sel <- factor(bioreptable_waddons[sample_sel, COI_sel], levels = COI_colref_list[[COI_sel]][["ORDER"]])
    
    rect_plot <- create_rect_plot(bioreptable_sel = bioreptable_sel, COI_colref_list = COI_colref_list, addnumbers = TRUE) 
    pdf(paste0(outfilepathmaster, "cohort_tabulations/", "rectplots/", omicstype_sel, "_", COI_sel, "_rectplot.pdf"), 
        useDingbats = FALSE, width = 40, height = 4)
    print(rect_plot)
    junk <- dev.off()
    
    rect_plot <- create_rect_plot(bioreptable_sel = bioreptable_sel, COI_colref_list = COI_colref_list, plotaspercents = TRUE, addnumbers = TRUE) 
    pdf(paste0(outfilepathmaster, "cohort_tabulations/", "rectplots/", omicstype_sel, "_", COI_sel, "_rectplot_w_percent.pdf"), 
        useDingbats = FALSE, width = 40, height = 4)
    print(rect_plot)
    junk <- dev.off()
}




# --------------------------------- clean volcano plots and gsea tables - RNA --------------------------------- 
dir.create(paste0(outfilepathmaster, "clean_dge_figures/"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(outfilepathmaster, "clean_dge_figures/dge_volcanoplots/"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(outfilepathmaster, "clean_dge_figures/dge_gseaplots/"), showWarnings = FALSE, recursive = TRUE)

# Read in DESeq and GSEA results
dge_comp_list <- c("comp_ischemia__Sev_v_MildNone", "comp_anatomy70__3v_v_1v")
deseq_table_list <- list()
deseq_path <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/deseq/"
gsea_table_list <- list()
gsea_path <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/gsea/"
# GOI_table <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run3_rmoutliers2/integrative_analyses/volcano_GOI_table.csv", sep = ",", header = TRUE, check.names = FALSE)
# POI_table <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run3_rmoutliers2/integrative_analyses/dge_POI_table.csv", sep = ",", header = TRUE, check.names = FALSE)
GOI_table <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/integrative_analyses/manual_annotation_tables/volcano_GOI_table_v3.csv", sep = ",", header = TRUE, check.names = FALSE)
POI_table <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/integrative_analyses/manual_annotation_tables/dge_POI_table_v2.csv", sep = ",", header = TRUE, check.names = FALSE, na.strings = c("", NA, "NA", " "))


for (dge_comp in dge_comp_list) {
    intable_deseq <- read.table(paste0(deseq_path, dge_comp, "/deseq_results_", dge_comp, ".csv"), 
                          sep = ",", header = TRUE, row.names = 1)
    deseq_table_list[[dge_comp]] <- intable_deseq
    GOI <- GOI_table[,dge_comp]
    if (sum(!is.na(GOI)) == 0) {GOI <- NULL}
    
    # set stat params:
    if (dge_comp == "comp_ischemia__Sev_v_MildNone") {
        pval_test = comp_ischemia__Sev_v_MildNone__PVALTEST
        pval_cutoff = comp_ischemia__Sev_v_MildNone__PVALCUTOFF
        log2fc_cutoff = comp_ischemia__Sev_v_MildNone__FCCUTOFF
    }
    if (dge_comp %in% c("comp_anatomy__3v_v_1v", "comp_anatomy70__3v_v_1v")) {
        pval_test = comp_anatomy70__3v_v_1v__PVALTEST
        pval_cutoff = comp_anatomy70__3v_v_1v__PVALCUTOFF
        log2fc_cutoff = comp_anatomy70__3v_v_1v__FCCUTOFF
    }
    
    volcanoplot_out <- clean_volcano_plot_function(deseqtable=intable_deseq, nameparam=dge_comp, labeledgenes = GOI,
                                                   pval_test = pval_test, pval_cutoff = pval_cutoff, log2fc_cutoff = log2fc_cutoff, returnplottab = TRUE)
    if (dge_comp == "comp_ischemia__Sev_v_MildNone") {
        volcanoplot_out[[1]] <- volcanoplot_out[[1]] + scale_y_continuous(expand = expansion(mult = c(0.1, 0.3), 
                                                                                      add = c(0, 0)))
    }
    pdf(paste0(outfilepathmaster, "clean_dge_figures/", "dge_volcanoplots/", dge_comp, "_volcano_plot.pdf"), width = 10, height = 10)
    print(volcanoplot_out[[1]])
    junk <- dev.off()
    write.table(volcanoplot_out[[2]], paste0(outfilepathmaster, "clean_dge_figures/", "dge_volcanoplots/", dge_comp, "_volcano_plottab.csv"),
                sep = ",", col.names = NA, row.names = TRUE)
    
    
    # pval_summary <- destatplot(intable_deseq[,c("log2FoldChange", "pvalue")], log2fclevels = c(0, 0.25, 0.5, 1.0))[[2]]
    pval_summary <- destatplot(intable_deseq[,c("log2FoldChange", "pvalue")], log2fclevels = c(0, 0.25, 0.5, 1.0), split_posneg = TRUE)[[2]]
    write.table(pval_summary, paste0(outfilepathmaster, "clean_dge_figures/", "dge_volcanoplots/", dge_comp, "_pval_summary_table.csv"), sep = ",", col.names = NA, row.names = TRUE)
    padj_summary <- destatplot(intable_deseq[,c("log2FoldChange", "padj")], log2fclevels = c(0, 0.25, 0.5, 1.0), split_posneg = TRUE)[[2]]
    write.table(padj_summary, paste0(outfilepathmaster, "clean_dge_figures/", "dge_volcanoplots/", dge_comp, "_padj_summary_table.csv"), sep = ",", col.names = NA, row.names = TRUE)
    
    
    intable_gsea <- read.table(paste0(gsea_path, dge_comp, "/", dge_comp, "_gsea_GO.csv"), 
                          sep = ",", header = TRUE, row.names = 1)
    gsea_table_list[[dge_comp]] <- intable_gsea
    write.table(intable_gsea, paste0(outfilepathmaster, "clean_dge_figures/", "dge_gseaplots/", dge_comp, "_gsea_table.csv"), sep = ",", col.names = NA, row.names = TRUE)
    
    ## Select POI and plot
    POI <- na.omit(POI_table[,dge_comp])
    gsea_out <- clean_gsea_plot_function(gseatable=intable_gsea, nameparam=dge_comp, pathwayselect = POI)
    pdf(paste0(outfilepathmaster, "clean_dge_figures/", "dge_gseaplots/", dge_comp, "_gsea_plot.pdf"), width = 10, height = 10)
    print(gsea_out)
    junk <- dev.off()
    write.table(intable_gsea[POI,], paste0(outfilepathmaster, "clean_dge_figures/", "dge_gseaplots/", dge_comp, "_gsea_plottable.csv"),
                sep = ",", col.names = TRUE, row.names = FALSE)
    
}





# --------------------------------- clean volcano plots and psea tables - Meth --------------------------------- 
dir.create(paste0(outfilepathmaster, "clean_dmp_figures/"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(outfilepathmaster, "clean_dmp_figures/dmp_volcanoplots/"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(outfilepathmaster, "clean_dmp_figures/dmp_pseaplots/"), showWarnings = FALSE, recursive = TRUE)

# Read in DESeq and PSEA results
dmp_comp_list <- c("DMP_ischemia__Sev_v_MildNone", "DMP_anatomy70__3v_v_1v")
dmp_table_list <- list()
dmp_path <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/meth_processing/run3/asr_corrected_DMP/"
psea_table_list <- list()
psea_path <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/meth_processing/run3/asr_corrected_DMP/PSEA/run2/c5/"
GOI_table <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/integrative_analyses/manual_annotation_tables/volcano_GOI_table_v3.csv", sep = ",", header = TRUE, check.names = FALSE)
POI_table <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/integrative_analyses/manual_annotation_tables/dge_POI_table_v2.csv", sep = ",", header = TRUE, check.names = FALSE, na.strings = c("", NA, "NA", " "))


for (dmp_comp in dmp_comp_list) {
    # set stat params:
    if (dmp_comp == "DMP_ischemia__Sev_v_MildNone") {
        pval_test = DMP_ischemia__Sev_v_MildNone__PVALTEST
        pval_cutoff = DMP_ischemia__Sev_v_MildNone__PVALCUTOFF
        log2fc_cutoff = DMP_ischemia__Sev_v_MildNone__FCCUTOFF
    }
    if (dmp_comp %in% c("DMP_anatomy70__3v_v_1v")) {
        pval_test = DMP_anatomy70__3v_v_1v__PVALTEST
        pval_cutoff = DMP_anatomy70__3v_v_1v__PVALCUTOFF
        log2fc_cutoff = DMP_anatomy70__3v_v_1v__FCCUTOFF
    }
    
    intable_dmp <- read.table(paste0(dmp_path, dmp_comp, "/DMP_analysis_and_figs/", dmp_comp, "_DMP_table.csv"), 
                                sep = ",", header = TRUE, row.names = 1)
    dmp_table_list[[dmp_comp]] <- intable_dmp
    GOI <- GOI_table[,dmp_comp]
    if (sum(!is.na(GOI)) == 0) {GOI <- NULL}
    
    # Do some formatting for plotting purposes:
    DMPtable_formatted <- intable_dmp[,c("gene", "deltaBeta", colnames(intable_dmp)[grepl("_avg", colnames(intable_dmp))], "pvalue", "padj")]
    colnames(DMPtable_formatted)[2] <- "log2FoldChange"
    DMPtable_formatted[,"cg_probe"] <- rownames(DMPtable_formatted)
    
    # Grab labeledgenes
    grab_sigprobes <- DMPtable_formatted[DMPtable_formatted[,"gene"] %in% GOI &
                                         DMPtable_formatted[,pval_test] < pval_cutoff & 
                                         abs(DMPtable_formatted[,"log2FoldChange"]) > log2fc_cutoff,]
    grab_sigprobes <- rownames(grab_sigprobes[!duplicated(grab_sigprobes[,"gene"]),])
    DMPtable_formatted[,"labeledgenes"] <- ifelse(rownames(DMPtable_formatted) %in% grab_sigprobes, as.character(DMPtable_formatted[,"gene"]),
                                                  rownames(DMPtable_formatted))
    rownames(DMPtable_formatted) <- make.unique(DMPtable_formatted[,"labeledgenes"])
    labeledgenes <- make.unique(DMPtable_formatted[,"labeledgenes"])[!grepl("^cg", make.unique(DMPtable_formatted[,"labeledgenes"]))]
    
    volcanoplot_out <- clean_volcano_plot_function(deseqtable=DMPtable_formatted, nameparam=dmp_comp, labeledgenes = labeledgenes,
                                                   pval_test = pval_test, pval_cutoff = pval_cutoff, log2fc_cutoff = log2fc_cutoff, returnplottab = TRUE)
    # if (dge_comp == "comp_ischemia__Sev_v_MildNone") {
    #     volcanoplot_out <- volcanoplot_out + scale_y_continuous(expand = expansion(mult = c(0.1, 0.3), 
    #                                                                                add = c(0, 0)))
    # }
    pdf(paste0(outfilepathmaster, "clean_dmp_figures/", "dmp_volcanoplots/", dmp_comp, "_volcano_plot.pdf"), width = 10, height = 10)
    print(volcanoplot_out[[1]])
    junk <- dev.off()
    
    # I need to reverse the hack here and then add the names back for ref - so all rownames need to be probes, but save the naming hack to a separate column
    write.table(cbind(volcanoplot_out[[2]], DMPtable_formatted[,c("cg_probe", "labeledgenes")]),
                paste0(outfilepathmaster, "clean_dmp_figures/", "dmp_volcanoplots/", dmp_comp, "_volcano_plottab.csv"),
                sep = ",", col.names = NA, row.names = TRUE)
    
    
    # pval sum tabs
    pval_summary <- destatplot(intable_dmp[,c("deltaBeta", "pvalue")], log2fclevels = c(0, 0.01, 0.03, 0.05, 0.1), split_posneg = TRUE)[[2]]
    write.table(pval_summary, paste0(outfilepathmaster, "clean_dmp_figures/", "dmp_volcanoplots/", dmp_comp, "_pval_summary_table.csv"), sep = ",", col.names = NA, row.names = TRUE)
    padj_summary <- destatplot(intable_dmp[,c("deltaBeta", "padj")], log2fclevels = c(0, 0.01, 0.03, 0.05, 0.1), split_posneg = TRUE)[[2]]
    write.table(padj_summary, paste0(outfilepathmaster, "clean_dmp_figures/", "dmp_volcanoplots/", dmp_comp, "_padj_summary_table.csv"), sep = ",", col.names = NA, row.names = TRUE)
    
    
    # Clean psea plots
    intable_psea <- read.table(paste0(psea_path, dmp_comp, "/PSEA/", dmp_comp, "_PSEA_CUSTOM_GO.csv"), 
                               sep = ",", header = TRUE, row.names = 1)
    psea_table_list[[dmp_comp]] <- intable_psea
    write.table(intable_psea, paste0(outfilepathmaster, "clean_dmp_figures/", "dmp_pseaplots/", dmp_comp, "_psea_table.csv"), sep = ",", col.names = NA, row.names = TRUE)
    
    POI <- na.omit(POI_table[,dmp_comp])
    psea_out <- clean_gsea_plot_function(gseatable=intable_psea, nameparam=dmp_comp, pathwayselect = POI)
    pdf(paste0(outfilepathmaster, "clean_dmp_figures/", "dmp_pseaplots/", dmp_comp, "_psea_plot.pdf"), width = 10, height = 10)
    print(psea_out)
    junk <- dev.off()
    write.table(intable_psea[POI,], paste0(outfilepathmaster, "clean_dmp_figures/", "dmp_pseaplots/", dmp_comp, "_psea_plottable.csv"),
                sep = ",", col.names = TRUE, row.names = FALSE)
    
}




# --------------------------------- Compare IMGDEGIS and CTNDV70 so we can look for things that are similar and unique --------------------------------- 
dir.create(paste0(outfilepathmaster, "IMGDEGIS_v_CTNDV70/"), showWarnings = FALSE, recursive = TRUE)
diffexp_list <- list(deseq_comp = deseq_table_list, gsea_comp = gsea_table_list,
                     dmp_comp = dmp_table_list, psea_comp = psea_table_list)

for (analysis_comp_num in seq_along(diffexp_list)) {
    # create sub folder out
    analysis_label_sel <- names(diffexp_list)[analysis_comp_num]
    suboutfolder <- paste0(outfilepathmaster, "IMGDEGIS_v_CTNDV70/", analysis_label_sel, "/")
    dir.create(suboutfolder, showWarnings = FALSE, recursive = TRUE)

    # grab all of the parameters from the run we are running
    if (analysis_label_sel %in% c("deseq_comp", "gsea_comp")) {
        IMGDEGIS_diffexp_table <- diffexp_list[[analysis_comp_num]][["comp_ischemia__Sev_v_MildNone"]]
        CTNDV70_diffexp_table <- diffexp_list[[analysis_comp_num]][["comp_anatomy70__3v_v_1v"]]
    }
    if (analysis_label_sel %in% c("dmp_comp", "psea_comp")) {
        IMGDEGIS_diffexp_table <- diffexp_list[[analysis_comp_num]][["DMP_ischemia__Sev_v_MildNone"]]
        CTNDV70_diffexp_table <- diffexp_list[[analysis_comp_num]][["DMP_anatomy70__3v_v_1v"]]
    }

    if (analysis_label_sel %in% "deseq_comp") {
        effectsize_label <- "log2FoldChange"
        IMGDEGIS_stattype <- comp_ischemia__Sev_v_MildNone__PVALTEST
        CTNDV70_stattype <- comp_anatomy70__3v_v_1v__PVALTEST
        IMGDEGIS_diffexp_table[,"sigfeat"] <- ifelse(abs(IMGDEGIS_diffexp_table[,effectsize_label]) > comp_ischemia__Sev_v_MildNone__FCCUTOFF &
                                                     IMGDEGIS_diffexp_table[,IMGDEGIS_stattype] < comp_ischemia__Sev_v_MildNone__PVALCUTOFF & 
                                                     !is.na(IMGDEGIS_diffexp_table[,IMGDEGIS_stattype]), "sigfeat", NA)
        CTNDV70_diffexp_table[,"sigfeat"] <- ifelse(abs(CTNDV70_diffexp_table[,effectsize_label]) > comp_anatomy70__3v_v_1v__FCCUTOFF &
                                                    CTNDV70_diffexp_table[,CTNDV70_stattype] < comp_anatomy70__3v_v_1v__PVALCUTOFF & 
                                                    !is.na(CTNDV70_diffexp_table[,CTNDV70_stattype]), "sigfeat", NA)
    }
    if (analysis_label_sel %in% "dmp_comp") {
        effectsize_label <- "deltaBeta"
        IMGDEGIS_stattype <- DMP_ischemia__Sev_v_MildNone__PVALTEST
        CTNDV70_stattype <- DMP_anatomy70__3v_v_1v__PVALTEST
        IMGDEGIS_diffexp_table[,"sigfeat"] <- ifelse(abs(IMGDEGIS_diffexp_table[,effectsize_label]) > DMP_ischemia__Sev_v_MildNone__FCCUTOFF &
                                                         IMGDEGIS_diffexp_table[,IMGDEGIS_stattype] < DMP_ischemia__Sev_v_MildNone__PVALCUTOFF & 
                                                         !is.na(IMGDEGIS_diffexp_table[,IMGDEGIS_stattype]), "sigfeat", NA)
        CTNDV70_diffexp_table[,"sigfeat"] <- ifelse(abs(CTNDV70_diffexp_table[,effectsize_label]) > DMP_anatomy70__3v_v_1v__FCCUTOFF &
                                                        CTNDV70_diffexp_table[,CTNDV70_stattype] < DMP_anatomy70__3v_v_1v__PVALCUTOFF & 
                                                        !is.na(CTNDV70_diffexp_table[,CTNDV70_stattype]), "sigfeat", NA)
        
        # Im gonna only take probes with genes here too - cause why not
        IMGDEGIS_diffexp_table <- IMGDEGIS_diffexp_table[IMGDEGIS_diffexp_table[,"gene"] != "",]
        CTNDV70_diffexp_table <- CTNDV70_diffexp_table[CTNDV70_diffexp_table[,"gene"] != "",]
        rownames(IMGDEGIS_diffexp_table) <- make.unique(IMGDEGIS_diffexp_table[,"gene"])
        rownames(CTNDV70_diffexp_table) <- make.unique(CTNDV70_diffexp_table[,"gene"])
        
    }
    if (analysis_label_sel %in% c("gsea_comp", "psea_comp")) {
        effectsize_label <- "NES"
        IMGDEGIS_stattype <- "p.adjust"
        CTNDV70_stattype <- "p.adjust"
        IMGDEGIS_diffexp_table[,"sigfeat"] <- ifelse(abs(IMGDEGIS_diffexp_table[,IMGDEGIS_stattype]) < 0.05, "sigfeat", NA)
        CTNDV70_diffexp_table[,"sigfeat"] <- ifelse(abs(CTNDV70_diffexp_table[,CTNDV70_stattype]) < 0.05, "sigfeat", NA)
    }
    
    
    # ## Process meth table for plotting
    feature_plottable <- suppressWarnings(merge(IMGDEGIS_diffexp_table[,c(effectsize_label, IMGDEGIS_stattype, "sigfeat")], 
                                                CTNDV70_diffexp_table[,c(effectsize_label, CTNDV70_stattype, "sigfeat")],
                                                by = "row.names", suffixes = c("", ""), all = TRUE))
    dimnames(feature_plottable) <- list(feature_plottable[,"Row.names"], c("feature", paste0(colnames(feature_plottable)[2:ncol(feature_plottable)], c(rep("_IMGDEGIS", 3), rep("_CTNDV70", 3)))))
    
    
    feature_plottable[,"featgroup"] <- ifelse(!is.na(feature_plottable[,"sigfeat_IMGDEGIS"]) & !is.na(feature_plottable[,"sigfeat_CTNDV70"]) & 
                                              abs(rowSums(sign(feature_plottable[,grepl(effectsize_label, colnames(feature_plottable))]))) == 2, "sig_BOTH_samedir",
                                       ifelse(!is.na(feature_plottable[,"sigfeat_IMGDEGIS"]) & !is.na(feature_plottable[,"sigfeat_CTNDV70"]) & 
                                              abs(rowSums(sign(feature_plottable[,grepl(effectsize_label, colnames(feature_plottable))]))) == 0, "sig_BOTH_diffdir",
                                       ifelse(!is.na(feature_plottable[,"sigfeat_IMGDEGIS"]) & feature_plottable[,paste0(CTNDV70_stattype, "_CTNDV70")] > 0.2, "sig_IMGDEGIS",
                                       ifelse(!is.na(feature_plottable[,"sigfeat_CTNDV70"]) & feature_plottable[,paste0(IMGDEGIS_stattype, "_IMGDEGIS")] > 0.2, "sig_CTNDV70",
                                       "notsig"))))
    
    # write out this table
    write.table(feature_plottable, paste0(suboutfolder, analysis_label_sel, "_feature_plottable.csv"), sep = ",", col.names = TRUE, row.names = FALSE)
    # # Use this to grab some genes to label per quadrant (5 max?)
    maxnumlabeledgenesquad <- 7
    datalabels <- unlist(unname(lapply(split(feature_plottable[!feature_plottable[,"featgroup"] %in% "notsig",],
                                             feature_plottable[!feature_plottable[,"featgroup"] %in% "notsig","featgroup"]),
                                       function(subdf) {
                                           ## in an extra few steps here for the pathway analysis - remove HP pathways (theyre dumb)
                                           subdf <- subdf[!grepl("HP_", subdf[,"feature"]),]
                                           featgroup_label <- gsub("sig_", "", unique(subdf[,"featgroup"]))
                                           rankcol <- paste0(effectsize_label, "_", featgroup_label)
                                           if (grepl("BOTH", rankcol)) {
                                               subdf[,"comboeffect"] <- rowSums(abs(subdf[,grepl(effectsize_label, colnames(subdf))]))
                                               rankcol <- "comboeffect"
                                           }
                                           ## Lets rank and grab in order actually as opposed to a random sampling
                                           subdf <- subdf[order(abs(subdf[,rankcol]), decreasing = TRUE),]
                                           # subdf[sample(1:nrow(subdf), min(maxnumlabeledgenesquad, nrow(subdf))), "feature"]
                                           subdf[1:min(maxnumlabeledgenesquad, nrow(subdf)), "feature"]
                                           })))
    
    plotting_color_vector <- c(
        "notsig" = "lightgrey", "sig_BOTH_samedir" = "darkgreen","sig_BOTH_diffdir" = "blue",
        "sig_IMGDEGIS" = colorguide[colorguide[,"category"] %in% "IMGDEGIS" & colorguide[,"feature"] %in% "Severe", "color"],
        "sig_CTNDV70" = colorguide[colorguide[,"category"] %in% "CTNDV70" & colorguide[,"feature"] %in% "3", "color"]
    )
    
    # Plot the gene comparison scatter
    pout <- scatter_plotter(indata = feature_plottable[,grepl(effectsize_label, colnames(feature_plottable))],
                            # datalabels = datalabels,
                            colorvar = feature_plottable[,"featgroup", drop=FALSE],
                            labsparam = list(title = paste0(analysis_label_sel, " ", effectsize_label, " comparison, labeled sig genes"),
                                             x = paste0("IMGDEGIS ", effectsize_label), y = paste0("CTNDV70 ", effectsize_label))
    )
    pout <- pout + geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 0, linetype = 2)
    pout <- pout + scale_color_manual(breaks = names(plotting_color_vector), values = unname(plotting_color_vector))
    # pdf(paste0(suboutfolder, analysis_label_sel, "_sig_gene_scatterplot.pdf"), 7, 7, useDingbats = FALSE)
    png(paste0(suboutfolder, analysis_label_sel, "_sig_gene_scatterplot.png"), 1400, 1400, res = 300)
    print(pout)
    junk <- dev.off()
    # 
    # 
    # ## Hypergeo tests of each quadrant:
    if (analysis_label_sel %in% c("deseq_comp", "dmp_comp")) {
        for (sigfeat_type in unique(feature_plottable[!feature_plottable[,"featgroup"] %in% "notsig","featgroup"])) {
            speciesparam <- "Homo sapiens"
            sigfeat_genes <- feature_plottable[feature_plottable[,"featgroup"] %in% sigfeat_type, "feature"]
            if (analysis_label_sel == "dmp_comp") {sigfeat_genes <- unique(gsub("\\.[^\\.]*$", "", sigfeat_genes))}
            hypergeo_genetest_out = data.frame(hypergeo_genetest(data.frame(sigfeat_genes), statcutoffparam = "DUMMY",
                                                                 genesetparam = c("C5"), speciesparam = speciesparam)$enricherUPout)
            write.table(hypergeo_genetest_out, paste0(suboutfolder, analysis_label_sel, sigfeat_type, "_hypergeo.csv"), sep = ",", col.names = TRUE, row.names = FALSE)
        }
    }
    
}


# Furthering this process by then also comparing to our RNA vs Meth comparison we are doing - just tryna find the perfect pathways to focus on.............
dir.create(paste0(outfilepathmaster, "IMGDEGIS_v_CTNDV70/", "overlap_pathway_lists/"), showWarnings = FALSE, recursive = TRUE)
IMGDEGIS_v_CTNDV70_gsea_comptable <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/integrative_analyses/run8/IMGDEGIS_v_CTNDV70/gsea_comp/gsea_comp_feature_plottable.csv", sep = ",", header = TRUE, row.names = 1)
IMGDEGIS_v_CTNDV70_psea_comptable <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/integrative_analyses/run8/IMGDEGIS_v_CTNDV70/psea_comp/psea_comp_feature_plottable.csv", sep = ",", header = TRUE, row.names = 1)
RNA_v_Meth_IMGDEGIS_comptable <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/integrative_analyses/run8/RNA_v_Meth/comp_ischemia__Sev_v_MildNone/gsa_comp/comp_ischemia__Sev_v_MildNone_gsacomp_plottable.csv", sep = ",", header = TRUE, row.names = 1)
RNA_v_Meth_CTNDV70_comptable <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/integrative_analyses/run8/RNA_v_Meth/comp_anatomy70__3v_v_1v/gsa_comp/comp_anatomy70__3v_v_1v_gsacomp_plottable.csv", sep = ",", header = TRUE, row.names = 1)

## overlap inlists
RNA_Meth_IMGDEGIS_CTNDV70_overlapinlists <- list(
    # sig IMGDEGIS only GSEA, # sig IMGDEGIS only PSEA, # UL corner in RNA vs Meth for IMGDEGIS
    overlapinlist1 = list(
        overlap_label = "IMGDEGIS_uniquesig_UL_quadrant",
        gsea_sig_IMGDEGIS = rownames(IMGDEGIS_v_CTNDV70_gsea_comptable[IMGDEGIS_v_CTNDV70_gsea_comptable[,"featgroup"] == "sig_IMGDEGIS",]),
        psea_sig_IMGDEGIS = rownames(IMGDEGIS_v_CTNDV70_psea_comptable[IMGDEGIS_v_CTNDV70_psea_comptable[,"featgroup"] == "sig_IMGDEGIS",]),
        RNA_v_Meth_IMGDEGIS_comptable = rownames(RNA_v_Meth_IMGDEGIS_comptable[RNA_v_Meth_IMGDEGIS_comptable[,"quadrant"] == "UL",])
    ),
    # sig CTNDV70 only GSEA, # sig CTNDV70 only PSEA, # UL corner in RNA vs Meth for CTNDV70
    overlapinlist2 = list(
        overlap_label = "CTNDV70_uniquesig_UL_quadrant",
        gsea_sig_CTNDV70 = rownames(IMGDEGIS_v_CTNDV70_gsea_comptable[IMGDEGIS_v_CTNDV70_gsea_comptable[,"featgroup"] == "sig_CTNDV70",]),
        psea_sig_CTNDV70 = rownames(IMGDEGIS_v_CTNDV70_psea_comptable[IMGDEGIS_v_CTNDV70_psea_comptable[,"featgroup"] == "sig_CTNDV70",]),
        RNA_v_Meth_CTNDV70_comptable = rownames(RNA_v_Meth_CTNDV70_comptable[RNA_v_Meth_CTNDV70_comptable[,"quadrant"] == "UL",])
    )
)

# inlist = RNA_Meth_IMGDEGIS_CTNDV70_overlapinlists[[1]]
for (inlist in RNA_Meth_IMGDEGIS_CTNDV70_overlapinlists) {
    suboutfolder <- paste0(outfilepathmaster, "IMGDEGIS_v_CTNDV70/", "overlap_pathway_lists/", inlist[["overlap_label"]], "/")
    dir.create(suboutfolder, showWarnings = FALSE, recursive = TRUE)

    overlapout <- overlap_finder(inlist[2:4])
    overlaptable <- overlapout$overlaptable
    vennplot <- overlapout$vennplot
    overlapgrouptab <- overlapout$overlapgrouptab
    overlapsummary <- overlapout$overlapsummary

    write.table(overlaptable, paste0(suboutfolder, "testoverlap_overlaptable.csv"), sep = ",", col.names = TRUE, row.names = FALSE)
    write.table(overlapgrouptab, paste0(suboutfolder, "testoverlap_overlapgrouptable.csv"), sep = ",", col.names = TRUE, row.names = FALSE)
    write.table(overlapsummary, paste0(suboutfolder, "testoverlap_overlapsummary.csv"), sep = ",", col.names = TRUE, row.names = TRUE)
    pdf(paste0(suboutfolder, "testoverlap_overlapvenn.pdf"))
    grid.draw(vennplot)
    junk <- dev.off()
}

# dir.create(paste0(compplot_outfilepath, "protect_v_risk/"), showWarnings = FALSE, recursive = TRUE)
# protect_genes <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run3_rmoutliers2/integrative_analyses/run6/subtype-combo/diff_comp_plots/rna_RS1_v_meth_MS2/UL_genes.csv", sep = ",")[,1]
# risk_genes <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run3_rmoutliers2/integrative_analyses/run6/subtype-combo/diff_comp_plots/rna_RS2_v_meth_MS3/UL_genes.csv", sep = ",")[,1]

# overlapinlist <- list(
#     protect_genes = protect_genes,
#     risk_genes = risk_genes
# )
# overlapout <- overlap_finder(overlapinlist)
# overlaptable <- overlapout$overlaptable
# vennplot <- overlapout$vennplot
# overlapgrouptab <- overlapout$overlapgrouptab
# colnames(overlapgrouptab)[1] <- "Row.names"
# overlapsummary <- overlapout$overlapsummary
# 
# dir.create(paste0(outfilepathmaster, "subtype-combo/", "protect_v_risk/", "gene_comp/"), showWarnings = FALSE, recursive = TRUE)
# write.table(overlaptable, paste0(outfilepathmaster, "subtype-combo/", "protect_v_risk/", "gene_comp/", "testoverlap_overlaptable.csv"),
#             sep = ",", col.names = TRUE, row.names = FALSE)
# write.table(overlapgrouptab, paste0(outfilepathmaster, "subtype-combo/", "protect_v_risk/", "gene_comp/", "testoverlap_overlapgrouptable.csv"),
#             sep = ",", col.names = TRUE, row.names = FALSE)
# write.table(overlapsummary, paste0(outfilepathmaster, "subtype-combo/", "protect_v_risk/", "gene_comp/", "testoverlap_overlapsummary.csv"),
#             sep = ",", col.names = TRUE, row.names = TRUE)
# pdf(paste0(outfilepathmaster, "subtype-combo/", "protect_v_risk/", "gene_comp/", "testoverlap_overlapvenn.pdf"))
# grid.draw(vennplot)
# junk <- dev.off()



# --------------------------------- RNAseq vs Meth DESeq v DMP and GSEA v PSEA --------------------------------- 
dir.create(paste0(outfilepathmaster, "RNA_v_Meth/"), showWarnings = FALSE, recursive = TRUE)

# Full pathway table
go_rna_genesettab = as.data.frame(msigdbr(species = "Homo sapiens", category = c("C5"))[,c("gs_name", "gene_symbol")])
go_meth_genesttab <- readRDS("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/meth_processing/PSEA_reference_tables/probe_set_C5_reftable.rds")
## Need this table for reference
DMP_table_in <- read.table(paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/meth_processing/run3/asr_corrected_DMP/",
                                  "DMP_ischemia__Sev_v_MildNone/DMP_analysis_and_figs/DMP_ischemia__Sev_v_MildNone_DMP_table.csv"), sep = ",", row.names = 1, header = TRUE)
methgenegrab <- unique(DMP_table_in[,"gene"])
rnagenes_w_probes_geneset_counttab <- aggregate(go_rna_genesettab[,"gene_symbol",drop=FALSE], 
                                                by = list(go_rna_genesettab[,"gs_name"]), function(x) {
                                                    length(unique(x)[unique(x) %in% methgenegrab])})
colnames(rnagenes_w_probes_geneset_counttab) <- c("gs_name", "genes_with_probe_count")

## Set analyses for comparison
Meth_path <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/meth_processing/run3/asr_corrected_DMP/"
RNAseq_path <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/"
analysis_comp_list <- list(
    run1 = list(analysis_label = "comp_ischemia__Sev_v_MildNone",
                DMP_table_path = paste0(Meth_path, "DMP_ischemia__Sev_v_MildNone/DMP_analysis_and_figs/DMP_ischemia__Sev_v_MildNone_DMP_table.csv"),
                DMP_cutoffs = c(meth_deltaBeta_cutoff = DMP_ischemia__Sev_v_MildNone__FCCUTOFF,
                                meth_pvalue_cutoff = DMP_ischemia__Sev_v_MildNone__PVALCUTOFF,
                                meth_pstattype = DMP_ischemia__Sev_v_MildNone__PVALTEST),
                RNA_table_path = paste0(RNAseq_path, "deseq/comp_ischemia__Sev_v_MildNone/deseq_results_comp_ischemia__Sev_v_MildNone.csv"),
                RNA_cutoffs = c(RNA_log2fc_cutoff = comp_ischemia__Sev_v_MildNone__FCCUTOFF,
                                RNA_pvalue_cutoff = comp_ischemia__Sev_v_MildNone__PVALCUTOFF,
                                RNA_pstattype = comp_ischemia__Sev_v_MildNone__PVALTEST),
                meth_annotation_label = "feature",
                DMP_PSEA_table_path = paste0(Meth_path, "PSEA/run2/c5/DMP_ischemia__Sev_v_MildNone/PSEA/DMP_ischemia__Sev_v_MildNone_PSEA_CUSTOM_GO.csv"),
                RNA_GSEA_table_path = paste0(RNAseq_path, "gsea/comp_ischemia__Sev_v_MildNone/comp_ischemia__Sev_v_MildNone_gsea_GO.csv"),
                GSOIvec = c("GOCC_IMMUNOGLOBULIN_COMPLEX", "GOCC_PLATELET_ALPHA_GRANULE", "GOBP_PLATELET_AGGREGATION",
                            "GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY", "GOBP_B_CELL_DIFFERENTIATION", "GOBP_REGULATION_OF_RESPONSE_TO_WOUNDING", "GOBP_RESPONSE_TO_INTERFERON_GAMMA")
    ),
    run2 = list(analysis_label = "comp_anatomy70__3v_v_1v",
                DMP_table_path = paste0(Meth_path, "DMP_anatomy70__3v_v_1v/DMP_analysis_and_figs/DMP_anatomy70__3v_v_1v_DMP_table.csv"), 
                DMP_cutoffs = c(meth_deltaBeta_cutoff = DMP_anatomy70__3v_v_1v__FCCUTOFF,
                                meth_pvalue_cutoff = DMP_anatomy70__3v_v_1v__PVALCUTOFF,
                                meth_pstattype = DMP_anatomy70__3v_v_1v__PVALTEST), 
                RNA_table_path = paste0(RNAseq_path, "deseq/comp_anatomy70__3v_v_1v/deseq_results_comp_anatomy70__3v_v_1v.csv"),
                RNA_cutoffs = c(RNA_log2fc_cutoff = comp_anatomy70__3v_v_1v__FCCUTOFF,
                                RNA_pvalue_cutoff = comp_anatomy70__3v_v_1v__PVALCUTOFF,
                                RNA_pstattype = comp_anatomy70__3v_v_1v__PVALTEST),
                meth_annotation_label = "feature",
                DMP_PSEA_table_path = paste0(Meth_path, "PSEA/run2/c5/DMP_anatomy70__3v_v_1v/PSEA/DMP_anatomy70__3v_v_1v_PSEA_CUSTOM_GO.csv"),
                RNA_GSEA_table_path = paste0(RNAseq_path, "gsea/comp_anatomy__3v_v_1v/comp_anatomy__3v_v_1v_gsea_GO.csv"),
                GSOIvec = c(
                    "GOMF_IMMUNOGLOBULIN_RECEPTOR_BINDING", "GOCC_MHC_CLASS_II_PROTEIN_COMPLEX", "GOCC_SPECIFIC_GRANULE_LUMEN", "GOBP_IMMUNE_EFFECTOR_PROCESS",
                            "GOBP_INNATE_IMMUNE_RESPONSE", "GOBP_LEUKOCYTE_MEDIATED_IMMUNITY", "GOCC_SPECIFIC_GRANULE_LUMEN",
                            "GOBP_ANTIGEN_RECEPTOR_MEDIATED_SIGNALING_PATHWAY", "GOBP_B_CELL_RECEPTOR_SIGNALING_PATHWAY")
    )
)

for (analysis_comp_num in seq_along(analysis_comp_list)) {
    # create sub folder out
    analysis_label_sel <- analysis_comp_list[[analysis_comp_num]][["analysis_label"]]
    suboutfolder <- paste0(outfilepathmaster, "RNA_v_Meth/", analysis_label_sel, "/")
    dir.create(paste0(suboutfolder, "diffexp_comp/"), showWarnings = FALSE, recursive = TRUE)
    
    # grab all of the parameters from the run we are running
    DMP_table_in <- read.table(analysis_comp_list[[analysis_comp_num]][["DMP_table_path"]], sep = ",", row.names = 1, header = TRUE)
    DMP_cutoffs_sel <- analysis_comp_list[[analysis_comp_num]][["DMP_cutoffs"]]
    RNA_table_in <- read.table(analysis_comp_list[[analysis_comp_num]][["RNA_table_path"]], sep = ",", row.names = 1, header = TRUE)
    RNA_cutoffs_sel <- analysis_comp_list[[analysis_comp_num]][["RNA_cutoffs"]]
    meth_annotation_label <- analysis_comp_list[[analysis_comp_num]][["meth_annotation_label"]]
    
    # First grab our selected DMP table
    DMP_sig_sel <- DMP_table_in[
        DMP_table_in[,DMP_cutoffs_sel["meth_pstattype"]] < as.numeric(DMP_cutoffs_sel["meth_pvalue_cutoff"]) &
            abs(DMP_table_in[,"deltaBeta"]) > as.numeric(DMP_cutoffs_sel["meth_deltaBeta_cutoff"]) & 
            !DMP_table_in[,"gene"] %in% c("NA", NA, ""), 
        c("gene", DMP_cutoffs_sel["meth_pstattype"], "deltaBeta", meth_annotation_label)]
    
    # Then grab selected RNA table
    DESEQ_sig_sel <- RNA_table_in[
        RNA_table_in[,RNA_cutoffs_sel["RNA_pstattype"]] < as.numeric(RNA_cutoffs_sel["RNA_pvalue_cutoff"]) &
            abs(RNA_table_in[,"log2FoldChange"]) > as.numeric(RNA_cutoffs_sel["RNA_log2fc_cutoff"]),
        c(RNA_cutoffs_sel["RNA_pstattype"], "log2FoldChange")]
    
    ## Process meth table for plotting
    feature_plottable <- merge(DMP_sig_sel, DESEQ_sig_sel, by.x = "gene", by.y = "row.names", all = TRUE, suffixes = c("_meth", "_rna"))
    feature_plottable[,c("deltaBeta", "log2FoldChange")][is.na(feature_plottable[,c("deltaBeta", "log2FoldChange")])] <- 0
    feature_plottable[,"gene"] <- make.unique(feature_plottable[,"gene"])
    row.names(feature_plottable) <- feature_plottable[,"gene"]
    feature_plottable[is.na(feature_plottable[,"feature"]),"feature"] <- "RNA"
    feature_plottable <- na.omit(feature_plottable)
    
    # Write out the table for checking later for convenience:
    grabstatcols <- c(paste0(DMP_cutoffs_sel["meth_pstattype"], "_meth"), paste0(RNA_cutoffs_sel["RNA_pstattype"], "_rna"))
    gene_outtable <- cbind(feature_plottable, 
                           quadrant = ifelse(rowSums(feature_plottable[,grabstatcols] < 0.05) == 2 & feature_plottable[,"deltaBeta"] > 0 & feature_plottable[,"log2FoldChange"] > 0 , "UR",
                                             ifelse(rowSums(feature_plottable[,grabstatcols] < 0.05) == 2 & feature_plottable[,"deltaBeta"] < 0 & feature_plottable[,"log2FoldChange"] > 0 , "UL",
                                                    ifelse(rowSums(feature_plottable[,grabstatcols] < 0.05) == 2 & feature_plottable[,"deltaBeta"] > 0 & feature_plottable[,"log2FoldChange"] < 0 , "LR",
                                                           ifelse(rowSums(feature_plottable[,grabstatcols] < 0.05) == 2 & feature_plottable[,"deltaBeta"] < 0 & feature_plottable[,"log2FoldChange"] < 0 , "LL",
                                                                  ""))))
    )
    write.table(gene_outtable, paste0(suboutfolder, "diffexp_comp/", analysis_label_sel, "_genecomp_plottable.csv"),
                sep = ",", col.names = TRUE, row.names = FALSE)
    
    # Use this to grab some genes to label per quadrant (5 max?)
    maxnumlabeledgenesquad <- 5
    datalabels <- do.call(rbind, unname(lapply(split(gene_outtable, gene_outtable$quadrant),
                                               function(subdf) subdf[sample(1:nrow(subdf), min(maxnumlabeledgenesquad, nrow(subdf))),]
    )))
    
    # Plot the gene comparison scatter
    pout <- scatter_plotter(indata = feature_plottable[,c("deltaBeta", "log2FoldChange")],
                            datalabels = rownames(datalabels),
                            labsparam = list(title = paste0(analysis_label_sel, " deltaBeta (meth) vs log2FoldChange (RNA)"), x = "deltaBeta (meth)", y = "log2FoldChange (RNA)")
    )
    pout <- pout + geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 0, linetype = 2)
    pdf(paste0(suboutfolder, "diffexp_comp/", analysis_label_sel, "_sig_gene_scatterplot.pdf"), 7, 7, useDingbats = FALSE)
    print(pout)
    junk <- dev.off()
    
    
    ## Hypergeo tests of each quadrant:
    speciesparam <- "Homo sapiens"
    ULgenes <- feature_plottable[feature_plottable[, "deltaBeta"] < 0 & feature_plottable[, "log2FoldChange"] > 0,"gene"]
    hypergeo_genetest_out_UL = data.frame(hypergeo_genetest(data.frame(ULgenes), statcutoffparam = "DUMMY", 
                                                            genesetparam = c("C5"), speciesparam = speciesparam)$enricherUPout)
    write.table(data.frame(ULgenes), paste0(suboutfolder, "diffexp_comp/", "UL_genes.csv"), sep = ",", col.names = FALSE, row.names = FALSE)
    write.table(hypergeo_genetest_out_UL, paste0(suboutfolder, "diffexp_comp/", "UL_hypergeo.csv"), sep = ",", col.names = TRUE, row.names = FALSE)
    
    URgenes <- feature_plottable[feature_plottable[, "deltaBeta"] > 0 & feature_plottable[, "log2FoldChange"] > 0,"gene"]
    hypergeo_genetest_out_UR = data.frame(hypergeo_genetest(data.frame(URgenes), statcutoffparam = "DUMMY", 
                                                            genesetparam = c("C5"), speciesparam = speciesparam)$enricherUPout)
    write.table(data.frame(URgenes), paste0(suboutfolder, "diffexp_comp/", "UR_genes.csv"), sep = ",", col.names = FALSE, row.names = FALSE)
    write.table(hypergeo_genetest_out_UR, paste0(suboutfolder, "diffexp_comp/", "UR_hypergeo.csv"), sep = ",", col.names = TRUE, row.names = FALSE)
    
    LLgenes <- feature_plottable[feature_plottable[, "deltaBeta"] < 0 & feature_plottable[, "log2FoldChange"] < 0,"gene"]
    hypergeo_genetest_out_LL = data.frame(hypergeo_genetest(data.frame(LLgenes), statcutoffparam = "DUMMY", 
                                                            genesetparam = c("C5"), speciesparam = speciesparam)$enricherUPout)
    write.table(data.frame(LLgenes), paste0(suboutfolder, "diffexp_comp/", "LL_genes.csv"), sep = ",", col.names = FALSE, row.names = FALSE)
    write.table(hypergeo_genetest_out_LL, paste0(suboutfolder, "diffexp_comp/", "LL_hypergeo.csv"), sep = ",", col.names = TRUE, row.names = FALSE)
    
    LRgenes <- feature_plottable[feature_plottable[, "deltaBeta"] > 0 & feature_plottable[, "log2FoldChange"] < 0,"gene"]
    hypergeo_genetest_out_LR = data.frame(hypergeo_genetest(data.frame(LRgenes), statcutoffparam = "DUMMY", 
                                                            genesetparam = c("C5"), speciesparam = speciesparam)$enricherUPout)
    write.table(data.frame(LRgenes), paste0(suboutfolder, "diffexp_comp/", "LR_genes.csv"), sep = ",", col.names = FALSE, row.names = FALSE)
    write.table(hypergeo_genetest_out_LR, paste0(suboutfolder, "diffexp_comp/", "LR_hypergeo.csv"), sep = ",", col.names = TRUE, row.names = FALSE)
    
    
    # Just do the GSEA vs PSEA - cleaner and easier
    dir.create(paste0(suboutfolder, "gsa_comp/"), showWarnings = FALSE, recursive = TRUE)
    meth_gsa_table <- read.table(analysis_comp_list[[analysis_comp_num]][["DMP_PSEA_table_path"]], sep = ",", header = TRUE, row.names = 1)
    rna_gsa_table <- read.table(analysis_comp_list[[analysis_comp_num]][["RNA_GSEA_table_path"]], sep = ",", header = TRUE, row.names = 1)
    
    gsa_combo_table <- merge(meth_gsa_table[,c("ID", "NES", "pvalue", "p.adjust", "setSize")],
                             rna_gsa_table[,c("ID", "NES", "pvalue", "p.adjust", "setSize")], by = "ID", suffixes = c("_meth", "_rna"))
    rownames(gsa_combo_table) <- gsa_combo_table[,"ID"]
    stattype <- "p.adjust"
    methmetric <- "NES_meth"
    rnametric <- "NES_rna"
    
    
    # Plot all
    plottable1 <- gsa_combo_table[,c(methmetric, rnametric)]
    pout1 <- scatter_plotter(indata = plottable1, labsparam = list(title = analysis_label_sel, x = methmetric, y = rnametric))
    pout1 <- pout1 + geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 0, linetype = 2)
    if (grepl("hypergeo", analysis_comp_list[[analysis_comp_num]][["DMP_PSEA_table_path"]])) {
        pout1 <- pout1 + geom_hline(yintercept = -log10(0.05), linetype = 2) + geom_vline(xintercept = -log10(0.05), linetype = 2)}
    pdf(paste0(suboutfolder, "gsa_comp/", analysis_label_sel, "_gsa_all_scatterplot", ".pdf"), 7, 7, useDingbats = FALSE)
    print(pout1)
    junk <- dev.off()
    
    # Plot only those that meet stat cutoff
    grabstatcols <- colnames(gsa_combo_table)[grepl(stattype, colnames(gsa_combo_table))]
    plottable2 <- gsa_combo_table[,c(methmetric, rnametric, grabstatcols)]
    plottable2[,"colorcol"] <- factor(ifelse(rowSums(plottable2[,grabstatcols] < 0.05) == 2, "black", "grey"), levels = c("grey", "black"))
    plottable2 <- plottable2[, c(methmetric, rnametric, "colorcol")]
    # pout2 <- scatter_plotter(indata = plottable2, labsparam = list(title = analysis_label_sel, x = methmetric, y = rnametric))
    pout2 <- scatter_plotter(indata = plottable2[,c(methmetric, rnametric)],
                             labsparam = list(title = analysis_label_sel, x = methmetric, y = rnametric), colorvar = plottable2[,"colorcol",drop=FALSE])
    pout2 <- pout2 + geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 0, linetype = 2)
    save_plot_limits <- unlist(get_plot_limits(pout1))
    pout2 <- pout2 + coord_cartesian(xlim = save_plot_limits[c("xmin", "xmax")], ylim = save_plot_limits[c("ymin", "ymax")], expand = FALSE)
    if (grepl("hypergeo", analysis_comp_list[[analysis_comp_num]][["DMP_PSEA_table_path"]])) {
        pout2 <- pout2 + geom_hline(yintercept = -log10(0.05), linetype = 2) + geom_vline(xintercept = -log10(0.05), linetype = 2)}
    pdf(paste0(suboutfolder, "gsa_comp/", analysis_label_sel, "_gsa_filtered_scatterplot", ".pdf"), 7, 7, useDingbats = FALSE)
    print(pout2)
    junk <- dev.off()
    
    # Write out the table for checking later for convenience:
    gs_outtable <- cbind(gsa_combo_table, 
                         quadrant = ifelse(rowSums(gsa_combo_table[,grabstatcols] < 0.05) == 2 & gsa_combo_table[,methmetric] > 0 & gsa_combo_table[,rnametric] > 0 , "UR",
                                           ifelse(rowSums(gsa_combo_table[,grabstatcols] < 0.05) == 2 & gsa_combo_table[,methmetric] < 0 & gsa_combo_table[,rnametric] > 0 , "UL",
                                                  ifelse(rowSums(gsa_combo_table[,grabstatcols] < 0.05) == 2 & gsa_combo_table[,methmetric] > 0 & gsa_combo_table[,rnametric] < 0 , "LR",
                                                         ifelse(rowSums(gsa_combo_table[,grabstatcols] < 0.05) == 2 & gsa_combo_table[,methmetric] < 0 & gsa_combo_table[,rnametric] < 0 , "LL",
                                                                ""))))
    )
    # Add number of genes and probes
    gs_outtable <- merge(gs_outtable, rnagenes_w_probes_geneset_counttab, by.x = "ID", by.y = "gs_name")
    write.table(gs_outtable, paste0(suboutfolder, "gsa_comp/", analysis_label_sel, "_gsacomp_plottable.csv"),
                sep = ",", col.names = TRUE, row.names = FALSE)
    
    
    # GSOI for plotting
    GSOIvec <- analysis_comp_list[[analysis_comp_num]][["GSOIvec"]]
    
    # Combine our diff exp tables
    fullfeature_plottable <- merge(cbind(DMP_table_in, probe = rownames(DMP_table_in)),
                                   RNA_table_in, by.x = "gene", by.y = "row.names", all = TRUE, suffixes = c("_meth", "_rna"))
    
    pdf(paste0(suboutfolder, "gsa_comp/", analysis_label_sel, "_GSOI_plots.pdf"), useDingbats = FALSE)
    GOI_sel_list <- list()
    for (GSOI in GSOIvec) {
        print(GSOI)
        # Grab the RNA genes and meth probes part of our GSOI
        grab_rna_genes <- unique(go_rna_genesettab[go_rna_genesettab[,"gs_name"] %in% GSOI, "gene_symbol"])[
            unique(go_rna_genesettab[go_rna_genesettab[,"gs_name"] %in% GSOI, "gene_symbol"]) %in% rownames(normcounttable)]
        grab_meth_probes <- unique(go_meth_genesttab[go_meth_genesttab[,"gs_name"] %in% GSOI, "gene_symbol"])[
            unique(go_meth_genesttab[go_meth_genesttab[,"gs_name"] %in% GSOI, "gene_symbol"]) %in% rownames(methcounttable)]
        grab_rna_counttable <- normcounttable[grab_rna_genes, intersect(colnames(normcounttable), colnames(methcounttable))]
        grab_meth_counttable <- methcounttable[grab_meth_probes, intersect(colnames(normcounttable), colnames(methcounttable))]
        
        # For each gene probe combo that is available - query our corr table and grab the corr value if we can
        GSOI_plottable <- data.frame(t(apply(fullfeature_plottable[fullfeature_plottable[,"gene"] %in% grab_rna_genes,
                                                                   c("gene", "probe", "deltaBeta", "log2FoldChange")], 1, function(x) {
                                                                       # Get individual gene correlations:
                                                                       if (!is.na(unname(x[["gene"]])) & !is.na(unname(x[["probe"]]))) {
                                                                           corr_out <- corr.test(t(grab_rna_counttable)[,unname(x[["gene"]])], t(grab_meth_counttable)[,unname(x[["probe"]])],
                                                                                                 method = "spearman", use = "pairwise.complete.obs")
                                                                           out1 <- c(x, gene_probe_rval = corr_out$r, gene_to_probe_pval = corr_out$p)
                                                                       } else {
                                                                           out1 <- c(x, gene_probe_rval = NA, gene_to_probe_pval = NA)
                                                                       }
                                                                       out1
                                                                   })))
        GSOI_plottable[,c("deltaBeta", "log2FoldChange", "gene_probe_rval", "gene_to_probe_pval")] <- apply(
            GSOI_plottable[,c("deltaBeta", "log2FoldChange", "gene_probe_rval", "gene_to_probe_pval")], 2, function(x) as.numeric(as.character(x)))
        
        
        ## Make the first plot
        pout <- scatter_plotter(GSOI_plottable[,c("deltaBeta", "log2FoldChange")], 
                                labsparam = list(title = paste0(GSOI, " feature plot"), x = "deltaBeta", y = "log2FoldChange"),
                                colorvar = GSOI_plottable[,"gene_probe_rval", drop = FALSE], sizevar = abs(GSOI_plottable[,"gene_probe_rval", drop = FALSE]))
        pout <- pout + scale_color_gradient2(low="darkblue", mid="grey", high="darkred", limits = c(-0.5, 0.5), oob = scales::squish)
        pout <- pout + scale_size_continuous(range = c(0,14), limits=c(0,0.8), breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8))
        pout <- pout + geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 0, linetype = 2)
        print(pout)
        
        ## Make a plot with a filtered corr value
        pout_corrfilt <- scatter_plotter(GSOI_plottable[GSOI_plottable[,"gene_to_probe_pval"] < 0.05, c("deltaBeta", "log2FoldChange")], 
                                         labsparam = list(title = paste0(GSOI, " feature plot CORR FILT"), x = "deltaBeta", y = "log2FoldChange"),
                                         colorvar = GSOI_plottable[GSOI_plottable[,"gene_to_probe_pval"] < 0.05, "gene_probe_rval", drop = FALSE],
                                         sizevar = abs(GSOI_plottable[GSOI_plottable[,"gene_to_probe_pval"] < 0.05, "gene_probe_rval", drop = FALSE]))
        pout_corrfilt <- pout_corrfilt + scale_color_gradient2(low="darkblue", mid="grey", high="darkred", limits = c(-0.5, 0.5), oob = scales::squish)
        pout_corrfilt <- pout_corrfilt + scale_size_continuous(range = c(0,10), limits=c(0,0.5), breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5))
        pout_corrfilt <- pout_corrfilt + geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 0, linetype = 2)
        print(pout_corrfilt)
        
        write.table(GSOI_plottable, paste0(suboutfolder, "gsa_comp/", analysis_label_sel, "__", GSOI, "_table.csv"), sep = ",", col.names = TRUE, row.names = FALSE)
        
        # GSOI plot table
        GSOI_plottable[,"quadrant"] <- ifelse(GSOI_plottable[,"deltaBeta"] < 0 & GSOI_plottable[,"log2FoldChange"] < 0, "LL", 
                                              ifelse(GSOI_plottable[,"deltaBeta"] < 0 & GSOI_plottable[,"log2FoldChange"] > 0, "UL", 
                                                     ifelse(GSOI_plottable[,"deltaBeta"] > 0 & GSOI_plottable[,"log2FoldChange"] < 0, "LR", 
                                                            ifelse(GSOI_plottable[,"deltaBeta"] > 0 & GSOI_plottable[,"log2FoldChange"] > 0, "UR", 
                                                                   NA))))
        GSOI_plottable[,"corrannot"] <- ifelse(GSOI_plottable[,"gene_probe_rval"] < 0 & GSOI_plottable[,"gene_to_probe_pval"] < 0.05, "signeg", 
                                               ifelse(GSOI_plottable[,"gene_probe_rval"] > 0 & GSOI_plottable[,"gene_to_probe_pval"] < 0.05, "sigpos", "notsig"))
        
        GOI_sel_list[[GSOI]] <- list(GSOI = GSOI, GOI = unique(GSOI_plottable[,"gene"]), GSOI_plottable = GSOI_plottable)
        
        # Helper table for this: counts of genes and correlation values and such
        GSOI_plottable[,"sigposcor"] <- ifelse(GSOI_plottable[,"gene_probe_rval"] > 0 & GSOI_plottable[,"gene_to_probe_pval"] < 0.05, 1, 0)
        GSOI_plottable[,"signegcor"] <- ifelse(GSOI_plottable[,"gene_probe_rval"] < 0 & GSOI_plottable[,"gene_to_probe_pval"] < 0.05, 1, 0)
        
        # Helper count table for GSOI
        GSOI_helper_counttable <- merge(data.frame(table(GSOI_plottable[GSOI_plottable[,"gene_to_probe_pval"] < 0.05, c("quadrant")])),
                                        aggregate(GSOI_plottable[,c("sigposcor", "signegcor")], by = list(GSOI_plottable[,"quadrant"]), sum),
                                        by.x = "Var1", by.y = "Group.1")
        colnames(GSOI_helper_counttable) <- c("quadrant", "sig_gene_probe_pairs", "sigposcor", "signegcor")
        write.table(GSOI_helper_counttable, paste0(suboutfolder, "gsa_comp/", analysis_label_sel, "__", GSOI, "_helper_counttable.csv"),
                    sep = ",", col.names = TRUE, row.names = FALSE)
        
    }
    junk <- dev.off()
    
    # Stripchart for some GOI - dont separately so we blow up the axes and fit all of the genes
    pdf(paste0(suboutfolder, "gsa_comp/", analysis_label_sel, "_GOI_for_GSOI_plots.pdf"), useDingbats = FALSE, 20, 20)
    for (GOI_sel in GOI_sel_list) {
        GSOI <- GOI_sel[["GSOI"]]
        GSOI_plottable <- GOI_sel[["GSOI_plottable"]]
        GOI <- GOI_sel[["GOI"]]
        GOI_plottable <- GSOI_plottable[GSOI_plottable[,"gene_to_probe_pval"] < 0.05 & GSOI_plottable[,"gene"] %in% unique(GSOI_plottable[,"gene"]),]
        pout_GOI <- scatter_plotter(GOI_plottable[, c("deltaBeta", "gene")],
                                    labsparam = list(title = paste0(GSOI, " feature plot CORR FILT"), x = "deltaBeta", y = "log2FoldChange"),
                                    colorvar = GOI_plottable[GOI_plottable[,"gene_to_probe_pval"] < 0.05, "gene_probe_rval", drop = FALSE],
                                    sizevar = abs(GOI_plottable[GOI_plottable[,"gene_to_probe_pval"] < 0.05, "gene_probe_rval", drop = FALSE]))
        pout_GOI <- pout_GOI + scale_color_gradient2(low="darkblue", mid="grey", high="darkred", limits = c(-0.5, 0.5), oob = scales::squish)
        pout_GOI <- pout_GOI + scale_size_continuous(range = c(0,10), limits=c(0,0.5), breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5))
        pout_GOI <- pout_GOI + theme(panel.grid.major.x = element_line(size = 1, linetype = 2))
        print(pout_GOI)
    }
    junk <- dev.off()
    
}




# --------------------------------- Redo ML Models for our diffexp features ---------------------------------
dir.create(paste0(outfilepathmaster, "diff_feat_signature_modeling/"), showWarnings = FALSE, recursive = TRUE)
## Need a params grid for train-test splits (0.4, 0.6, 0.75, 1.0)
rna_ischemia__Sev_v_MildNone_MLfeature_table <- read.table(paste0(outfilepathmaster, "/clean_dge_figures/dge_volcanoplots/comp_ischemia__Sev_v_MildNone_volcano_plottab.csv"), sep = ",", header = TRUE, row.names = 1)
rna_ischemia__Sev_v_MildNone_MLfeatures <- rownames(rna_ischemia__Sev_v_MildNone_MLfeature_table[rna_ischemia__Sev_v_MildNone_MLfeature_table[,"color"] %in% c("sigbothup", "sigbothdown"),])
rna_anatomy70__3v_v_1v_MLfeature_table <- read.table(paste0(outfilepathmaster, "/clean_dge_figures/dge_volcanoplots/comp_anatomy70__3v_v_1v_volcano_plottab.csv"), sep = ",", header = TRUE, row.names = 1)
rna_anatomy70__3v_v_1v_MLfeatures <- rownames(rna_anatomy70__3v_v_1v_MLfeature_table[rna_anatomy70__3v_v_1v_MLfeature_table[,"color"] %in% c("sigbothup", "sigbothdown"),])

meth_ischemia__Sev_v_MildNone_MLfeature_table <- read.table(paste0(outfilepathmaster, "/clean_dmp_figures/dmp_volcanoplots/DMP_ischemia__Sev_v_MildNone_volcano_plottab.csv"), sep = ",", header = TRUE, row.names = 1)
meth_ischemia__Sev_v_MildNone_MLfeatures <- meth_ischemia__Sev_v_MildNone_MLfeature_table[meth_ischemia__Sev_v_MildNone_MLfeature_table[,"color"] %in% c("sigbothup", "sigbothdown"), "cg_probe"]
meth_anatomy70__3v_v_1v_MLfeature_table <- read.table(paste0(outfilepathmaster, "/clean_dmp_figures/dmp_volcanoplots/DMP_anatomy70__3v_v_1v_volcano_plottab.csv"), sep = ",", header = TRUE, row.names = 1)
meth_anatomy70__3v_v_1v_MLfeatures <- meth_anatomy70__3v_v_1v_MLfeature_table[meth_anatomy70__3v_v_1v_MLfeature_table[,"color"] %in% c("sigbothup", "sigbothdown"), "cg_probe"]



# I think I need a custom ML metatable - only cause i need it for the meth data really...
# RNA metatable
metatable_file <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/rna_processing/metatable_in.csv"
metatable <- read.table(metatable_file, sep = ",", header = TRUE, row.names = 1)
rna_ML_metatable <- cbind(Sample_Patnum = rownames(metatable), metatable[,c("comp_ischemia__Sev_v_MildNone", "comp_anatomy70__3v_v_1v")])
rna_ML_metatable[,c("comp_ischemia__Sev_v_MildNone")] <- ifelse(rna_ML_metatable[,c("comp_ischemia__Sev_v_MildNone")] == 1, "Severe",
                                                         ifelse(rna_ML_metatable[,c("comp_ischemia__Sev_v_MildNone")] == 0, "MildNone", NA))
rna_ML_metatable[,c("comp_anatomy70__3v_v_1v")] <- ifelse(rna_ML_metatable[,c("comp_anatomy70__3v_v_1v")] == 1, "3v",
                                                   ifelse(rna_ML_metatable[,c("comp_anatomy70__3v_v_1v")] == 0, "1v", NA))
colnames(rna_ML_metatable) <- c("PATNUM", "rna_comp_ischemia__Sev_v_MildNone", "rna_comp_anatomy70__3v_v_1v")
# Meth metatable
meth_imgdegis_metatable_temp <- readRDS("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/final_meth_data_v1_20230207/DMP_with_asr_correction/metadata_DMP_on_IMGDEGIS_control.rds")[,c("Sample_Patnum", "IMGDEGIS")]
meth_imgdegis_metatable_temp[,"meth_comp_ischemia__Sev_v_MildNone"] <- ifelse(meth_imgdegis_metatable_temp[,"IMGDEGIS"] == "Severe", "Severe", 
                                                                              ifelse(meth_imgdegis_metatable_temp[,"IMGDEGIS"] == "No_mild", "No_mild", NA))
meth_ctndv70_metatable_temp <- readRDS("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/final_meth_data_v1_20230207/DMP_with_asr_correction/metadata_DMP_on_CTNDV70_control.rds")[,c("Sample_Patnum", "CTNDV70")]
meth_ctndv70_metatable_temp[,"meth_comp_anatomy70__3v_v_1v"] <- ifelse(meth_ctndv70_metatable_temp[,"CTNDV70"] == "3", "3v",
                                                                       ifelse(meth_ctndv70_metatable_temp[,"CTNDV70"] == "1", "1v", NA))
meth_ML_metatable <- merge(meth_imgdegis_metatable_temp[,c("Sample_Patnum", "meth_comp_ischemia__Sev_v_MildNone")],
                           meth_ctndv70_metatable_temp[,c("Sample_Patnum", "meth_comp_anatomy70__3v_v_1v")], by = "Sample_Patnum", all = TRUE)
## Remove "048003-005" as a duplicate in this metatable as well
meth_ML_metatable <- meth_ML_metatable[!meth_ML_metatable[,"Sample_Patnum"] %in% "048003-005",]

# Combine for ML
ML_metatable <- merge(rna_ML_metatable, meth_ML_metatable, by.x = "PATNUM", by.y = "Sample_Patnum", all = TRUE)
rownames(ML_metatable) <- ML_metatable[,"PATNUM"]
ML_metatable <- ML_metatable[,!colnames(ML_metatable) %in% "PATNUM"]

ML_feature_list <- list(
    rna_ischemia__Sev_v_MildNone_MLfeatures = rna_ischemia__Sev_v_MildNone_MLfeatures,
    rna_anatomy70__3v_v_1v_MLfeatures = rna_anatomy70__3v_v_1v_MLfeatures,
    meth_ischemia__Sev_v_MildNone_MLfeatures = meth_ischemia__Sev_v_MildNone_MLfeatures,
    meth_anatomy70__3v_v_1v_MLfeatures = meth_anatomy70__3v_v_1v_MLfeatures
)

ML_params_grid <- expand.grid(ML_features = names(ML_feature_list), trainsplit = c(0.4, 0.6, 0.75, 1.0))
ML_param_inlist <- split(ML_params_grid, seq(nrow(ML_params_grid)))

for (MLrun in ML_param_inlist) {
# for (MLrun in ML_param_inlist[14:16]) {
    ML_features = MLrun[["ML_features"]]
    ML_trainsplit = MLrun[["trainsplit"]]
    complabel <- gsub("rna_|meth_|_MLfeatures", "", ML_features)
    suboutfolder <- paste0(outfilepathmaster, "diff_feat_signature_modeling/", ML_features, "/", paste0("trainsplit", gsub("\\.", "", ML_trainsplit)), "/")
    dir.create(suboutfolder, showWarnings = FALSE, recursive = TRUE)
    
    # Grab feature table and outcome table
    if (grepl("rna", ML_features)) { ## RNA run
        # First derive outcome table
        comptable <- ML_metatable[,grepl("rna_", colnames(ML_metatable)) & grepl(complabel, colnames(ML_metatable)),drop=FALSE]
        outcometable <- na.omit(comptable)
        outcomelevels <- unique(outcometable[,1])
        # Then derive featuretable
        grab_features <- ML_feature_list[[ML_features]]
        featuretable <- t(normcounttable[grab_features, rownames(outcometable)])
    }
    if (grepl("meth", ML_features)) { ## RNA run
        # First derive outcome table
        comptable <- ML_metatable[,grepl("meth_", colnames(ML_metatable)) & grepl(complabel, colnames(ML_metatable)),drop=FALSE]
        outcometable <- na.omit(comptable)
        outcomelevels <- unique(outcometable[,1])
        # Then derive featuretable
        grab_features <- ML_feature_list[[ML_features]]
        featuretable <- t(methcounttable[grab_features, rownames(outcometable)])
    }
    
    # Run modeling
    cohort_multi_classification_ml_analysis_out<- cohort_multi_classification_ml_analysis(
        featuretable = featuretable, outcometable = outcometable, outcomelevels = outcomelevels, seedparam=11111,
        presplitdata = NULL, train_partition_percent = ML_trainsplit, OHE_featuretable = FALSE,
        subsampleparam = "down", models_to_run = c("glm", "xgbTree", "svmRadial", "glmnet"))
        # subsampleparam = "down", models_to_run = c("glm", "glmnet"))
    print("modeling done")
    saveRDS(cohort_multi_classification_ml_analysis_out, paste0(suboutfolder, "cohort_multi_classification_ml_analysis_out.rds"))
    # cohort_multi_classification_ml_analysis_out <- readRDS(paste0(suboutfolder, "cohort_multi_classification_ml_analysis_out.rds"))
    
    # Write out all of the things from this run
    for (outtablenum in seq_along(cohort_multi_classification_ml_analysis_out)) {
        filename <- names(cohort_multi_classification_ml_analysis_out)[outtablenum]
        if (filename %in% (c("training", "testing", "train_modelstats_outtable", "train_overallstat_outtable",
                             "test_modelstats_outtable", "test_overallstat_outtable", "train_roc_auc_with_ci_res", "test_roc_auc_with_ci_res", "varImp_outtable"))) {
            write.table(cohort_multi_classification_ml_analysis_out[[outtablenum]], paste0(suboutfolder, filename, ".csv"), 
                        sep = ",", col.names = TRUE, row.names = FALSE)
        } else if (filename %in% c("train_rocplot", "test_rocplot")) {
            pdf(paste0(suboutfolder, filename, ".pdf"), useDingbats = FALSE)
            print(cohort_multi_classification_ml_analysis_out[[filename]])
            junk <- dev.off()
        } else {
            saveRDS(cohort_multi_classification_ml_analysis_out, paste0(suboutfolder, filename, ".rds"))
        }
    }
    
    
    
    # best_model_stat
    MOI <- c("Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value", "Precision", "Recall", "F1", 
             "Prevalence", "Detection Rate", "Detection Prevalence", "Balanced Accuracy")
    best_model_stat <- "AUC"
    ## TEST
    test_modelstats_outtable <- cohort_multi_classification_ml_analysis_out[["test_modelstats_outtable"]]
    test_roc_auc_with_ci_res <- cohort_multi_classification_ml_analysis_out[["test_roc_auc_with_ci_res"]]
    
    ## Orders should be the same, but I should double check...
    test_modelstats_clean <- data.frame(test_modelstats_outtable[,2:ncol(test_modelstats_outtable)], row.names = test_modelstats_outtable[,"metric"])
    colnames(test_modelstats_clean) <- paste0(c("c1_", "c2_"), unlist(lapply(strsplit(colnames(test_modelstats_clean), split = "__"), 
                                                                             function(x) x[length(x)])))
    test_auc_stats_clean_temp <- data.frame(t(data.frame(test_roc_auc_with_ci_res, 
                                                         row.names = unlist(lapply(strsplit(rownames(test_roc_auc_with_ci_res), split = "_"), function(x) x[1])))[,c(2:4)]), row.names = c("AUC", "AUClowerCI", "AUCupperCI"))
    test_auc_stats_clean <- test_auc_stats_clean_temp[,rep(1:ncol(test_auc_stats_clean_temp), each = 2)]
    colnames(test_auc_stats_clean) <- apply(expand.grid(c("c1_", "c2_"), colnames(test_auc_stats_clean_temp)), 1, paste, collapse="")
    test_allstat_table <- rbind(test_modelstats_clean, test_auc_stats_clean)
    test_allstat_table[] <- apply(test_allstat_table, 2, as.numeric)
    
    # Now with the combined table - find the model with the best stat
    maxmodelstat <- which(test_allstat_table[best_model_stat, ] == max(test_allstat_table[best_model_stat, ]))
    # maxmodelstat <- which.max(test_allstat_table[best_model_stat, ])
    if (length(maxmodelstat) > 1) {  ## If we get a max stat with more than one model, always default to AUC > balanced acc > random
        maxmodelstat <- tiebreak_maxstat_1 <- maxmodelstat[which(test_allstat_table["AUC", maxmodelstat]== max(test_allstat_table["AUC", maxmodelstat]))]
        if (length(tiebreak_maxstat_1) > 1) {
            maxmodelstat <- tiebreak_maxstat_2 <- maxmodelstat[which(test_allstat_table["Balanced Accuracy", maxmodelstat]== max(test_allstat_table["Balanced Accuracy", maxmodelstat]))]
            if (length(tiebreak_maxstat_2) > 1) {
                maxmodelstat <- tiebreak_maxstat_3 <- sample(tiebreak_maxstat_2, size = 1)
            }
        }
    }
    # maxmodelstat
    bestmodel_label <- gsub("c1_|c2_", "", colnames(test_allstat_table)[maxmodelstat])
    
    ## TRAIN
    train_modelstats_outtable <- cohort_multi_classification_ml_analysis_out[["train_modelstats_outtable"]]
    train_roc_auc_with_ci_res <- cohort_multi_classification_ml_analysis_out[["train_roc_auc_with_ci_res"]]
    train_modelstats_clean <- data.frame(train_modelstats_outtable[,2:ncol(train_modelstats_outtable)], row.names = train_modelstats_outtable[,"metric"])
    colnames(train_modelstats_clean) <- paste0(c("c1_", "c2_"), unlist(lapply(strsplit(colnames(train_modelstats_clean), split = "__"), 
                                                                              function(x) x[length(x)])))
    train_auc_stats_clean_temp <- data.frame(t(data.frame(train_roc_auc_with_ci_res, 
                                                          row.names = unlist(lapply(strsplit(rownames(train_roc_auc_with_ci_res), split = "_"), function(x) x[1])))[,c(2:4)]), row.names = c("AUC", "AUClowerCI", "AUCupperCI"))
    train_auc_stats_clean <- train_auc_stats_clean_temp[,rep(1:ncol(train_auc_stats_clean_temp), each = 2)]
    colnames(train_auc_stats_clean) <- apply(expand.grid(c("c1_", "c2_"), colnames(train_auc_stats_clean_temp)), 1, paste, collapse="")
    train_allstat_table <- rbind(train_modelstats_clean, train_auc_stats_clean)
    train_allstat_table[] <- apply(train_allstat_table, 2, as.numeric)
    
    
    ## COMBINE
    savestat_table <- cbind(test_allstat_table[,grepl(paste0(bestmodel_label, "$"), colnames(test_allstat_table))],
                            train_allstat_table[,grepl(paste0(bestmodel_label, "$"), colnames(train_allstat_table))])
    colnames(savestat_table) <- paste0(c(rep("test_", 2), rep("train_", 2)), gsub("c1_|c2_", "", colnames(savestat_table)))
    write.table(savestat_table, paste0(suboutfolder, complabel, "_bestmodel_stattable.csv"), sep = ",", col.names = NA, row.names = TRUE)
    
    # Also output the bestmodel ROC curves in a cleaner way:
    bestmodel_out <- cohort_multi_classification_ml_analysis_out[["modellist"]][[names(cohort_multi_classification_ml_analysis_out[["modellist"]])[grepl(paste0(bestmodel_label, "_"), names(cohort_multi_classification_ml_analysis_out[["modellist"]]))]]]
    training = cohort_multi_classification_ml_analysis_out[["training"]]
    testing = cohort_multi_classification_ml_analysis_out[["testing"]]
    
    # Replot the AUC curve
    ROC_model_data_plot_list <- list(combo1 = list(ROC_model = bestmodel_out, ROC_data = testing, plot_number = 1, 
                                                   ROC_name = paste0(bestmodel_label, "_test")),
                                     combo2 = list(ROC_model = bestmodel_out, ROC_data = training, plot_number = 1,
                                                   ROC_name = paste0(bestmodel_label, "_train")))
    ROC_custom_out <- ROC_custom_plotter(ROC_model_data_plot_list = ROC_model_data_plot_list, outcome_label = colnames(outcometable), plot_ci = TRUE)
    pdf(paste0(suboutfolder, complabel, "_bestmodel_ROCplots.pdf"), 7, 7, useDingbats = FALSE)
    # lapply(ROC_custom_out[["ROC_plots"]], print)
    for (grab_ROC_plot in ROC_custom_out[["ROC_plots"]]) {print(grab_ROC_plot)}
    junk <- dev.off()
    write.table(do.call(rbind, ROC_custom_out[["ROC_AUC_tables"]]), paste0(suboutfolder, complabel, "_bestmodel_AUCtable.csv"),
                sep = ",", col.names = NA, row.names = TRUE)

}

# Summarize all runs by splits:
summary_outfolder <- paste0(outfilepathmaster, "diff_feat_signature_modeling/", "modeling_runs_summary/")
dir.create(summary_outfolder, recursive = TRUE, showWarnings = FALSE)
summary_model_stats <- c("Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value", "Precision", "Recall", "F1", 
                         "Prevalence", "Detection Rate", "Detection Prevalence", "Balanced Accuracy", "AUC", "AUClowerCI", "AUCupperCI")
best_model_stat <- "AUC"
allstat_table_list <- list()
for (MLrun in ML_param_inlist) {
    ML_features = MLrun[["ML_features"]]
    ML_trainsplit = MLrun[["trainsplit"]]
    complabel <- gsub("_MLfeatures", "", ML_features)
    dir.create(suboutfolder, showWarnings = FALSE, recursive = TRUE)
    subinfolder <- paste0(outfilepathmaster, "diff_feat_signature_modeling/", ML_features, "/", paste0("trainsplit", gsub("\\.", "", ML_trainsplit)), "/")
    
    # Grab the model and test results for ROC plotting based off of AUC
    modeling_output_sel <- readRDS(paste0(subinfolder, "cohort_multi_classification_ml_analysis_out.rds"))
    test_modelstats_outtable <- modeling_output_sel[["test_modelstats_outtable"]]
    test_roc_auc_with_ci_res <- modeling_output_sel[["test_roc_auc_with_ci_res"]]
    allstat_table_temp1 <- combine_mgc_modelstats_and_auc_table(test_modelstats_outtable, test_roc_auc_with_ci_res)[summary_model_stats,]
    allstat_table_temp2 <- data.frame(t(allstat_table_temp1[!duplicated(t(allstat_table_temp1))]))
    allstat_table <- cbind(complabel = complabel, ML_trainsplit = ML_trainsplit, model = gsub("c1_", "", rownames(allstat_table_temp2)), allstat_table_temp2)
    allstat_table[which.max(allstat_table[,best_model_stat]),"bestmodel_help"] <- "best"
    allstat_table_list[[paste0(MLrun[1,1], MLrun[1,2])]] <- allstat_table
}

allstat_table <- do.call(rbind, allstat_table_list)
write.table(allstat_table, paste0(outfilepathmaster, "diff_feat_signature_modeling/", "modeling_runs_summary/", "allstat_summary_table.csv"),
            sep = ",", col.names = TRUE, row.names = FALSE)




# --------------------------------- ML Modeling followup and plotting --------------------------------- 
# Now that we have done everything - lets go back and grab from our models and runs of interest what we want
## Which for now - is the glmnet test ROC curve at a 0.6 split for each model
dir.create(paste0(outfilepathmaster, "diff_feat_signature_modeling/post_modeling_results/"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(outfilepathmaster, "diff_feat_signature_modeling/post_modeling_results/", "glmnet_test_ROC_plots/"), showWarnings = FALSE, recursive = TRUE)
modeling_object_list <-Sys.glob(paste0(outfilepathmaster, "diff_feat_signature_modeling/", "*", "/trainsplit06/cohort_multi_classification_ml_analysis_out.rds"))

MLrun_sel <- ML_param_inlist[unlist(lapply(ML_param_inlist, function(x) x[["trainsplit"]] == 0.6))]

for (MLrun in MLrun_sel) {
    ML_features = MLrun[["ML_features"]]
    ML_trainsplit = MLrun[["trainsplit"]]
    complabel <- gsub("_MLfeatures", "", ML_features)
    subinfolder <- paste0(outfilepathmaster, "diff_feat_signature_modeling/", ML_features, "/", paste0("trainsplit", gsub("\\.", "", ML_trainsplit)), "/")
    
    # Grab the things we need
    modeling_output_sel <- readRDS(paste0(subinfolder, "cohort_multi_classification_ml_analysis_out.rds"))
    testing <- modeling_output_sel[["testing"]]
    training <- modeling_output_sel[["training"]]
    test_modelstats_outtable <- modeling_output_sel[["test_modelstats_outtable"]]
    
    # Grab the model and test results for ROC plotting based off of AUC
    suboutfolder <- paste0(outfilepathmaster, "diff_feat_signature_modeling/post_modeling_results/", "glmnet_test_ROC_plots/")
    dir.create(suboutfolder, showWarnings = FALSE, recursive = TRUE)
    glmnet_model <- modeling_output_sel[["modellist"]][["glmnet_modfit"]] ## "xgbTree_modfit"
    ROC_model_data_plot_list <- list(combo1 = list(ROC_model = glmnet_model, ROC_data = testing, plot_number = 1, 
                                                   ROC_name = paste0("glmnet_test_", complabel)))
    ROC_custom_out <- ROC_custom_plotter(ROC_model_data_plot_list = ROC_model_data_plot_list, outcome_label = colnames(testing)[1], plot_ci = TRUE)
    pdf(paste0(suboutfolder, complabel, "_glmnet_test_ROCplot.pdf"), 7, 7, useDingbats = FALSE)
    for (grab_ROC_plot in ROC_custom_out[["ROC_plots"]]) {print(grab_ROC_plot)}
    junk <- dev.off()

    glmnet_cm <- test_modelstats_outtable[grepl("predicted_class", test_modelstats_outtable[,"metric"]),
                                          c("metric", colnames(test_modelstats_outtable)[grepl("glmnet", colnames(test_modelstats_outtable))])]
    write.table(glmnet_cm, paste0(suboutfolder, complabel, "_glmnet_cmtable.csv"), sep = ",", col.names = TRUE, row.names = FALSE)
    write.table(rbind(training = table(training[,1]), testing = table(testing[,1]), total = table(rbind(training, testing)[,1])),
                paste0(suboutfolder, complabel, "_glmnet_outcome_totals.csv"), sep = ",", col.names = TRUE, row.names = TRUE)
    
    
    # Grab the model and test results for ROC plotting based off of AUC
    suboutfolder <- paste0(outfilepathmaster, "diff_feat_signature_modeling/post_modeling_results/", "xgbTree_test_ROC_plots/")
    dir.create(suboutfolder, showWarnings = FALSE, recursive = TRUE)
    xgbTree_model <- modeling_output_sel[["modellist"]][["xgbTree_modfit"]]
    ROC_model_data_plot_list <- list(combo1 = list(ROC_model = xgbTree_model, ROC_data = testing, plot_number = 1, 
                                                   ROC_name = paste0("xgbTree_test_", complabel)))
    ROC_custom_out <- ROC_custom_plotter(ROC_model_data_plot_list = ROC_model_data_plot_list, outcome_label = colnames(testing)[1], plot_ci = TRUE)
    pdf(paste0(suboutfolder, complabel, "_xgbTree_test_ROCplot.pdf"), 7, 7, useDingbats = FALSE)
    for (grab_ROC_plot in ROC_custom_out[["ROC_plots"]]) {print(grab_ROC_plot)}
    junk <- dev.off()

    xgbTree_cm <- test_modelstats_outtable[grepl("predicted_class", test_modelstats_outtable[,"metric"]),
                                          c("metric", colnames(test_modelstats_outtable)[grepl("xgbTree", colnames(test_modelstats_outtable))])]
    write.table(xgbTree_cm, paste0(suboutfolder, complabel, "_xgbTree_cmtable.csv"), sep = ",", col.names = TRUE, row.names = FALSE)
    write.table(rbind(training = table(training[,1]), testing = table(testing[,1]), total = table(rbind(training, testing)[,1])),
                paste0(suboutfolder, complabel, "_xgbTree_outcome_totals.csv"), sep = ",", col.names = TRUE, row.names = TRUE)
    
}


# --------------------------------- DANGER-LONG-RUN-TIME cluster_eigen_meta_summary_heatmap --------------------------------- 
methcounttable_file <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/final_meth_data_v1_20230207/limma_noCovariate.rds"
methcounttable <- readRDS(methcounttable_file)
colnames(methcounttable) <- methmetatable[match(methmetatable[,"Sample_Name"], colnames(methcounttable)), "Sample_Patnum"]
methcounttable <- methcounttable[,-578] # Need to remove the duplicate of this sample: "048003-005", columns 204 and 578

meth_eigengenecounttable <- read.table(paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/integrative_analyses/run7/", "meth_WGCNA/", "meth_eigengene_nogrey_table.csv"),sep = ",", header = TRUE, row.names = 1)
meth_wgcna_genestocolors <- read.table(paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/integrative_analyses/run7/", "meth_WGCNA/", "meth_wgcna_genestocolors.csv"), sep = ",", header = TRUE, row.names = 1)


# I manually curated and created module annotation tables
rna_wgcna_chosen_annotation_table <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/WGCNA/rna_wgcna_notation_table.csv",
                                                sep = ",", header = TRUE)
rna_wgcna_full_C5_annotation_table <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/WGCNA/WGCNA_power14_size30/C5_wgcna_module_annotations.csv",
                                                 sep = ",", header = TRUE, row.names = 1)

meth_wgcna_chosen_annotation_table_temp <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/integrative_analyses/manual_annotation_tables/meth_WGCNA_custom_selected_pathways_v1.csv", sep = ",", header = TRUE)
meth_wgcna_chosen_annotation_table <- data.frame(cbind(eigengene = meth_wgcna_chosen_annotation_table_temp[,"module"], goterm = meth_wgcna_chosen_annotation_table_temp[,"ID"]))
meth_wgcna_full_C5_annotation_table <- read.table(paste0(outfilepathmaster, "meth_WGCNA/", "GO_meth_wgcna_module_annotations.csv"), sep = ",", header = TRUE, row.names = 1)


## Summary plots for each clustering
group_columns_for_testing <- c("rna_4cluster_w3AB", "rna_4cluster", "meth_3cluster")
for (group_column_selected in group_columns_for_testing) {
    # ordered_combined_clustermembership_table <- combined_clustermembership_table[order(combined_clustermembership_table[,group_column_selected]),]
    dir.create(paste0(outfilepathmaster, "nmf_wgcna_integration/", group_column_selected, "/"), showWarnings = FALSE, recursive = TRUE)
    
    if (grepl("rna", group_column_selected)) {
        incounttable <- normcounttable
        ingenestocolortable <- genestocolorstab
        eigengeneorder <- c("magenta", "salmon", "purple", "midnightblue", "turquoise", 
                                     "yellow", "pink", "red", 
                                     "blue", "tan", "brown",
                                     "greenyellow", "black", "cyan", "green")
        ineigencounttable <- eigengenecounttable
                                     
        module_chosen_annotation_table <- rna_wgcna_chosen_annotation_table
        module_full_annotation_table <- rna_wgcna_full_C5_annotation_table   
                                     
                                     
    } else {
        incounttable <- methcounttable
        ingenestocolortable <- meth_wgcna_genestocolors
        eigengeneorder <- c("yellow", "cyan", "tan", "turquoise", "black", 
                            "magenta", "pink", "brown",
                            "blue", "red", "salmon", "green", "purple")
        ## 2023-03-28 - removing the "dark" ones - too small and no annotation
        # eigengeneorder <- c("yellow", "cyan", "darkgreen", "tan", "turquoise", "black", 
        #                     "magenta", "pink", "brown", "darkturquoise",
        #                     "blue", "red", "salmon", "green", "purple")
        genesselect <- rownames(ingenestocolortable[ingenestocolortable[,"moduleColors"] %in% eigengeneorder,])

        # Downsample to only 30% of probes - better for plotting
        downsample_gene_selection <- cbind(probe = rownames(ingenestocolortable[ingenestocolortable[,"moduleColors"] %in% eigengeneorder,]), 
                                           ingenestocolortable[ingenestocolortable[,"moduleColors"] %in% eigengeneorder,])
        set.seed(1234)
        genesselect <- do.call(rbind, lapply(split(downsample_gene_selection, downsample_gene_selection[,"moduleColors"]), 
                                           function(x) x[sample(nrow(x), floor(0.3 * nrow(x))),]))[,"probe"]
        incounttable <- incounttable[genesselect,]
        ineigencounttable <- meth_eigengenecounttable
                                 
        module_chosen_annotation_table <- meth_wgcna_chosen_annotation_table
        module_full_annotation_table <- meth_wgcna_full_C5_annotation_table   
        
        # Addinig in this analysis here - what KIND of probes are in each of our modules?
        dir.create(paste0(outfilepathmaster, "nmf_wgcna_integration/", group_column_selected, "/", "probe_feature_label_plots/"), showWarnings = FALSE, recursive = TRUE)
        ## Need this table for reference
        DMP_table_in <- read.table(paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/meth_processing/run3/asr_corrected_DMP/",
                                          "DMP_ischemia__Sev_v_MildNone/DMP_analysis_and_figs/DMP_ischemia__Sev_v_MildNone_DMP_table.csv"), sep = ",", row.names = 1, header = TRUE)
        probe_feature_labels <- c("feature", "cgi", "feat.cgi", "CHR")
        module_probetype_table <- cbind(ingenestocolortable[,"moduleColors",drop=FALSE], DMP_table_in[rownames(ingenestocolortable), probe_feature_labels])
        module_probetype_table <- module_probetype_table[!module_probetype_table[,"moduleColors"] %in% "grey",]
        # tt1 <- aggregate(module_probetype_table[,c("feature", "cgi")], by = list(module_probetype_table[,"moduleColors"]), table)
        
        probe_feature_label_params_grid <- expand.grid(probe_feature_labels, c("allmodules", "selectedmodules", "cluster"))
        colnames(probe_feature_label_params_grid) <- c("probe_feature_label", "grouping_label")
        for (feature_label_grouping_sel in seq_len(nrow(probe_feature_label_params_grid))) {
            # Select params
            probe_feature_label <- as.character(probe_feature_label_params_grid[feature_label_grouping_sel, "probe_feature_label"])
            grouping_label <- as.character(probe_feature_label_params_grid[feature_label_grouping_sel, "grouping_label"])
            # Create feature tabulation table
            featuretype_plottable <- do.call(rbind, lapply(split(module_probetype_table, module_probetype_table[,"moduleColors"]), function(x) {
                feature_type_table <- data.frame(table(factor(x[,probe_feature_label], 
                                                              levels = unique(module_probetype_table[,probe_feature_label]))), descriptor = probe_feature_label)
                feature_type_table[,"feature_percent"] <- round(feature_type_table[,2] / sum(feature_type_table[,2]), 3)
                feature_type_table
            }))
            if (probe_feature_label == "CHR") {featuretype_plottable[,"Var1"] <- factor(featuretype_plottable[,"Var1"], seq(1,22))}
            featuretype_plottable[,"module"] <- gsub("\\.[^.]*$", "", rownames(featuretype_plottable))
            featuretype_plottable[,"cluster"] <- ifelse(featuretype_plottable[,"module"] %in% c("yellow", "cyan", "darkgreen", "tan", "turquoise", "black"), "MS1",
                                                 ifelse(featuretype_plottable[,"module"] %in% c("magenta", "pink", "brown", "darkturquoise"), "MS2", 
                                                 ifelse(featuretype_plottable[,"module"] %in% c("blue", "red", "salmon", "green", "purple"), "MS3", NA)))
            
            # Here we select which plot we are making based on the grouping:
            if (grouping_label %in% c("selectedmodules", "cluster")) {
                # select only our modules of interest
                featuretype_plottable <- featuretype_plottable[featuretype_plottable[,"module"] %in% eigengeneorder,]
                featuretype_plottable[,"module"] <- factor(featuretype_plottable[,"module"], levels = eigengeneorder)
                if (grouping_label %in% "cluster") {
                    featuretype_plottable[,"module"] <- factor(featuretype_plottable[,"cluster"], levels = c("MS1", "MS2", "MS3"))
                    featuretype_plottable <- aggregate(featuretype_plottable[,c("Freq", "feature_percent"),drop=FALSE],
                                                       by = featuretype_plottable[,c("Var1", "descriptor", "module")], sum)
                }
            }
            featuretype_plottable <- featuretype_plottable[order(
                factor(featuretype_plottable[,"module"], levels = levels(featuretype_plottable[,"module"])),
                factor(featuretype_plottable[,"Var1"], levels = levels(featuretype_plottable[,"Var1"])),
                -featuretype_plottable[,"feature_percent"], decreasing = FALSE
            ),]

            
            pout <- ggplot(featuretype_plottable, aes(x = module, y = Freq, fill = Var1))
            pout <- pout + geom_bar(position="fill", stat="identity", color = "black")
            pout <- pout + theme_bw()
            pdf(paste0(outfilepathmaster, "nmf_wgcna_integration/", group_column_selected, "/", "probe_feature_label_plots/", 
                       probe_feature_label, "__", grouping_label, "_stackedbar.pdf"), 15, 7, useDingbats = FALSE)
            print(pout)
            junk <- dev.off()
        }
                                 
    }
    
    # rewrite this
    # input: genestocolortab, counttable, biorep_waddons, the cluster label we want, the order of eigengenes
    # for loop over annotations - grab data  type and counttable and genes to colortab by cluster label
    # plottable is features by samples
    # col annotation is the cluster label and th emetadata we want to know
    # row annotation is the gene module, and the module rank (and annotation?)
    
    
    # ideally, this is all done, and then you simply aggregate by module to get the avg summary hm
    
    summary_heatmap <- cluster_module_meta_summary_heatmap_function(ingenestocolortable = ingenestocolortable,
                                                                    incounttable = incounttable,
                                                                    bioreptable_waddons,
                                                                    selected_clusterlabel = group_column_selected,
                                                                    eigengeneorder = eigengeneorder,
                                                                    metacols_sel = c("SEX", "CAGE_RND", "IMGDEGIS_01_2_3", "CTNDV70_combo0noneval"),
                                                                    colorguide)
    pdf(paste0(outfilepathmaster, "nmf_wgcna_integration/", group_column_selected, "/", group_column_selected, "_eigen_meta_summary_heatmap.pdf"),
        12, 10, useDingbats = FALSE)
    draw(summary_heatmap)
    junk <- dev.off()
    
    # I also want the aggregate plot - wheres its just samples vs eigengenes
    # Create the side annotation table that annotates our eigengenes
    # sel_module_annot_table_temp <- module_full_annotation_table[module_full_annotation_table[,"ID"] %in% module_chosen_annotation_table[,"goterm"], c("module", "p.adjust", "ID")]
    dotplotout <- WGCNA_annotation_dotplot(WGCNA_annotation_table = module_full_annotation_table,
                                           padjcutoffparam = NULL, pathwayselect = module_chosen_annotation_table[,"goterm"],
                                           cleannames = TRUE)
    writeout_dotplot <- dotplotout[[1]] + scale_x_discrete(limits=rev(eigengeneorder)) + coord_flip()
    pdf(paste0(outfilepathmaster, "nmf_wgcna_integration/", group_column_selected, "/", group_column_selected, "_eigen_annot_dotplot.pdf"), 15, 10, useDingbats = FALSE)
    print(writeout_dotplot)
    junk <- dev.off()
    
    # Write out the avg eigengene with ranking annotation
    eigen_avg_heatmap <- eigen_v_cluster_heatmap(eigengenecounttable = ineigencounttable,
                                                 selected_clusterlabel = group_column_selected, eigengeneorder =  eigengeneorder, colorguide =  colorguide)
    pdf(paste0(outfilepathmaster, "nmf_wgcna_integration/", group_column_selected, "/", group_column_selected, "_eigen_avg_heatmap.pdf"), 12, 10, useDingbats = FALSE)
    draw(eigen_avg_heatmap[[1]])
    junk <- dev.off()
    write.table(eigen_avg_heatmap[[2]], paste0(outfilepathmaster, "nmf_wgcna_integration/", group_column_selected, "/", group_column_selected, "_eigen_avg_table.csv"), sep = ",", col.names = NA, row.names = TRUE)
    write.table(eigen_avg_heatmap[["cluster_eigenrank_table_value"]], paste0(outfilepathmaster, "nmf_wgcna_integration/", group_column_selected, "/", group_column_selected, "_cluster_eigenrank_table_value.csv"), sep = ",", col.names = NA, row.names = TRUE)
    cluster_eigenrank_table_value[,paste0("ME", eigengeneorder)]
    
    # Also want the boxplots across subtypes
    ## Need to factorize the annotations with the levels you want for proper ordering in the boxplot
    boxplot_annotationtable <- bioreptable_waddons[,group_column_selected, drop = FALSE]
    boxplot_annotationtable[,1] <- factor(boxplot_annotationtable[,1], levels = sort(unique(boxplot_annotationtable[,1])))
    bp_outfilepath <- paste0(outfilepathmaster, "nmf_wgcna_integration/", group_column_selected, "/", "eigen_bps_by_cluster/")
    
    bp_color_ref <- custom_annotation_list_from_colorguide(COI=group_column_selected, colorguide)
    eigen_cluster_bp_out <- boxplots_from_counttable_by_annotation(counttable = t(ineigencounttable),
                                                                   boxplot_annotationtable = boxplot_annotationtable,
                                                                   outfilepath = bp_outfilepath, calculate_corrvalues = FALSE,
                                                                   bp_color_ref = bp_color_ref  # Adding this on for this specific analysis so I can colorthe bps correctly
    )
    # bp_outfilepath_logscale <- paste0(outfilepathmaster, "nmf_wgcna_integration/", group_column_selected, "/", "eigen_bps_by_cluster_logscale/")
    # eigen_cluster_bp_out_logscale <- boxplots_from_counttable_by_annotation(
    #     # counttable = t(log2(ineigencounttable * 1000)),
    #     counttable = t(ineigencounttable),
    #     boxplot_annotationtable = boxplot_annotationtable,
    #     outfilepath = bp_outfilepath_logscale, calculate_corrvalues = FALSE,
    #     bp_color_ref = bp_color_ref  # Adding this on for this specific analysis so I can colorthe bps correctly
    # )
    
    write.table(eigen_cluster_bp_out[["wilcox_summarytable"]], paste0(outfilepathmaster, "nmf_wgcna_integration/", group_column_selected, "/", "eigen_bps_by_cluster/", group_column_selected, "bp_wilcox_stat_summary_table.csv"), sep = ",", col.names = TRUE, row.names = FALSE)
    write.table(eigen_cluster_bp_out[["spearman_summarytable"]], paste0(outfilepathmaster, "nmf_wgcna_integration/", group_column_selected, "/", "eigen_bps_by_cluster/", group_column_selected, "bp_spearman_stat_summary_table.csv"), sep = ",", col.names = TRUE, row.names = FALSE)
    
    
}



# --------------------------------- In group vs outgroup enrichment and summary --------------------------------- 
## In group vs outgroup analysis enrichment of characteristics
COI <- c("SEX__Female", "AGE_RAND__NUMERIC",  "ETHNIC__Hispanic or Latino", 
         "RACE__American Indian or Alaska Native",  "RACE__Asian",  "RACE__Black or African American",  "RACE__White", "RACE__Multiple Races", 
         "BMI__NUMERIC", "HEMOGLOB__NUMERIC",  "HDLC__NUMERIC",  "LDLC__NUMERIC",  "TRIGLYC__NUMERIC",  "TOTCHOL__NUMERIC", "EGFR__NUMERIC",
         "RVBPDIA__NUMERIC", "RVBPSYS__NUMERIC", "LVEFCONT__NUMERIC",
         "MESTATIN__Yes",  "MEAPASP__Yes", 
         "SMOKSTAT__Current Smoker", "SMOKSTAT__Former Smoker", "SMOKSTAT__Never Smoked",
         "CKD__Yes",  "HYPTENSE__Yes",  "DIABETES__Yes", 
         "IMGDEGIS__None",  "IMGDEGIS__Mild",  "IMGDEGIS__Moderate",  "IMGDEGIS__Severe", 
         "CUSTOM_IMGDEGIS_NONEMILD__NoneMild",  "CUSTOM_IMGDEGIS_NONEMILD__Moderate", "CUSTOM_IMGDEGIS_NONEMILD__Severe",
         "IMGDEGIS_01_2_3__NoneMild", "IMGDEGIS_01_2_3__Moderate", "IMGDEGIS_01_2_3__Severe", 
         "CTNDV50__1",  "CTNDV50__2",  "CTNDV50__3",  "CTMULT50__Yes",
         "CTNDV70__0",  "CTNDV70__1",  "CTNDV70__2",  "CTNDV70__3",  "CTMULT70__Yes",
         "CTNDV70_combo0noneval__0orNoneval", "CTNDV70_combo0noneval__1", "CTNDV70_combo0noneval__2", "CTNDV70_combo0noneval__3", 
         "DUKESCORE__1 Vessel with at least Moderate (>=50%) Plaque", 
         "DUKESCORE__2 Vessels with at least Moderate (>=50%) Plaque or 1 Severe (>=70%) Plaque", 
         "DUKESCORE__3 Vessels with at least Moderate (>=50%) Plaque or 2 Severe (>=70%) Plaque or Severe (>=70%) Proximal LAD Plaque", 
         "DUKESCORE__3 Vessels with Severe (>=70%) Plaque or 2 Severe (>=70%) Plaque Including Proximal LAD",  "DUKESCORE__Left Main >=50%", 
         "C_PRIMARY__0",  "C_PRIMARY__1", 
         "C_CVDMI_P__0",  "C_CVDMI_P__1",
         
         "SAQ7_AF__NUMERIC", "SAQ7_PL__NUMERIC", "SAQ7_QL__NUMERIC", "SAQ7_SS__NUMERIC",
         "TREATMNT__INV",
         
         
         "rna_4cluster_w3AB__nmf_cluster_1", "rna_4cluster_w3AB__nmf_cluster_2", "rna_4cluster_w3AB__nmf_cluster_3A", 
         "rna_4cluster_w3AB__nmf_cluster_3B", "rna_4cluster_w3AB__nmf_cluster_4",
         "rna_4cluster__nmf_cluster_1", "rna_4cluster__nmf_cluster_2", "rna_4cluster__nmf_cluster_3", "rna_4cluster__nmf_cluster_4",
         "meth_4cluster__meth_4cluster_1", "meth_4cluster__meth_4cluster_2", "meth_4cluster__meth_4cluster_3", "meth_4cluster__meth_4cluster_4", 
         "meth_3cluster__meth_3cluster_1", "meth_3cluster__meth_3cluster_2", "meth_3cluster__meth_3cluster_3"                     
)

tabulation_table <- read.table(paste0(outfilepathmaster, "cohort_tabulations/", "tabulation_table.csv"), sep =",", row.names = 1, header = TRUE)

group_columns_for_testing <- c("rna_4cluster_w3AB", "rna_4cluster", "meth_3cluster")
for (group_column_selected in group_columns_for_testing) {
    COI_tabulation_subpath <- paste0(outfilepathmaster, "cohort_tabulations/", group_column_selected, "/")
    dir.create(COI_tabulation_subpath, showWarnings = FALSE, recursive = TRUE)
    out1 <- ingroup_vs_outgroup_cohort_enrichment_tests(tabulation_table, group_column = group_column_selected,
                                                        cohort_tabulations_outpath = COI_tabulation_subpath)
    out1_ratiotable <- out1[[1]]
    out1_citable <- out1[[2]]
    
    # Add Ns to the col labels
    N_annot_table <- table(tabulation_table[,group_column_selected])[colnames(out1_ratiotable)]
    colnames_w_N <- paste0(names(N_annot_table), paste0("__N=", N_annot_table))
    
    # heatmapcolorparam <- colorRamp2(breaks = c(-0.05, -0.000000001, 0.000000001, 0.05), c("white", "darkblue", "darkred", "white"))
    # heatmapcolorparam <- colorRamp2(breaks = c(-5, -0.5, 0, 0.5, 5), c("darkblue", "blue", "white", "red", "darkred"))
    heatmapcolorparam <- colorRamp2(breaks = c(0, 0.99999, 1, 1.00001, 5), c("#000099", "#ccccff", "white", "#ffb2b2", "#b20000"))
    rowmetatable <- data.frame(category = unlist(lapply(strsplit(rownames(out1_ratiotable), split = "__"), function(x) x[1])),
                               row.names = rownames(out1_ratiotable))[COI[COI %in% rownames(out1_ratiotable)],,drop=FALSE]
    rowannotationlist <- annotationlist_builder(rowmetatable) 
    plottable <- out1_ratiotable[COI[COI %in% rownames(out1_ratiotable)],]
    colnames(plottable) <- colnames_w_N
    hm1 <- create_heatmap(
        counttab = plottable,
        # counttab = out1_ratiotable[COI[COI %in% rownames(out1_ratiotable)],],
        subsetnum = FALSE, scale_data = FALSE,
        rowmetatable = rowmetatable,
        # rowmetatable = rowmetatable[COI[COI %in% rownames(out1_ratiotable)],,drop=FALSE], 
        rowannotationlist = rowannotationlist,
        separate_legend = TRUE, heatmapcolorparam = heatmapcolorparam, addborders = TRUE)
    pdf(paste0(COI_tabulation_subpath, "signed_filtered_sigfeature_cluster_hm.pdf"), 10, 10, useDingbats = FALSE)
    draw(hm1[[1]])
    junk <- dev.off()
    
    write.table(out1_ratiotable, paste0(COI_tabulation_subpath, "signed_filtered_sigfeature_cluster_ratio_table.csv"), col.names = NA, row.names = TRUE, sep = ",")
    write.table(out1_citable, paste0(COI_tabulation_subpath, "signed_filtered_sigfeature_cluster_ci_table.csv"), col.names = NA, row.names = TRUE, sep = ",")
}


# --------------------------------- cluster vs event summary heatmaps --------------------------------- 
# First we are going to add a combo here:
bioreptable_waddons[,"rna_4cluster__AND__meth_3cluster"] <- apply(bioreptable_waddons, 1, function(x) {
    out <- ifelse((!is.na(x["rna_4cluster"]) & !is.na(x["meth_3cluster"])), 
                  paste(c(x["rna_4cluster"],  x["meth_3cluster"]), collapse = "__AND__"), NA)
    out
})
# What happens if I add on the intervetion variable to this..........
bioreptable_waddons[,paste0(c("rna_4cluster_w3AB", "rna_4cluster", "meth_3cluster", "meth_4cluster"), "__WithINV")] <- apply(bioreptable_waddons[,c("rna_4cluster_w3AB", "rna_4cluster", "meth_3cluster", "meth_4cluster")], 2, function(x) {
    out1 <- paste0(x, "__", bioreptable_waddons[,"INV"])
    out1 <- ifelse(is.na(x), NA, out1)
    out1
})
## In group vs outgroup analysis enrichment of characteristics
group_columns_for_testing <- c("rna_4cluster_w3AB", "rna_4cluster",
                               "meth_4cluster", "meth_3cluster", "rna_4cluster__AND__meth_3cluster",
                               "rna_4cluster_w3AB__WithINV", "rna_4cluster__WithINV", 
                               "meth_3cluster__WithINV", "meth_4cluster__WithINV"
)
# Want to do a manual addition just for this where we do a 3A vs 3B comparison
group_columns_for_testing <- c(group_columns_for_testing, "rna_3A_v_3B_only")

for (group_column_selected in group_columns_for_testing) {
    
    # Want to do a manual addition just for this where we do a 3A vs 3B comparison
    if (group_column_selected == "rna_3A_v_3B_only") {
        selected_cluster_table <- combined_clustermembership_table[, c("PATNUM", "rna_4cluster_w3AB")]
        selected_cluster_table[,"rna_4cluster_w3AB"] <- ifelse(selected_cluster_table[,"rna_4cluster_w3AB"] %in% c("nmf_cluster_3A", "nmf_cluster_3B"),
                                                               selected_cluster_table[,"rna_4cluster_w3AB"], NA)
        colnames(selected_cluster_table) <- c("PATNUM", "rna_3A_v_3B_only")
    } else {
        # selected_cluster_table <- combined_clustermembership_table[, c("PATNUM", group_column_selected)]
        selected_cluster_table <- bioreptable_waddons[, c("PATNUM", group_column_selected)]
    }
    selected_cluster_table <- selected_cluster_table[order(selected_cluster_table[,group_column_selected]), ]
    
    # need to be only pairwise for the summary to work correctly for now...
    highlowage <- continuous_to_named_quantile(bioreptable_waddons[,"AGE_RAND",drop=FALSE], 3)[[1]]
    highlowage[,1] <- ifelse(highlowage[,1] == "Q3_High", "AGE_RAND_highQ3", "AGE_RAND_low")
    
    COIref_tablist <- list(
        # IMGDEGIS = bioreptable_waddons[,"IMGDEGIS", drop=FALSE],
        IMGDEGIS_01_3 = bioreptable_waddons[,"IMGDEGIS_01_3", drop=FALSE],
        # CTNDV50 = bioreptable_waddons[,"CTNDV50", drop=FALSE],
        # CTNDV50_13_clean = bioreptable_waddons[,"CTNDV50_13_clean", drop=FALSE],
        CTNDV70_13_clean = bioreptable_waddons[,"CTNDV70_13_clean", drop=FALSE],
        CKD = bioreptable_waddons[,"CKD",drop=FALSE],
        highlowage = highlowage
    )
    
    plotCOI <- c(names(COIref_tablist),
                 paste0(group_column_selected, "_",
                        unique(na.omit(selected_cluster_table[,group_column_selected]))[order(unique(na.omit(selected_cluster_table[,group_column_selected])))])
    )
    
    ## Preselect our EOI? nah...
    # EOI = c("")
    
    survivalplot_outfilepath <- paste0(outfilepathmaster, "cluster_v_event_survival/", group_column_selected, "/")
    dir.create(survivalplot_outfilepath, showWarnings = FALSE, recursive = TRUE)
    summary_heatmap <- cluster_vs_event_kmanalysis_summary_heatmap(selected_cluster_table = selected_cluster_table,
                                                                   COIref_tablist = COIref_tablist,
                                                                   survivalplot_outfilepath = survivalplot_outfilepath, 
                                                                   bioreptable_waddons = bioreptable_waddons,
                                                                   eventpvalcutoff = 0.05, return_output_tables = FALSE, plotCOI)
    pdf(paste0(outfilepathmaster, "cluster_v_event_survival/", group_column_selected, "/", group_column_selected, "_event_summary_heatmap.pdf"),
        12, 10, useDingbats = FALSE)
    draw(summary_heatmap)
    junk <- dev.off()
    
    
    # I also want to do the KM analysis with all of the groups together for each group
    if (group_column_selected == "rna_3A_v_3B_only") { next }
    
    survivalplot_outfilepath <- paste0(outfilepathmaster, "cluster_v_event_survival_NONBINARY/", group_column_selected, "/")
    dir.create(survivalplot_outfilepath, showWarnings = FALSE, recursive = TRUE)
    COIref_tablist <- list(
        CUSTOM_IMGDEGIS_NONEMILD = bioreptable_waddons[,"CUSTOM_IMGDEGIS_NONEMILD", drop=FALSE],
        CTNDV70 = bioreptable_waddons[,"CTNDV70", drop=FALSE],
        CKD = bioreptable_waddons[,"CKD",drop=FALSE]
        # age_cat = bioreptable_waddons[,"CAGE_RND",drop=FALSE]
    )
    plotCOI <- c(names(COIref_tablist), group_column_selected)
    
    coxph_control_table <- factorize_metatable(bioreptable_waddons[,c("CAGE_RND", "SEX", "RACE")])
    
    summary_heatmap <- cluster_NONBINARY_vs_event_kmanalysis_summary_heatmap(selected_cluster_table = selected_cluster_table,
                                                                             COIref_tablist = COIref_tablist,
                                                                             survivalplot_outfilepath = survivalplot_outfilepath, 
                                                                             bioreptable_waddons = bioreptable_waddons,
                                                                             eventpvalcutoff = 0.05,
                                                                             return_pairwise_coxph_values = TRUE, coxph_control_table = coxph_control_table,
                                                                             plotCOI)
    pdf(paste0(outfilepathmaster, "cluster_v_event_survival_NONBINARY/", group_column_selected, "/", group_column_selected, "_event_summary_heatmap.pdf"),
        12, 10, useDingbats = FALSE)
    draw(summary_heatmap)
    junk <- dev.off()
    
}



# --------------------------------- cluster vs event summary heatmaps - ONLY INV AND CON --------------------------------- 
# First we are going to add a combo here:
bioreptable_waddons[,"rna_4cluster__AND__meth_3cluster"] <- apply(bioreptable_waddons, 1, function(x) {
    out <- ifelse((!is.na(x["rna_4cluster"]) & !is.na(x["meth_3cluster"])),
                  paste(c(x["rna_4cluster"],  x["meth_3cluster"]), collapse = "__AND__"), NA)
    out
})
# What happens if I add on the intervetion variable to this..........
bioreptable_waddons[,paste0(c("rna_4cluster_w3AB", "rna_4cluster", "meth_3cluster", "meth_4cluster"), "__INVONLY")] <- apply(bioreptable_waddons[,c("rna_4cluster_w3AB", "rna_4cluster", "meth_3cluster", "meth_4cluster")], 2, function(x) {
    out1 <- paste0(x, "__", bioreptable_waddons[,"INV"])
    out1 <- ifelse(is.na(x) | bioreptable_waddons[,"INV"] == "CON", NA, out1)
    out1
})
bioreptable_waddons[,paste0(c("rna_4cluster_w3AB", "rna_4cluster", "meth_3cluster", "meth_4cluster"), "__CONONLY")] <- apply(bioreptable_waddons[,c("rna_4cluster_w3AB", "rna_4cluster", "meth_3cluster", "meth_4cluster")], 2, function(x) {
    out1 <- paste0(x, "__", bioreptable_waddons[,"INV"])
    out1 <- ifelse(is.na(x) | bioreptable_waddons[,"INV"] == "INV", NA, out1)
    out1
})
## In group vs outgroup analysis enrichment of characteristics
group_columns_for_testing <- c("rna_4cluster_w3AB__INVONLY", "rna_4cluster__INVONLY", "meth_3cluster__INVONLY", "meth_4cluster__INVONLY",
                               "rna_4cluster_w3AB__CONONLY", "rna_4cluster__CONONLY", "meth_3cluster__CONONLY", "meth_4cluster__CONONLY"
)
for (group_column_selected in group_columns_for_testing) {
    
    # Want to do a manual addition just for this where we do a 3A vs 3B comparison
    if (group_column_selected == "rna_3A_v_3B_only") {
        selected_cluster_table <- combined_clustermembership_table[, c("PATNUM", "rna_4cluster_w3AB")]
        selected_cluster_table[,"rna_4cluster_w3AB"] <- ifelse(selected_cluster_table[,"rna_4cluster_w3AB"] %in% c("nmf_cluster_3A", "nmf_cluster_3B"),
                                                               selected_cluster_table[,"rna_4cluster_w3AB"], NA)
        colnames(selected_cluster_table) <- c("PATNUM", "rna_3A_v_3B_only")
    } else {
        # selected_cluster_table <- combined_clustermembership_table[, c("PATNUM", group_column_selected)]
        selected_cluster_table <- bioreptable_waddons[, c("PATNUM", group_column_selected)]
    }
    selected_cluster_table <- selected_cluster_table[order(selected_cluster_table[,group_column_selected]), ]
    
    # need to be only pairwise for the summary to work correctly for now...
    highlowage <- continuous_to_named_quantile(bioreptable_waddons[,"AGE_RAND",drop=FALSE], 3)[[1]]
    highlowage[,1] <- ifelse(highlowage[,1] == "Q3_High", "AGE_RAND_highQ3", "AGE_RAND_low")
    
    COIref_tablist <- list(
        # IMGDEGIS = bioreptable_waddons[,"IMGDEGIS", drop=FALSE],
        IMGDEGIS_01_3 = bioreptable_waddons[,"IMGDEGIS_01_3", drop=FALSE],
        # CTNDV50 = bioreptable_waddons[,"CTNDV50", drop=FALSE],
        # CTNDV50_13_clean = bioreptable_waddons[,"CTNDV50_13_clean", drop=FALSE],
        CTNDV70_13_clean = bioreptable_waddons[,"CTNDV70_13_clean", drop=FALSE],
        CKD = bioreptable_waddons[,"CKD",drop=FALSE],
        highlowage = highlowage
    )
    
    plotCOI <- c(names(COIref_tablist),
                 paste0(group_column_selected, "_",
                        unique(na.omit(selected_cluster_table[,group_column_selected]))[order(unique(na.omit(selected_cluster_table[,group_column_selected])))])
    )
    
    survivalplot_outfilepath <- paste0(outfilepathmaster, "cluster_v_event_survival/", group_column_selected, "/")
    dir.create(survivalplot_outfilepath, showWarnings = FALSE, recursive = TRUE)
    summary_heatmap <- cluster_vs_event_kmanalysis_summary_heatmap(selected_cluster_table = selected_cluster_table,
                                                                   COIref_tablist = COIref_tablist,
                                                                   survivalplot_outfilepath = survivalplot_outfilepath,
                                                                   bioreptable_waddons = bioreptable_waddons,
                                                                   eventpvalcutoff = 0.05, return_output_tables = FALSE, plotCOI)
    pdf(paste0(outfilepathmaster, "cluster_v_event_survival/", group_column_selected, "/", group_column_selected, "_event_summary_heatmap.pdf"),
        12, 10, useDingbats = FALSE)
    draw(summary_heatmap)
    junk <- dev.off()
    
    
    # I also want to do the KM analysis with all of the groups together for each group
    if (group_column_selected == "rna_3A_v_3B_only") { next }
    
    survivalplot_outfilepath <- paste0(outfilepathmaster, "cluster_v_event_survival_NONBINARY/", group_column_selected, "/")
    dir.create(survivalplot_outfilepath, showWarnings = FALSE, recursive = TRUE)
    COIref_tablist <- list(
        CUSTOM_IMGDEGIS_NONEMILD = bioreptable_waddons[,"CUSTOM_IMGDEGIS_NONEMILD", drop=FALSE],
        CTNDV70 = bioreptable_waddons[,"CTNDV70", drop=FALSE],
        CKD = bioreptable_waddons[,"CKD",drop=FALSE]
        # age_cat = bioreptable_waddons[,"CAGE_RND",drop=FALSE]
    )
    plotCOI <- c(names(COIref_tablist), group_column_selected)
    
    coxph_control_table <- factorize_metatable(bioreptable_waddons[,c("CAGE_RND", "SEX", "RACE")])
    
    summary_heatmap <- cluster_NONBINARY_vs_event_kmanalysis_summary_heatmap(selected_cluster_table = selected_cluster_table,
                                                                             COIref_tablist = COIref_tablist,
                                                                             survivalplot_outfilepath = survivalplot_outfilepath,
                                                                             bioreptable_waddons = bioreptable_waddons,
                                                                             eventpvalcutoff = 0.05,
                                                                             return_pairwise_coxph_values = TRUE, coxph_control_table = coxph_control_table,
                                                                             plotCOI)
    pdf(paste0(outfilepathmaster, "cluster_v_event_survival_NONBINARY/", group_column_selected, "/", group_column_selected, "_event_summary_heatmap.pdf"),
        12, 10, useDingbats = FALSE)
    draw(summary_heatmap)
    junk <- dev.off()
    
}



# --------------------------------- cluster deseq summary --------------------------------- 
# Read in DESeq and GSEA results
dir.create(paste0(outfilepathmaster, "subtype_deseq_figures/"), showWarnings = FALSE, recursive = TRUE)
subtype_list <- c("RNAtype1_vs_NOTRNAtype1", "RNAtype2_vs_NOTRNAtype2", "RNAtype3_vs_NOTRNAtype3", "RNAtype4_vs_NOTRNAtype4",
                  "RNAtype2_vs_RNAtype1", "RNAtype3A_vs_NOTRNAtype3A", "RNAtype3B_vs_NOTRNAtype3B", "RNAtype3B_vs_RNAtype3A")
deseq_table_list <- list()
deseq_path <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/deseq/"
for (subtype in subtype_list) {
    intable <- read.table(paste0(deseq_path, "comp_RNAsubtype_", subtype, "/deseq_results_comp_RNAsubtype_", subtype, ".csv"), 
                          sep = ",", header = TRUE, row.names = 1)
    deseq_table_list[[subtype]] <- intable
}

gsea_table_list <- list()
gsea_path <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/gsea/"
for (subtype in subtype_list) {
    intable <- read.table(paste0(gsea_path, "comp_RNAsubtype_", subtype, "/comp_RNAsubtype_", subtype, "_gsea_GO.csv"), 
                          sep = ",", header = TRUE, row.names = 1)
    gsea_table_list[[subtype]] <- intable
}


# Total summary heatmap
# selected_deseq_tables <- c("RNAtype1_vs_NOTRNAtype1", "RNAtype2_vs_NOTRNAtype2", "RNAtype3A_vs_NOTRNAtype3A", "RNAtype3B_vs_NOTRNAtype3B", "RNAtype4_vs_NOTRNAtype4")
selected_deseq_tables <- c("RNAtype1_vs_NOTRNAtype1", "RNAtype2_vs_NOTRNAtype2", "RNAtype3_vs_NOTRNAtype3", "RNAtype4_vs_NOTRNAtype4")
gene_selection_vector <- c("all", "positive", "negative")

param_intable <- expand.grid(selected_genes = c("all", "positive", "negative"), rowclusterparam = c(TRUE, FALSE))
param_inlist <- split(param_intable, seq(nrow(param_intable)))

for (selected_parameters in param_inlist) {
    outhm <- subtype_deseq_summary_heatmap(selected_deseq_table_list = deseq_table_list[selected_deseq_tables],
                                           whichgenes = selected_parameters[,"selected_genes"], 
                                           normcounttable=normcounttable[,SOIrna], 
                                           selected_rna_clustermembership_table=rna_clustermembership_table[,"rna_4cluster_w3AB",drop=FALSE],
                                           rowclusterparam = selected_parameters[,"rowclusterparam"])
    pdf(paste0(outfilepathmaster, "subtype_deseq_figures/", "nmfcluster_", selected_parameters[,"selected_genes"], "_diffexpgenes_",
               selected_parameters[,"rowclusterparam"], "cluster_heatmap.pdf"), useDingbats = FALSE, height = 7, width = 13)
    draw(outhm)
    junk <- dev.off()
}


# Individual published volcano plots
dir.create(paste0(outfilepathmaster, "subtype_deseq_figures/", "nmfcluster_dge_volcanoplots/"), showWarnings = FALSE, recursive = TRUE)
for (subtype in c("RNAtype1_vs_NOTRNAtype1", "RNAtype2_vs_NOTRNAtype2", "RNAtype3_vs_NOTRNAtype3", "RNAtype4_vs_NOTRNAtype4",
                  "RNAtype3A_vs_NOTRNAtype3A", "RNAtype3B_vs_NOTRNAtype3B", "RNAtype3B_vs_RNAtype3A")) {
    table_select <- deseq_table_list[[subtype]]
    name_select <- strsplit(subtype, split = "_")[[1]][1]
    if (subtype == "RNAtype3B_vs_RNAtype3A") {name_select <- "RNAtype3Bv3A"}
    
    GOI_table <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run3_rmoutliers2/integrative_analyses/run4/subtype_deseq_figures/subtype_GOI_table.csv", sep = ",", header = TRUE, check.names = FALSE)
    GOI <- GOI_table[,subtype]
    if (sum(!is.na(GOI)) == 0) {GOI <- NULL}
    
    volcanoplot_out <- clean_volcano_plot_function(deseqtable=table_select, nameparam=name_select, labeledgenes = GOI)
    pdf(paste0(outfilepathmaster, "subtype_deseq_figures/", "nmfcluster_dge_volcanoplots/", name_select, "_volcano_plot.pdf"), width =8, height = 8)
    # jpeg(paste0(outfilepathmaster, "subtype_deseq_figures/", "nmfcluster_dge_volcanoplots/", name_select, "_volcano_plot.jpeg"), width =768, height = 768)
    print(volcanoplot_out)
    junk <- dev.off()
}



# Individual GSEA summary outputs
dir.create(paste0(outfilepathmaster, "subtype_deseq_figures/", "nmfcluster_dge_gsea/"), showWarnings = FALSE, recursive = TRUE)
for (subtype in c("RNAtype1_vs_NOTRNAtype1", "RNAtype2_vs_NOTRNAtype2", "RNAtype3_vs_NOTRNAtype3", "RNAtype4_vs_NOTRNAtype4",
                  "RNAtype3A_vs_NOTRNAtype3A", "RNAtype3B_vs_NOTRNAtype3B", "RNAtype3B_vs_RNAtype3A")) {
    table_select <- gsea_table_list[[subtype]]
    name_select <- strsplit(subtype, split = "_")[[1]][1]
    if (subtype == "RNAtype3B_vs_RNAtype3A") {name_select <- "RNAtype3Bv3A"}
    
    POI_table <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run3_rmoutliers2/integrative_analyses/run4/subtype_deseq_figures/subtype_POI_table.csv", sep = ",", header = TRUE, check.names = FALSE)
    POI <- POI_table[,subtype]
    
    gsea_out <- clean_gsea_plot_function(gseatable=table_select, nameparam=subtype, pathwayselect = POI)
    pdf(paste0(outfilepathmaster, "subtype_deseq_figures/", "nmfcluster_dge_gsea/", name_select, "_gsea_plot.pdf"), width =8, height = 8)
    # jpeg(paste0(outfilepathmaster, "subtype_deseq_figures/", "nmfcluster_dge_gsea/", name_select, "_gsea_plot.jpeg"), width =768, height = 768)
    print(gsea_out)
    junk <- dev.off()
}






# --------------------------------- cohort signature scores ---------------------------------
# Read in meth count table and metatable
methmetatable_file <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/final_meth_data_v1_20230207/kpp_metadata.rds"
methmetatable <- readRDS(methmetatable_file)
methcounttable_file <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/final_meth_data_v1_20230207/limma_noCovariate.rds"
methcounttable <- readRDS(methcounttable_file)
colnames(methcounttable) <- methmetatable[match(methmetatable[,"Sample_Name"], colnames(methcounttable)), "Sample_Patnum"]
methcounttable <- methcounttable[,-578] # Need to remove the duplicate of this sample: "048003-005", columns 204 and 578
methmetatable <- methmetatable[-578,] # Need to remove the duplicate of this sample: "048003-005", columns 204 and 578

# Need to grab this metatable unfortunately
metatable_file <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run3b_rmoutliers2_addoncomps/rna_processing/metatable_in.csv"
metatable_forsigscore <- read.table(metatable_file, sep = ",", header = TRUE, row.names = 1)
## Need to add on our meth comparisons
meth_metatable_forsigscore <- bioreptable_waddons[SOImeth, c("DEGRISCH", "CTNDV70", "meth_3cluster")]
meth_metatable_forsigscore[,"DMP_on_degrisch_sev_nomild"] <- ifelse(meth_metatable_forsigscore[,"DEGRISCH"] %in% "Severe", 1,
                                                             ifelse(meth_metatable_forsigscore[,"DEGRISCH"] %in% c("None", "Mild"), 0, NA))
meth_metatable_forsigscore[,"DMP_on_CTNDV70_3vs1"] <- ifelse(meth_metatable_forsigscore[,"CTNDV70"] %in% "3", 1,
                                                      ifelse(meth_metatable_forsigscore[,"CTNDV70"] %in% c("1"), 0, NA))
meth_metatable_forsigscore[,"DMP_on_krcpp_1_not1"] <- ifelse(meth_metatable_forsigscore[,"meth_3cluster"] %in% "meth_3cluster_1", 1, 0)
meth_metatable_forsigscore[,"DMP_on_krcpp_2_not2"] <- ifelse(meth_metatable_forsigscore[,"meth_3cluster"] %in% "meth_3cluster_2", 1, 0)
meth_metatable_forsigscore[,"DMP_on_krcpp_3_not3"] <- ifelse(meth_metatable_forsigscore[,"meth_3cluster"] %in% "meth_3cluster_3", 1, 0)

metatable_forsigscore <- merge(metatable_forsigscore, meth_metatable_forsigscore, by = "row.names", all = TRUE)
rownames(metatable_forsigscore) <- metatable_forsigscore[,"Row.names"]
metatable_forsigscore <- metatable_forsigscore[,!grepl("Row.names", colnames(metatable_forsigscore))]


# Start with RNA - just easier for me to conceptualize
dir.create(paste0(outfilepathmaster, "signature_scoring/"), showWarnings = FALSE, recursive = TRUE)

# Set paths for where the DGE files are
deseq_path <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/deseq/"
# DMP_path <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/methylation_data_20220426/methylation_output_Tosh/"
# DMP_path <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/methylation_data_20220608/DMP/"
DMP_path <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/final_meth_data_v1_20230207/DMP_with_asr_correction/"
# Create list of dges to score
rna_dge_cutoffs <- c(pval_test = "padj", pval_cutoff = 0.05, log2fc_cutoff = 0.25)
meth_dmp_cutoffs <- c(pval_test = "adj.P.Val", pval_cutoff = 0.05, log2fc_cutoff = 0.03)
dge_comp_list <- list(
    rna_imgdegis_score1 = c(complabel = "comp_ischemia__Sev_v_MildNone", pval_test = comp_ischemia__Sev_v_MildNone__PVALTEST,
                            pval_cutoff = comp_ischemia__Sev_v_MildNone__PVALCUTOFF, log2fc_cutoff = comp_ischemia__Sev_v_MildNone__FCCUTOFF, COI = "IMGDEGIS"),
    rna_ctndv70_score1 = c(complabel = "comp_anatomy70__3v_v_1v", pval_test = comp_anatomy70__3v_v_1v__PVALTEST,
                           pval_cutoff = comp_anatomy70__3v_v_1v__PVALCUTOFF, log2fc_cutoff = comp_anatomy70__3v_v_1v__FCCUTOFF, COI = "CTNDV70"),
    rna_subtype1 = c(complabel = "comp_RNAsubtype_RNAtype1_vs_NOTRNAtype1", rna_dge_cutoffs, COI = "rna_4cluster_w3AB"),
    rna_subtype2 = c(complabel = "comp_RNAsubtype_RNAtype2_vs_NOTRNAtype2", rna_dge_cutoffs, COI = "rna_4cluster_w3AB"),
    rna_subtype3 = c(complabel = "comp_RNAsubtype_RNAtype3_vs_NOTRNAtype3", rna_dge_cutoffs, COI = "rna_4cluster_w3AB"),
    rna_subtype3A = c(complabel = "comp_RNAsubtype_RNAtype3A_vs_NOTRNAtype3A", rna_dge_cutoffs, COI = "rna_4cluster_w3AB"),
    rna_subtype3B = c(complabel = "comp_RNAsubtype_RNAtype3B_vs_NOTRNAtype3B", rna_dge_cutoffs, COI = "rna_4cluster_w3AB"),
    rna_subtype4 = c(complabel = "comp_RNAsubtype_RNAtype4_vs_NOTRNAtype4", rna_dge_cutoffs, COI = "rna_4cluster_w3AB"),
    rna_subtype3Bv3A = c(complabel = "comp_RNAsubtype_RNAtype3B_vs_RNAtype3A", rna_dge_cutoffs, COI = "rna_4cluster_w3AB"),
    
    meth_imgdegis_score1 = c(complabel = "DMP_on_degrisch_sev_nomild", meth_dmp_cutoffs, COI = "DEGRISCH"),
    meth_ctndv70_score1 = c(complabel = "DMP_on_CTNDV70_3vs1", c(pval_test = "P.Value", pval_cutoff = 0.05, log2fc_cutoff = 0.03), COI = "CTNDV70"),
    # Note that for all of these, we have to grab the same object and then select the actual one we want
    # meth_subtype1 = c(complabel = "DMP_on_krcpp_3_xVSothers", c(pval_test = "adj.P.Val", pval_cutoff = 0.05, log2fc_cutoff = 0.1), COI = "meth_3cluster"),
    # meth_subtype2 = c(complabel = "DMP_on_krcpp_3_xVSothers", c(pval_test = "adj.P.Val", pval_cutoff = 0.05, log2fc_cutoff = 0.1), COI = "meth_3cluster"),
    # meth_subtype3 = c(complabel = "DMP_on_krcpp_3_xVSothers", c(pval_test = "adj.P.Val", pval_cutoff = 0.05, log2fc_cutoff = 0.1), COI = "meth_3cluster"),
    meth_subtype1 = c(complabel = "DMP_on_krcpp_1_not1", c(pval_test = "adj.P.Val", pval_cutoff = 0.05, log2fc_cutoff = 0.03), COI = "meth_3cluster"),
    meth_subtype2 = c(complabel = "DMP_on_krcpp_2_not2", c(pval_test = "adj.P.Val", pval_cutoff = 0.05, log2fc_cutoff = 0.03), COI = "meth_3cluster"),
    meth_subtype3 = c(complabel = "DMP_on_krcpp_3_not3", c(pval_test = "adj.P.Val", pval_cutoff = 0.05, log2fc_cutoff = 0.03), COI = "meth_3cluster")
)


# Ok I just need to clean this object, cause theres a couple things that need to be adjusted
# meth_subtype_DMP_object <- readRDS(paste0(DMP_path, "DMP_k3_limma_noCov_oneVSother.rds"))
# OLDTEMP <- readRDS(paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/methylation_data_20220608/DMP/", "DMP_k3_limma_noCov_oneVSother.rds"))
meth_subtype_DMP_object <- readRDS(paste0(DMP_path, "limma_BCcorrected_DMP_on_kpp_oneVSother_control.rds"))
# cleaned_meth_subtype_DMP_object <- clean_meth_subtype_DMP_object_function(meth_subtype_DMP_object)


deseq_table_list <- list()
for (dge_comp_num in seq_len(length(dge_comp_list))) {
    # Grab the comp we are running
    dge_comp_sel <- dge_comp_list[[dge_comp_num]]
    dge_comp_label <- names(dge_comp_list)[dge_comp_num]
    complabel_sel <- dge_comp_sel["complabel"]
    pval_test_sel <- dge_comp_sel["pval_test"]
    pval_cutoff_sel <- as.numeric(dge_comp_sel["pval_cutoff"])
    log2fc_cutoff_sel <- as.numeric(dge_comp_sel["log2fc_cutoff"])
    COI_sel <- dge_comp_sel["COI"]
    # Also grab the data type - which i am just going to standardize for myself as appedning to the beginngin of the label
    datatype_sel <- strsplit(dge_comp_label[[1]], split = "_")[[1]][1]
    
    # Create the subfolder output
    signaturescore_outfilepath <- paste0(outfilepathmaster, "signature_scoring/", dge_comp_label, "/")
    dir.create(signaturescore_outfilepath, showWarnings = FALSE, recursive = TRUE)
    
    # Grab our dge table and counttable
    if (datatype_sel == "rna") {
        intable_deseq <- read.table(paste0(deseq_path, complabel_sel, "/deseq_results_", complabel_sel, ".csv"), 
                                    sep = ",", header = TRUE, row.names = 1)
        counttable_sel <- normcounttable
        run_WGCNAscore_param <- TRUE
    }
    if (datatype_sel == "meth") {
        
        counttable_sel <- methcounttable
        # the meth data is a little more complicated - so we need some specific rules for other things
        if (dge_comp_label %in% c("meth_subtype1", "meth_subtype2", "meth_subtype3")) {
            grab_subtype_number <- gsub("meth_subtype", "", dge_comp_label)
            # meth_subtype_DMP_object <- readRDS(paste0(DMP_path, "DMP_on_krcpp_3_xVSothers.rds"))
            # meth_subtype_DMP_object <- readRDS(paste0(DMP_path, "DMP_k3_limma_noCov_oneVSother.rds"))
            # intable_deseq <- meth_subtype_DMP_object[[grab_subtype_number]][[1]]
            
            # OMG - for the DMP objects, they arent all the same direction (cmon Ze......) cleaned and fixed now
            # cleaned_meth_subtype_DMP_object
            intable_deseq <- cleaned_meth_subtype_DMP_object[[grep(grab_subtype_number, names(cleaned_meth_subtype_DMP_object))]]

        } else {
            intable_deseq <- readRDS(paste0(DMP_path, complabel_sel, ".rds"))
            colnames(intable_deseq)[colnames(intable_deseq) == "logFC"] <- "signed_deltaBeta"
        }
        colnames(intable_deseq)[colnames(intable_deseq) == "signed_deltaBeta"] <- "log2FoldChange"
        
        run_WGCNAscore_param <- FALSE
    }
    
    # grab the metadata we need
    comp_metatable <- merge(bioreptable_waddons[,COI_sel,drop = FALSE], metatable_forsigscore[,complabel_sel,drop=FALSE], by = "row.names")
    rownames(comp_metatable) <- comp_metatable[,"Row.names"]
    comp_metatable <- comp_metatable[colnames(counttable_sel), !grepl("Row.names", colnames(comp_metatable))]
    
    # Now i want a wrapper funciton, that i give it a DESeq2 output with cutoffs, and it gives me tables, plots, etc.
    sample_signature_analysis(intable_deseq, dge_comp_label, pval_test_sel, pval_cutoff_sel, log2fc_cutoff_sel,
                              comp_metatable, counttable_sel=counttable_sel, signaturescore_outfilepath,
                              run_singscore = TRUE, run_WGCNAscore = run_WGCNAscore_param,
                              scoring_genesetsize = c("all", 1000, 500, 100, 50))
    
}

# Test these signatures against events
# Function for taking events with continuous time to event, and creating new columns filtered by day cutoffs - censoring appropriately
EOI_grab <- c("PRIMARY", "CVDMI_S", "ACDMI_S", "CVD", "ACD", "MIPRIM")
toe_filter_inlist <- lapply(EOI_grab, function(x) {bioreptable_waddons[,c(paste0("T_", x), paste0("C_", x))]})
time_filters <- c(182.625, 365.250, 547.875, 730.500, 1095.750)
toe_table_outlist <- list()
for (toe_table in toe_filter_inlist) {
    filtered_toe_table <- filter_time_to_event_table(time_to_event_table=toe_table, time_filters=time_filters)
    toe_label <- gsub("T_", "", colnames(toe_table)[1])
    toe_table_outlist[[toe_label]] <- filtered_toe_table
}
EOI_table <- do.call(cbind.data.frame, toe_table_outlist)
colnames(EOI_table) <- unlist(lapply(toe_table_outlist, colnames))


# I want some boxplots now for our events based on these types........ Like each score for event (like in validation)
singscore_path <- paste0(outfilepathmaster, "signature_scoring/")
singscore_list <- c("rna_subtype1", "rna_subtype2", "rna_subtype3", "rna_subtype4"
                    # "meth_subtype1", "meth_subtype2", "meth_subtype3"
                    )
singscore_table_list <- list()
singscore_path <- paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/integrative_analyses/run7/",
                         "signature_scoring/")
scoring_genesetsize = c("all", 1000, 500, 100, 50)
singscore_grid <- expand.grid(singscore_list, scoring_genesetsize)
for (singscore_param_num in seq_len(nrow(singscore_grid))) {
    singscore_result <- singscore_grid[singscore_param_num,"Var1"]
    genesetsize <- singscore_grid[singscore_param_num,"Var2"]
    outlabel <- paste0(singscore_result, "__", genesetsize)
    singscore_basepath <- paste0(singscore_path, singscore_result, "/singscore/")
    singscore_subpath <- list.files(singscore_basepath)[grepl(paste0(genesetsize, "UP"), list.files(singscore_basepath)) | grepl(paste0(genesetsize, "DOWN"), list.files(singscore_basepath))]
    scoretable_filename <- list.files(paste0(singscore_basepath, singscore_subpath, "/"))[
        grepl("singscore_outtable", list.files(paste0(singscore_basepath, singscore_subpath, "/")))]
    singscore_table_list[[outlabel]] <- singscore_table_grab <- read.table(paste0(singscore_basepath, singscore_subpath, "/", scoretable_filename),
                                                                               sep = ",", header = TRUE)
}

combined_singscore_table <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Row.names", all = TRUE, sort = FALSE), 
    lapply(singscore_table_list, function(x) x[,c("Row.names", "TotalScore", "UpScore", "DownScore")]))
dimnames(combined_singscore_table) <- list(combined_singscore_table[,"Row.names"], 
                                                c("Row.names", apply(expand.grid(paste0(c("TotalScore", "UpScore", "DownScore"), "__"),
                                                                                 names(singscore_table_list)), 1, paste0, collapse = "")))
combined_singscore_table <- combined_singscore_table[,!grepl("Row.names", colnames(combined_singscore_table))]


## grab for test
tt1 <- singscore_table_list[["rna_subtype2__50"]]
tt2 <- merge(tt1, EOI_table[rownames(combined_singscore_table),!grepl("T_|C_", colnames(EOI_table))], by.x = "Row.names", by.y = "row.names")
aggregate(tt2[,"TotalScore"], by = list(tt2[,"filtered1095.75_PRIMARY"]), summary)


## Add "median centered" scores as well

# param_table <- setNames(data.frame(t(expand.grid(c("singscore_table", "medcentered_singscore_table"), "EOI_table"))), c("counttable_sel", "eventtable_sel"))
param_table <- expand.grid("counttable_sel" = c("singscore_table", "medcentered_singscore_table"), "eventtable_sel" = "EOI_table")
param_table[,"outpathlabel"] <- c("sstable_allevents", "medsstable_allevents")
for (annot_param_num in seq_len(nrow(param_table))) {
    
    ## Grab counttable
    if(param_table[annot_param_num, "counttable_sel"] == "singscore_table") {counttable_sel <- combined_singscore_table}
    if(param_table[annot_param_num, "counttable_sel"] == "medcentered_singscore_table") {
        counttable_sel <- apply(combined_singscore_table, 2, function(x) x/median(x, na.rm = TRUE))
    }
    ## Grab annot table
    if(param_table[annot_param_num, "eventtable_sel"] == "EOI_table") {annottable_sel <- EOI_table}
    ## Create outfolder label
    outpathlabel_sel <- param_table[annot_param_num, "outpathlabel"]
    
    ## Create our annotation table from what we have
    boxplot_annotationtable <- annottable_sel[rownames(counttable_sel),!grepl("T_|C_", colnames(annottable_sel))]
    boxplot_annotationtable[boxplot_annotationtable == 2] <- NA
    outfilepath <- paste0(outfilepathmaster, "signature_scoring/", "event_comparisons/", outpathlabel_sel, "/")
    dir.create(outfilepath, showWarnings = FALSE, recursive = TRUE)
    # write.table(boxplot_annotationtable, paste0(outfilepathmaster, "eigen_v_annot/", boxplot_label, "/", boxplot_label, "_boxplot_annotationtable.csv"), sep = ",", col.names = NA, row.names = TRUE)
    
    ## Failsafe here:
    boxplot_annotationtable <- boxplot_annotationtable[,!colSums(is.na(boxplot_annotationtable)) == nrow(boxplot_annotationtable)]
    boxplots_from_counttable_by_annotation_out <- boxplots_from_counttable_by_annotation(
        counttable = t(counttable_sel),
        boxplot_annotationtable = boxplot_annotationtable,
        outfilepath = outfilepath, 
        calculate_corrvalues = FALSE, bp_color_ref = NULL)
    bp_stat_table <- data.frame(boxplots_from_counttable_by_annotation_out[[1]])
    bp_stat_table[,"signed_wilcox"] <- as.numeric(bp_stat_table[,"wilcox_pval"]) * apply(bp_stat_table[,c("Feature1_mean", "Feature2_mean")], 1, function(x) ifelse(as.numeric(x[1]) > as.numeric(x[2]), 1, -1))
    write.table(bp_stat_table, paste0(outfilepath, "EOI_table", "_stattable.csv"), sep = ",", col.names = TRUE, row.names = FALSE)
    
    # Now I want a table to make a heatmap viz out of
    wilcox_pval_table <- reshape2::dcast(bp_stat_table, Cohort ~ Feature, value.var = "signed_wilcox")
    rownames(wilcox_pval_table) <- wilcox_pval_table[,"Cohort"]
    wilcox_pstat_table <- -log10(abs(wilcox_pval_table[colnames(boxplot_annotationtable), !colnames(wilcox_pval_table) %in% "Cohort"]))
    ## FLIPPING THE SIGN HERE TO MAKE THE HEATMAP MORE READABLE
    wilcox_sign_table <- -sign(wilcox_pval_table[colnames(boxplot_annotationtable), !colnames(wilcox_pval_table) %in% "Cohort"])
    ## FLIPPING THE SIGN HERE TO MAKE THE HEATMAP MORE READABLE
    wilcox_plottable <- data.frame(t(wilcox_pstat_table * 
                                         wilcox_sign_table))[, colnames(boxplot_annotationtable), drop = FALSE]
    
    write.table(wilcox_plottable, paste0(outfilepath, "EOI_table", "_plottable.csv"), sep = ",", col.names = NA, row.names = TRUE)
    
    # Plot out a heatmap - and add a reference for the eigengenes as being "risk" ones of "protective" ones
    rowmetatable <- eigenegene_to_cluster_table <- data.frame(eigengene = rownames(wilcox_plottable), 
        rna_4cluster = ifelse(grepl("rna_subtype1", rownames(wilcox_plottable)), "nmf_cluster_1",
                       ifelse(grepl("rna_subtype2", rownames(wilcox_plottable)), "nmf_cluster_2", 
                       ifelse(grepl("rna_subtype3", rownames(wilcox_plottable)), "nmf_cluster_3", 
                       ifelse(grepl("rna_subtype4", rownames(wilcox_plottable)), "nmf_cluster_4", NA)))), 
                                                              row.names = rownames(wilcox_plottable))
    rowmetatable <- rowmetatable[order(rowmetatable[,2]),,drop=FALSE]
    rowannotationlist <- custom_annotation_list_from_colorguide("rna_4cluster", colorguide = colorguide)
    
    heatmapcolorparam <- colorRamp2(breaks = c(-3, log10(0.05), 0, -log10(0.05), 3), c("blue", "white", "white", "white", "red"))
    hmout <- create_heatmap(counttab = wilcox_plottable[rownames(rowmetatable),,drop=FALSE], subsetnum = FALSE, scale_data = FALSE,
                            colmetatable = NULL, colannotationlist = NULL,
                            rowmetatable = rowmetatable, rowannotationlist = rowannotationlist,
                            addborders = TRUE, heatmapcolorparam = heatmapcolorparam
    )
    pdf(paste0(outfilepath, "EOI_table", "_heatmap.pdf"), useDingbats = FALSE, width = 15, height = 10)
    draw(hmout[[1]])
    junk <- dev.off()
}





# --------------------------------- cohort singscore summaries ---------------------------------
# For each score metric (combo of score metric ands gene size)
## TOSH WORKING HERE - FLIP DISPEERISON TO MIN INSTEAD OF MAX, ADD BACK IN ML, RERUN, SORT OUT THE END TABLE
## I want the 3 heatmaps I am making - (1) full cluster, (2) RNA split, (3) Meth split
## And then I want a basic prediction - if we take the top score for each sample, how well does that predict the assigned cluster
dir.create(paste0(outfilepathmaster, "signature_scoring/singscore_summaries/"), showWarnings = FALSE, recursive = TRUE)
scoring_ss_metric <- c("TotalScore", "TotalDispersion", "UpScore", "UpDispersion", "DownScore", "DownDispersion")
scoring_genesetsize <- rev(c("all", 1000, 500, 100, 50))
# scoring_ss_metric <- c("TotalScore", "UpScore")
# scoring_genesetsize <- rev(c("all", 500))
singscore_grid <- expand.grid(scoring_ss_metric, scoring_genesetsize)
pred_out_stat_list <- list()
for (singscore_param_num in seq_len(nrow(singscore_grid))) {
    ss_metric <- as.character(singscore_grid[singscore_param_num,"Var1"])
    genesetsize <- as.character(singscore_grid[singscore_param_num,"Var2"])
    outlabel <- paste0(ss_metric, "__", genesetsize)
    dir.create(paste0(outfilepathmaster, "signature_scoring/singscore_summaries/", outlabel, "/"), showWarnings = FALSE, recursive = TRUE)

    singscore_path <- paste0(outfilepathmaster, "signature_scoring/")
    singscore_list <- c("rna_subtype1", "rna_subtype2", "rna_subtype3", "rna_subtype4", "meth_subtype1", "meth_subtype2", "meth_subtype3")
    ss_table_list <- list()
    for (subtype_label in singscore_list) {
        singscore_basepath <- paste0(singscore_path, subtype_label, "/singscore/")
        subpath <- list.files(singscore_basepath)[grepl(paste0(genesetsize, "UP"), list.files(singscore_basepath)) | grepl(paste0(genesetsize, "DOWN"), list.files(singscore_basepath))]
        ss_table_list[[subtype_label]] <- read.table(paste0(singscore_basepath, subpath, "/", subtype_label, "_singscore_outtable.csv"),
                                                     sep =",", header = TRUE)[,c("Row.names", ss_metric)]
    }
    combined_singscore_table <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Row.names", all = TRUE, sort = FALSE), 
                                       lapply(ss_table_list, function(x) x[,c("Row.names", ss_metric)]))
    dimnames(combined_singscore_table) <- list(combined_singscore_table[,"Row.names"], c("Row.names", names(ss_table_list)))
    combined_singscore_table <- combined_singscore_table[,!grepl("Row.names", colnames(combined_singscore_table))]
    
    
    ## Make heatmaps
    dir.create(paste0(outfilepathmaster, "signature_scoring/singscore_summaries/", outlabel, "/heatmaps/"), showWarnings = FALSE, recursive = TRUE)
    summary_hm_plottab <- t(na.omit(combined_singscore_table))
    singscore_sum_hm_metatable <- bioreptable_waddons[colnames(summary_hm_plottab), c("IMGDEGIS_01_2_3", "CTNDV70", "rna_4cluster", "meth_3cluster")]
    singscore_sum_hm_annottable <- custom_annotation_list_from_colorguide(
        COI = c("IMGDEGIS_01_2_3", "CTNDV70", "rna_4cluster", "meth_3cluster"), colorguide)
    
    singscore_summary_hm <- create_heatmap(counttab = summary_hm_plottab, scale_data = FALSE,
                                           colmetatable = singscore_sum_hm_metatable, colannotationlist = singscore_sum_hm_annottable,
                                           colclusterparam = TRUE, rowclusterparam = FALSE)
    pdf(paste0(outfilepathmaster, "signature_scoring/singscore_summaries/", outlabel, "/heatmaps/", "singscore_summary_hm.pdf"), 10, 8, useDingbats = FALSE)
    draw(singscore_summary_hm[[1]])
    junk <- dev.off()
    
    singscore_summary_hm_rnasplit <- create_heatmap(counttab = summary_hm_plottab, scale_data = TRUE,
                                           colmetatable = singscore_sum_hm_metatable, colannotationlist = singscore_sum_hm_annottable,
                                           colclusterparam = TRUE, rowclusterparam = FALSE,
                                           columnsplitparam = singscore_sum_hm_metatable[,"rna_4cluster"],
                                           heatmapcolorparam = colorRamp2(breaks = c(-2, 0, 2), c("darkblue", "white", "darkred")))
    pdf(paste0(outfilepathmaster, "signature_scoring/singscore_summaries/", outlabel, "/heatmaps/", "singscore_summary_hm_rnasplit.pdf"), 10, 8, useDingbats = FALSE)
    draw(singscore_summary_hm_rnasplit[[1]])
    junk <- dev.off()
    
    singscore_summary_hm_methsplit <- create_heatmap(counttab = summary_hm_plottab, scale_data = TRUE,
                                           colmetatable = singscore_sum_hm_metatable, colannotationlist = singscore_sum_hm_annottable,
                                           colclusterparam = TRUE, rowclusterparam = FALSE,
                                           columnsplitparam = singscore_sum_hm_metatable[,"meth_3cluster"],
                                           heatmapcolorparam = colorRamp2(breaks = c(-2, 0, 2), c("darkblue", "white", "darkred")))
    pdf(paste0(outfilepathmaster, "signature_scoring/singscore_summaries/", outlabel, "/heatmaps/", "singscore_summary_hm_methsplit.pdf"), 10, 8, useDingbats = FALSE)
    draw(singscore_summary_hm_methsplit[[1]])
    junk <- dev.off()
    
    
    ## NEED A LOT OF FUNCTIONALIZATION AND LOOPING HERE
    dir.create(paste0(outfilepathmaster, "signature_scoring/singscore_summaries/", outlabel, "/prediction/"), showWarnings = FALSE, recursive = TRUE)
    ss_intable_list <- list(
        list(label = "rna_4cluster", true_ref_col = "rna_4cluster", ss_intable = combined_singscore_table[,c("rna_subtype1", "rna_subtype2", "rna_subtype3", "rna_subtype4")]),
        list(label = "meth_3cluster", true_ref_col = "meth_3cluster", ss_intable = combined_singscore_table[,c("meth_subtype1", "meth_subtype2", "meth_subtype3")])
    )
    for (ss_intable_list_sel in ss_intable_list) {
        # Grab VOI and create outfile path
        ss_test_label = ss_intable_list_sel[["label"]]
        ss_intable = ss_intable_list_sel[["ss_intable"]]
        ## Rename our ss_intables to match our true tables - actually cleaner in the end (and makes more sense)
        if (ss_test_label == "rna_4cluster") {colnames(ss_intable) <- c("nmf_cluster_1", "nmf_cluster_2", "nmf_cluster_3", "nmf_cluster_4")}
        if (ss_test_label == "meth_3cluster") {colnames(ss_intable) <- c("meth_3cluster_1", "meth_3cluster_2", "meth_3cluster_3")}
        subanalysis_outfilepath <- paste0(outfilepathmaster, "signature_scoring/singscore_summaries/", outlabel, "/prediction/", ss_test_label, "/")
        dir.create(subanalysis_outfilepath, showWarnings = FALSE, recursive = TRUE)
        
        # Set Truth
        true_label_table <- na.omit(bioreptable_waddons[,ss_intable_list_sel[["true_ref_col"]],drop=FALSE])
        true_label_table[,1] <- factorize_metatable(true_label_table[,1,drop = FALSE])
        prediction_factors <- levels(true_label_table[,1])
        
        # First test - raw value prediction:
        dir.create(paste0(subanalysis_outfilepath, "rawval_pred/"), showWarnings = F, recursive = T)
        if (grepl("Dispersion", ss_metric)) {
            raw_predict_table <- data.frame(raw_predict_class = factor(unlist(apply(ss_intable, 1, function(x) names(x)[which.min(x)])), 
                                                                       levels = prediction_factors))
        } else {
            raw_predict_table <- data.frame(raw_predict_class = factor(unlist(apply(ss_intable, 1, function(x) names(x)[which.max(x)])), 
                                                                       levels = prediction_factors))
        }
        raw_pred_stats_out <- return_pred_stats(true_labels = true_label_table[,1], test_labels = raw_predict_table[rownames(true_label_table),])
        write.table(raw_pred_stats_out[["confstatsout"]], paste0(subanalysis_outfilepath, "rawval_pred/", "rawval_confstatsout.csv"), sep = ",", row.names = FALSE, col.names = TRUE)
        write.table(raw_pred_stats_out[["overall_statout"]], paste0(subanalysis_outfilepath, "rawval_pred/", "rawval_overall_statout.csv"), sep = ",", row.names = FALSE, col.names = TRUE)
        
        # Second test - zscored values
        dir.create(paste0(subanalysis_outfilepath, "zval_pred/"), showWarnings = F, recursive = T)
        zscore_rna_ss_table <- apply(ss_intable, 2, zscore)
        if (grepl("Dispersion", ss_metric)) {
            zscore_predict_table <- data.frame(raw_predict_class = factor(unlist(apply(zscore_rna_ss_table, 1, function(x) names(x)[which.min(x)])), 
                                                                          levels = prediction_factors))
        } else {
            zscore_predict_table <- data.frame(raw_predict_class = factor(unlist(apply(zscore_rna_ss_table, 1, function(x) names(x)[which.max(x)])), 
                                                                          levels = prediction_factors))
        }
        zscore_pred_stats_out <- return_pred_stats(true_labels = true_label_table[,1], test_labels = zscore_predict_table[rownames(true_label_table),])
        write.table(zscore_pred_stats_out[["confstatsout"]], paste0(subanalysis_outfilepath, "zval_pred/", "zval_confstatsout.csv"), sep = ",", row.names = FALSE, col.names = TRUE)
        write.table(zscore_pred_stats_out[["overall_statout"]], paste0(subanalysis_outfilepath, "zval_pred/", "zval_overall_statout.csv"), sep = ",", row.names = FALSE, col.names = TRUE)
        
        # Third test - ML
        dir.create(paste0(subanalysis_outfilepath, "ml_pred/"), showWarnings = F, recursive = T)
        cohort_multi_classification_ml_analysis_out <- cohort_multi_classification_ml_analysis(
            featuretable = ss_intable, outcometable = true_label_table, outcomelevels = prediction_factors,
            seedparam=11111, presplitdata = NULL, train_partition_percent = 0.6, subsampleparam = NULL, models_to_run = c("glmnet"))
        saveRDS(cohort_multi_classification_ml_analysis_out, paste0(subanalysis_outfilepath, "ml_pred/", "cohort_multi_classification_ml_analysis_out.RDS"))
        write.table(cohort_multi_classification_ml_analysis_out[["test_modelstats_outtable"]], paste0(subanalysis_outfilepath, "ml_pred/", "ml_test_confstatsout.csv"), sep = ",", row.names = FALSE, col.names = TRUE)
        write.table(cohort_multi_classification_ml_analysis_out[["test_overallstat_outtable"]], paste0(subanalysis_outfilepath, "ml_pred/", "ml_test_overall_statout.csv"), sep = ",", row.names = FALSE, col.names = TRUE)
        
        # I want every accuracy from every combo of: Score metric (TotalScore), genesetsize (50), prediction_label (rna_4cluster) and predction_method ("zval")
        pred_out_stat_list[[paste(ss_metric, genesetsize, ss_test_label,sep = "__")]] <- cbind(ss_metric = ss_metric, genesetsize = genesetsize, ss_test_label = ss_test_label,
            pred_test = c("rawval", "zval",
                          "ml"
                          ),
            Accuracy = c(raw_pred_stats_out[["overall_statout"]]["Accuracy",2], 
                         zscore_pred_stats_out[["overall_statout"]]["Accuracy",2],
                         cohort_multi_classification_ml_analysis_out[["test_overallstat_outtable"]]["Accuracy",2]
              )
        )
    }
}
pred_out_stat_table <- do.call(rbind, pred_out_stat_list)
pred_out_stat_table <- pred_out_stat_table[order(pred_out_stat_table[,"ss_test_label"], pred_out_stat_table[,"pred_test"], pred_out_stat_table[,"Accuracy"], decreasing = TRUE),]
write.table(pred_out_stat_table, paste0(outfilepathmaster, "signature_scoring/singscore_summaries/", "pred_stat_full_table.csv"),
            sep = ",", col.names = TRUE, row.names = FALSE)

# --------------------------------- cohort signature score vs events - confusion matrix optimization ---------------------------------
# dir.create(paste0(outfilepathmaster, "score_v_event_analysis/"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(outfilepathmaster, "score_v_event_analysis/CM_custom_cutoffs/"), showWarnings = FALSE, recursive = TRUE)

# Read in singscore results results
singscore_list <- c("rna_imgdegis_score1", "rna_ctndv70_score1", "rna_subtype1", "rna_subtype2", "rna_subtype3", "rna_subtype4",
                    "meth_imgdegis_score1", "meth_ctndv70_score1", "meth_subtype1", "meth_subtype2", "meth_subtype3")
singscore_table_list <- list()
singscore_path <- paste0(outfilepathmaster, "signature_scoring/")
for (singscore_result in singscore_list) {
    intable_singscore <- read.table(paste0(singscore_path, singscore_result, "/singscore/", singscore_result, "_singscore_outtable.csv"),
                                    sep = ",", header = TRUE)
    singscore_table_list[[singscore_result]] <- intable_singscore
}
combined_singscore_table <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Row.names", all = TRUE, sort = FALSE), 
                                   lapply(singscore_table_list, function(x) x[,c("Row.names", "TotalScore")]))
dimnames(combined_singscore_table) <- list(combined_singscore_table[,"Row.names"], c("Row.names", names(singscore_table_list)))
combined_singscore_table <- combined_singscore_table[,!grepl("Row.names", colnames(combined_singscore_table))]

# Create filtered event table
# Function for taking events with continuous time to event, and creating new columns filtered by day cutoffs - censoring appropriately
toe_filter_inlist <- list(table1 = bioreptable_waddons[,c("T_PRIMARY", "C_PRIMARY")],
                          table2 = bioreptable_waddons[,c("T_CVDMI_P", "C_CVDMI_P")])
time_filters <- c(182, 365, 548, 731, 913, 1096)
toe_table_outlist <- list()
for (toe_table in toe_filter_inlist) {
    filtered_toe_table <- filter_time_to_event_table(time_to_event_table=toe_table, time_filters=time_filters)
    toe_label <- gsub("T_", "", colnames(toe_table)[1])
    print(toe_label)
    toe_table_outlist[[toe_label]] <- filtered_toe_table
}
EOI_table <- do.call(cbind.data.frame, toe_table_outlist)
colnames(EOI_table) <- unlist(lapply(toe_table_outlist, colnames))


# Now combined the singscore table and the EOI_table
singscore_testtable <- merge(combined_singscore_table, EOI_table, by = "row.names", all = TRUE)
rownames(singscore_testtable) <- singscore_testtable[,"Row.names"]
param_table <- expand.grid(rna_score = c("rna_imgdegis_score1", "rna_ctndv70_score1", "rna_subtype1", "rna_subtype2", "rna_subtype3", "rna_subtype4"),
                           meth_score = c("meth_imgdegis_score1", "meth_ctndv70_score1", "meth_subtype1", "meth_subtype2", "meth_subtype3"),
                           # events = c(colnames(EOI_table[grepl("filtered", colnames(EOI_table))]))
                           events = c(colnames(EOI_table[grepl("filtered", colnames(EOI_table))]))
)
param_list <- split(param_table, seq(nrow(param_table)))

# Have to test for which go the OPPOSITE way in scoring - DO I HAVE TO FLIP ALL OF THE METH SCORES??!?!?!
out1 <- ggpairs(combined_singscore_table)

# Ok so for each pair of scores, we want to plot it out, and then see which have the best confusion matrix (and stats) at each diagonal cut of the plot (cutting at the Qs)
number_of_quantiles <- c(2,3,4,5,10)

allstatoutlist <- list()
for (quantile_sel in number_of_quantiles) {
    for (param_sel in param_list) {
        # select parameters
        rna_score_sel <- as.character(param_sel[["rna_score"]])
        meth_score_sel <- as.character(param_sel[["meth_score"]])
        events_sel <- as.character(param_sel[["events"]])
        test_label <- paste(as.character(unlist(param_sel)), collapse = "__")
        
        # create a plottable and plot the 2 scores with event labels
        plottable <- singscore_testtable[,c(rna_score_sel, meth_score_sel, events_sel)]
        plottable <- plottable[rowSums(is.na(plottable)) == 0 & plottable[,events_sel] != 2,]
        # midpoint axes calculation:
        midpoint_param <- apply(plottable[,c(1,2)], 2, function(x) {
            out1 = summary(x)[c("Min.", "Max.")]
            out1["Mid"] = c(out1["Max."] + out1["Min."]) / 2
            out1["Mid"]
        })
        pout <- scatter_plotter(indata = plottable[,c(rna_score_sel, meth_score_sel)], colorvar = plottable[,events_sel, drop=FALSE],
                                labsparam = list(title = test_label,
                                                 x = rna_score_sel, y = meth_score_sel), plotstats = FALSE)
        pout <- pout + geom_vline(xintercept = midpoint_param[1], linetype = 2) + geom_hline(yintercept = midpoint_param[2], linetype = 2)
        
        
        # Now use this table to cut into quantiles and create CM
        # continuous value to named quantile:
        
        # continuous_dataframe <- quartile_testtable[,"rna_subtype2",drop=FALSE]
        rna_quantile_table_out <- continuous_to_named_quantile(continuous_dataframe=plottable[,rna_score_sel,drop=FALSE],
                                                           number_of_quantiles=quantile_sel)
        rna_quantile_table <- rna_quantile_table_out[[1]]
        meth_quantile_table_out <- continuous_to_named_quantile(continuous_dataframe=plottable[,meth_score_sel,drop=FALSE],
                                                            number_of_quantiles=quantile_sel)
        meth_quantile_table <- meth_quantile_table_out[[1]]
        quantile_test_table <- cbind(rna_quantile_table, meth_quantile_table, plottable[rownames(rna_quantile_table),events_sel])
        
        
        apply(data.frame(as.character(rna_quantile_table_out[[3]])), 1, function(x) {
            out1 <- strsplit(gsub("\\[|\\]|\\)|\\(", "", x), split = ",")[[1]]
            out2 <- str_trim(out1[2], side = "both")
            out2
        })
        
        # pp_noevent <- table(quantile_test_table[quantile_test_table[,3] == 0,])
        # pp_event <- table(quantile_test_table[quantile_test_table[,3] == 1,])
        # pp_event / (pp_noevent + pp_event)
        
        # NEed to make a tablet to grab from - RNA on x axis like the plot
        pp_noevent <- as.data.frame.matrix(table(quantile_test_table[quantile_test_table[,3] == 0, c(meth_score_sel, rna_score_sel)]))
        pp_event <- as.data.frame.matrix(table(quantile_test_table[quantile_test_table[,3] == 1, c(meth_score_sel, rna_score_sel)]))
        
        # TOSH WORKING HERE - creating CM and then fisher tests or whatever to evaluate different cuts for different events based off of the scores
        quantile_selection_table <- expand.grid(seq_len(quantile_sel), seq_len(quantile_sel))
        quantile_selection_table[,"sum"] <- rowSums(quantile_selection_table)
        quantile_selection_table[,c(1,2)] <- apply(quantile_selection_table[,c(1,2)], 2, function(x) {
            out1 <- paste0("Q", x)
            out1 <- ifelse(out1 == "Q1", "Q1_Low", out1)
            out1 <- ifelse(out1 == paste0("Q", quantile_sel), paste0("Q", quantile_sel, "_High"), out1)
        })
        colnames(quantile_selection_table) <- c(colnames(quantile_test_table)[1:2], "sum")
        
        # Diagonal test
        statoutlist <- list()
        for (num_sel in seq_len(quantile_sel * 2)) {
            number_parameter <- (quantile_sel * 2) - num_sel + 1
            selected_rows_ingroup <- quantile_selection_table[quantile_selection_table[,"sum"] >= number_parameter, ]
            selected_rows_outgroup <- quantile_selection_table[quantile_selection_table[,"sum"] < number_parameter, ]
            
            if ((nrow(selected_rows_outgroup) == 0) | (nrow(selected_rows_ingroup) == 0)) {next}
            
            # Sum all the events and non events for the ingroup
            ingroup_event <- sum(apply(selected_rows_ingroup, 1, function(x){pp_event[unname(x[1]), unname(x[2])]}))
            ingroup_noevent <- sum(apply(selected_rows_ingroup, 1, function(x){pp_noevent[unname(x[1]), unname(x[2])]}))
            
            outgroup_event <- sum(apply(selected_rows_outgroup, 1, function(x){pp_event[unname(x[1]), unname(x[2])]}))
            outgroup_noevent <- sum(apply(selected_rows_outgroup, 1, function(x){pp_noevent[unname(x[1]), unname(x[2])]}))
            
            # Now with this - create a fisher table/confusion matrix and calculate stats
            # CM_table <- matrix(c(outgroup_noevent, ingroup_noevent, outgroup_event, ingroup_event), nrow=2, dimnames = list(c("outgroup", "ingroup"), c("noevent", "event")))
            CM_table <- matrix(c(outgroup_noevent, ingroup_noevent, outgroup_event, ingroup_event), nrow=2)
            CM_statout <- writeout_confmat(confusionMatrix(CM_table))
            rownames(CM_statout)[(nrow(CM_statout)-3):nrow(CM_statout)] <- c("outgroup_noevent", "ingroup_noevent", "outgroup_event", "ingroup_event")
            fisher_statout <- fisher.test(CM_table)
            fisher_writeout <- data.frame(val = c(fisher_pval = fisher_statout[["p.value"]], fisher_OR = unname(fisher_statout[["estimate"]]),
                                                  fisher_CI = fisher_statout[["conf.int"]]))
            
            # Descriptive table
            testrun_description <- data.frame(val = c(score1 = rna_score_sel, score2 = meth_score_sel, event = events_sel,
                                              cutpoint = paste0("diagonal_", number_parameter),
                                              quantile = paste0("quantile_", quantile_sel)))
            
            # Combine and write out the whole thing
            fullstatwriteout <- t(rbind(testrun_description, CM_statout, fisher_writeout))
            statoutlist[[num_sel]] <- fullstatwriteout
        }
        statouttable <- do.call(rbind, statoutlist)
        allstatoutlist[[test_label]] <- statouttable
        
        allstatouttable <- do.call(rbind, allstatoutlist)
        write.table(allstatouttable, 
                    paste0(outfilepathmaster, "score_v_event_analysis/CM_custom_cutoffs/", "quantile", quantile_sel, "_tempout.csv"),
                    sep = ",", col.names = TRUE, row.names = FALSE)
    }
}

# apply(data.frame(paste0(outfilepathmaster, "score_v_event_analysis/CM_custom_cutoffs/", list.files(paste0(outfilepathmaster, "score_v_event_analysis/CM_custom_cutoffs/")))), 1, function(x) read.table(x, header = TRUE, sep = ","))


# ## TOSH WORKING HERE
# # Plot out the one we care about:
# quantile_sel <- 4
# rna_score_sel <- "rna_subtype4"
# meth_score_sel <- "meth_subtype2"
# events_sel <- "filtered548_PRIMARY"
# test_label <- "test"
# 
# # rna_subtype4	meth_subtype2	filtered548_PRIMARY
# # create a plottable and plot the 2 scores with event labels
# plottable <- singscore_testtable[,c(rna_score_sel, meth_score_sel, events_sel)]
# plottable <- plottable[rowSums(is.na(plottable)) == 0 & plottable[,events_sel] != 2,]
# # midpoint axes calculation:
# midpoint_param <- apply(plottable[,c(1,2)], 2, function(x) {
#     out1 = summary(x)[c("Min.", "Max.")]
#     out1["Mid"] = c(out1["Max."] + out1["Min."]) / 2
#     out1["Q1"] = c(out1["Mid"] + out1["Min."]) / 2
#     out1["Q3"] = c(out1["Mid"] + out1["Max."]) / 2
#     out1
# })
# pout <- scatter_plotter(indata = plottable[,c(rna_score_sel, meth_score_sel)], colorvar = plottable[,events_sel, drop=FALSE],
#                         labsparam = list(title = test_label,
#                                          x = rna_score_sel, y = meth_score_sel), plotstats = FALSE)
# pout <- pout + geom_vline(xintercept = midpoint_param[c("Q1", "Mid", "Q3"),1], linetype = 2) + 
#                geom_hline(yintercept = midpoint_param[c("Q1", "Mid", "Q3"),2], linetype = 2)
# 
# ycoords <- c(midpoint_param[c("Q3", "Mid"),"meth_subtype2"])
# xcoords <- c(midpoint_param[c("Mid", "Q3"),"rna_subtype4"])
# slope <- diff(ycoords)/diff(xcoords)
# intercept <- ycoords[1]-slope*xcoords[1]
# 
# pout + geom_abline(slope,intercept)
# pout     
#       

      
# --------------------------------- alluvial and subtype-combo analyses ---------------------------------
dir.create(paste0(outfilepathmaster, "subtype-combo/"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(outfilepathmaster, "subtype-combo/", "alluvials/"), showWarnings = FALSE, recursive = TRUE)

# Create alluvial table
meth_rna_crosstab <- as.data.frame.matrix(table(combined_clustermembership_table[,c("rna_4cluster", "meth_3cluster")]))
## Chi sq value for independence test of the MS and RS variables
chisq.test(meth_rna_crosstab)


alluvial_intable <- data.frame(melt(cbind(rownames(meth_rna_crosstab), meth_rna_crosstab)))
colnames(alluvial_intable) <- c("rna_cluster", "meth_cluster", "value")

library(ggalluvial)

is_alluvia_form(alluvial_intable, axes = 1:2, silent = TRUE)

pout <- ggplot(alluvial_intable, aes(y = value, axis1 = rna_cluster, axis2 = meth_cluster))
pout <- pout + geom_alluvium(aes(fill = rna_cluster), width = 1/4)
pout <- pout + geom_stratum(aes(fill = rna_cluster), width = 1/4, color = "grey")
pout <- pout + geom_label(stat = "stratum", aes(label = after_stat(stratum)))
pout <- pout + scale_x_discrete(limits = c("RNA cluster", "Meth cluster"), expand = c(.05, .05))
pout



# try the lode form - which is more adjustable and can allow for more speicifc control
alluvial_long <- to_lodes_form(alluvial_intable, key = "Cluster", value = "Group", id = "Cohort", axes = 1:2)

pout <- ggplot(data = alluvial_long, aes(x = Cluster, stratum = Group, alluvium = Cohort, y = value))
pout <- pout + geom_alluvium(aes(fill=Group), width = 1/4)
pout <- pout + geom_stratum(aes(fill=Group), width = 1/4)
pout <- pout + geom_text(stat = "stratum", aes(label = Group))
pout <- pout + theme_minimal()
# pout <- pout + scale_fill_manual(values = c())#### TOSH ADD COLOR GUIDE

# manual coloring
alluvial_colors <- unlist(custom_annotation_list_from_colorguide(COI = c("rna_4cluster", "meth_3cluster"), colorguide))
names(alluvial_colors) <- gsub(".*\\.", "", names(alluvial_colors))
pout <- pout + scale_fill_manual(values = c(alluvial_colors), breaks = names(alluvial_colors))#### TOSH ADD COLOR GUIDE

pdf(paste0(outfilepathmaster, "subtype-combo/", "alluvials/", "cohort_alluvial_plot.pdf"), 15, 8, useDingbats = FALSE)
print(pout)
junk <- dev.off()


# Ok - now what about a simple TOTAL event analyses, so not with a time variable, but just overall
# meth_rna_crosstab - this is our denominator
EOI_vec <- c("CVDMI_P", "CVDMI_S", "PRIMARY")
EOI_table <- bioreptable_waddons[,grepl(paste(EOI_vec, collapse = "|"), colnames(bioreptable_waddons))]
timeparam <- c(182, 365, 548, 731)
for (EOI in EOI_vec) {
    for (timevar in timeparam) {
        EOI_table[,paste0("filtered", timevar, "_", EOI)] <- ifelse(EOI_table[, paste0("T_", EOI)] >= timevar & 
                                                                        EOI_table[, paste0("C_", EOI)] == 1, 2, # filter out too late event
                                                                    ifelse(EOI_table[, paste0("T_", EOI)] < timevar & 
                                                                               EOI_table[, paste0("C_", EOI)] == 0, 2, # filter out too soon non-event
                                                                           EOI_table[, paste0("C_", EOI)]))
    }
}
EOI_table <- EOI_table[,!grepl("^DT_|^T_", colnames(EOI_table))]
# Now for each event - do a ratio (and this will probably be a fisher test eventually...)
# Simple ratio plot first
for (EOI_sel in colnames(EOI_table)) {
    event_combotype_analysis_table <- merge(combined_clustermembership_table[,c("rna_4cluster", "meth_3cluster")],
                                            EOI_table[,EOI_sel,drop=FALSE], by = "row.names", all = TRUE)
    rownames(event_combotype_analysis_table) <- event_combotype_analysis_table[,"Row.names"]
    event_combotype_analysis_table <- na.omit(event_combotype_analysis_table[event_combotype_analysis_table[,EOI_sel] %in% c(0,1),
                                                                             !grepl("Row.names", colnames(event_combotype_analysis_table))])

    total_table <- dcast(event_combotype_analysis_table, rna_4cluster ~ meth_3cluster, value.var = EOI_sel, fun.aggregate = length)
    noevent_table <- dcast(event_combotype_analysis_table[event_combotype_analysis_table[,EOI_sel] == 0,],
                           rna_4cluster ~ meth_3cluster, value.var = EOI_sel, fun.aggregate = length)
    event_table <- dcast(event_combotype_analysis_table[event_combotype_analysis_table[,EOI_sel] == 1,],
                         rna_4cluster ~ meth_3cluster, value.var = EOI_sel, fun.aggregate = length)
    rownames(total_table) <- rownames(noevent_table) <- rownames(event_table) <- total_table[,1]
    

    event_ratio_table <- event_table[,!grepl("rna_4cluster", colnames(event_table))] / total_table[,!grepl("rna_4cluster", colnames(total_table))]
    noevent_ratio_table <- noevent_table[,!grepl("rna_4cluster", colnames(noevent_table))] / total_table[,!grepl("rna_4cluster", colnames(total_table))]
}


# TOSH - need to make a column with each of the combos for labels, and then do a fisher test 
# "rna_4cluster_w3AB", "rna_4cluster", "meth_4cluster", "meth_3cluster"
# tabulation_table <- bioreptable_waddons[, c("PATNUM", metaCOI)]
# out1 <- ingroup_vs_outgroup_cohort_enrichment_tests(tabulation_table, group_column = group_column_selected,
#                                                     cohort_tabulations_outpath = COI_tabulation_subpath)



# --------------------------------- modules with subtype combos ---------------------------------
dir.create(paste0(outfilepathmaster, "subtype-combo/", "eigen_heatmaps/"), showWarnings = FALSE, recursive = TRUE)

# RNA - WGCNA
eigengenecounttable_file <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run3_rmoutliers2/WGCNA/WGCNA_power14_size30/wgcna_eigengenes.csv"
eigengenecounttable <- read.table(eigengenecounttable_file, sep = ",", header = TRUE, row.names = 1)
rna_eigengeneOI <- paste0("ME", c("magenta", "salmon", "purple", "midnightblue", "turquoise", 
                                  "yellow", "pink", "red", 
                                  "blue", "tan", "brown",
                                  "greenyellow", "black", "cyan", "green"))
meth_eigengenecounttable <- read.table(paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/integrative_analyses/run7/", "meth_WGCNA/", "meth_eigengene_nogrey_table.csv"),sep = ",", header = TRUE, row.names = 1)
meth_eigengeneOI <- paste0("ME", c("yellow", "cyan", "tan", "turquoise", "black", 
                                   "magenta", "pink", "brown",
                                   "blue", "red", "salmon", "green", "purple"))
combo_subtypelabel_table <- data.frame(combo_subtype = apply(na.omit(combined_clustermembership_table[,c("rna_4cluster", "meth_3cluster")]), 1, 
                                                             function(x) paste0(x[1], "__", x[2])))

meth_eigen_seltable <- meth_eigengenecounttable[rownames(combo_subtypelabel_table), meth_eigengeneOI]
colnames(meth_eigen_seltable) <- paste0("meth_", colnames(meth_eigen_seltable))
rna_eigen_seltable <- eigengenecounttable[rownames(combo_subtypelabel_table), rna_eigengeneOI]
colnames(rna_eigen_seltable) <- paste0("rna_", colnames(rna_eigen_seltable))

combo_eigentable <- cbind(combo_subtypelabel_table, rna_eigen_seltable, meth_eigen_seltable)
combo_agg_eigentable <- aggregate(combo_eigentable[,!colnames(combo_eigentable) %in% "combo_subtype"], by = list(combo_eigentable[,"combo_subtype"]), mean)

# I want a few heatmaps
## (1) RNA eigengenes with both subtypes as col annots (RNA override, then subbed by meth)
## (1A) - side annotation with full ranking of the module across all RNAxMeth subtyps
## (1B) - side annotation where I only do ranking within the RNAsubtype

## (2) Meth eigengenes with both subtypes as col annots (Meth override, then subbed by RNA)
## (2A) - side annotation with full ranking of the module across all RNAxMeth subtyps
## (2B) - side annotation where I only do ranking within the Methsubtype

# Ok - really the hard part here is building the side annotation ... the rest is really not so bad
combo_eigen_plot_inlist <- list(
    run1 = list(eigen_sel = colnames(combo_eigentable)[grepl("rna_", colnames(combo_eigentable))], 
                primary_col_annot = "rna_4cluster", secondary_col_annot = "meth_3cluster", row_annot_type = "fullrank",
                outlabel = "rnaprim_rnaeigen_fullrankannot"),
    run2 = list(eigen_sel = colnames(combo_eigentable)[grepl("rna_", colnames(combo_eigentable))], 
                primary_col_annot = "rna_4cluster", secondary_col_annot = "meth_3cluster", row_annot_type = "subrank",
                outlabel = "rnaprim_rnaeigen_subrankannot"),
    run3 = list(eigen_sel = colnames(combo_eigentable)[grepl("meth_", colnames(combo_eigentable))], 
                primary_col_annot = "meth_3cluster", secondary_col_annot = "rna_4cluster", row_annot_type = "fullrank",
                outlabel = "methprim_metheigen_fullrankannot"),
    run4 = list(eigen_sel = colnames(combo_eigentable)[grepl("meth_", colnames(combo_eigentable))], 
                primary_col_annot = "meth_3cluster", secondary_col_annot = "rna_4cluster", row_annot_type = "subrank",
                outlabel = "methprim_metheigen_subrankannot"),
    run5 = list(eigen_sel = colnames(combo_eigentable)[grepl("rna_|meth_", colnames(combo_eigentable))], 
                primary_col_annot = "rna_4cluster", secondary_col_annot = "meth_3cluster", row_annot_type = "fullrank",
                outlabel = "rnaprim_alleigen_fullrankannot")
)

for (combo_eigen_plotparam in combo_eigen_plot_inlist) {
    ## Grab parameters
    eigen_sel = combo_eigen_plotparam[["eigen_sel"]]
    primary_col_annot = combo_eigen_plotparam[["primary_col_annot"]]
    secondary_col_annot = combo_eigen_plotparam[["secondary_col_annot"]]
    row_annot_type = combo_eigen_plotparam[["row_annot_type"]]
    outlabel_sel = combo_eigen_plotparam[["outlabel"]]
    
    ## Create outstring and outpath:
    ### Should be something like rnaprim_rnaeigen_fullrankannot
    outpath_sel <- paste0(outfilepathmaster, "subtype-combo/", "eigen_heatmaps/", outlabel_sel, "/")
    dir.create(outpath_sel, showWarnings = FALSE, recursive = TRUE)
    
    ## Create plottable
    hmplottable <- t(combo_eigentable[,eigen_sel])
    plotSOI <- colnames(hmplottable)
    
    ## Based off the primary and secondary col annotation - we will want different orders of the row annot
    if (primary_col_annot == "rna_4cluster") {
        combo_subtype_order = c("nmf_cluster_1__meth_3cluster_1", "nmf_cluster_1__meth_3cluster_2", "nmf_cluster_1__meth_3cluster_3",
                                "nmf_cluster_2__meth_3cluster_1", "nmf_cluster_2__meth_3cluster_2", "nmf_cluster_2__meth_3cluster_3", 
                                "nmf_cluster_3__meth_3cluster_1", "nmf_cluster_3__meth_3cluster_2", "nmf_cluster_3__meth_3cluster_3", 
                                "nmf_cluster_4__meth_3cluster_1", "nmf_cluster_4__meth_3cluster_2", "nmf_cluster_4__meth_3cluster_3")
    }
    if (primary_col_annot == "meth_3cluster") {
        combo_subtype_order = c("nmf_cluster_1__meth_3cluster_1", "nmf_cluster_2__meth_3cluster_1", "nmf_cluster_3__meth_3cluster_1", "nmf_cluster_4__meth_3cluster_1",
                                "nmf_cluster_1__meth_3cluster_2", "nmf_cluster_2__meth_3cluster_2", "nmf_cluster_3__meth_3cluster_2", "nmf_cluster_4__meth_3cluster_2",
                                "nmf_cluster_1__meth_3cluster_3", "nmf_cluster_2__meth_3cluster_3", "nmf_cluster_3__meth_3cluster_3", "nmf_cluster_4__meth_3cluster_3")
    }
    
    ## col annotation is just the cluster labels (easy)
    colmetatable <- bioreptable_waddons[plotSOI, c(primary_col_annot, secondary_col_annot), drop = FALSE]
    colsplitparam <- data.frame(factor(combo_eigentable[plotSOI, "combo_subtype"], levels = combo_subtype_order))
    dimnames(colsplitparam) <- list(plotSOI, "combo_subtype")
    colannotationlist <- custom_annotation_list_from_colorguide(COI = c(primary_col_annot, secondary_col_annot), colorguide)
    
    ## row annotation is much harder:
    # But i also want the rank table by cluster, so i need a box plot style melted table, and then go from there
    cluster_eigenrank_table_value <- aggregate(combo_eigentable[,eigen_sel], by = list(combo_eigentable[,"combo_subtype"]), mean)
    rownames(cluster_eigenrank_table_value) <- cluster_eigenrank_table_value[,"Group.1"]
    cluster_eigenrank_table_value <- cluster_eigenrank_table_value[,!grepl("Group.1", colnames(cluster_eigenrank_table_value))]
    
    ## Now with this - what kind of rowannot do we want?
    if (row_annot_type == "fullrank") {
        cluster_eigenrank_table_rank <- data.frame(t(data.frame(apply(cluster_eigenrank_table_value, 2, function(x) rev(rank(x))))[combo_subtype_order,]))
    }
    if (row_annot_type == "subrank") {
        cluster_eigenrank_table_rank <- data.frame(t(data.frame(apply(cluster_eigenrank_table_value, 2, function(x) rev(rank(x))))[combo_subtype_order,]))
        
        cluster_eigenrank_table_rank[,"associated_subtype_for_eigengene"] <- 
            ifelse(rownames(cluster_eigenrank_table_rank) %in% paste0("rna_ME", c("magenta", "salmon", "purple", "midnightblue", "turquoise")), "nmf_cluster_1",
            ifelse(rownames(cluster_eigenrank_table_rank) %in% paste0("rna_ME", c("yellow", "pink", "red")), "nmf_cluster_2",
            ifelse(rownames(cluster_eigenrank_table_rank) %in% paste0("rna_ME", c("blue", "tan", "brown")), "nmf_cluster_3",
            ifelse(rownames(cluster_eigenrank_table_rank) %in% paste0("rna_ME", c("greenyellow", "black", "cyan", "green")), "nmf_cluster_4",
            ifelse(rownames(cluster_eigenrank_table_rank) %in% paste0("meth_ME", c("yellow", "cyan", "tan", "turquoise", "black")), "meth_3cluster_1",
            ifelse(rownames(cluster_eigenrank_table_rank) %in% paste0("meth_ME", c("magenta", "pink", "brown")), "meth_3cluster_2",
            ifelse(rownames(cluster_eigenrank_table_rank) %in% paste0("meth_ME", c("blue", "red", "salmon", "green", "purple")), "meth_3cluster_3",
                   NA)))))))
    
        cluster_eigenrank_table_rank_temp <- data.frame(t(apply(cluster_eigenrank_table_rank, 1, function(x) {
            # x = cluster_eigenrank_table_rank[1,]
            associated_subtype_for_eigengene <- x["associated_subtype_for_eigengene"]
            x[!grepl(associated_subtype_for_eigengene, names(x))] <- NA
            rank(as.numeric(x), na.last = "keep")
        })))
        colnames(cluster_eigenrank_table_rank_temp) <- colnames(cluster_eigenrank_table_rank)
        cluster_eigenrank_table_rank <- cluster_eigenrank_table_rank_temp
    }
    rowmetatable <- cbind(eigengene_color = gsub("rna_ME|meth_ME", "", rownames(cluster_eigenrank_table_rank)), cluster_eigenrank_table_rank[, combo_subtype_order])
    
    
    ## grab the table for tbe stacked bar plot analysis
    grabbed_table_for_barplot <- cluster_eigenrank_table_rank
    # grabbed_table_for_barplot <- data.frame(t(cluster_eigenrank_table_value))
    
    ## We need the fullrank table - but make the others NA for subrank
    if (row_annot_type == "subrank") {
        grabbed_table_for_barplot <- data.frame(t(data.frame(apply(cluster_eigenrank_table_value, 2, function(x) rev(rank(x))))[combo_subtype_order,]))
        grabbed_table_for_barplot[is.na(cluster_eigenrank_table_rank[,colnames(grabbed_table_for_barplot)])] <- NA
    }
    ## We need the fullrank table, but make the others NA for subrank
    
    
    ## With this - create a stacked bar for each column - aka each pairing of subtypes
    barplot_intable <- melt(cbind(
        eigengene = rownames(grabbed_table_for_barplot[,!colnames(grabbed_table_for_barplot) %in% "associated_subtype_for_eigengene"]),
        grabbed_table_for_barplot[,!colnames(grabbed_table_for_barplot) %in% "associated_subtype_for_eigengene"]), 
        id.vars = "eigengene")
    barplot_intable[,"eigengene"] <- factor(barplot_intable[,"eigengene"], levels = rownames(grabbed_table_for_barplot))
    barplot_intable[,"associated_subtype_for_eigengene"] <- 
        ifelse(rownames(grabbed_table_for_barplot) %in% paste0("rna_ME", c("magenta", "salmon", "purple", "midnightblue", "turquoise")), "nmf_cluster_1",
        ifelse(rownames(grabbed_table_for_barplot) %in% paste0("rna_ME", c("yellow", "pink", "red")), "nmf_cluster_2",
        ifelse(rownames(grabbed_table_for_barplot) %in% paste0("rna_ME", c("blue", "tan", "brown")), "nmf_cluster_3",
        ifelse(rownames(grabbed_table_for_barplot) %in% paste0("rna_ME", c("greenyellow", "black", "cyan", "green")), "nmf_cluster_4",
        ifelse(rownames(grabbed_table_for_barplot) %in% paste0("meth_ME", c("yellow", "cyan", "tan", "turquoise", "black")), "meth_3cluster_1",
        ifelse(rownames(grabbed_table_for_barplot) %in% paste0("meth_ME", c("magenta", "pink", "brown")), "meth_3cluster_2",
        ifelse(rownames(grabbed_table_for_barplot) %in% paste0("meth_ME", c("blue", "red", "salmon", "green", "purple")), "meth_3cluster_3",
                                                         NA)))))))
    
    # I have to scale this maually here for it to work
    grab_maxcolval <- max(aggregate(barplot_intable[,"value"], by = list(barplot_intable[,"variable"]), function(x) sum(x, na.rm = TRUE))[,2], na.rm = TRUE)
    barplot_intable <- do.call(rbind, lapply(split(barplot_intable, barplot_intable[,"variable"]), function(x) {
        scalefactor <- grab_maxcolval/sum(x[,"value"], na.rm = TRUE)
        x[,"scaled_value"] <- x[,"value"] * scalefactor
        x
    }))
    
    # If full rank - scale up, if subrank, leave as raw
    if (row_annot_type == "fullrank") {
        barplot_intable <- do.call(rbind, lapply(split(barplot_intable, barplot_intable[,"variable"]), function(x) {
            scalefactor <- grab_maxcolval/sum(x[,"value"], na.rm = TRUE)
            x[,"scaled_value"] <- x[,"value"] * scalefactor
            x
        }))
    }
    if (row_annot_type == "subrank") {
        barplot_intable[,"scaled_value"] <- barplot_intable[,"value"]
    }
    
    # Add in the color
    barplot_intable[,"eigengene_color"] <- factor(gsub("meth_ME|rna_ME", "", barplot_intable[,"eigengene"]),
                                                  levels = unique(gsub("meth_ME|rna_ME", "", rownames(grabbed_table_for_barplot))))
    
    # Create a separate barplot for outlines
    barplot_intable_subtypeagg <- aggregate(barplot_intable[,"scaled_value",drop=FALSE],
                                            by = barplot_intable[,c("variable", "associated_subtype_for_eigengene")], sum)
    barplot_intable_subtypeagg[,"eigengene_color"] <- NA
    barplot_intable_subtypeagg[,"eigengene"] <- NA

    pout <- ggplot(barplot_intable, aes(x = variable, y = scaled_value, fill = eigengene_color))
    pout <- pout + geom_bar(position = "stack", stat = "identity", alpha = 0.3)
    pout <- pout + scale_fill_identity()
    pout <- pout + geom_bar(data = barplot_intable_subtypeagg, mapping = aes(x = variable, y = scaled_value, color = associated_subtype_for_eigengene),
                            position = "stack", stat = "identity", alpha = 1, size = 2)
    
    
    # pout <- ggplot(barplot_intable, aes(x = variable, y = scaled_value, fill = eigengene))
    # pout <- pout + geom_bar(position = "stack", stat = "identity", alpha = 0.3)
    # pout <- pout + scale_fill_manual(breaks = as.character(barplot_intable[barplot_intable[,"variable"] == unique(barplot_intable[,"variable"])[1],"eigengene"]),
    #                                  values = as.character(barplot_intable[barplot_intable[,"variable"] == unique(barplot_intable[,"variable"])[1],"eigengene_color"]))
    # pout <- pout + geom_bar(data = barplot_intable_subtypeagg, mapping = aes(x = variable, y = scaled_value, color = associated_subtype_for_eigengene),
    #                         position = "stack", stat = "identity", alpha = 1, size = 2)
    # pout
    
    ## Add outlining colors
    outlinecolorvector <- c()
    if (sum(grepl("nmf_cluster", unique(barplot_intable_subtypeagg[,"associated_subtype_for_eigengene"]))) > 0) {
        outlinecolorvector <- custom_annotation_list_from_colorguide(COI = "rna_4cluster", colorguide)[[1]]
    }
    if (sum(grepl("meth_3cluster", unique(barplot_intable_subtypeagg[,"associated_subtype_for_eigengene"]))) > 0) {
        outlinecolorvector <- c(outlinecolorvector, custom_annotation_list_from_colorguide(COI = "meth_3cluster", colorguide)[[1]])
    }
    pout <- pout + scale_color_manual(breaks = names(outlinecolorvector), values = unname(outlinecolorvector))
    ## Add outlining colors
    
    pout <- pout + labs(title = outlabel_sel, x = "module_pair", y = "module_rank_weight") + theme_pubr(x.text.angle = 90)
    pout
    
    pdf(paste0(outpath_sel, outlabel_sel, "_stacked_module_plot.pdf"), 15, 10)
    print(pout)
    junk <- dev.off()
    
    
    ## Create a row annootation from this:
    blackwhite_color_range <- colorRampPalette(c("#e5e5e5", "#191919"))
    rowcustomcolorlist <- c(
        list(setNames(rowmetatable[,"eigengene_color"], rowmetatable[,"eigengene_color"])),
        rep(list(colorRamp2(
                    breaks = seq(1, max(unique(unlist(rowmetatable[, combo_subtype_order])), na.rm = TRUE)),
                    colors = blackwhite_color_range(max(unique(unlist(rowmetatable[, combo_subtype_order])), na.rm = TRUE))
            )), length(combo_subtype_order))
    )
    names(rowcustomcolorlist) <- colnames(rowmetatable)
    rowannotationlist <- annotationlist_builder(rowmetatable, customcolorlist = rowcustomcolorlist)
    
    ## Now with all of thse pieces - make our heatmap
    heatmapcolorparam <- colorRamp2(c(-0.1, 0, 0.1), c("blue", "white", "red"))
    outhm_wsplit <- create_heatmap(counttab = hmplottable, scale_data = FALSE, separate_legend = TRUE,
                                colmetatable = colmetatable, colannotationlist = colannotationlist,
                                rowmetatable = rowmetatable, rowannotationlist = rowannotationlist,
                                colclusterparam = TRUE, rowclusterparam = FALSE,
                                heatmapcolorparam = heatmapcolorparam,
                                columnsplitparam = colsplitparam
    )
    pdf(paste0(outpath_sel, outlabel_sel, "_hm.pdf"), 15, 10)
    draw(outhm_wsplit[[1]])
    junk <- dev.off()
    write.table(hmplottable, paste0(outpath_sel, outlabel_sel, "hm_table.csv"), sep = ",", col.names = NA, row.names = TRUE)
    write.table(rowmetatable, paste0(outpath_sel, outlabel_sel, "hm_eigenrank_table.csv"), sep = ",", col.names = NA, row.names = TRUE)
    
    # Unclustered heatmap
    outhm_nosplit <- create_heatmap(counttab = hmplottable, scale_data = FALSE, separate_legend = TRUE,
                                colmetatable = colmetatable, colannotationlist = colannotationlist,
                                rowmetatable = rowmetatable, rowannotationlist = rowannotationlist,
                                colclusterparam = TRUE, rowclusterparam = FALSE,
                                heatmapcolorparam = heatmapcolorparam,
    )
    pdf(paste0(outpath_sel, outlabel_sel, "_unsplit_hm.pdf"), 15, 10)
    draw(outhm_nosplit[[1]])
    junk <- dev.off()
    
    outhm_colclust <- as.hclust(column_dend(outhm_nosplit[[1]]))
    multi_subtype_clusterlabel_table <- data.frame(multi_subtype_clusterlabel = cutree(outhm_colclust, k = seq(2,10)))
    
    ## Add in the subtypes and do some tests to look for enrichment:
    subtype_label_table <- cbind(do.call(rbind, strsplit(combo_eigentable[,"combo_subtype"], split = "__")), combo_eigentable[,"combo_subtype"])
    dimnames(subtype_label_table) <- list(rownames(combo_eigentable), c("rna_4cluster", "meth_3cluster", "combo_subtype"))
    
    outtable_list <- list()
    for (k in 2:10) {
        dend_tab_table <- cbind(PATNUM = rownames(multi_subtype_clusterlabel_table), 
                                multi_subtype_clusterlabel_table[,grepl(k, colnames(multi_subtype_clusterlabel_table)),drop=FALSE], subtype_label_table)
        out1 <- summarize_table(dend_tab_table,
                                groupvar = colnames(multi_subtype_clusterlabel_table[,grepl(k, colnames(multi_subtype_clusterlabel_table)), drop=FALSE]),
                                outfile = paste0(outpath_sel, outlabel_sel, "_statdummy.csv"), calc_stats = TRUE)
        outtable_list[[k]] <- do.call(rbind, out1[2:length(out1)])
    }
    ## TOSH FINISH THIS THOUGHT HERE FOR WORK!
    multisubtype_per_subtype_stattable <- do.call(cbind, Filter(Negate(is.null), outtable_list))
    write.table(multisubtype_per_subtype_stattable, paste0(outpath_sel, outlabel_sel, "hm_table.csv"), sep = ",", col.names = NA, row.names = TRUE)
    
    ## column clustering
    # pdf(paste0(outpath_sel, outlabel_sel, "_col_dendrogram.pdf"), 15, 10)
    # plot(outhm_colclust)
    # junk <- dev.off()
    write.table(multi_subtype_clusterlabel_table, paste0(outpath_sel, outlabel_sel, "multi_subtype_clusterlabel_table.csv"),
                sep = ",", col.names = NA, row.names = TRUE)
    pdf(paste0(outpath_sel, outlabel_sel, "_col_dendrogram.pdf"), 15, 10)
    for (k in 1:10) {
        dend_out <- color_branches(outhm_colclust, k = k)
        plot(dend_out)
    }
    junk <- dev.off()
    
}



# --------------------------------- subtype-combo analyses scatters ---------------------------------
dir.create(paste0(outfilepathmaster, "subtype-combo/", "diff_comp_plots/"), showWarnings = FALSE, recursive = TRUE)

# Read in the deseq tables that we need
deseq_rna_CTNDV70_3v_v_1v_table <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run3b_rmoutliers2_addoncomps/deseq/comp_anatomy70__3v_v_1v/deseq_results_comp_anatomy70__3v_v_1v.csv", sep = ",", header = TRUE, row.names = 1)
deseq_rna_IMGDEGIS_Sev_v_MildNone_table <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run3b_rmoutliers2_addoncomps/deseq/comp_ischemia__Sev_v_MildNone/deseq_results_comp_ischemia__Sev_v_MildNone.csv", sep = ",", header = TRUE, row.names = 1)
deseq_rna_RS1_v_notRS1_table <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run3_rmoutliers2/deseq/comp_RNAsubtype_RNAtype1_vs_NOTRNAtype1/deseq_results_comp_RNAsubtype_RNAtype1_vs_NOTRNAtype1.csv", sep = ",", header = TRUE, row.names = 1)
deseq_rna_RS2_v_notRS2_table <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run3_rmoutliers2/deseq/comp_RNAsubtype_RNAtype2_vs_NOTRNAtype2/deseq_results_comp_RNAsubtype_RNAtype2_vs_NOTRNAtype2.csv", sep = ",", header = TRUE, row.names = 1)
deseq_rna_RS3_v_notRS3_table <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run3_rmoutliers2/deseq/comp_RNAsubtype_RNAtype3_vs_NOTRNAtype3/deseq_results_comp_RNAsubtype_RNAtype3_vs_NOTRNAtype3.csv", sep = ",", header = TRUE, row.names = 1)
deseq_rna_RS4_v_notRS4_table <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run3_rmoutliers2/deseq/comp_RNAsubtype_RNAtype4_vs_NOTRNAtype4/deseq_results_comp_RNAsubtype_RNAtype4_vs_NOTRNAtype4.csv", sep = ",", header = TRUE, row.names = 1)


DMP_meth_CTNDV70_3v_v_1v_table <- readRDS("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/methylation_data_20220608/DMP/DMP_on_CTNDV70_3vs1.rds")
DMP_meth_degrisch_sev_nomild_table <- readRDS("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/methylation_data_20220608/DMP/DMP_on_degrisch_sev_nomild.rds")
colnames(DMP_meth_CTNDV70_3v_v_1v_table)[colnames(DMP_meth_CTNDV70_3v_v_1v_table) == "logFC"] <- "signed_deltaBeta"
colnames(DMP_meth_degrisch_sev_nomild_table)[colnames(DMP_meth_degrisch_sev_nomild_table) == "logFC"] <- "signed_deltaBeta"
meth_subtype_DMP_object <- readRDS(paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/methylation_data_20220608/DMP/DMP_k3_limma_noCov_oneVSother.rds"))
cleaned_meth_subtype_DMP_object <- clean_meth_subtype_DMP_object_function(meth_subtype_DMP_object)
DMP_meth_MS1_v_notMS1_table <- cleaned_meth_subtype_DMP_object[["not_to_1"]]    
DMP_meth_MS2_v_notMS2_table <- cleaned_meth_subtype_DMP_object[["2_to_not"]]    
DMP_meth_MS3_v_notMS3_table <- cleaned_meth_subtype_DMP_object[["not_to_3"]]
## signed_deltaBeta - THIS IS THE EFFECT SIZE VARIABLE THAT WE NEED


multiomic_comp_list <- list(comp1 = list(deseq_rna_CTNDV70_3v_v_1v_table = deseq_rna_CTNDV70_3v_v_1v_table,
                                         DMP_meth_CTNDV70_3v_v_1v_table = DMP_meth_CTNDV70_3v_v_1v_table),
                            comp2 = list(deseq_rna_IMGDEGIS_Sev_v_MildNone_table = deseq_rna_IMGDEGIS_Sev_v_MildNone_table,
                                         DMP_meth_degrisch_sev_nomild_table = DMP_meth_degrisch_sev_nomild_table),
                            comp3 = list(deseq_rna_RS1_v_notRS1_table = deseq_rna_RS1_v_notRS1_table,
                                         DMP_meth_MS2_v_notMS2_table = DMP_meth_MS2_v_notMS2_table),
                            comp4 = list(deseq_rna_RS2_v_notRS2_table = deseq_rna_RS2_v_notRS2_table,
                                         DMP_meth_MS3_v_notMS3_table = DMP_meth_MS3_v_notMS3_table)
                            )
rna_dge_cutoffs <- c(pval_test = "padj", pval_cutoff = 0.05, log2fc_cutoff = 0.25)
meth_dmp_cutoffs <- c(pval_test = "adj.P.Val", pval_cutoff = 0.05, log2fc_cutoff = 0.03)

for (multiomic_comp_num in seq_along(multiomic_comp_list)) {
    comp_sel <- multiomic_comp_list[[multiomic_comp_num]]
    deseq_table_1 <- comp_sel[[1]]
    deseq_table_1_label <- names(comp_sel)[1]
    deseq_table_2 <- comp_sel[[2]]
    deseq_table_2_label <- names(comp_sel)[2]
    comp_label <- paste0(paste(strsplit(deseq_table_1_label, split = "_")[[1]][c(2,3)], collapse = "_"), "_v_", 
                         paste(strsplit(deseq_table_2_label, split = "_")[[1]][c(2,3)], collapse = "_"))
    compplot_outfilepath <- paste0(outfilepathmaster, "subtype-combo/", "diff_comp_plots/", comp_label, "/")
    dir.create(compplot_outfilepath, showWarnings = FALSE, recursive = TRUE)
    
    deseq_table_1_filt <- deseq_table_1[deseq_table_1[,"pvalue"] < 0.05,]
    deseq_table_2_filt <- deseq_table_2[deseq_table_2[,"P.Value"] < 0.05,]
    # & abs(deseq_table_1[,"log2FoldChange"]) > 0.25
    # & abs(deseq_table_2[,"signed_deltaBeta"]) > 0.03
    
    comptable <- merge(deseq_table_1_filt[,c("log2FoldChange", "pvalue")], deseq_table_2_filt[,c("gene", "signed_deltaBeta", "P.Value")],
                       by.x = "row.names", by.y = "gene", all = TRUE)
    comptable2 <- na.omit(comptable[!comptable[,1] %in% "",])
    dim(comptable2)
    
    agg_probe_counttable <- aggregate(comptable2[,c("signed_deltaBeta", "log2FoldChange")], by = list(comptable2[,"Row.names"]), mean)
    agg_probe_counttable <- agg_probe_counttable[abs(agg_probe_counttable[,"signed_deltaBeta"]) > 0.02 &
                                                 abs(agg_probe_counttable[,"log2FoldChange"]) > 0.25,]
    
    pout <- scatter_plotter(agg_probe_counttable[,c("signed_deltaBeta", "log2FoldChange")],
                    labsparam = list(title = comp_label, x = "METH - signed_deltaBeta", y = "RNA - log2FoldChange"))
    pout <- pout + geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, linetype = 2)
    pdf(paste0(compplot_outfilepath, "rna_v_meth_logfc_scatter.pdf"))
    print(pout)
    junk <- dev.off()
    
    # We will inevitably have to save out gsea from our group of interest (top left, decrease meth -> increase expression)
    speciesparam <- "Homo sapiens"
    ULgenes <- agg_probe_counttable[agg_probe_counttable[, "signed_deltaBeta"] < 0 & agg_probe_counttable[, "log2FoldChange"] > 0,1]
    hypergeo_genetest_out_UL = data.frame(hypergeo_genetest(data.frame(ULgenes), statcutoffparam = "DUMMY", 
                                                            genesetparam = c("C5"), speciesparam = speciesparam)$enricherUPout)
    write.table(data.frame(ULgenes), paste0(compplot_outfilepath, "UL_genes.csv"), sep = ",", col.names = FALSE, row.names = FALSE)
    write.table(hypergeo_genetest_out_UL, paste0(compplot_outfilepath, "UL_hypergeo.csv"), sep = ",", col.names = TRUE, row.names = FALSE)
    
}


# How is this possible - the genes are pathways are the same that are up in our protective and vulnerable groups... wtf
dir.create(paste0(compplot_outfilepath, "protect_v_risk/"), showWarnings = FALSE, recursive = TRUE)
protect_genes <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run3_rmoutliers2/integrative_analyses/run6/subtype-combo/diff_comp_plots/rna_RS1_v_meth_MS2/UL_genes.csv", sep = ",")[,1]
risk_genes <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run3_rmoutliers2/integrative_analyses/run6/subtype-combo/diff_comp_plots/rna_RS2_v_meth_MS3/UL_genes.csv", sep = ",")[,1]

overlapinlist <- list(
    protect_genes = protect_genes,
    risk_genes = risk_genes
)
overlapout <- overlap_finder(overlapinlist)
overlaptable <- overlapout$overlaptable
vennplot <- overlapout$vennplot
overlapgrouptab <- overlapout$overlapgrouptab
colnames(overlapgrouptab)[1] <- "Row.names"
overlapsummary <- overlapout$overlapsummary

dir.create(paste0(outfilepathmaster, "subtype-combo/", "protect_v_risk/", "gene_comp/"), showWarnings = FALSE, recursive = TRUE)
write.table(overlaptable, paste0(outfilepathmaster, "subtype-combo/", "protect_v_risk/", "gene_comp/", "testoverlap_overlaptable.csv"),
            sep = ",", col.names = TRUE, row.names = FALSE)
write.table(overlapgrouptab, paste0(outfilepathmaster, "subtype-combo/", "protect_v_risk/", "gene_comp/", "testoverlap_overlapgrouptable.csv"),
            sep = ",", col.names = TRUE, row.names = FALSE)
write.table(overlapsummary, paste0(outfilepathmaster, "subtype-combo/", "protect_v_risk/", "gene_comp/", "testoverlap_overlapsummary.csv"),
            sep = ",", col.names = TRUE, row.names = TRUE)
pdf(paste0(outfilepathmaster, "subtype-combo/", "protect_v_risk/", "gene_comp/", "testoverlap_overlapvenn.pdf"))
grid.draw(vennplot)
junk <- dev.off()


# What about pathways
dir.create(paste0(outfilepathmaster, "subtype-combo/", "protect_v_risk/", "gsea_comp/"), showWarnings = FALSE, recursive = TRUE)
protect_gsea <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run3_rmoutliers2/integrative_analyses/run6/subtype-combo/diff_comp_plots/rna_RS1_v_meth_MS2/UL_hypergeo.csv", sep = ",", header = TRUE, row.names = 1)
risk_gsea <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run3_rmoutliers2/integrative_analyses/run6/subtype-combo/diff_comp_plots/rna_RS2_v_meth_MS3/UL_hypergeo.csv", sep = ",", header = TRUE, row.names = 1)

gsea_comp_table <- merge(protect_gsea[,"p.adjust",drop=FALSE], risk_gsea[,"p.adjust",drop=FALSE],
                         by = "row.names", suffixes = c("_protect", "_risk"), all = TRUE)
gsea_comp_table[is.na(gsea_comp_table)] <- 1  # recode NA values to a pvalue of 1
gsea_comp_table[,c("pstat_protect", "pstat_risk")] <- -log10(gsea_comp_table[,c("p.adjust_protect", "p.adjust_risk")])
gsea_comp_table <- gsea_comp_table[(gsea_comp_table[,"p.adjust_protect"] < 0.05) | (gsea_comp_table[,"p.adjust_risk"] < 0.05), ]

# Add some helper columns
gsea_comp_table[,"quadrant"] <- ifelse(gsea_comp_table[,"p.adjust_protect"] < 0.05 & gsea_comp_table[,"p.adjust_risk"] < 0.05, "both", 
                                ifelse(gsea_comp_table[,"p.adjust_protect"] < 0.05 & gsea_comp_table[,"p.adjust_risk"] > 0.05, "protect",
                                ifelse(gsea_comp_table[,"p.adjust_protect"] > 0.05 & gsea_comp_table[,"p.adjust_risk"] < 0.05, "risk", NA)))

gsea_comp_table[,"doublesig"] <- ifelse(gsea_comp_table[,"p.adjust_protect"] < 0.05 & gsea_comp_table[,"p.adjust_risk"] < 0.05, "both", 
                                 ifelse(gsea_comp_table[,"p.adjust_protect"] < 0.05 & gsea_comp_table[,"p.adjust_risk"] > 0.2, "protect",
                                 ifelse(gsea_comp_table[,"p.adjust_protect"] > 0.2 & gsea_comp_table[,"p.adjust_risk"] < 0.05, "risk", NA)))

write.table(gsea_comp_table, paste0(outfilepathmaster, "subtype-combo/", "protect_v_risk/", "gsea_comp/", "gsea_comp_table.csv"),
            sep = ",", col.names = TRUE, row.names = FALSE)

pout <- scatter_plotter(gsea_comp_table[, c("pstat_protect", "pstat_risk")], labsparam = list(title = "protect gsea vs risk gsea", x = "pstat_protect", y = "pstat_risk"))
pout <- pout + geom_hline(yintercept = -log10(0.05), linetype = 2) + geom_vline(xintercept = -log10(0.05), linetype = 2)
pdf(paste0(outfilepathmaster, "subtype-combo/", "protect_v_risk/", "gsea_comp/", "gsea_stat_scatter.pdf"))
print(pout)
junk <- dev.off()





# Clinical enrichment for comobs
COI_tabulation_subpath <- paste0(outfilepathmaster, "subtype-combo/", "clinical_enrichment_plots/")
dir.create(COI_tabulation_subpath, showWarnings = FALSE, recursive = TRUE)

rna_4cluster__meth_3cluster_COMBO_table <- data.frame(rna_4cluster__meth_3cluster_COMBO = apply(
    na.omit(bioreptable_waddons[,c("rna_4cluster", "meth_3cluster")]), 1, function(x) paste(x, collapse = "__")))
tabulation_table <- merge(bioreptable_waddons[, c("PATNUM", metaCOI)], rna_4cluster__meth_3cluster_COMBO_table,
                          by.x = "PATNUM", by.y = "row.names", all.x = TRUE)
rownames(tabulation_table) <- tabulation_table[,"PATNUM"]

out1 <- ingroup_vs_outgroup_cohort_enrichment_tests(tabulation_table, group_column = "rna_4cluster__meth_3cluster_COMBO",
                                                    cohort_tabulations_outpath = COI_tabulation_subpath)
out1_ratiotable <- out1[[1]]
out1_citable <- out1[[2]]


heatmapcolorparam <- colorRamp2(breaks = c(0, 0.99999, 1, 1.00001, 5), c("#000099", "#ccccff", "white", "#ffb2b2", "#b20000"))
rowmetatable <- data.frame(category = unlist(lapply(strsplit(rownames(out1_ratiotable), split = "__"), function(x) x[1])),
                           row.names = rownames(out1_ratiotable))[COI[COI %in% rownames(out1_ratiotable)],,drop=FALSE]
rowannotationlist <- annotationlist_builder(rowmetatable) 

# Add Ns to the col labels
N_annot_table <- table(tabulation_table[,"rna_4cluster__meth_3cluster_COMBO"])[colnames(out1_ratiotable)]
colnames_w_N <- paste0(names(N_annot_table), paste0("__N=", N_annot_table))
plottable <- out1_ratiotable[COI[COI %in% rownames(out1_ratiotable)],]
colnames(plottable) <- colnames_w_N

hm1 <- create_heatmap(counttab = plottable, subsetnum = FALSE, scale_data = FALSE,
                      rowmetatable = rowmetatable, rowannotationlist = rowannotationlist,
                      separate_legend = TRUE, heatmapcolorparam = heatmapcolorparam, addborders = TRUE)
pdf(paste0(COI_tabulation_subpath, "signed_filtered_sigfeature_hm.pdf"), 10, 10, useDingbats = FALSE)
draw(hm1[[1]])
junk <- dev.off()

write.table(out1_ratiotable, paste0(COI_tabulation_subpath, "signed_filtered_sigfeature_ratio_table.csv"), col.names = NA, row.names = TRUE, sep = ",")
write.table(out1_citable, paste0(COI_tabulation_subpath, "signed_filtered_sigfeature_ci_table.csv"), col.names = NA, row.names = TRUE, sep = ",")


gene_probe_corrtable <- readRDS("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/ze_analysis/rna_meth_cor/r_value_output/integrative_correlation_all_group_1_50000.rds")




# --------------------------------- cluster vs event summary heatmaps - combined omics type --------------------------------- 
combo_cluster_event_subpath <- paste0(outfilepathmaster, "subtype-combo/", "combocluster_survival/")
dir.create(combo_cluster_event_subpath, showWarnings = FALSE, recursive = TRUE)
# First we are going to add a combo here:
combined_clustermembership_table[,"rna_4cluster__AND__meth_3cluster"] <- apply(combined_clustermembership_table, 1, function(x) {
    out <- ifelse((!is.na(x["rna_4cluster"]) & !is.na(x["meth_3cluster"])), paste(c(x["rna_4cluster"],  x["meth_3cluster"]), collapse = "__AND__"), NA)
    out
})

subtype_combo_cluster_table <- combined_clustermembership_table[,c("PATNUM", "rna_4cluster__AND__meth_3cluster")]
subtype_combo_cluster_table[,"rna_4cluster_nmf1_AND_allmeth"] <- ifelse(grepl("nmf_cluster_1", combined_clustermembership_table[,"rna_4cluster__AND__meth_3cluster"]), combined_clustermembership_table[,"rna_4cluster__AND__meth_3cluster"], NA)
subtype_combo_cluster_table[,"rna_4cluster_nmf2_AND_allmeth"] <- ifelse(grepl("nmf_cluster_2", combined_clustermembership_table[,"rna_4cluster__AND__meth_3cluster"]), combined_clustermembership_table[,"rna_4cluster__AND__meth_3cluster"], NA)
subtype_combo_cluster_table[,"rna_4cluster_nmf3_AND_allmeth"] <- ifelse(grepl("nmf_cluster_3", combined_clustermembership_table[,"rna_4cluster__AND__meth_3cluster"]), combined_clustermembership_table[,"rna_4cluster__AND__meth_3cluster"], NA)
subtype_combo_cluster_table[,"rna_4cluster_nmf4_AND_allmeth"] <- ifelse(grepl("nmf_cluster_4", combined_clustermembership_table[,"rna_4cluster__AND__meth_3cluster"]), combined_clustermembership_table[,"rna_4cluster__AND__meth_3cluster"], NA)
subtype_combo_cluster_table[,"meth_3cluster_meth1_AND_allnmf"] <- ifelse(grepl("meth_3cluster_1", combined_clustermembership_table[,"rna_4cluster__AND__meth_3cluster"]), combined_clustermembership_table[,"rna_4cluster__AND__meth_3cluster"], NA)
subtype_combo_cluster_table[,"meth_3cluster_meth2_AND_allnmf"] <- ifelse(grepl("meth_3cluster_2", combined_clustermembership_table[,"rna_4cluster__AND__meth_3cluster"]), combined_clustermembership_table[,"rna_4cluster__AND__meth_3cluster"], NA)
subtype_combo_cluster_table[,"meth_3cluster_meth3_AND_allnmf"] <- ifelse(grepl("meth_3cluster_3", combined_clustermembership_table[,"rna_4cluster__AND__meth_3cluster"]), combined_clustermembership_table[,"rna_4cluster__AND__meth_3cluster"], NA)


bioreptable_waddons <- merge(bioreptable_waddons, combined_clustermembership_table[,c("PATNUM", "rna_4cluster__AND__meth_3cluster"),drop=FALSE], by = "PATNUM", all.x = TRUE)
rownames(bioreptable_waddons) <- bioreptable_waddons[,"PATNUM"]

bioreptable_waddons <- merge(bioreptable_waddons, subtype_combo_cluster_table[c(1,3:ncol(subtype_combo_cluster_table))],
                             by = "PATNUM", all.x = TRUE)
rownames(bioreptable_waddons) <- bioreptable_waddons[,"PATNUM"]
# So we want RNA1 M1 vs not RNA1 M1 IN RNA1, and then just R1 across all Ms
# Then RNA M2 vs not R1 M2, and then R1 across all Ms
# THEN the vice versa - of start with M and break down by R


# Adding in here - our hierarchical clustering results with both eigengenes:
eigengene_combo_hierarchical_label_table <- read.table(paste0(outfilepathmaster, "subtype-combo/eigen_heatmaps/rnaprim_alleigen_fullrankannot/rnaprim_alleigen_fullrankannothm_table.csv"), sep = ",", header = TRUE, row.names = 1)
eigengene_combo_hierarchical_label_table_sel <- cbind(
    PATNUM = rownames(eigengene_combo_hierarchical_label_table), 
    comboeigen_hier4 = paste0("comboeigen_hier4_", eigengene_combo_hierarchical_label_table[,"multi_subtype_clusterlabel.4"]),
    comboeigen_hier5 = paste0("comboeigen_hier5_", eigengene_combo_hierarchical_label_table[,"multi_subtype_clusterlabel.5"]),
    comboeigen_hier7 = paste0("comboeigen_hier7_", eigengene_combo_hierarchical_label_table[,"multi_subtype_clusterlabel.7"])
)

## Im gonna remove the outliers here for survival analysis
eigengene_combo_hierarchical_label_table_sel[,2:ncol(eigengene_combo_hierarchical_label_table_sel)] <- apply(
    eigengene_combo_hierarchical_label_table_sel[,2:ncol(eigengene_combo_hierarchical_label_table_sel)], 2, function(x) {
        subtypecount <- table(x)
        blankout_labels <- names(subtypecount[subtypecount < 5])
        x[x %in% blankout_labels] <- NA
        x
    })
## Im gonna remove the outliers here for survival analysis

subtype_combo_cluster_table <- merge(subtype_combo_cluster_table, eigengene_combo_hierarchical_label_table_sel, by = "PATNUM", all = TRUE)
rownames(subtype_combo_cluster_table) <- subtype_combo_cluster_table[,"PATNUM"]


## In group vs outgroup analysis enrichment of characteristics
group_columns_for_testing <- colnames(subtype_combo_cluster_table)[2:ncol(subtype_combo_cluster_table)]
# Want to do a manual addition just for this where we do a 3A vs 3B comparison
# group_columns_for_testing <- c(group_columns_for_testing, "rna_3A_v_3B_only")


for (group_column_selected in group_columns_for_testing) {
    
    # Want to do a manual addition just for this where we do a 3A vs 3B comparison
    # if (group_column_selected == "rna_3A_v_3B_only") {
    #     selected_cluster_table <- combined_clustermembership_table[, c("PATNUM", "rna_4cluster_w3AB")]
    #     selected_cluster_table[,"rna_4cluster_w3AB"] <- ifelse(selected_cluster_table[,"rna_4cluster_w3AB"] %in% c("nmf_cluster_3A", "nmf_cluster_3B"),
    #                                                            selected_cluster_table[,"rna_4cluster_w3AB"], NA)
    #     colnames(selected_cluster_table) <- c("PATNUM", "rna_3A_v_3B_only")
    # } else {
    selected_cluster_table <- subtype_combo_cluster_table[, c("PATNUM", group_column_selected)]
    # }
    selected_cluster_table <- selected_cluster_table[order(selected_cluster_table[,group_column_selected]), ]
    
    # need to be only pairwise for the summary to work correctly for now...
    highlowage <- continuous_to_named_quantile(bioreptable_waddons[,"AGE_RAND",drop=FALSE], 3)[[1]]
    highlowage[,1] <- ifelse(highlowage[,1] == "Q3_High", "AGE_RAND_highQ3", "AGE_RAND_low")
    
    COIref_tablist <- list(
        # IMGDEGIS = bioreptable_waddons[,"IMGDEGIS", drop=FALSE],
        # IMGDEGIS_01_3 = bioreptable_waddons[,"IMGDEGIS_01_3", drop=FALSE],
        # CTNDV50 = bioreptable_waddons[,"CTNDV50", drop=FALSE],
        # CTNDV50_13_clean = bioreptable_waddons[,"CTNDV50_13_clean", drop=FALSE],
        CTNDV70_13_clean = bioreptable_waddons[,"CTNDV70_13_clean", drop=FALSE]
        # CKD = bioreptable_waddons[,"CKD",drop=FALSE],
        # highlowage = highlowage
    )
    
    plotCOI <- c(names(COIref_tablist),
                 paste0(group_column_selected, "_",
                        unique(na.omit(selected_cluster_table[,group_column_selected]))[order(unique(na.omit(selected_cluster_table[,group_column_selected])))])
    )
    
    survivalplot_outfilepath <- paste0(combo_cluster_event_subpath, group_column_selected, "/")
    dir.create(survivalplot_outfilepath, showWarnings = FALSE, recursive = TRUE)
    # summary_heatmap <- cluster_vs_event_kmanalysis_summary_heatmap(selected_cluster_table = selected_cluster_table,
    #                                                                COIref_tablist = COIref_tablist,
    #                                                                survivalplot_outfilepath = survivalplot_outfilepath, 
    #                                                                bioreptable_waddons = bioreptable_waddons,
    #                                                                eventpvalcutoff = 0.05, plotCOI)
    summary_heatmap <- cluster_vs_event_kmanalysis_summary_heatmap(selected_cluster_table = selected_cluster_table,
                                                                   COIref_tablist = COIref_tablist,
                                                                   survivalplot_outfilepath = survivalplot_outfilepath, 
                                                                   bioreptable_waddons = bioreptable_waddons,
                                                                   eventpvalcutoff = 0.05, return_output_tables = FALSE, plotCOI)
    pdf(paste0(combo_cluster_event_subpath, group_column_selected, "/", group_column_selected, "_event_summary_heatmap.pdf"),
        12, 10, useDingbats = FALSE)
    draw(summary_heatmap)
    junk <- dev.off()
    
    
    # I also want to do the KM analysis with all of the groups together for each group
    survivalplot_outfilepath <- paste0(combo_cluster_event_subpath, group_column_selected, "_NONBINARY/")
    dir.create(survivalplot_outfilepath, showWarnings = FALSE, recursive = TRUE)
    COIref_tablist <- list(
        # CUSTOM_IMGDEGIS_NONEMILD = bioreptable_waddons[,"CUSTOM_IMGDEGIS_NONEMILD", drop=FALSE],
        CTNDV70 = bioreptable_waddons[,"CTNDV70", drop=FALSE]
        # CKD = bioreptable_waddons[,"CKD",drop=FALSE],
        # age_cat = bioreptable_waddons[,"CAGE_RND",drop=FALSE]
    )
    plotCOI <- c(names(COIref_tablist), group_column_selected)
    coxph_control_table <- factorize_metatable(bioreptable_waddons[,c("CAGE_RND", "SEX", "RACE")])
    # summary_heatmap <- cluster_NONBINARY_vs_event_kmanalysis_summary_heatmap(group_column_selected = group_column_selected,
    #                                                                          COIref_tablist = COIref_tablist,
    #                                                                          survivalplot_outfilepath = survivalplot_outfilepath, 
    #                                                                          bioreptable_waddons = bioreptable_waddons_special,
    #                                                                          eventpvalcutoff = 0.05, plotCOI)
    summary_heatmap <- cluster_NONBINARY_vs_event_kmanalysis_summary_heatmap(selected_cluster_table = selected_cluster_table,
                                                                             COIref_tablist = COIref_tablist,
                                                                             survivalplot_outfilepath = survivalplot_outfilepath,
                                                                             bioreptable_waddons = bioreptable_waddons,
                                                                             eventpvalcutoff = 0.05,
                                                                             return_pairwise_coxph_values = TRUE, coxph_control_table = coxph_control_table,
                                                                             plotCOI)
    pdf(paste0(combo_cluster_event_subpath, group_column_selected, "_NONBINARY/", group_column_selected, "_event_summary_heatmap.pdf"),
        12, 10, useDingbats = FALSE)
    draw(summary_heatmap)
    junk <- dev.off()
    
}


# --------------------------------- Deconvolution analysis ---------------------------------
dir.create(paste0(outfilepathmaster, "deconvolution_analysis/"), showWarnings = FALSE, recursive = TRUE)

deconv_alg <- c("mcp_counter", "epic", "quantiseq", "xcell", "cibersort", "cibersort_abs")
deconv_outtable_list <- deconvolution_analysis(normcounttable, deconv_alg = deconv_alg)

cleantable_all_outlist <- list()
# Now we want immune composition differences for some groups of interest:
compcols <- c("IMGDEGIS_01_2_3", "CTNDV70_0123_clean", "rna_4cluster", "meth_3cluster")
cleantable_comp_outlist <- list()
for (comp_sel in compcols) {
# for (deconv_table_num in seq_len(length(deconv_outtable_list))) {

    # dir.create(paste0(outfilepathmaster, "deconvolution_analysis/", comp_sel, "/"), showWarnings = FALSE, recursive = TRUE)
    # # Now we want immune composition differences for some groups of interest:
    # compcols <- c("IMGDEGIS_01_2_3", "CTNDV70_0123_clean", "rna_4cluster", "meth_3cluster")
    # cleantable_comp_outlist <- list()
    # for (comp_sel in compcols) {
    for (deconv_table_num in seq_len(length(deconv_outtable_list))) {

        deconv_table_sel <- deconv_outtable_list[[deconv_table_num]]
        alg_sel <- names(deconv_outtable_list)[deconv_table_num]
        
        dir.create(paste0(outfilepathmaster, "deconvolution_analysis/", comp_sel, "/", alg_sel, "/"), showWarnings = FALSE, recursive = TRUE)
        
        # Lets get some by group comparisons:
        deconv_comptable <- data.frame(merge(bioreptable_waddons, deconv_table_sel, by.x = "PATNUM", by.y = "row.names"))
        write.table(deconv_table_sel, paste0(outfilepathmaster, "deconvolution_analysis/", comp_sel, "/", alg_sel, "/", alg_sel, "_outtable.csv"),
                    sep = ",", col.names = NA, row.names = TRUE)
        write.table(deconv_comptable, paste0(outfilepathmaster, "deconvolution_analysis/", comp_sel, "/", alg_sel, "/", alg_sel, "_outtable_wbiorep.csv"),
                    sep = ",", col.names = NA, row.names = TRUE)
        
        tabulation_table <- deconv_comptable[,c("PATNUM", comp_sel, colnames(deconv_table_sel))]
        # Run summary to get the anova across all groups for the immune composition
        dir.create(paste0(outfilepathmaster, "deconvolution_analysis/", comp_sel, "/", alg_sel, "/", alg_sel, "_", comp_sel),
                   showWarnings = FALSE, recursive = TRUE)
        outfile <- paste0(outfilepathmaster, "deconvolution_analysis/", comp_sel, "/", alg_sel, "/", alg_sel, "_", comp_sel, "_statstable.csv")
        summary_stattable_all <- summarize_table(tabulation_table, groupvar = comp_sel, outfile = outfile,
                                                 calc_stats = TRUE, calc_ci = FALSE)
        deconv_cleantable_all <- clean_summarize_table_output(
            sumstattable_input = summary_stattable_all, addpercents = "vertical", 
            contsummary = c("mean", "sd"), roundpvaldigits = 3)
        write.table(deconv_cleantable_all, paste0(outfilepathmaster, "deconvolution_analysis/", 
                                                  comp_sel, "/", alg_sel, "/", alg_sel, "_", comp_sel, "_cleanstatstable.csv"),
                    col.names = TRUE, row.names = FALSE, sep = ",")
        
        summarize_table_figures(intable = tabulation_table, groupvar = comp_sel, outfilepath = 
                                    paste0(outfilepathmaster, "deconvolution_analysis/", comp_sel, "/", alg_sel, "/", alg_sel, "_", comp_sel))
        
        cleantable_comp_outlist[[comp_sel]] <- melt(deconv_cleantable_all, id.vars = c("category", "feature"),
                                                    measure.vars = colnames(deconv_cleantable_all)[3:ncol(deconv_cleantable_all)])
    }
}
full_deconv_outtable <-do.call(cbind, deconv_outtable_list)
write.table(full_deconv_outtable, paste0(outfilepathmaster, "deconvolution_analysis/", "full_deconv_outtable.csv"),
            sep = ",", col.names = NA, row.names = TRUE)

outfile <- paste0(outfilepathmaster, "deconvolution_analysis/", "full_deconv_statstable.csv")
tabulation_table <- 
summary_stattable_all <- summarize_table(tabulation_table, groupvar = comp_sel, outfile = outfile,
                                         calc_stats = TRUE, calc_ci = FALSE)
deconv_cleantable_all <- clean_summarize_table_output(
    sumstattable_input = summary_stattable_all, addpercents = "vertical", 
    contsummary = c("mean", "sd"), roundpvaldigits = 3)
write.table(deconv_cleantable_all, paste0(outfilepathmaster, "deconvolution_analysis/", 
                                          alg_sel, "/", alg_sel, "_", comp_sel, "/", comp_sel, "_cleanstatstable.csv"),
            col.names = TRUE, row.names = FALSE, sep = ",")

# Should output this heatmap with the methlation and RNAseq annotations on columns, and the decon alg on the rows
deconv_plottable <- data.frame(t(full_deconv_outtable), check.names = FALSE)

colmetatable <- bioreptable_waddons[SOIrna, c("rna_4cluster", "meth_3cluster")]
colannotationlist <- custom_annotation_list_from_colorguide(COI = c("rna_4cluster", "meth_3cluster"), colorguide)

rowmetatable <- data.frame(deconv_alg = rownames(deconv_plottable), row.names = rownames(deconv_plottable))
rowmetatable[,"deconv_alg"] <- do.call(rbind, lapply(rownames(deconv_plottable), function(x) strsplit(x, split = "\\.")[[1]][1]))
rowannotationlist <- annotationlist_builder(metatable = rowmetatable)

heatmapcolorparam = colorRamp2(breaks = c(-5,0,5), colors = c("blue", "white", "red"))
outhm <- create_heatmap(counttab = deconv_plottable, scale_data = TRUE, heatmapcolorparam = heatmapcolorparam,
                        colmetatable = colmetatable, colannotationlist = colannotationlist,
                        rowmetatable = rowmetatable, rowannotationlist = rowannotationlist,
                        colclusterparam = TRUE, rowclusterparam = TRUE,
                        separate_legend = TRUE,
                        columnsplitparam = colmetatable[,"rna_4cluster"],
                        # rowsplitparam = rowmetatable[,"deconv_alg"]
                        )




# --------------------------------- Eigengene vs Event Analysis ---------------------------------
dir.create(paste0(outfilepathmaster, "eigengene_vs_event/"), showWarnings = FALSE, recursive = TRUE)
bioreptable_waddons <- read.table(paste0(outfilepathmaster, "bioreptable_waddons.csv"), sep = ",", header = TRUE, row.names = 1)

## We want RNA and methylation Eigengenes, can test all at once and just adjust population as necessary
meth_eigengenecounttable <- read.table(paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/integrative_analyses/run7/", "meth_WGCNA/", "meth_eigengene_nogrey_table.csv"), sep = ",", header = TRUE, row.names = 1)
meth_eigenGOI <- paste0("ME", c("tan", "cyan", "darkgreen", "yellow", "turquoise", "black",
  "brown", "pink", "magenta", "darkturquoise",
  "purple", "salmon", "blue", "green", "red"))
meth_eigengenecounttable <- meth_eigengenecounttable[,meth_eigenGOI]


eigengenecounttable <- read.table("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run3_rmoutliers2/WGCNA/WGCNA_power14_size30/wgcna_eigengenes.csv", sep = ",", header = TRUE, row.names = 1)


combo_eigengenecounttable <- merge(eigengenecounttable, meth_eigengenecounttable, all = TRUE, by = "row.names")
colnames(combo_eigengenecounttable) <- c("gene", paste0("rna_", colnames(eigengenecounttable)), paste0("meth_", colnames(meth_eigengenecounttable)))
rownames(combo_eigengenecounttable) <- combo_eigengenecounttable[,"gene"]

## (1) ML Eigengene vs event - with events at various time cut offs
dir.create(paste0(outfilepathmaster, "eigengene_vs_event/", "eigen_v_binaryevent_ML/"), showWarnings = FALSE, recursive = TRUE)

# Create filtered event table
# Function for taking events with continuous time to event, and creating new columns filtered by day cutoffs - censoring appropriately
toe_filter_inlist <- list(table1 = bioreptable_waddons[,c("T_PRIMARY", "C_PRIMARY")],
                          table2 = bioreptable_waddons[,c("T_CVDMI_P", "C_CVDMI_P")])
time_filters <- c(182, 365, 548, 731, 913, 1096)
toe_table_outlist <- list()
for (toe_table in toe_filter_inlist) {
    filtered_toe_table <- filter_time_to_event_table(time_to_event_table=toe_table, time_filters=time_filters)
    toe_label <- gsub("T_", "", colnames(toe_table)[1])
    print(toe_label)
    toe_table_outlist[[toe_label]] <- filtered_toe_table
}
EOI_table <- do.call(cbind.data.frame, toe_table_outlist)
colnames(EOI_table) <- unlist(lapply(toe_table_outlist, colnames))

# Now for each event - test each of our eigengenes separately, then all RNA, all METH, then ALL eigengenes
param_table <- expand.grid(eigengene = c(colnames(combo_eigengenecounttable)[!colnames(combo_eigengenecounttable) %in% "gene"], 
                           "rna_eigen", "meth_eigen", "all_eigen"),
                           events = c(colnames(EOI_table[grepl("filtered", colnames(EOI_table))]))
)
param_list <- split(param_table, seq(nrow(param_table)))
outstatcounter <- 1
outstat_list <- outstat_table_list <- list()
pdf(paste0(outfilepathmaster, "eigengene_vs_event/", "eigen_v_binaryevent_ML/", "eigen_v_binaryevent_ML_ROC_plots.pdf"))
for (param_sel in param_list) {
    # Select the data that we need
    eigen_sel <- as.character(param_sel[["eigengene"]])
    if ("rna_eigen" %in% eigen_sel) { eigen_sel <- colnames(combo_eigengenecounttable)[grepl("rna_", colnames(combo_eigengenecounttable))] }
    if ("meth_eigen" %in% eigen_sel) { eigen_sel <- colnames(combo_eigengenecounttable)[grepl("meth_", colnames(combo_eigengenecounttable))] }
    if ("all_eigen" %in% eigen_sel) { eigen_sel <- colnames(combo_eigengenecounttable)[!colnames(combo_eigengenecounttable) %in% "gene"] }
    event_sel <- as.character(param_sel[["events"]])
    
    ## Run our ML with the given input and ouput
    featuretable <- combo_eigengenecounttable[,eigen_sel,drop=FALSE]
    outcometable <- ifelse(EOI_table[,event_sel,drop=FALSE] == 1, event_sel, 
                    ifelse(EOI_table[,event_sel,drop=FALSE] == 0, "noevent", NA))
    outcomelevels <- c(event_sel, "noevent")

    # Then run our modeling with these genes and labels
    out1 <- cohort_classification_ml_analysis(featuretable=featuretable, outcometable=outcometable, outcomelevels=outcomelevels,
                                              seedparam=11111, subsampleparam = "down",
                                              # models_to_run = c("glm", "xgbTree", "svmRadial", "rf", "glmnet")
                                              models_to_run = c("glm", "svmRadial", "rf", "glmnet")
    )
    
    ## I really dont want to write out 400 stat tables... But I can store them and grab them as needed?
    outstat_table_list[[outstatcounter]] <- out1[["modelstats_outtable"]]
    
    ## Save the stats we really want
    outstat <- cbind( param_sel, out1[["modelstats_outtable"]][out1[["modelstats_outtable"]][,1] %in% "AUC_w_CI",])
    outstat_list[[outstatcounter]] <- outstat
    outstatcounter <- outstatcounter + 1
    
    print(out1[["rocplot"]] + ggtitle(paste0(as.character(param_sel[["eigengene"]]), " vs ", event_sel)))
    
}
junk <- dev.off()
outstat_table <- do.call(rbind, outstat_list)
write.table(outstat_table, paste0(outfilepathmaster, "eigengene_vs_event/", "eigen_v_binaryevent_ML/", "eigen_v_binaryevent_ML_stattable.csv"),
            row.names = FALSE, col.names = TRUE, sep = ",")
save.RDS(outstat_table_list, paste0(outfilepathmaster, "eigengene_vs_event/", "eigen_v_binaryevent_ML/", "eigen_v_binaryevent_ML_fullstattables.rds"))

outstat_table <- read.table(paste0(outfilepathmaster, "eigengene_vs_event/", "eigen_v_binaryevent_ML/", "eigen_v_binaryevent_ML_stattable.csv"), 
                            header = TRUE, sep = ",")

# Plotting out some stats
# plottable1 <- t(apply(outstat_table[outstat_table[,"events"] %in% "filtered182_PRIMARY", "rf_out",drop=FALSE], 1, function(x) strsplit(unlist(x), split = " ", fixed = TRUE)[[1]]))
# plottable2 <- cbind.data.frame(plottable1[,2], t(apply(plottable1[,3,drop=FALSE], 1, function(x) strsplit(gsub("\\(|\\)", "", x), split = "-")[[1]])))
# colnames(plottable2) <- c("mean", "ci_lower", "ci_upper")
# plottable2 <- apply(plottable2, 2, as.numeric)
# rownames(plottable2) <- outstat_table[outstat_table[,"events"] %in% "filtered182_PRIMARY","eigengene"]
# 
# plottable3 <- cbind.data.frame(eigengene = rownames(plottable2), plottable2[,1,drop=FALSE], errorbarvar = data.frame(error = c(plottable2[,1] - plottable2[,2])))
# plottable3 <- plottable3[order(plottable3[,2], decreasing = TRUE),]
# plottable3[,1] <- factor(plottable3[,1])
# 
# plot_barchart(plottable3[,c(1,2)], errorbarvar = plottable3[,3,drop=FALSE], 
#               labsparam = list(title = "rf stats for eigen event prediction", x = "eigengene", y = "AUC"))



## (2) COXPH for eigengenes vs event - all at once? indivudally with covariates? + 
## (3) Quantiled Eigengene vs event KM analysis - WE HAVE TO QUANTILE FOR COXPH ANYWAYS - SO WHY NOT?
dir.create(paste0(outfilepathmaster, "eigengene_vs_event/", "quanteigen_v_event_km_coxph/"), showWarnings = FALSE, recursive = TRUE)
# For every event, and for every eigengene - run a coxph, probably with some other variables?
EOI <- gsub("^C_", "", colnames(bioreptable_waddons)[grepl("^C_", colnames(bioreptable_waddons))])
EOIcols <- apply(expand.grid(c("T_", "C_"), EOI, stringsAsFactors = FALSE), 1, paste, collapse = "")

# Now for each event - test each of our eigengenes separately, then all RNA, all METH, then ALL eigengenes
param_table <- expand.grid(
    eigengene = c(colnames(combo_eigengenecounttable)[!colnames(combo_eigengenecounttable) %in% "gene"]),
    events = EOI,
    quantiles = c(2,3,4,5,10)
)

param_list <- split(param_table, seq(nrow(param_table)))
outstatcounter <- 1
outstat_list <- outstat_table_list <- list()
coxph_controltable <- factorize_metatable(bioreptable_waddons[,c("SEX", "RACE", "CAGE_RND")])
coxph_controltable[,"CAGE_RND"] <- factor(ifelse(as.character(coxph_controltable[,"CAGE_RND"]) %in% c("<= 34", "35 - 44", "45 - 54"),
                                          "<= 54", as.character(coxph_controltable[,"CAGE_RND"])),
                                          levels = c("<= 54", "55 - 64", "65 - 74", ">= 75"))
coxph_controltable[,"RACE"] <- factor(ifelse(as.character(coxph_controltable[,"RACE"]) %in% c("American Indian or Alaska Native", "Multiple Races", "Native Hawaiian or Other Pacific Islander", ""), "Other", as.character(coxph_controltable[,"RACE"])),
                                          levels = c("White", "Black or African American", "Asian", "Other"))
# pdf(paste0(outfilepathmaster, "eigengene_vs_event/", "eigen_v_binaryevent_ML/", "eigen_v_binaryevent_ML_ROC_plots.pdf"))
coxph_outstat_list <- km_outstat_list <- list()
outstatcounter <- 1
for (param_sel in param_list) {
    # Select the data that we need
    eigen_sel <- as.character(param_sel[["eigengene"]])
    # if ("rna_eigen" %in% eigen_sel) { eigen_sel <- colnames(combo_eigengenecounttable)[grepl("rna_", colnames(combo_eigengenecounttable))] }
    # if ("meth_eigen" %in% eigen_sel) { eigen_sel <- colnames(combo_eigengenecounttable)[grepl("meth_", colnames(combo_eigengenecounttable))] }
    # if ("all_eigen" %in% eigen_sel) { eigen_sel <- colnames(combo_eigengenecounttable)[!colnames(combo_eigengenecounttable) %in% "gene"] }
    event_sel <- as.character(param_sel[["events"]])
    quantile_sel <- param_sel[["quantiles"]]
    
    ## Run our ML with the given input and ouput
    quantiled_eigen_out <- continuous_to_named_quantile(combo_eigengenecounttable[,eigen_sel,drop=FALSE], number_of_quantiles = quantile_sel)
    eigentable <- quantiled_eigen_out$quantile_dataframe
    # eigentable <- cbind.data.frame(quantiled_eigen_out$quantile_dataframe, combo_eigengenecounttable[,eigen_sel,drop=FALSE])
    # colnames(eigentable) <- c(paste0(eigen_sel, "_quantile", quantile_sel), eigen_sel)
    
    outcometable <- bioreptable_waddons[,paste0(c("T_", "C_"), event_sel)]
    outcometable[,2] <- as.numeric(gsub(2, NA, outcometable[,2]))
    
    # coxph_controltable
    shared_samples <- Reduce(intersect, list(rownames(eigentable), rownames(outcometable), rownames(coxph_controltable)))
    survivaldata <- na.omit(cbind.data.frame(outcometable[shared_samples,], eigentable[shared_samples,,drop=FALSE],
                                           coxph_controltable[shared_samples,]))
    try(coxph_out <- coxph_analysis(survivaldata)) {}
    coxph_table <- coxph_out$coxph_outtable
    coxph_plot <- coxph_out$coxph_outplot
    
    ## Run the KM analysis
    survivalanalysisout <- create_survival_plot(survivaldata = survivaldata[,c(1:3)], timebreakparam = NULL, ylimitparam = c(0.25,1))
    outsurvtable <- survivalanalysisout$outsurvtable
    outsurvplot <- survivalanalysisout$outsurvplot
    outsurvpvalue <- survivalanalysisout$outsurvpvalue
    
    # Ok......... so now what do we save out?
    ## I think just the COXph and KM pvals...
    coxph_outstat_list[[outstatcounter]] <- cbind(event = event_sel, eigengene = eigen_sel, quantile = quantile_sel,
                                                  coxph_table[grepl(eigen_sel, coxph_table[,1]), c("Row.names", "Pr(>|z|)")])
    km_outstat_list[[outstatcounter]] <- cbind(event = event_sel, eigengene = eigen_sel, quantile = quantile_sel,
                                               km_pval = outsurvpvalue[1,"pval"])
    outstatcounter <- outstatcounter + 1
                                              
}
coxph_outstat_table <- do.call(rbind, coxph_outstat_list)
km_outstat_table <- do.call(rbind, km_outstat_list)

write.table(coxph_outstat_table, paste0(outfilepathmaster, "eigengene_vs_event/", "quanteigen_v_event_km_coxph/", "coxph_outstat_table.csv"), sep = ",", col.names = TRUE, row.names = FALSE)
write.table(km_outstat_table, paste0(outfilepathmaster, "eigengene_vs_event/", "quanteigen_v_event_km_coxph/", "km_outstat_table.csv"), sep = ",", col.names = TRUE, row.names = FALSE)






# --------------------------------- END ---------------------------------


# --------------------------------- cluster performance tests - DONT RUN --------------------------------- 
## Basically we want the confusion matrix for the classification of our samples based on their DGE
# so input is the DGE and the class labels, and output is the confusion matrix... pretty straightforward
# dir.create(paste0(outfilepathmaster, "cluster_performance_modeling"), showWarnings = FALSE, recursive = TRUE)
# dir.create(paste0(outfilepathmaster, "cluster_performance_modeling_downsample/"), showWarnings = FALSE, recursive = TRUE)
dge_comp_list <- c("comp_ischemia__Sev_v_MildNone", "comp_anatomy__3v_v_1v")
deseq_table_list <- list()
deseq_path <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run3_rmoutliers2/deseq/"
for (dge_comp in dge_comp_list) {
    intable_deseq <- read.table(paste0(deseq_path, dge_comp, "/deseq_results_", dge_comp, ".csv"), 
                                sep = ",", header = TRUE, row.names = 1)
    deseq_table_list[[dge_comp]] <- intable_deseq
    dir.create(paste0(outfilepathmaster, "cluster_performance_modeling_down/", dge_comp, "/"), showWarnings = FALSE, recursive = TRUE)
    
    # set stat params:
    if (dge_comp == "comp_ischemia__Sev_v_MildNone") {
        pval_test = "padj"
        pval_cutoff = 0.05
        log2fc_cutoff = 0.25
    }
    if (dge_comp == "comp_anatomy__3v_v_1v") {
        pval_test = "pvalue"
        pval_cutoff = 0.05
        log2fc_cutoff = 0.25
    }
    
    GOI <- rownames(intable_deseq[intable_deseq[,pval_test] < pval_cutoff & !is.na(intable_deseq[,pval_test]) & 
                                      abs(intable_deseq[,"log2FoldChange"]) > log2fc_cutoff,])
    write.table(data.frame(GOI), paste0(outfilepathmaster, "cluster_performance_modeling_down/", dge_comp, "/", dge_comp, "_modeling_features.txt"),
                sep = "\t", col.names = FALSE, row.names = FALSE)
    
    # Set up out featuretable and outcometable
    outcomelabel <- dge_comp
    strsplit_temp <- strsplit(dge_comp, split = "_")[[1]]
    outcomelevels <- strsplit_temp[c(length(strsplit_temp) - 2, length(strsplit_temp))]
    featurelabels <- GOI
    outcometable <- metatable[,outcomelabel,drop=FALSE]
    outcometable[,1] <- ifelse(outcometable[,1] == 1, outcomelevels[1], ifelse(outcometable[,1] == 0, outcomelevels[2], NA))
    featuretable <- t(normcounttable[featurelabels,rownames(metatable)])
    
    # Then run our modeling with these genes and labels
    out1 <- cohort_classification_ml_analysis(featuretable=featuretable, outcometable=outcometable, outcomelevels=outcomelevels,
                                              seedparam=11112, subsampleparam = "down",
                                              models_to_run = c("svmRadial", "rf", "glmnet")
    )
    
    # Write out the results
    write.table(out1$modelstats_outtable, paste0(outfilepathmaster, "cluster_performance_modeling_down/", dge_comp, "/", "model_stats_outtable.csv"),
                sep = ",", col.names = TRUE, row.names = FALSE)
    pdf(paste0(outfilepathmaster, "cluster_performance_modeling_down/", dge_comp, "/", "ROCcurves.pdf"), 10, 10, useDingbats = FALSE)
    print(out1$rocplot)
    junk <- dev.off()
    
    ## NEED TO LOOK INTO COHORT MODELING - CAUSE GETTING PRETTY GOOD AUCS, BUT ONLY CAUSE ITS CALL ALL ONE WAY
    # https://topepo.github.io/caret/subsampling-for-class-imbalances.html
    
}






# --------------------------------- cohort signature score vs events - ML (FAILED) ---------------------------------
dir.create(paste0(outfilepathmaster, "score_v_event_analysis/"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(outfilepathmaster, "score_v_event_analysis/tempout/"), showWarnings = FALSE, recursive = TRUE)
# Read in singscore results results
singscore_list <- c("rna_imgdegis_score1", "rna_ctndv70_score1", "rna_subtype1", "rna_subtype2", "rna_subtype3", "rna_subtype4",
                    "meth_imgdegis_score1", "meth_ctndv70_score1", "meth_subtype1", "meth_subtype2", "meth_subtype3")
singscore_table_list <- list()
singscore_path <- paste0(outfilepathmaster, "signature_scoring/")
for (singscore_result in singscore_list) {
    intable_singscore <- read.table(paste0(singscore_path, singscore_result, "/singscore/", singscore_result, "_singscore_outtable.csv"),
                                    sep = ",", header = TRUE)
    singscore_table_list[[singscore_result]] <- intable_singscore
}
combined_singscore_table <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Row.names", all = TRUE, sort = FALSE), 
                                   lapply(singscore_table_list, function(x) x[,c("Row.names", "TotalScore")]))
dimnames(combined_singscore_table) <- list(combined_singscore_table[,"Row.names"], c("Row.names", names(singscore_table_list)))
combined_singscore_table <- combined_singscore_table[,!grepl("Row.names", colnames(combined_singscore_table))]

EOI_vec <- c("CVDMI_P", "PRIMARY")
EOI_table <- bioreptable_waddons[,grepl(paste(EOI_vec, collapse = "|"), colnames(bioreptable_waddons))]
timeparam <- c(182, 365, 548, 731)
for (EOI in EOI_vec) {
    for (timevar in timeparam) {
        EOI_table[,paste0("filtered", timevar, "_", EOI)] <- ifelse(EOI_table[, paste0("T_", EOI)] >= timevar & 
                                                                        EOI_table[, paste0("C_", EOI)] == 1, 2, # filter out too late event
                                                                    ifelse(EOI_table[, paste0("T_", EOI)] < timevar & 
                                                                               EOI_table[, paste0("C_", EOI)] == 0, 2, # filter out too soon non-event
                                                                           EOI_table[, paste0("C_", EOI)]))
    }
}
# apply(EOI_table[,grepl("filtered", colnames(EOI_table))], 2, table)



singscore_testtable <- merge(combined_singscore_table, EOI_table, by = "row.names", all = TRUE)
rownames(singscore_testtable) <- singscore_testtable[,"Row.names"]
param_table <- expand.grid(rna_score = c("rna_imgdegis_score1", "rna_ctndv70_score1", "rna_subtype1", "rna_subtype2", "rna_subtype3", "rna_subtype4"),
                           meth_score = c("meth_imgdegis_score1", "meth_ctndv70_score1", "meth_subtype1", "meth_subtype2", "meth_subtype3"),
                           events = c(colnames(EOI_table[grepl("filtered", colnames(EOI_table))]))
                           # events = c(colnames(EOI_table[grepl("filtered", colnames(EOI_table))]), "T_PRIMARY", "T_CVDMI_P")
)
param_list <- split(param_table, seq(nrow(param_table)))

pairsout <- ggpairs(singscore_testtable[,c(singscore_list, "filtered731_PRIMARY")], aes(color = filtered731_PRIMARY))

# manually add in all of our scores.........?  This was still really not very good........
# manual_add <- list(addon1 = c(rna_score = c("rna_imgdegis_score1", "rna_ctndv70_score1", "rna_subtype1", "rna_subtype2", "rna_subtype3", "rna_subtype4"),
#                    meth_score = c("meth_imgdegis_score1", "meth_ctndv70_score1", "meth_subtype1", "meth_subtype2", "meth_subtype3"),
#                    events = "filtered548_PRIMARY"))

outstat_list <- list()
outstatcounter <- 1
for (param_sel in param_list[1:length(param_list)]) {
    
    
    event_sel <- as.character(param_sel[["events"]])
    outcomelevels <- paste0(event_sel, c("_event", "_noevent"))
    selected_testtable <- singscore_testtable[, as.character(unlist(param_sel))]
    # ## MANUAL TESTING
    # selected_testtable <- singscore_testtable[,c("rna_subtype1", "meth_subtype1", "T_PRIMARY")]
    # event_sel <- "T_PRIMARY"
    # ## MANUAL TESTING
    if (!grepl("^T_", event_sel)) {
        selected_testtable <- selected_testtable[rowSums(is.na(selected_testtable)) == 0 & selected_testtable[,event_sel] != 2,]
        selected_testtable[,event_sel] <- ifelse(selected_testtable[,event_sel] == 1, paste0(event_sel, "_event"), paste0(event_sel, "_noevent"))
    } else {
        selected_testtable <- na.omit(selected_testtable)
    }
    
    featuretable <- selected_testtable[,!colnames(selected_testtable) %in% event_sel]
    outcometable <- selected_testtable[,event_sel,drop=FALSE]
    
    # Set up out featuretable and outcometable
    # strsplit_temp <- strsplit(dge_comp, split = "_")[[1]]
    # outcomelevels <- strsplit_temp[c(length(strsplit_temp) - 2, length(strsplit_temp))]
    # featurelabels <- GOI
    # outcometable <- metatable[,outcomelabel,drop=FALSE]
    # outcometable[,1] <- ifelse(outcometable[,1] == 1, outcomelevels[1], ifelse(outcometable[,1] == 0, outcomelevels[2], NA))
    # featuretable <- t(normcounttable[featurelabels,rownames(metatable)])
    
    # modeling_outfilepath = paste0(outfilepathmaster, "cluster_performance_modeling_down/", dge_comp, "/"),
    
    # Then run our modeling with these genes and labels
    out1 <- cohort_classification_ml_analysis(featuretable=featuretable, outcometable=outcometable, outcomelevels=outcomelevels,
                                              seedparam=11111, subsampleparam = "down",
                                              # models_to_run = c("glm", "xgbTree", "svmRadial", "rf", "glmnet")
                                              models_to_run = c("glm", "svmRadial", "rf", "glmnet")
    )
    outstat <- cbind( param_sel, out1[[1]][out1[[1]][,1] %in% "AUC_w_CI",])
    
    outstat_list[[outstatcounter]] <- outstat
    
    # temp writing this out
    write.table(outstat, paste0(outfilepathmaster, "score_v_event_analysis/tempout/outstat", outstatcounter, ".csv"), sep = ",", col.names = TRUE, row.names = FALSE)
    # For this first run - i really dont need the table and plot for each analysis - thats wasteful, lets just grab AUCs
    # Write out the results
    # write.table(out1$modelstats_outtable, paste0(outfilepathmaster, "cluster_performance_modeling_down/", dge_comp, "/", "model_stats_outtable.csv"),
    #             sep = ",", col.names = TRUE, row.names = FALSE)
    # pdf(paste0(outfilepathmaster, "cluster_performance_modeling_down/", dge_comp, "/", "ROCcurves.pdf"), 10, 10, useDingbats = FALSE)
    # print(out1$rocplot)
    # junk <- dev.off()
    outstatcounter <- outstatcounter + 1
    
}
outstattable <- do.call(rbind, outstat_list)
write.table(outstattable, paste0(outfilepathmaster, "score_v_event_analysis/tempout/outstatfinal.csv"), sep = ",", col.names = TRUE, row.names = FALSE)



# Output this class table
event_sel <- "filtered182_PRIMARY"
outcomelevels <- paste0(event_sel, c("_event", "_noevent"))
selected_testtable <- singscore_testtable[, c("rna_subtype1", "meth_subtype1", "filtered182_PRIMARY")]
selected_testtable <- selected_testtable[rowSums(is.na(selected_testtable)) == 0 & selected_testtable[,event_sel] != 2,]
selected_testtable[,event_sel] <- ifelse(selected_testtable[,event_sel] == 1, paste0(event_sel, "_event"), paste0(event_sel, "_noevent"))
featuretable <- selected_testtable[,!colnames(selected_testtable) %in% event_sel]
outcometable <- selected_testtable[,event_sel,drop=FALSE]

out1 <- cohort_classification_ml_analysis(featuretable=featuretable, outcometable=outcometable, outcomelevels=outcomelevels,
                                          seedparam=11111, subsampleparam = "down",
                                          # models_to_run = c("glm", "xgbTree", "svmRadial", "rf", "glmnet")
                                          models_to_run = c("rf", "nnet")
)






## TRYING A TIME TO EVENT ANALYSIS INSTEAD!
# aggregate(EOI_table[,"T_PRIMARY"], by = list(EOI_table[,"C_PRIMARY"]), summary)
# table(EOI_table[,"C_PRIMARY"])
# 
# singscore_testtable <- merge(combined_singscore_table, EOI_table, by = "row.names", all = TRUE)
# rownames(singscore_testtable) <- singscore_testtable[,"Row.names"]
# param_table <- expand.grid(rna_score = c("rna_imgdegis_score1", "rna_ctndv70_score1", "rna_subtype1", "rna_subtype2", "rna_subtype3", "rna_subtype4"),
#                            meth_score = c("meth_imgdegis_score1", "meth_ctndv70_score1", "meth_subtype1", "meth_subtype2", "meth_subtype3"),
#                            # events = c(colnames(EOI_table[grepl("filtered", colnames(EOI_table))]))
#                            events = c("T_PRIMARY", "T_CVDMI_P")
# )
# param_list <- split(param_table, seq(nrow(param_table)))
# outstat_list <- list()
# outstatcounter <- 1
# for (param_sel in param_list[1:length(param_list)]) {
#     
#     event_sel <- as.character(param_sel[["events"]])
#     outcomelevels <- paste0(event_sel, c("_event", "_noevent"))
#     selected_testtable <- singscore_testtable[, as.character(unlist(param_sel))]
#     selected_testtable <- na.omit(selected_testtable)
#     featuretable <- selected_testtable[,!colnames(selected_testtable) %in% event_sel]
#     outcometable <- selected_testtable[,event_sel,drop=FALSE]
#     
#     # Then run our modeling with these genes and labels
#     out1 <- regression_ml_analysis(featuretable=featuretable, outcometable=outcometable, outcomelevels=outcomelevels,
#                                               seedparam=11111, 
#                                             # subsampleparam = "down",
#                                               # models_to_run = c("glm", "xgbTree", "svmRadial", "rf", "glmnet")
#                                               models_to_run = c("glm", "svmRadial", "rf", "glmnet")
#     )
#     outstat <- cbind( param_sel, out1[[1]][out1[[1]][,1] %in% "AUC_w_CI",])
#     
#     outstat_list[[outstatcounter]] <- outstat
#     
#     # temp writing this out
#     write.table(outstat, paste0(outfilepathmaster, "score_v_event_analysis/tempout/outstat", outstatcounter, ".csv"), sep = ",", col.names = TRUE, row.names = FALSE)
#     # For this first run - i really dont need the table and plot for each analysis - thats wasteful, lets just grab AUCs
#     # Write out the results
#     # write.table(out1$modelstats_outtable, paste0(outfilepathmaster, "cluster_performance_modeling_down/", dge_comp, "/", "model_stats_outtable.csv"),
#     #             sep = ",", col.names = TRUE, row.names = FALSE)
#     # pdf(paste0(outfilepathmaster, "cluster_performance_modeling_down/", dge_comp, "/", "ROCcurves.pdf"), 10, 10, useDingbats = FALSE)
#     # print(out1$rocplot)
#     # junk <- dev.off()
#     outstatcounter <- outstatcounter + 1
#     
# }
# outstattable <- do.call(rbind, outstat_list)
# write.table(outstattable, paste0(outfilepathmaster, "score_v_event_analysis/tempout/outstatfinal.csv"), sep = ",", col.names = TRUE, row.names = FALSE)





## Ok what about a regression INTO class model - TOSH WORKING HERE
aggregate(EOI_table[,"T_PRIMARY"], by = list(EOI_table[,"C_PRIMARY"]), summary)
table(EOI_table[,"C_PRIMARY"])

singscore_testtable <- merge(combined_singscore_table, EOI_table, by = "row.names", all = TRUE)
rownames(singscore_testtable) <- singscore_testtable[,"Row.names"]
param_table <- expand.grid(rna_score = c("rna_imgdegis_score1", "rna_ctndv70_score1", "rna_subtype1", "rna_subtype2", "rna_subtype3", "rna_subtype4"),
                           meth_score = c("meth_imgdegis_score1", "meth_ctndv70_score1", "meth_subtype1", "meth_subtype2", "meth_subtype3"),
                           # events = c(colnames(EOI_table[grepl("filtered", colnames(EOI_table))]))
                           events = c("PRIMARY", "CVDMI_P")
)
param_list <- split(param_table, seq(nrow(param_table)))
outstat_list <- list()
outstatcounter <- 1
# for (param_sel in param_list[1:length(param_list)]) {
for (param_sel in param_list[31:length(param_list)]) {
    
    event_sel <- as.character(param_sel[["events"]])
    outcomelevels <- paste0("C_", event_sel, c("_event", "_noevent"))
    selected_testtable <- singscore_testtable[, c(as.character(unlist(param_sel)[c(1,2)]), paste0(c("T_", "C_"), event_sel))]
    selected_testtable <- selected_testtable[rowSums(is.na(selected_testtable)) == 0 & selected_testtable[,paste0("C_", event_sel)] != 2,]
    selected_testtable[,paste0("C_", event_sel)] <- ifelse(selected_testtable[,paste0("C_", event_sel)] == 1,
                                                           paste0("C_", event_sel, "_event"), paste0(paste0("C_", event_sel), "_noevent"))
    
    featuretable <- selected_testtable[,!grepl(event_sel, colnames(selected_testtable))]
    outcometable <- selected_testtable[,grepl(event_sel, colnames(selected_testtable)),drop=FALSE]
    
    # Then run our modeling with these genes and labels
    out1 <- regression_to_classification_ml_analysis(
        featuretable=featuretable, outcometable=outcometable, outcomelevels=outcomelevels,
        seedparam=11111, subsampleparam = "down",
        # models_to_run = c("glm", "xgbTree", "svmRadial", "rf", "glmnet")
        models_to_run = c("glm", "svmRadial", "rf", "glmnet")
    )
    outstat <- cbind( param_sel, out1[out1[,1] %in% "AUC_w_CI",])
    outstat_list[[outstatcounter]] <- outstat
    
    # temp writing this out
    write.table(outstat, paste0(outfilepathmaster, "score_v_event_analysis/tempout_toe/outstat", outstatcounter, ".csv"), sep = ",", col.names = TRUE, row.names = FALSE)
    # For this first run - i really dont need the table and plot for each analysis - thats wasteful, lets just grab AUCs
    # Write out the results
    # write.table(out1$modelstats_outtable, paste0(outfilepathmaster, "cluster_performance_modeling_down/", dge_comp, "/", "model_stats_outtable.csv"),
    #             sep = ",", col.names = TRUE, row.names = FALSE)
    # pdf(paste0(outfilepathmaster, "cluster_performance_modeling_down/", dge_comp, "/", "ROCcurves.pdf"), 10, 10, useDingbats = FALSE)
    # print(out1$rocplot)
    # junk <- dev.off()
    outstatcounter <- outstatcounter + 1
    
}
outstattable <- do.call(rbind, outstat_list)
write.table(outstattable, paste0(outfilepathmaster, "score_v_event_analysis/tempout_toe/outstatfinal.csv"), sep = ",", col.names = TRUE, row.names = FALSE)





# --------------------------------- cohort signature score vs events - COXPH ---------------------------------
# dir.create(paste0(outfilepathmaster, "score_v_event_analysis/"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(outfilepathmaster, "score_v_event_analysis/coxph/"), showWarnings = FALSE, recursive = TRUE)
# Read in singscore results results
singscore_list <- c("rna_imgdegis_score1", "rna_ctndv70_score1", "rna_subtype1", "rna_subtype2", "rna_subtype3", "rna_subtype4",
                    "meth_imgdegis_score1", "meth_ctndv70_score1", "meth_subtype1", "meth_subtype2", "meth_subtype3")
singscore_table_list <- list()
singscore_path <- paste0(outfilepathmaster, "signature_scoring/")
for (singscore_result in singscore_list) {
    intable_singscore <- read.table(paste0(singscore_path, singscore_result, "/singscore/", singscore_result, "_singscore_outtable.csv"),
                                    sep = ",", header = TRUE)
    singscore_table_list[[singscore_result]] <- intable_singscore
}
combined_singscore_table <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Row.names", all = TRUE, sort = FALSE), 
                                   lapply(singscore_table_list, function(x) x[,c("Row.names", "TotalScore")]))
dimnames(combined_singscore_table) <- list(combined_singscore_table[,"Row.names"], c("Row.names", names(singscore_table_list)))
combined_singscore_table <- combined_singscore_table[,!grepl("Row.names", colnames(combined_singscore_table))]

# Now for the COXPH we need some baseline characteristics, and then our events of interest
# coxph_COI <- c(
#             # "CUSTOM_maxfollow", 
#              "AGE_RAND", "SEX",
#              # "RACE",
#              "CUSTOM_RACE_AGGR", "ETHNIC",
#              # "BMI",
#              "CUSTOM_BMI_obese", "HYPTENSE", "DIABETES",
#              # "DIABTRT",
#              # "SMOKSTAT",
#              "CUSTOM_SMOKSTAT_currentsmoker", 
#              # "FAMHXPHD", "PRIORMI", "PRIPCIYN", "HXCABG", "CUSTOM_PCIorCABG", "XPRIFALYN", "PRICEREB", "PRIPVDYN",
#              # "LVEFCONT", "CUSTOM_LVEFCONT_LVEFund55", "EFLT45", "EGFR", "CUSTOM_EGFR_EGFRund60", "DIALYSIS",
#              # "HDLC", "TRIGLYC", "TOTCHOL", "RVBPDIA",  "HEMOGLOB",
#              # "LDLC", "CUSTOM_LDLC_LDLCund70", "RVBPSYS", "CUSTOM_RVBPSYS_RVBPSYSund140", "HEMOGA1C",
#              # "MESTATIN", "CUSTOM_RVBPSYS_RVBPSYSund140", "CUSTOM_LDLCund70andMESTATIN", "MEAPAC",
#              # "DEGRISCH", "CUSTOM_DEGRISCH_NONEMILD", "IMGDEGIS",
#              "CUSTOM_IMGDEGIS_NONEMILD",
#              "CTNDV70", "CTMULT70", "CTNDV50", "CTMULT50"
#              # "DUKESCORE"
#              # "TRT", "CKD", 
#              # "CUSTOM_CKD_TRT",
#              # "MEAPASP", "MEAPCLOP", "HIGHSTAT", "MELLEZET", "MEAHACE", "MEAHARB", "CUSTOM_ACEIandARB", "INSULDBT"
#              # "C_PRIMARY", "C_CVDMI_P", "C_ACD", "C_CVD", "C_ACDMI_P", "C_MIPRIM", "C_RCA", "C_HOSP_UA", "C_HOSP_HF",
#              # colnames(biomarkertable)[grepl("_clean$", colnames(biomarkertable))],
#              # "rna_4cluster_w3AB", "rna_4cluster", "meth_4cluster", "meth_3cluster"
# )

# "AGE_RAND", "SEX","CUSTOM_RACE_AGGR", "ETHNIC","CUSTOM_BMI_obese", "HYPTENSE", "DIABETES","SMOKSTAT", "CUSTOM_IMGDEGIS_NONEMILD", "CTNDV70", "CTMULT70", "CTNDV50", "CTMULT50"
EOI_vec <- c("PRIMARY", "CVDMI_P")
# coxph_baseline_table <- bioreptable_waddons[,coxph_COI]
# coxph_analysis_table <- merge(coxph_baseline_table, combined_singscore_table, by = "row.names", all = TRUE)

# coxph_analysis_table <- merge(bioreptable_waddons[,"PATNUM",drop=FALSE], combined_singscore_table,
#                               by.x = "PATNUM", by.y = "row.names", all = TRUE)
# rownames(coxph_analysis_table) <- coxph_analysis_table[,"PATNUM"]

coxph_baseline_table <- bioreptable_waddons[,"PATNUM",drop=FALSE]
coxph_baseline_table[,"AGE_RAND"] <- bioreptable_waddons[,"AGE_RAND"]
coxph_baseline_table[,"SEX"] <- factor(bioreptable_waddons[,"SEX"], levels = c("Female", "Male"))
coxph_baseline_table[,"CUSTOM_RACE_AGGR"] <- factor(bioreptable_waddons[,"CUSTOM_RACE_AGGR"], levels = c("White", "Black or African American", "Asian", "Other or multiple ethnic groups"))
coxph_baseline_table[,"ETHNIC"] <- factor(bioreptable_waddons[,"ETHNIC"], levels = c("Not Hispanic or Latino", "Hispanic or Latino"))
coxph_baseline_table[,"CUSTOM_BMI_obese"] <- factor(bioreptable_waddons[,"CUSTOM_BMI_obese"], levels = c("NotObese", "Obese"))
coxph_baseline_table[,"HYPTENSE"] <- factor(bioreptable_waddons[,"HYPTENSE"], levels = c("No", "Yes"))
coxph_baseline_table[,"SMOKSTAT"] <- factor(bioreptable_waddons[,"SMOKSTAT"], levels = c("Never Smoked", "Former Smoker", "Current Smoker"))
# coxph_analysis_table[,"CUSTOM_IMGDEGIS_NONEMILD"] <- factor(bioreptable_waddons[,"CUSTOM_IMGDEGIS_NONEMILD"], levels = c("NoneMild", "Moderate", "Severe"))
# coxph_analysis_table[,"CTNDV70"] <- factor(bioreptable_waddons[,"CTNDV70"], levels = c("0", "1", "2", "3", "Non-evaluable"))
# coxph_analysis_table[,"CTMULT70"] <- factor(bioreptable_waddons[,"CTMULT70"], levels = c("No", "Yes", "Not evaluable"))
# coxph_analysis_table[,"CTNDV50"] <- factor(bioreptable_waddons[,"CTNDV50"], levels = c("1", "2", "3", "Non-evaluable"))
# coxph_analysis_table[,"CTMULT50"] <- factor(bioreptable_waddons[,"CTMULT50"], levels = c("No", "Yes", "Not evaluable"))

coxph_baseline_table <- coxph_baseline_table[rownames(bioreptable_waddons), !grepl("PATNUM", colnames(coxph_baseline_table))]



# For each sigscore - we want a coxph of all of our variables + the sigscore for each event, also for each event just have everything
for (EOI in EOI_vec) {
    event_table_sel <- bioreptable_waddons[,grepl(EOI, colnames(bioreptable_waddons))]
    timevariable <- colnames(event_table_sel)[grepl("^T_", colnames(event_table_sel))]
    eventvariable <- colnames(event_table_sel)[grepl("^C_", colnames(event_table_sel))]
    
    ## Format to (1) time to, (2) event, (3) everything else
    survivaldata_covar <- merge(event_table_sel[,c(timevariable, eventvariable)], coxph_baseline_table,
                                by = "row.names", all = TRUE)
    rownames(survivaldata_covar) <- survivaldata_covar[,1]
    survivaldata_covar <- survivaldata_covar[,!grepl("Row.names", colnames(survivaldata_covar))]
    survivaldata_covar[,eventvariable] <- as.numeric(survivaldata_covar[,eventvariable])
    # Remove competing events
    survivaldata_covar <- survivaldata_covar[survivaldata_covar[,eventvariable] %in% c(1,0),]
    
    for (score_metric in singscore_list) {
        score_sel <- na.omit(combined_singscore_table[,score_metric,drop=FALSE])
        survivaldata_selscore <- na.omit(cbind.data.frame(survivaldata_covar[rownames(score_sel),], score_sel))
        
        # Preprocess our data?
        preprocvalues <- preProcess(survivaldata_selscore[,!colnames(survivaldata_selscore) %in% c(timevariable, eventvariable)],
                                    method = c("center", "scale"))
        survivaldata_selscore[,!colnames(survivaldata_selscore) %in% c(timevariable, eventvariable)] <- predict(
            preprocvalues, survivaldata_selscore[,!colnames(survivaldata_selscore) %in% c(timevariable, eventvariable)])
        
        
        coxph_out <- coxph_analysis(survivaldata = survivaldata_selscore)
        coxph_outtable <- coxph_out$coxph_outtable
        coxph_outplot <- coxph_out$coxph_outplot
        coxph_outplot
        
        write.table(coxph_outtable, paste0(outfilepathmaster, "score_v_event_analysis/coxph/", score_metric, "_coxph_table.csv"), sep = ",",
                    col.names = TRUE, row.names = FALSE)
        pdf(paste0(outfilepathmaster, "score_v_event_analysis/coxph/", score_metric, "_coxph_plot.pdf"), 11.5, 15)
        print(coxph_outplot)
        junk <- dev.off()
    }
    
    # First lets do an analysis with ALL of our scores in there - doesnt really work, so lets try each score individually?
    # survivaldata_allscores <- na.omit(survivaldata_full)
    # coxph_out <- coxph_analysis(survivaldata = survivaldata_allscores)
    # coxph_outtable <- coxph_out$coxph_outtable
    # coxph_outplot <- coxph_out$coxph_outplot
    # coxph_outplot
    # 
    # write.table(coxph_outtable, paste0(outfilepathmaster, "multivariate_analysis/", outcome_sel, "_coxph_table.csv"), sep = ",",
    #             col.names = TRUE, row.names = FALSE)
    # pdf(paste0(outfilepathmaster, "multivariate_analysis/", outcome_sel, "_coxph_plot.pdf"), 11.5, 15)
    # print(coxph_outplot)
    # junk <- dev.off()
    
    
    
}





EOI_table <- bioreptable_waddons[,grepl(paste(EOI_vec, collapse = "|"), colnames(bioreptable_waddons))]
timeparam <- c(182, 365, 548, 731)
for (EOI in EOI_vec) {
    for (timevar in timeparam) {
        EOI_table[,paste0("filtered", timevar, "_", EOI)] <- ifelse(EOI_table[, paste0("T_", EOI)] >= timevar & 
                                                                        EOI_table[, paste0("C_", EOI)] == 1, 2, # filter out too late event
                                                                    ifelse(EOI_table[, paste0("T_", EOI)] < timevar & 
                                                                               EOI_table[, paste0("C_", EOI)] == 0, 2, # filter out too soon non-event
                                                                           EOI_table[, paste0("C_", EOI)]))
    }
}
# apply(EOI_table[,grepl("filtered", colnames(EOI_table))], 2, table)



singscore_testtable <- merge(combined_singscore_table, EOI_table, by = "row.names", all = TRUE)
rownames(singscore_testtable) <- singscore_testtable[,"Row.names"]
param_table <- expand.grid(rna_score = c("rna_imgdegis_score1", "rna_ctndv70_score1", "rna_subtype1", "rna_subtype2", "rna_subtype3", "rna_subtype4"),
                           meth_score = c("meth_imgdegis_score1", "meth_ctndv70_score1", "meth_subtype1", "meth_subtype2", "meth_subtype3"),
                           # events = c(colnames(EOI_table[grepl("filtered", colnames(EOI_table))]))
                           events = c(colnames(EOI_table[grepl("filtered", colnames(EOI_table))]), "T_PRIMARY", "T_CVDMI_P")
)
param_list <- split(param_table, seq(nrow(param_table)))

# manually add in all of our scores.........?  This was still really not very good........
# manual_add <- list(addon1 = c(rna_score = c("rna_imgdegis_score1", "rna_ctndv70_score1", "rna_subtype1", "rna_subtype2", "rna_subtype3", "rna_subtype4"),
#                    meth_score = c("meth_imgdegis_score1", "meth_ctndv70_score1", "meth_subtype1", "meth_subtype2", "meth_subtype3"),
#                    events = "filtered548_PRIMARY"))





# # find the time variable and add it on
# timevariable <- paste0("time_to_", gsub(pattern = "censor_", replacement = "", outcome_sel))
# coxph_intable_temp <- na.omit(cbind(mv_analysis_intable_sel, 
#                                     fullmetatable[rownames(mv_analysis_intable_sel), timevariable, drop=FALSE]))
# coxph_intable_temp[,outcome_sel] <- ifelse(coxph_intable_temp[,outcome_sel] %in% c("Yes"), 1, 0)
# 
# ## Turn age into tertiles
# coxph_intable_temp[,"age"] <- cut2(t(coxph_intable_temp[,"age",drop=FALSE]),g=3)
# # levels(coxph_intable_temp[,"age"])[match(levels(coxph_intable_temp[,"age"])[1],levels(coxph_intable_temp[,"TESscore"]))] <- "Q1_Low"
# # levels(coxph_intable_temp[,"age"])[match(levels(coxph_intable_temp[,"age"])[2],levels(coxph_intable_temp[,"TESscore"]))] <- "Q2"
# # levels(coxph_intable_temp[,"age"])[match(levels(coxph_intable_temp[,"age"])[3],levels(coxph_intable_temp[,"TESscore"]))] <- "Q3_High"
# 
# ## Turn TES into tertiles
# # coxph_intable_temp[,"TESscore"] <- cut2(t(coxph_intable_temp[,"TESscore",drop=FALSE]),g=3)
# # levels(coxph_intable_temp[,"TESscore"])[match(levels(coxph_intable_temp[,"TESscore"])[1],levels(coxph_intable_temp[,"TESscore"]))] <- "Q1_Low"
# # levels(coxph_intable_temp[,"TESscore"])[match(levels(coxph_intable_temp[,"TESscore"])[2],levels(coxph_intable_temp[,"TESscore"]))] <- "Q2"
# # levels(coxph_intable_temp[,"TESscore"])[match(levels(coxph_intable_temp[,"TESscore"])[3],levels(coxph_intable_temp[,"TESscore"]))] <- "Q3"
# # levels(coxph_intable_temp[,"TESscore"])[match(levels(coxph_intable_temp[,"TESscore"])[4],levels(coxph_intable_temp[,"TESscore"]))] <- "Q4_High"
# 
# ## Turn TES into tertiles with SD and mean
# # meanval <- mean(coxph_intable_temp[,"TESscore"])
# # sdval <- sd(coxph_intable_temp[,"TESscore"])
# # coxph_intable_temp[,"TESscore"] <- factor(ifelse(coxph_intable_temp[,"TESscore"] < meanval-0.5*sdval, "Low", 
# #                                           ifelse(coxph_intable_temp[,"TESscore"] < meanval+0.5*sdval, "Mid", "High")),
# #                                           levels = c("Mid", "Low", "High"))
# coxph_intable_temp[,"TESscore"] <- scale(coxph_intable_temp[,"TESscore"])
# 
# ## Format to (1) time to, (2) event, (3) everything else
# survivaldata <- coxph_intable_temp[,c(timevariable, outcome_sel, 
#                                       colnames(mv_analysis_intable_sel)[!colnames(mv_analysis_intable_sel) %in% c(timevariable, outcome_sel)])]
# coxph_out <- coxph_analysis(survivaldata, covardata = NULL)
# coxph_outtable <- coxph_out$coxph_outtable
# coxph_outplot <- coxph_out$coxph_outplot
# coxph_outplot
# 
# write.table(coxph_outtable, paste0(outfilepathmaster, "multivariate_analysis/", outcome_sel, "_coxph_table.csv"), sep = ",",
#             col.names = TRUE, row.names = FALSE)
# pdf(paste0(outfilepathmaster, "multivariate_analysis/", outcome_sel, "_coxph_plot.pdf"), 11.5, 15)
# print(coxph_outplot)
# junk <- dev.off()




