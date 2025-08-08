################################
# Name: MacIntosh Cornwell
# Email: mcornwell1957@gmail.com
################################
## Read in and format methylation data

## Load in Libraries
# options(scipen=999)
packagelist = c()
junk <- lapply(packagelist, function(xxx) suppressMessages(
    require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))
source("/Users/tosh/Desktop/Ruggles_Lab/code/mgc_plotting_functions.R")
source("/Users/tosh/Desktop/Ruggles_Lab/code/rnaseq_scripts/deseq_functions.R")
source("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/code/ischemia2021_exploratory_analysis_functions.R")


## Outfilepath
outfilepathmaster <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/meth_processing/run3/asr_corrected_DMP/"
# outfilepathmaster <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/meth_processing/run3/uncorrected_DMP/"
dir.create(outfilepathmaster, showWarnings = FALSE, recursive = TRUE)


# --------------------------------- Read in Methylation Count Table and metadata table ---------------------------------
# Read in meth count table and metatable
methmetatable_file <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/final_meth_data_v1_20230207/kpp_metadata.rds"
methmetatable <- readRDS(methmetatable_file)
methcounttable_file <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/final_meth_data_v1_20230207/limma_noCovariate.rds"
methcounttable <- readRDS(methcounttable_file)
colnames(methcounttable) <- methmetatable[match(methmetatable[,"Sample_Name"], colnames(methcounttable)), "Sample_Patnum"]
methcounttable <- methcounttable[,-578] # Need to remove the duplicate of this sample: "048003-005", columns 204 and 578
methmetatable <- methmetatable[-578,] # Need to remove the duplicate of this sample: "048003-005", columns 204 and 578

# Colorguide
colorguidefile <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/ischemia2021_colorguide.csv"
colorguide <- read.table(colorguidefile, sep = ",", header = TRUE, comment.char = "", colClasses = c("character", "character", "character"))

# Meth cluster membership table
meth_clustermembership_table_file <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/MULTIOMICS_SUBTYPE_LABELS/meth_cluster_membership_table_4_nodup_20220606.csv"
meth_clustermembership_table <- read.table(meth_clustermembership_table_file, sep = ",", header = TRUE, row.names = 1)


# --------------------------------- Custom functions for this script ---------------------------------
clean_gsea_hypergeo_plot_out <- function(gsea_table, direction, number_pathways, titleparam) {
    plottable <- gsea_table[1:number_pathways,c("ID", "pvalue", "p.adjust")]
    plottable[,"pstat"] <- -log10(plottable[,"pvalue"]) * ifelse(direction == "up", 1, -1)
    plottable[,"colorvar"] <- factor(ifelse(plottable[,"pstat"] > 0 & plottable[,"p.adjust"] < 0.05, "darkred", 
                                     ifelse(plottable[,"pstat"] > 0 & plottable[,"pvalue"] < 0.05, "red", 
                                     ifelse(plottable[,"pstat"] < 0 & plottable[,"p.adjust"] < 0.05, "darkblue", 
                                     ifelse(plottable[,"pstat"] < 0 & plottable[,"pvalue"] < 0.05, "blue", "grey")))),
                                     levels = c("darkred", "red", "grey", "blue", "darkblue"))
    plottable[,"ID"] <- gsub("GOMF_|GOBP_|GOCC_|HALLMARK_", "", plottable[,"ID"])
    plottable[,"ID"] <- sapply(tolower(gsub("_", " ", plottable[,"ID"])), simpleCap)
    plottable <- plottable[order(plottable[,"pstat"], decreasing = direction != "up"),]
    
    GSEA_pout <- plot_barchart(bartab = plottable[,c("ID", "pstat")], labsparam = list(title = titleparam, x = "Pathway", y = "signed -log10(p.value)"),
                               colorvar = plottable[,"colorvar",drop=FALSE], strwrap_characternum = 32)
    # GSEA_pout <- plot_barchart(bartab = plottable[,c("ID", "pstat")], labsparam = list(title = titleparam, x = "Pathway", y = "signed -log10(p.value)"),
    #                            strwrap_characternum = 32)
    GSEA_pout <- GSEA_pout + coord_flip() + theme_pubr(base_size = 20, legend = "right")
    GSEA_pout <- GSEA_pout + guides(fill = guide_legend(title = "Sig"))
    suppressMessages(if (direction == "up") {
        GSEA_pout <- GSEA_pout + scale_fill_identity(guide = "legend", breaks = c("darkred", "red", "grey"), 
                                                     labels = c("adj. pval < 0.05", "pval < 0.05", "pval > 0.05"))
    } else {
        GSEA_pout <- GSEA_pout + scale_fill_identity(guide = "legend", breaks = c("darkblue", "blue", "grey"), 
                                                     labels = c("adj. pval < 0.05", "pval < 0.05", "pval > 0.05"))
    })
    return(GSEA_pout)
}

# --------------------------------- Make some custom annotations for plotting ---------------------------------
meth_annottable <- merge(methmetatable, meth_clustermembership_table[,"meth_3cluster",drop=FALSE], by.x = "Sample_Patnum", by.y = "row.names")
rownames(meth_annottable) <- meth_annottable[,"Sample_Patnum"]

## Create custom annotation for AFIB
meth_annottable[,"AFIB_CUSTOM"] <- ifelse(meth_annottable[,"ECGAFIB"] == "Yes" | meth_annottable[,"XAFIBYN"] == "Yes", "AFIB", "noAFIB")
meth_annottable[,"AFIBECG_CUSTOM"] <- ifelse(meth_annottable[,"ECGAFIB"] == "Yes", "AFIBECG", "noAFIB")
meth_annottable[meth_annottable[,"XAFIBYN"] == "Yes" & meth_annottable[,"ECGAFIB"] == "No","AFIBECG_CUSTOM"] <- NA

## Create comp columns
meth_annottable[,"DMP_afib__afib_v_noafib"] <- ifelse(meth_annottable[,"ECGAFIB"] == "Yes" | meth_annottable[,"XAFIBYN"] == "Yes", "AFIB", "noAFIB")
meth_annottable[,"DMP_afibecg__afibecg_v_noafib"] <- ifelse(meth_annottable[,"ECGAFIB"] == "Yes", "AFIBECG", "noAFIB")
meth_annottable[meth_annottable[,"XAFIBYN"] == "Yes" & meth_annottable[,"ECGAFIB"] == "No","DMP_afibecg__afibecg_v_noafib"] <- NA

meth_annottable[,"DMP_anatomy70__3v_v_1v"] <- ifelse(meth_annottable[,"CTNDV70"] == "3", "3", 
                                              ifelse(meth_annottable[,"CTNDV70"] == "1", "1", NA))
meth_annottable[,"DMP_ischemia__Sev_v_MildNone"] <- ifelse(meth_annottable[,"IMGDEGIS"] == "Severe", "Severe", 
                                                    ifelse(meth_annottable[,"IMGDEGIS"] %in% c("None", "Mild"), "NoneMild", NA))

meth_annottable[,"DMP_methsubtype__methtype2_v_methtype1"] <- ifelse(meth_annottable[,"meth_3cluster"] == "meth_3cluster_2", "MS2", 
                                                              ifelse(meth_annottable[,"meth_3cluster"] == "meth_3cluster_1", "MS1", NA))
meth_annottable[,"DMP_methsubtype__methtype2_v_methtype3"] <- ifelse(meth_annottable[,"meth_3cluster"] == "meth_3cluster_2", "MS2", 
                                                              ifelse(meth_annottable[,"meth_3cluster"] == "meth_3cluster_3", "MS3", NA))
meth_annottable[,"DMP_methsubtype__methtype1_v_methtype3"] <- ifelse(meth_annottable[,"meth_3cluster"] == "meth_3cluster_1", "MS1", 
                                                              ifelse(meth_annottable[,"meth_3cluster"] == "meth_3cluster_3", "MS3", NA))

meth_annottable[,"DMP_methsubtype__methtype1_v_NOTmethtype1"] <- ifelse(meth_annottable[,"meth_3cluster"] == "meth_3cluster_1", "MS1", "notMS1")
meth_annottable[,"DMP_methsubtype__methtype2_v_NOTmethtype2"] <- ifelse(meth_annottable[,"meth_3cluster"] == "meth_3cluster_1", "MS2", "notMS2")
meth_annottable[,"DMP_methsubtype__methtype3_v_NOTmethtype3"] <- ifelse(meth_annottable[,"meth_3cluster"] == "meth_3cluster_1", "MS3", "notMS3")

write.table(meth_annottable, paste0(outfilepathmaster, "custom_meth_metatable.csv"), sep = ",", col.names = NA, row.names = TRUE)

## Doublechecking the metadata from the DMP analyses:
# infilepath <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/methylation_data_20221006/with_asr_correction/"
# metadata_file_vec <- list.files(infilepath)[grepl("metadata", list.files(infilepath))]
# names(metadata_file_vec) <- c("AFIB_CUSTOM", "AFIBECG_CUSTOM", "CTNDV70", "IMGDEGIS", "kpp_cluster", "kpp_cluster")
# metadata_check_list <- list()
# for (metadata_num in seq_along(metadata_file_vec)) {
#     print(metadata_file_vec[metadata_num])
#     metadata_in <- readRDS(paste0(infilepath, metadata_file_vec[metadata_num]))
#     print(dim(metadata_in))
#     metadata_check_list[[metadata_num]] <- metadata_in[,c("Sample_Patnum", names(metadata_file_vec[metadata_num]))]
# }
# metadata_check_table <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Sample_Patnum", all = TRUE, sort = FALSE), metadata_check_list)
# metadata_check_table <- metadata_check_table[metadata_check_table[,"Sample_Patnum"] %in% rownames(meth_annottable),]
# 
# metadata_check_table[metadata_check_table[,"Sample_Patnum"] %in% "048003-005",]
# "048003-005"

# --------------------------------- Read in DMP File ---------------------------------
## Read in the DMP file
## If it is multiple comparisons (like with the subtypes) - then split and analyze

DMP_path <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/final_meth_data_v1_20230207/DMP_with_asr_correction/"
dmp_comp_list <- list(
    meth_AFIB_CUSTOM = c(filelabel = "limma_BCcorrected_DMP_on_AFIB_CUSTOM_control.rds", complabel = "limma_BCcorrected_DMP_on_AFIB_CUSTOM_control",
        c(pval_test = "P.Value", pval_cutoff = 0.05, log2fc_cutoff = 0.03), COI = "AFIB_CUSTOM", outlabel = "DMP_afib__afib_v_noafib"),
    meth_AFIBECG_CUSTOM = c(filelabel = "limma_BCcorrected_DMP_on_AFIBECG_CUSTOM_control.rds", complabel = "limma_BCcorrected_DMP_on_AFIBECG_CUSTOM_control",
        c(pval_test = "P.Value", pval_cutoff = 0.05, log2fc_cutoff = 0.03), COI = "AFIBECG_CUSTOM", outlabel = "DMP_afibecg__afibecg_v_noafib"),
    meth_IMGDEGIS = c(filelabel = "limma_BCcorrected_DMP_on_IMGDEGIS_control.rds", complabel = "limma_BCcorrected_DMP_on_IMGDEGIS_control",
        c(pval_test = "adj.P.Val", pval_cutoff = 0.05, log2fc_cutoff = 0.03), COI = "IMGDEGIS", outlabel = "DMP_ischemia__Sev_v_MildNone"),
    meth_CTNDV70 = c(filelabel = "limma_BCcorrected_DMP_on_CTNDV70_control.rds", complabel = "limma_BCcorrected_DMP_on_CTNDV70_CUSTOM_control",
        c(pval_test = "P.Value", pval_cutoff = 0.05, log2fc_cutoff = 0.03), COI = "CTNDV70", outlabel = "DMP_anatomy70__3v_v_1v"),
    meth_methtype2_v_methtype1 = c(filelabel = "limma_BCcorrected_DMP_on_kpp_oneVSone_control.rds", complabel = "2VS1",
        c(pval_test = "adj.P.Val", pval_cutoff = 0.05, log2fc_cutoff = 0.03), COI = "meth_3cluster", outlabel = "DMP_methsubtype__methtype2_v_methtype1"),
    meth_methtype2_v_methtype3 = c(filelabel = "limma_BCcorrected_DMP_on_kpp_oneVSone_control.rds", complabel = "2VS3",
        c(pval_test = "adj.P.Val", pval_cutoff = 0.05, log2fc_cutoff = 0.03), COI = "meth_3cluster", outlabel = "DMP_methsubtype__methtype2_v_methtype3"),
    meth_methtype1_v_methtype3 = c(filelabel = "limma_BCcorrected_DMP_on_kpp_oneVSone_control.rds", complabel = "1VS3",
        c(pval_test = "adj.P.Val", pval_cutoff = 0.05, log2fc_cutoff = 0.03), COI = "meth_3cluster", outlabel = "DMP_methsubtype__methtype1_v_methtype3"),
    meth_methtype1_v_NOTmethtype1 = c(filelabel = "limma_BCcorrected_DMP_on_kpp_oneVSother_control.rds", complabel = "1VSother",
        c(pval_test = "adj.P.Val", pval_cutoff = 0.05, log2fc_cutoff = 0.03), COI = "meth_3cluster", outlabel = "DMP_methsubtype__methtype1_v_NOTmethtype1"),
    meth_methtype2_v_NOTmethtype2 = c(filelabel = "limma_BCcorrected_DMP_on_kpp_oneVSother_control.rds", complabel = "2VSother",
        c(pval_test = "adj.P.Val", pval_cutoff = 0.05, log2fc_cutoff = 0.03), COI = "meth_3cluster", outlabel = "DMP_methsubtype__methtype2_v_NOTmethtype2"),
    meth_methtype3_v_NOTmethtype3 = c(filelabel = "limma_BCcorrected_DMP_on_kpp_oneVSother_control.rds", complabel = "3VSother",
        c(pval_test = "adj.P.Val", pval_cutoff = 0.05, log2fc_cutoff = 0.03), COI = "meth_3cluster", outlabel = "DMP_methsubtype__methtype3_v_NOTmethtype3")

)


# DMP_path <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/final_meth_data_v1_20230207/DMP_nocovar/"
# dmp_comp_list <- list(
#     meth_AFIB_CUSTOM = c(filelabel = "limma_DMP_on_AFIB_CUSTOM_control.rds", complabel = "noAFIB_to_AFIB", 
#         c(pval_test = "P.Value", pval_cutoff = 0.05, log2fc_cutoff = 0.03), COI = "AFIB_CUSTOM", outlabel = "DMP_afib__afib_v_noafib"),
#     meth_AFIBECG_CUSTOM = c(filelabel = "limma_DMP_on_AFIBECG_CUSTOM_control.rds", complabel = "noAFIB_to_AFIB", 
#         c(pval_test = "P.Value", pval_cutoff = 0.05, log2fc_cutoff = 0.03), COI = "AFIBECG_CUSTOM", outlabel = "DMP_afibecg__afibecg_v_noafib"),
#     meth_IMGDEGIS = c(filelabel = "limma_DMP_on_IMGDEGIS_control.rds", complabel = "Severe_to_No_mild", 
#         c(pval_test = "adj.P.Val", pval_cutoff = 0.05, log2fc_cutoff = 0.03), COI = "IMGDEGIS", outlabel = "DMP_ischemia__Sev_v_MildNone"),
#     meth_CTNDV70 = c(filelabel = "limma_DMP_on_CTNDV70_control.rds", complabel = "1_to_3", ## This one is multicomp - so we have to dig into this
#         c(pval_test = "P.Value", pval_cutoff = 0.05, log2fc_cutoff = 0.03), COI = "CTNDV70", outlabel = "DMP_anatomy70__3v_v_1v")
# )


# Read in the DMPfile
# If there are multipl DMP tables in the file - split them out
# Then for each DMP table - run a summary analysis
## Volcano plot, GSEA, heatmap
for (dmp_comp in dmp_comp_list) {
    subanalysis_outfilepathmaster <- paste0(outfilepathmaster, dmp_comp[["outlabel"]], "/")
    dir.create(subanalysis_outfilepathmaster, showWarnings = FALSE, recursive = TRUE)
    DMPtable <- readRDS(paste0(DMP_path, dmp_comp[["filelabel"]]))
    ## If we have multi-table, then use the complabel to grab the proper table
    if (is.null(nrow(DMPtable))) {DMPtable <- DMPtable[[dmp_comp[["complabel"]]]]}
    
    
    # (1) DMP, volcano, heatmap
    dir.create(paste0(subanalysis_outfilepathmaster, "DMP_analysis_and_figs/"), showWarnings = FALSE, recursive = TRUE)
    
    DMPtable_formatted <- DMPtable[,c("gene", "logFC", colnames(DMPtable)[grepl("_avg", colnames(DMPtable))], "P.Value", "adj.P.Val")]
    colnames(DMPtable_formatted) <- c("gene", "log2FoldChange", colnames(DMPtable)[grepl("_avg", colnames(DMPtable))], "pvalue", "padj")
    DMPtable_formatted[,"cg_probe"] <- rownames(DMPtable_formatted)
    pval_test_param <- ifelse(dmp_comp[["pval_test"]] == "P.Value", "pvalue", "padj")

    
    # Grab labeledgenes
    genefinder_table <- DMPtable_formatted[,c("gene", "log2FoldChange", pval_test_param)]
    genefinder_table[,"sig_logfc"] <- ifelse(abs(genefinder_table[,"log2FoldChange"]) > as.numeric(dmp_comp[["log2fc_cutoff"]]), 1, 0)
    genefinder_table[,"sig_pstat"] <- ifelse(abs(genefinder_table[,pval_test_param]) < as.numeric(dmp_comp[["pval_cutoff"]]), 1, 0)
    genefinder_table[,"has_genename"] <- ifelse(genefinder_table[,"gene"] != "", 1, 0)
    genefinder_table[,"rankstat"] <- genefinder_table[,"log2FoldChange"] * -log10(genefinder_table[,pval_test_param])
    sig_genefinder_table <- genefinder_table[rowSums(genefinder_table[,c(4,5,6)]) == 3,]
    volc_labeled_genenumber = 10
    grab_sigprobes <- rownames(sig_genefinder_table[order(sig_genefinder_table[,"rankstat"], decreasing = TRUE),][c(1:volc_labeled_genenumber,(nrow(sig_genefinder_table)-volc_labeled_genenumber+1):nrow(sig_genefinder_table)),])
    
    DMPtable_formatted[,"labeledgenes"] <- ifelse(rownames(DMPtable_formatted) %in% grab_sigprobes, as.character(DMPtable_formatted[,"gene"]),
                                                  rownames(DMPtable_formatted))
    rownames(DMPtable_formatted) <- make.unique(DMPtable_formatted[,"labeledgenes"])
    labeledgenes <- make.unique(DMPtable_formatted[,"labeledgenes"])[!grepl("^cg", make.unique(DMPtable_formatted[,"labeledgenes"]))]
    
    ## Volcano plot
    vol_pout <- clean_volcano_plot_function(deseqtable = DMPtable_formatted, nameparam=dmp_comp[["complabel"]], labeledgenes = labeledgenes,
        pval_test = pval_test_param, pval_cutoff = as.numeric(dmp_comp[["pval_cutoff"]]), log2fc_cutoff = as.numeric(dmp_comp[["log2fc_cutoff"]]))
    jpeg(paste0(subanalysis_outfilepathmaster, "DMP_analysis_and_figs/", dmp_comp[["outlabel"]], "_volcano_plot.jpeg"), width =768, height = 768)
    print(vol_pout)
    junk <- dev.off()
    # pdf(paste0(subanalysis_outfilepathmaster, "DMP_analysis_and_figs/", dmp_comp[["outlabel"]], "_volcano_plot.pdf"), width =9, height = 9)
    # print(vol_pout)
    # junk <- dev.off()
    
    # Final format and write out the DMP table
    DMP_outtable <- DMPtable_formatted[,1:6]
    rownames(DMP_outtable) <- DMPtable_formatted[,"cg_probe"]
    ## Add back the other useful info
    DMP_outtable <- cbind(DMP_outtable, DMPtable[rownames(DMP_outtable), c("feature", "cgi", "feat.cgi", "CHR", "Strand", "Type")])
    colnames(DMP_outtable)[2] <- "deltaBeta"
    write.table(DMP_outtable, paste0(subanalysis_outfilepathmaster, "DMP_analysis_and_figs/", dmp_comp[["outlabel"]], "_DMP_table.csv"),
                sep = ",", col.names = NA, row.names = TRUE)
    
    pval_summary <- destatplot(DMPtable_formatted[,c("log2FoldChange", "pvalue")], log2fclevels = c(0, 0.01, 0.03, 0.05, 0.1), split_posneg = TRUE)[[2]]
    write.table(pval_summary, paste0(subanalysis_outfilepathmaster, "DMP_analysis_and_figs/", dmp_comp[["outlabel"]], "_pval_summary_table.csv"),
                sep = ",", col.names = NA, row.names = TRUE)
    padj_summary <- destatplot(DMPtable_formatted[,c("log2FoldChange", "padj")], log2fclevels = c(0, 0.01, 0.03, 0.05, 0.1), split_posneg = TRUE)[[2]]
    write.table(padj_summary, paste0(subanalysis_outfilepathmaster, "DMP_analysis_and_figs/", dmp_comp[["outlabel"]], "_padj_summary_table.csv"),
                sep = ",", col.names = NA, row.names = TRUE)
    
    ## Heatmap
    SOI <- rownames(na.omit(meth_annottable[,dmp_comp[["outlabel"]],drop=FALSE]))
    colmetatable <- meth_annottable[SOI, dmp_comp[["COI"]],drop = FALSE]
    colannotationlist <- custom_annotation_list_from_colorguide(COI = c("AFIB", "IMGDEGIS", "CTNDV70", "meth_3cluster"), colorguide)[dmp_comp[["COI"]]]
    orderref <- create_color_reference(COI = c("AFIB", "IMGDEGIS_01_2_3", "CTNDV70", "meth_3cluster"), colorguide)[dmp_comp[["COI"]]][[1]][["ORDER"]]
    
    # There is none for AFIB - so doing a custom insert here
    if (dmp_comp[["COI"]] == "AFIB_CUSTOM") {
        colannotationlist <- annotationlist_builder(colmetatable)
        orderref <- c("AFIB", "noAFIB")
    }
    if (dmp_comp[["COI"]] == "AFIBECG_CUSTOM") {
        colannotationlist <- annotationlist_builder(colmetatable)
        orderref <- c("AFIBECG", "noAFIB")
    }
    
    colmetatable <- colmetatable[order(match(colmetatable[,1], orderref)),,drop=FALSE]
    SOI <- rownames(colmetatable)
    sel_counttable <- methcounttable[grab_sigprobes, SOI]
    
    # rename counttable for heatmap
    sigprobe_to_gene <- data.frame(cbind(unique_genename = rownames(DMPtable_formatted[DMPtable_formatted[,"cg_probe"] %in% grab_sigprobes,]), 
                              cg_probe = DMPtable_formatted[DMPtable_formatted[,"cg_probe"] %in% grab_sigprobes, c("cg_probe")]),
                              row.names = DMPtable_formatted[DMPtable_formatted[,"cg_probe"] %in% grab_sigprobes, c("cg_probe")])[grab_sigprobes,]
    rownames(sel_counttable) <- sigprobe_to_gene[grab_sigprobes,"unique_genename"]
    
    
    
    # Sort samples in case we dont want to cluster
    heatmapcolorparam <- colorRamp2(breaks = c(-2, 0, 2), c("blue", "white", "red"))
    hmout_clust <- create_heatmap(counttab = sel_counttable[,SOI], scale_data = TRUE, heatmapcolorparam = heatmapcolorparam,
                                    colmetatable = colmetatable[SOI,,drop=FALSE], colannotationlist = colannotationlist,
                                    colclusterparam = TRUE, rowclusterparam = TRUE)
    hmout_noclust <- create_heatmap(counttab = sel_counttable[,SOI], scale_data = TRUE, heatmapcolorparam = heatmapcolorparam,
                            colmetatable = colmetatable[SOI,,drop=FALSE], colannotationlist = colannotationlist,
                            colclusterparam = TRUE, rowclusterparam = TRUE, columnsplitparam = colmetatable[SOI,,drop=FALSE])
    pdf(paste0(subanalysis_outfilepathmaster, "DMP_analysis_and_figs/", dmp_comp[["outlabel"]], "_heatmaps.pdf"), 10, 10)
    draw(hmout_clust[[1]])
    draw(hmout_noclust[[1]])
    junk <- dev.off()
    
    # Write out some other helper tables:
    write.table(colmetatable[SOI,,drop=FALSE], paste0(subanalysis_outfilepathmaster, "DMP_analysis_and_figs/", dmp_comp[["outlabel"]], "_colmetatable.csv"),
                sep = ",", col.names = NA, row.names = TRUE)
    write.table(data.frame(table(colmetatable[SOI,])), paste0(subanalysis_outfilepathmaster, "DMP_analysis_and_figs/", dmp_comp[["outlabel"]], "_cohort_counttable.csv"),
                sep = ",", col.names = NA, row.names = TRUE)
    
    
    # (2) Do a simple (unsophisticated) GSEA hypergeo of any genes associated with our significant probes
    dir.create(paste0(subanalysis_outfilepathmaster, "GSEA_on_sigprobes/"), showWarnings = FALSE, recursive = TRUE)
    
    GOIup <- sig_genefinder_table[sig_genefinder_table[,"log2FoldChange"] > 0,"gene"]
    GOIdown <- sig_genefinder_table[sig_genefinder_table[,"log2FoldChange"] < 0,"gene"]
    
    statcutoffparamlist = c("stattype" = "pvalue", "pstatcutoff" = 0.01, "log2fccutoff" = 0) ## Dummy
    speciesparam = "Homo sapiens"
    # Hallmark GSEA
    hypergeo_genetest_out_HALL_UP = hypergeo_genetest(data.frame(GOIup), statcutoffparam = statcutoffparamlist, 
                                                   genesetparam = c("H"), speciesparam = speciesparam)
    write.table(hypergeo_genetest_out_HALL_UP$enricherUPout, 
                file = paste0(subanalysis_outfilepathmaster, "GSEA_on_sigprobes/", dmp_comp[["outlabel"]], "_hypergeo_gene_UP_HALL.csv"), 
                sep = ",", row.names = TRUE, col.names = NA)
    hypergeo_genetest_out_HALL_DOWN = hypergeo_genetest(data.frame(GOIdown), statcutoffparam = statcutoffparamlist, 
                                                   genesetparam = c("H"), speciesparam = speciesparam)
    write.table(hypergeo_genetest_out_HALL_DOWN$enricherUPout, 
                file = paste0(subanalysis_outfilepathmaster, "GSEA_on_sigprobes/", dmp_comp[["outlabel"]], "_hypergeo_gene_DOWN_HALL.csv"), 
                sep = ",", row.names = TRUE, col.names = NA)
    # GO GSEA
    hypergeo_genetest_out_GO_UP = hypergeo_genetest(data.frame(GOIup), statcutoffparam = statcutoffparamlist, 
                                                      genesetparam = c("C5"), speciesparam = speciesparam)
    write.table(hypergeo_genetest_out_GO_UP$enricherUPout, 
                file = paste0(subanalysis_outfilepathmaster, "GSEA_on_sigprobes/", dmp_comp[["outlabel"]], "_hypergeo_gene_UP_GO.csv"), 
                sep = ",", row.names = TRUE, col.names = NA)
    hypergeo_genetest_out_GO_DOWN = hypergeo_genetest(data.frame(GOIdown), statcutoffparam = statcutoffparamlist, 
                                                        genesetparam = c("C5"), speciesparam = speciesparam)
    write.table(hypergeo_genetest_out_GO_DOWN$enricherUPout, 
                file = paste0(subanalysis_outfilepathmaster, "GSEA_on_sigprobes/", dmp_comp[["outlabel"]], "_hypergeo_gene_DOWN_GO.csv"), 
                sep = ",", row.names = TRUE, col.names = NA)
    
    
    # Simple plot function to grab top pos 10 and bot 10, then plot them all with a signed -log10(pval), and color only if p.adjust is sig or p.val is sig
    ## NOTE - preselect OUT the HP genesets - I dont want to write that functionality in here
    H_UP_gsea_pout <- clean_gsea_hypergeo_plot_out(gsea_table = hypergeo_genetest_out_HALL_UP$enricherUPout, direction = "up",
                                 number_pathways = 10, titleparam = paste0("HALLMARK GSEA UP ", dmp_comp[["outlabel"]]))
    H_DOWN_gsea_pout <- clean_gsea_hypergeo_plot_out(gsea_table = hypergeo_genetest_out_HALL_DOWN$enricherUPout, direction = "down",
                                 number_pathways = 10, titleparam = paste0("HALLMARK GSEA DOWN ", dmp_comp[["outlabel"]]))
    GO_UP_gsea_pout <- clean_gsea_hypergeo_plot_out(gsea_table = hypergeo_genetest_out_GO_UP$enricherUPout, direction = "up",
                                 number_pathways = 10, titleparam = paste0("GO GSEA UP ", dmp_comp[["outlabel"]]))
    GO_DOWN_gsea_pout <- clean_gsea_hypergeo_plot_out(gsea_table = hypergeo_genetest_out_GO_DOWN$enricherUPout, direction = "down",
                                 number_pathways = 10, titleparam = paste0("GO GSEA DOWN ", dmp_comp[["outlabel"]]))
    
    pdf(paste0(subanalysis_outfilepathmaster, "GSEA_on_sigprobes/", dmp_comp[["outlabel"]], "_GSEA_plots.pdf"), 13, 10)
    suppressWarnings(print(H_UP_gsea_pout))
    suppressWarnings(print(H_DOWN_gsea_pout))
    suppressWarnings(print(GO_UP_gsea_pout))
    suppressWarnings(print(GO_DOWN_gsea_pout))
    junk <- dev.off()
    
    
    
    # (2) PSEA Analysis - Too computationally intensive to do locally - has go be done on BP
        

}



# --------------------------------- Create Probe sets from Gene sets ---------------------------------

## So can we run GSEA, but with probes aligned to genes
### TOSH - THIS IS ALL DONE AND DONT NEED TO DO AGAIN
# Create the probe set reference table:
# PSEA_reference_table_outfilepath <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/meth_processing/PSEA_reference_tables/"
# dir.create(PSEA_reference_table_outfilepath, showWarnings = FALSE, recursive = TRUE)
# probe_set_out <- create_probe_set_reference(DMPtable, speciesparam = "Homo sapiens", genesetparam = "C5")
# saveRDS(probe_set_out[["probe_set_table"]], paste0(PSEA_reference_table_outfilepath, "probe_set_C5_reftable.rds"))

# probe_set_out <- create_probe_set_reference(DMPtable = DMPtable, speciesparam = "Homo sapiens", genesetparam = "H",
#                                        return_probe_info_table = c("Strand", "Type", "gene", "feature", "cgi", "CHR"))
# saveRDS(probe_set_out[["probe_set_table"]], paste0(PSEA_reference_table_outfilepath, "probe_set_H_reftable.rds"))
# saveRDS(probe_set_out[["probe_info_table"]], paste0(PSEA_reference_table_outfilepath, "probe_info_table.rds"))


# probeset_enrichment_analysis(DMPtable = DMPtable, rankmetric = "logFC", pvalcutoffparam = 1, pvaltype = "P.Value",
#                              existing_probe_ref_file = "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/meth_processing/PSEA_reference_tables/probe_set_H_reftable.rds", probe_type_selection_list = NULL,
#                              genesetparam= "dummy", speciesparam = "dummy",
#                              seedparam = 12345, customprobeset = NULL)

# PSEA_out <- probeset_enrichment_analysis(DMPtable = DMPtable, rankmetric = "logFC", pvalcutoffparam = 1, pvaltype = "P.Value",
#                              existing_probe_ref_file = "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/meth_processing/PSEA_reference_tables/probe_set_H_reftable.rds", probe_type_selection_list = probe_type_selection_list,
#                              genesetparam= "dummy", speciesparam = "dummy",
#                              seedparam = 12345, customprobeset = NULL)





# --------------------------------- Read and process the test train ischemia data ---------------------------------


sig_performace_check_fullsize_output_control_read_and_process <- function(rdsfile, outfilepath) {
    dir.create(outfilepath, showWarnings = FALSE, recursive = TRUE)
    rds_inlist <- readRDS(rdsfile)
    train_splits <- names(rds_inlist)
    for (train_split in train_splits) {
        dir.create(paste0(outfilepath, "trainsplit_", train_split, "/"), showWarnings = FALSE, recursive = TRUE)
        train_sel_inlist <- rds_inlist[[train_split]]
        # Has the df stat table for that run - "sig_performance_validation_ratio0.4_lasso_table.csv"
        sig_performance_validation_ratio_lasso_table <- train_sel_inlist[["df"]]
        write.table(sig_performance_validation_ratio_lasso_table, 
            paste0(outfilepath, "trainsplit_", train_split, "/", "sig_performance_validation_ratio", train_split, "_lasso_table.csv"), 
                    sep = ",", row.names = TRUE, col.names = NA)
        # And the output object with details per analysis
        DMPsets <- names(train_sel_inlist[["output"]])
        for (DMPset in DMPsets) {
            dir.create(paste0(outfilepath, "trainsplit_", train_split, "/", DMPset, "/"), showWarnings = FALSE, recursive = TRUE)
            DMP_sel_inlist <- train_sel_inlist[["output"]][[DMPset]]
            # I already have all of the high level stats, so for now the only thing I need to regen from this is the ROC plot
            ## This is done with the following:
            roc_plottab <- data.frame(cbind(FPR = slot(DMP_sel_inlist[["performance"]], "x.values")[[1]], 
                                            TPR = slot(DMP_sel_inlist[["performance"]], "y.values")[[1]]))
            pout <- ggplot(roc_plottab, aes(x = FPR, y = TPR))
            pout <- pout + geom_line() + geom_abline(intercept = 0, slope = 1, linetype = 2)
            pout <- pout + theme_pubr(border = TRUE) + ggtitle(paste0(train_split, " ", DMPset))
            pdf(paste0(outfilepath, "trainsplit_", train_split, "/", DMPset, "/", train_split, " ", DMPset, "_ROCplot.pdf"),
                10, 10, useDingbats = FALSE)
            print(pout)
            junk <- dev.off()
        }
    }
}

dir.create(paste0(outfilepathmaster, "signature_ml_output/", "testresults/"), showWarnings = FALSE, recursive = TRUE)
rdsfile <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/methylation_data_20221006/ml_output/sig_performace_check_fullsize_output_control.rds"
sig_performace_check_fullsize_output_control_read_and_process(rdsfile,
                                                              outfilepath = paste0(outfilepathmaster, "signature_ml_output/", "testresults/"))

dir.create(paste0(outfilepathmaster, "signature_ml_output/", "trainresults/"), showWarnings = FALSE, recursive = TRUE)
rdsfile <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/methylation_data_20221006/ml_output/sig_performace_selfcheck_fullsize_output_control.rds"
sig_performace_check_fullsize_output_control_read_and_process(rdsfile,
                                                              outfilepath = paste0(outfilepathmaster, "signature_ml_output/", "trainresults/"))

# # tt2 <- readRDS()
# # tt2a <- tt2[["0.6"]] ## This is the train split number (0.4, 0.6, 0.75)
# # tt2aa <- tt2a[["df"]] ## This is the stat table
# # tt2ab <- tt2a[["output"]] ## This has an output per run ("sig_meth_sevnomild_644"    "sig_rna_sevnomild_1934"    "sig_meth_ctndv70_3v1_1575" "sig_rna_ctndv70_3v1_181"  )
# # tt2aba <- tt2ab[["sig_meth_sevnomild_644"]] 
# ## "cm"                  "auc"                 "predicition"         "performance"         "seed_num"            "feature_coefficient" "feature_importance" 
# tt2aba[["cm"]] ## CM values (and class stats)
# tt2aba[["auc"]] ## AUC for run
# pred1 <- tt2aba[["predicition"]] ## prediction object?
# tt2aba[["performance"]] ## no idea
# tt2aba[["seed_num"]] ## seed
# tt2aba[["feature_coefficient"]] ## coeff per CPG (214 cpgs)
# tt2aba[["feature_importance"]] ## varImp output for lasso
# 
# 
# # library(ROCR)
# # plot(tt2aba[["predicition"]])
# # plot(tt2aba[["performance"]]) + ggtitle("test")
# # plot(tt2aba[["performance"]], main = "test")
# # abline(h = seq(0,1,0.2), lty = "dashed", col = "gray30")
# 
# 
# library(ggplot2)
# library(pROC)
# #define object to plot
# # rocobj <- roc(slot(pred1, "labels")[[1]], tt2aba[["predicition"]])
# # # rocobj <- roc(test$default, predicted)
# # 
# # #create ROC plot
# # ggroc(rocobj)
# roc_plottab <- data.frame(cbind(FPR = slot(tt2aba[["performance"]], "x.values")[[1]], 
#                                 TPR = slot(tt2aba[["performance"]], "y.values")[[1]]))
# pout <- ggplot(roc_plottab, aes(x = FPR, y = TPR))
# pout <- pout + geom_line()
# pout <- pout + theme_pubclean()
# pout




# --------------------------------- Have to recalc the p.adjust for PSEA analyses ---------------------------------
PSEA_filepath <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run4_rmoutliers2_asr_control/meth_processing/run3/asr_corrected_DMP/PSEA/run2"
PSEA_files <- list.files(path = PSEA_filepath, pattern = "*_PSEA_*", recursive = TRUE)
PSEA_files <- PSEA_files[!grepl("_CUSTOM_", PSEA_files)]
for (PSEA_file in PSEA_files) {
    intable <- read.table(paste0(PSEA_filepath, "/", PSEA_file), sep = ",", header = TRUE, row.names = 1)
    intable[,"p.adjust"] <- p.adjust(intable[,"pvalue"], method = "fdr")
    write.table(intable, gsub("_PSEA_", "_PSEA_CUSTOM_", paste0(PSEA_filepath, "/", PSEA_file)), sep = ",", col.names = NA, row.names = TRUE)
}

# --------------------------------- END ---------------------------------




