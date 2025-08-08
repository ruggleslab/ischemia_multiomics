################################
# Name: MacIntosh Cornwell
# Email: mcornwell1957@gmail.com
################################
## Traditional Metric comparisons - comapring DGE and GSEA across imaging, anatomy, and duke score - seeing if traditional metrics can detect useful changes.

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
outfilepathdgecomp = paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run2b_rmoutliers2_controlage/dge_comparison/")
dir.create(outfilepathdgecomp, recursive = TRUE, showWarnings = FALSE)

## Infiles
inmetafile <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run2b_rmoutliers2_controlage/rna_processing/metatable_in.csv"
inmetatable <- read.table(inmetafile, sep = ",", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
# SOI <- rownames(na.omit(inmetatable[,"comp_ischemia__Sev_v_MildNone",drop=FALSE]))
SOI <- rownames(inmetatable)

## Grab the matched samples
incountfile <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run2b_rmoutliers2_controlage/rna_processing/normcounttab.txt"
normcounttable <- read.table(incountfile, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

## Full metafile read in - to grab some more info...
inbioreptablefile <- "/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/data/biorep_10_27.csv"
inbioreptable <- read.table(inbioreptablefile, sep = ",", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("NA", "", NA))
# extraCOI <- c("PATNUM",
#               "HEMOGLOB", "PLATELET", "WBC",
#               "IMGDEGIS",
#               "CTNDV50"
# )
# addonmetatable <- inbioreptable[inbioreptable[,"PATNUM"] %in% colnames(normcounttable),colnames(inbioreptable) %in% extraCOI]
# rownames(addonmetatable) <- addonmetatable[,"PATNUM"]
# addonmetatable <- addonmetatable[,!grepl("PATNUM", colnames(addonmetatable))]

## Read in all of the DGE Files
deseq_imaging_file <-"/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run2b_rmoutliers2_controlage/deseq/comp_ischemia__Sev_v_MildNone/deseq_results_comp_ischemia__Sev_v_MildNone.csv"
deseq_imaging_table <- read.table(deseq_imaging_file, sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
deseq_anatomy_file <-"/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run2b_rmoutliers2_controlage/deseq/comp_anatomy__3v_v_1v/deseq_results_comp_anatomy__3v_v_1v.csv"
deseq_anatomy_table <- read.table(deseq_anatomy_file, sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
deseq_duke_file <-"/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run2b_rmoutliers2_controlage/deseq/comp_duke__67_v_3/deseq_results_comp_duke__67_v_3.csv"
deseq_duke_table <- read.table(deseq_duke_file, sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

gseaHALL_imaging_file <-"/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run2b_rmoutliers2_controlage/gsea/comp_ischemia__Sev_v_MildNone/comp_ischemia__Sev_v_MildNone_gsea_HALL.csv"
gseaHALL_imaging_table <- read.table(gseaHALL_imaging_file, sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
gseaHALL_anatomy_file <-"/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run2b_rmoutliers2_controlage/gsea/comp_anatomy__3v_v_1v/comp_anatomy__3v_v_1v_gsea_HALL.csv"
gseaHALL_anatomy_table <- read.table(gseaHALL_anatomy_file, sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
gseaHALL_duke_file <-"/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/run2b_rmoutliers2_controlage/gsea/comp_duke__67_v_3/comp_duke__67_v_3_gsea_HALL.csv"
gseaHALL_duke_table <- read.table(gseaHALL_duke_file, sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE)




## I would like to get an idea on classifiers - so simple heatmap of the labels for imaging, anatomy, and duke
# "Left Main >=50%" ## SCORE = 7
# "3 Vessels with Severe (>=70%) Plaque or 2 Severe (>=70%) Plaque Including Proximal LAD"  ## SCORE = 6                     
# "3 Vessels with at least Moderate (>=50%) Plaque or 2 Severe (>=70%) Plaque or Severe (>=70%) Proximal LAD Plaque" ## SCORE = 5
# "2 Vessels with at least Moderate (>=50%) Plaque or 1 Severe (>=70%) Plaque" ## SCORE = 4
# "1 Vessel with at least Moderate (>=50%) Plaque" ## SCORE = 3  
# Grab the data we want, make a numeric form for sorting and other purposes.. (?)
labeltable_numeric <- inmetatable[,c("IMGDEGIS", "CTNDV50", "DUKESCORE")]
labeltable_numeric[,"IMGDEGIS"] <- unlist(lapply(labeltable_numeric[,"IMGDEGIS"], function(x) ifelse(x == "None", 0, 
                                                                            ifelse(x == "Mild", 1, 
                                                                            ifelse(x == "Moderate", 2,
                                                                            ifelse(x == "Severe", 3))))))
labeltable_numeric[,"CTNDV50"] <- unlist(lapply(labeltable_numeric[,"CTNDV50"], function(x) ifelse(x == "Non-evaluable", 0, 
                                                                            ifelse(x == "1", 1, 
                                                                            ifelse(x == "2", 2,
                                                                            ifelse(x == "3", 3))))))
labeltable_numeric[,"DUKESCORE"] <- unlist(lapply(labeltable_numeric[,"DUKESCORE"], function(x) 
  ifelse(x == "1 Vessel with at least Moderate (>=50%) Plaque", 3, 
  ifelse(x == "2 Vessels with at least Moderate (>=50%) Plaque or 1 Severe (>=70%) Plaque", 4, 
  ifelse(x == "3 Vessels with at least Moderate (>=50%) Plaque or 2 Severe (>=70%) Plaque or Severe (>=70%) Proximal LAD Plaque", 5,
  ifelse(x == "3 Vessels with Severe (>=70%) Plaque or 2 Severe (>=70%) Plaque Including Proximal LAD", 6,
  ifelse(x == "Left Main >=50%", 7)))))))
labeltable_numeric <- labeltable_numeric[do.call(order, labeltable_numeric),]

# Grab the original labels into a table, and sort by the sorted numeric table
labeltable <- inmetatable[rownames(labeltable_numeric), c("IMGDEGIS", "CTNDV50", "DUKESCORE")]

## Make the annotation (the main plot we want)
IMGDEGIS_color_f <- colorRampPalette(colors = c("#e4b2b2", "#A70000"))
CTNDV50_color_f <- colorRampPalette(colors = c("#cf99d9", "#8700A1"))
DUKESCORE_color_f <- colorRampPalette(colors = c("#99c7a5", "#005115"))
annotationlist1 = annotationlist_builder(labeltable, 
                                         customcolorlist = list(IMGDEGIS = IMGDEGIS_color_f(4),
                                                                CTNDV50 = CTNDV50_color_f(4),
                                                                DUKESCORE = DUKESCORE_color_f(5)
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
pdf(paste0(outfilepathdgecomp, "ischemia_test_labels.pdf"), 11.5, 8, useDingbats = FALSE)
draw(labelplot[[1]])
junk <- dev.off()


## Ok, so our labels dont really match up.. But do the differential expression results compare at all (probably not!)

## Read in all of the DGE Files
deseq_comptab <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "PATNUM", all = TRUE, sort = FALSE), 
  list(cbind.data.frame(PATNUM = rownames(deseq_imaging_table), IMGDEGIS_log2fc = deseq_imaging_table[,c("log2FoldChange")],
                        IMGDEGIS_padj = deseq_imaging_table[,c("padj")]),
       cbind.data.frame(PATNUM = rownames(deseq_anatomy_table), CTNDV50_log2fc = deseq_anatomy_table[,c("log2FoldChange")],
                        CTNDV50_padj = deseq_anatomy_table[,c("padj")]),
       cbind.data.frame(PATNUM = rownames(deseq_duke_table), DUKE_log2fc = deseq_duke_table[,c("log2FoldChange")],
                        DUKE_padj = deseq_duke_table[,c("padj")])
       ))
deseq_complot <- ggpairs(deseq_comptab[,grepl("log2fc", colnames(deseq_comptab))], diag = "blank")
deseq_complot <- deseq_complot + theme_pubr()
deseq_complot <- deseq_complot + geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 0, linetype = 2)

write.table(deseq_comptab, paste0(outfilepathdgecomp, "deseq_comptable.csv"), sep = ",", col.names = TRUE, row.names = FALSE)
pdf(paste0(outfilepathdgecomp, "deseq_corrplot.pdf"), 11.5, 8, useDingbats = FALSE)
print(deseq_complot)
junk <- dev.off()

# deseq_complot <- custom_corrplot(deseqcomptab[2:4])
# deseq_complot <- deseq_complot + geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 0, linetype = 2)

## Now what about the GSEA comps
gsea_comptab <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "GENESET", all = TRUE, sort = FALSE), 
                       list(cbind.data.frame(GENESET = rownames(gseaHALL_imaging_table), IMGDEGIS_NES = gseaHALL_imaging_table[,c("NES")], 
                                             IMGDEGIS_padj = gseaHALL_imaging_table[,c("p.adjust")]),
                            cbind.data.frame(GENESET = rownames(gseaHALL_anatomy_table), CTNDV50_NES = gseaHALL_anatomy_table[,c("NES")],
                                             CTNDV50_padj = gseaHALL_anatomy_table[,c("p.adjust")]),
                            cbind.data.frame(GENESET = rownames(gseaHALL_duke_table), DUKE_NES = gseaHALL_duke_table[,c("NES")],
                                             DUKE_padj = gseaHALL_duke_table[,c("p.adjust")])
                            ))
gseapvalcutoff <- 0.05
gseaHALLuplist <- list(
  IMGDEGIS_up = gsea_comptab[gsea_comptab[,"IMGDEGIS_NES"] > 0 & gsea_comptab[,"IMGDEGIS_padj"] < gseapvalcutoff,1],
  CTNDV50_up = gsea_comptab[gsea_comptab[,"CTNDV50_NES"] > 0 & gsea_comptab[,"CTNDV50_padj"] < gseapvalcutoff,1],
  DUKE_up = gsea_comptab[gsea_comptab[,"DUKE_NES"] > 0 & gsea_comptab[,"DUKE_padj"] < gseapvalcutoff,1]
)
gseaHALLdownlist <- list(
  IMGDEGIS_down = gsea_comptab[gsea_comptab[,"IMGDEGIS_NES"] < 0 & gsea_comptab[,"IMGDEGIS_padj"] < gseapvalcutoff,1],
  CTNDV50_down = gsea_comptab[gsea_comptab[,"CTNDV50_NES"] < 0 & gsea_comptab[,"CTNDV50_padj"] < gseapvalcutoff,1],
  DUKE_down = gsea_comptab[gsea_comptab[,"DUKE_NES"] < 0 & gsea_comptab[,"DUKE_padj"] < gseapvalcutoff,1]
)
gseaHALLup_overlapout <- overlap_finder(overlap_inlist = gseaHALLuplist)
write.table(gseaHALLup_overlapout[[1]], paste0(outfilepathdgecomp, "gsea_HALLup_comptable.csv"), sep = ",", col.names = TRUE, row.names = FALSE)
write.table(gseaHALLup_overlapout[[4]], paste0(outfilepathdgecomp, "gsea_HALLup_summary.csv"), sep = ",", col.names = NA, row.names = TRUE)
pdf(paste0(outfilepathdgecomp, "gsea_HALLup_vennplot.pdf"), 7, 7, useDingbats = FALSE)
grid.draw(gseaHALLup_overlapout[[2]])
junk <- dev.off()

gseaHALLdown_overlapout <- overlap_finder(overlap_inlist = gseaHALLdownlist)
write.table(gseaHALLdown_overlapout[[1]], paste0(outfilepathdgecomp, "gsea_HALLdown_comptable.csv"), sep = ",", col.names = TRUE, row.names = FALSE)
write.table(gseaHALLdown_overlapout[[4]], paste0(outfilepathdgecomp, "gsea_HALLdown_summary.csv"), sep = ",", col.names = NA, row.names = TRUE)
pdf(paste0(outfilepathdgecomp, "gsea_HALLdown_vennplot.pdf"), 7, 7, useDingbats = FALSE)
grid.draw(gseaHALLdown_overlapout[[2]])
junk <- dev.off()






## I need to see if the DESeq for ischemia holds up to a non-parametric test, I think its drive by outliers in each gene, and theres just more samples, so you just get more outliers...
group1samples <- rownames(inmetatable[inmetatable[,"comp_ischemia__Sev_v_MildNone"] %in% 1,,drop=FALSE])
group2samples <- rownames(inmetatable[inmetatable[,"comp_ischemia__Sev_v_MildNone"] %in% 0,,drop=FALSE])
outlist <- list()
for (genenum in seq_len(nrow(normcounttable))) {
   genegroup1 <- as.vector(unlist(normcounttable[genenum,group1samples]))
   genegroup2 <- as.vector(unlist(normcounttable[genenum,group2samples]))
   outlist[[genenum]] <- c(gene = rownames(normcounttable)[genenum], pval = wilcox.test(genegroup1, genegroup2)$p.value, 
                           group1mean = mean(genegroup1), group2mean = mean(genegroup2), log2fc = log2(mean(genegroup1)/mean(genegroup2)))
}
wilcoxtesttab <- do.call(rbind, outlist)
wilcoxtesttab[,2:5] <- apply(wilcoxtesttab[,2:5], 2, as.numeric)
wilcoxtesttab <- data.frame(wilcoxtesttab, stringsAsFactors = FALSE)
wilcoxtesttab[,"wilcox_padj"] <- p.adjust(wilcoxtesttab[,"pval"], method = "fdr")
# write.table(wilcoxtesttab, paste0(outfilepathdgecomp, "wilcoxtesttab.csv"), sep = ",", col.names = TRUE, row.names = FALSE)

# Plot out some of the genes that are MOST different
GOItab <- deseq_imaging_table[order(deseq_imaging_table[,"stat"], decreasing = TRUE),]
GOI <- c(rownames(GOItab)[1:20], "RSAD2", "SERPING1")
topGOIboxplots_outpath <- paste0(outfilepathdgecomp, "topGOIboxplots/IMGDEGIS/")
dir.create(topGOIboxplots_outpath, recursive = TRUE, showWarnings = FALSE)
for (genenum in seq_len(length(GOI))) {
  genesel <- GOI[genenum]
  bptab <- merge(inmetatable[,"IMGDEGIS",drop=FALSE], t(normcounttable[genesel,]), by = "row.names")
  bptab[,1] <- "IMGDEGIS"
  
  bpout <- boxplot_plotter(boxplottable = bptab, xsplit = "feature", 
                           labsparam = list(title = paste0(genesel, " values for samples by ", "IMGDEGIS"), x = "IMGDEGIS", y = genesel, 
                                            featorder = "IMGDEGIS", catorder = c("None", "Mild", "Moderate", "Severe")), 
                           plotstats = "intra", testtypeparam = "wilcox.test"
  )
  pdf(paste0(topGOIboxplots_outpath, genesel, "_boxplot.pdf"), useDingbats = FALSE)
  print(bpout)
  junk <- dev.off()
}



## For each of the DGE genes (for IMGDEGIS severe vs mildnone) we want to see if the population with the top quartile for each gene is the same?
## Or can a collection of genes determine the severe population (with some genes for each severe)
# Final result is a HM with genes by samples, and then for each square is the quartile.
GOI <- rownames(deseq_imaging_table[deseq_imaging_table[,"padj"] < 0.01 & abs(deseq_imaging_table[,"log2FoldChange"]) > 0.25, ])
# GOI <- rownames(deseq_imaging_table[order(deseq_imaging_table[,"stat"], decreasing = TRUE), ][1:500,])
SOI <- rownames(na.omit(inmetatable[,"IMGDEGIS",drop=FALSE]))
subtab <- normcounttable[GOI,SOI]

numberquantiles = 2
gene_quantile_table <- t(apply(subtab, 1, function(x) {
  temptab <- cut2(x, g = numberquantiles)
  for (quantilenum in seq_len(numberquantiles)) {
    levels(temptab)[match(levels(temptab)[quantilenum],levels(temptab))] <- quantilenum
  }
  return(as.numeric(as.character(temptab)))
}))
colnames(gene_quantile_table) <- colnames(subtab)

## Plot it out as a hm
heatmapcolorparam <- colorRamp2(seq_len(numberquantiles), c("blue", rep("grey", numberquantiles - 2), "red"))
rowmetatable <- inmetatable[SOI,"IMGDEGIS",drop=FALSE]
rowmetatable <- rowmetatable[order(match(rowmetatable[,1], c("Severe", "Moderate", "Mild", "None"))),,drop=FALSE]
rowannotationlist <- annotationlist_builder(metatable = rowmetatable, 
                                            customcolorlist = list(IMGDEGIS = c(Severe = "darkred", Moderate = "red", Mild = "blue", None = "darkblue")))

hmplottab <- t(gene_quantile_table[,rownames(rowmetatable)])

# out1 <- create_heatmap(counttab = hmplottab, subsetnum = FALSE, scale_data = FALSE,
#                        rowmetatable = rowmetatable, rowannotationlist = rowannotationlist,
#                        colclusterparam = TRUE, rowclusterparam = FALSE, heatmapcolorparam = heatmapcolorparam)
out1 <- Heatmap(matrix = hmplottab,
        col = heatmapcolorparam,
        row_split = rowmetatable[,1],cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = FALSE, show_column_names = FALSE)
pdf(paste0(outfilepathdgecomp, "quantile_heatmap.pdf"))
draw(out1)
junk <- dev.off()

## So really what we want is for each Sample, we want to see which genes it has an upper quantile of expression


## Quantiles doesnt make any sense........ lets just take the GOI table and scale each genes expression...
scaled_subtab <- apply(subtab, 1, zscore)[rownames(rowmetatable),] # already transposes it in the operation
heatmapcolorparam <- colorRamp2(breaks = c(-3,0,3), c("blue", "white", "red"))
out1 <- Heatmap(matrix = scaled_subtab,
                col = heatmapcolorparam,
                row_split = rowmetatable[,1],cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = FALSE, show_column_names = FALSE)
pdf(paste0(outfilepathdgecomp, "scaled_heatmap.pdf"))
draw(out1)
junk <- dev.off()

