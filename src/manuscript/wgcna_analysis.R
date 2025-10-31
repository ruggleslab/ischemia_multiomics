###########################################################################
#
#                       ISCHEMIA WGCNA Analysis
#
###########################################################################
# Author: ISCHEMIA Study Team
# Date: Updated 2025-08-08
# Description: Weighted Gene Co-expression Network Analysis (WGCNA) for
#              ISCHEMIA whole blood RNA-seq data

# Load configuration and utilities
source("config.R")
source("utils.R")

# Load required packages for WGCNA
load_packages(c("WGCNA", "flashClust", "dynamicTreeCut"))

# Input files
inmetafile <- file.path(DATA_DIR, "metatable_in.csv")
inmetatable <- read.table(inmetafile, sep = ",", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
SOI <- rownames(inmetatable)

# Grab the matched samples
incountfile <- file.path(DATA_DIR, "normcounttab.txt")
normcounttable <- read.table(incountfile, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

# Create output directory for WGCNA
outfilepathwgcna <- create_output_dir(paste0("WGCNA_power", WGCNA_POWER))
powerparam = 14
modulesizeparam = 30
# FIXME: Update path - outfilepathwgcna = paste0("# PATH_UPDATED: output/run3_rmoutliers2/WGCNA/WGCNA_power", powerparam, "_size", modulesizeparam,"/")
dir.create(outfilepathwgcna, recursive = TRUE, showWarnings = FALSE)

## Step 0 - add metadata on that I want, turning some characters into numbers
# "Left Main >=50%" ## SCORE = 7
# "3 Vessels with Severe (>=70%) Plaque or 2 Severe (>=70%) Plaque Including Proximal LAD"  ## SCORE = 6                     
# "3 Vessels with at least Moderate (>=50%) Plaque or 2 Severe (>=70%) Plaque or Severe (>=70%) Proximal LAD Plaque" ## SCORE = 5
# "2 Vessels with at least Moderate (>=50%) Plaque or 1 Severe (>=70%) Plaque" ## SCORE = 4
# "1 Vessel with at least Moderate (>=50%) Plaque" ## SCORE = 3   
align_metatable_reflist <- list(
    "DEGRISCH" = list("0" = "None", "1" = "Mild", "2" = "Moderate", "3" = "Severe", "NA" = "Uninterpretable"),
    "IMGDEGIS" = list("0" = "None", "1" = "Mild", "2" = "Moderate", "3" = "Severe", "NA" = "Uninterpretable"),
    "CTNDV50" = list("1" = "1", "2" = "2", "3" = "3", "NA" = c(NA, "Non-evaluable")),
    "CTMULT50" = list("0" = "No", "1" = "Yes", "NA" = c("Not evaluable", NA)),
    "DUKESCORE" = list("7" = "Left Main >=50%", 
                       "6" = "3 Vessels with Severe (>=70%) Plaque or 2 Severe (>=70%) Plaque Including Proximal LAD", 
                       "5" = "3 Vessels with at least Moderate (>=50%) Plaque or 2 Severe (>=70%) Plaque or Severe (>=70%) Proximal LAD Plaque", 
                       "4" = "2 Vessels with at least Moderate (>=50%) Plaque or 1 Severe (>=70%) Plaque", 
                       "3" = "1 Vessel with at least Moderate (>=50%) Plaque"))
new_numeric_meta <- apply(align_metatable(intable = inmetatable[,c("DEGRISCH", "IMGDEGIS", "CTNDV50", "CTMULT50", "DUKESCORE")], align_metatable_reflist), 
                          2, as.numeric)
inmetatable[,c("DEGRISCH", "IMGDEGIS", "CTNDV50", "CTMULT50", "DUKESCORE")] <- new_numeric_meta


# Step 1 - find modules
# WGCNAcounttable <- t(normcounttable)
# WGCNAmetatable <- inmetatable[,c("Cohort", "sex", "race", "age", "bmi", "smoking", "diabetes", "hypertension", "hyperlipidemia",
#                                  colnames(inmetatable)[grepl("censor_", colnames(inmetatable))])]
WGCNAcounttable <- t(normcounttable[,SOI])
WGCNAmetatable <- inmetatable[SOI,c(colnames(inmetatable))]

WGCNA_out1 <- create_WGCNA_modules(counttable = WGCNAcounttable, 
                                   powerparam = powerparam, 
                                   minModuleSizeParam = modulesizeparam,
                                   savedendroplotfile = paste0(outfilepathwgcna, "dendrogram_color_plot.pdf"))
powerplot <- WGCNA_out1[[1]]
genestocolors <- WGCNA_out1[[2]]
moduleColors <- genestocolors[,"moduleColors"]

pdf(paste0(outfilepathwgcna, "wgcna_powerplot.pdf"))
print(powerplot)
junk <- dev.off()
write.table(genestocolors, paste0(outfilepathwgcna, "wgcna_genestocolors.csv"), sep = ",", col.names = NA, row.names = TRUE)

# Step 2 - annotate modules
WGCNA_out2_C5 <- annotate_WGCNA_modules(genestocolors, genesetparam = "C5")
write.table(WGCNA_out2_C5, paste0(outfilepathwgcna, "C5_wgcna_module_annotations.csv"), sep = ",", col.names = NA, row.names = TRUE)
WGCNA_out2_H <- annotate_WGCNA_modules(genestocolors, genesetparam = "H")
write.table(WGCNA_out2_H, paste0(outfilepathwgcna, "H_wgcna_module_annotations.csv"), sep = ",", col.names = NA, row.names = TRUE)

# Step 3 calc corr with metadata
WGCNA_out3 <- WGNCA_module_feature_correlation(counttable = WGCNAcounttable, moduleColors = moduleColors, metatable = WGCNAmetatable)
moduleTraitCor <- WGCNA_out3$moduleTraitCor
moduleTraitPvalue <- WGCNA_out3$moduleTraitPvalue
write.table(moduleTraitCor, paste0(outfilepathwgcna, "wgcna_moduleTraitCor.csv"), sep = ",", col.names = NA, row.names = TRUE)
write.table(moduleTraitPvalue, paste0(outfilepathwgcna, "wgcna_moduleTraitPvalue.csv"), sep = ",", col.names = NA, row.names = TRUE)


# Step 4 - create custom heatmap from output
WGNCA_hm <- WGCNA_custom_heatmap(t(moduleTraitPvalue), t(moduleTraitCor))
pdf(paste0(outfilepathwgcna, "wgcna_corr_hm.pdf"), 15, 40)
draw(WGNCA_hm)
junk <- dev.off()

## Now lets look at each of these as boxplots - more informative than the correlatoin values
eigengenes = orderMEs(moduleEigengenes(WGCNAcounttable, moduleColors)$eigengenes)
write.table(eigengenes, paste0(outfilepathwgcna, "wgcna_eigengenes.csv"), sep = ",", col.names = NA, row.names = TRUE)



## Now lets look at each of these as boxplots - more informative than the correlatoin values
eigengenes <- read.table(paste0(outfilepathwgcna, "wgcna_eigengenes.csv"), sep = ",", header = TRUE, row.names =1)
# eigengenebpplot_path <- paste0(outfilepathwgcna, "eigengene_boxplots/")
eigengenebpplot_path <- paste0(outfilepathwgcna, "eigengene_boxplots_noNAs/")
dir.create(eigengenebpplot_path, showWarnings = FALSE, recursive = TRUE)

## I want to do the bptab for each comparison and each module... So need to add a forloop around this with our categories of interest...
# FIXME: Update path - # inmetafile <- "# PATH_UPDATED: output/run1c_rmoutliers2/rna_processing/metatable_in.csv"
# inmetatable <- read.table(inmetafile, sep = ",", header = TRUE, stringsAsFactors = FALSE, row.names = 1)

## I want to add some more categories here - namely combining some metrics
inmetatable[,"IMGDEGIS_01_2_3"] <- inmetatable[,"IMGDEGIS",drop=FALSE]
inmetatable[inmetatable[,"IMGDEGIS_01_2_3"] %in% c("0", "1"),"IMGDEGIS_01_2_3"] <- "01"
inmetatable[,"IMGDEGIS_01_23"] <- inmetatable[,"IMGDEGIS",drop=FALSE]
inmetatable[inmetatable[,"IMGDEGIS_01_23"] %in% c("0", "1", "2", "3"),"IMGDEGIS_01_23"] <- ifelse(
    inmetatable[inmetatable[,"IMGDEGIS_01_23"] %in% c("0", "1", "2", "3"),"IMGDEGIS_01_23"] %in% c("0", "1"), "01", "23")

# COI <- list("IMGDEGIS" = c(NA, "None", "Mild", "Moderate", "Severe"))
# COI <- list("IMGDEGIS" = c("NA", "0", "1", "2", "3"),
#             "IMGDEGIS_01_2_3" = c("NA", "01", "2", "3"),
#             "IMGDEGIS_01_23" = c("NA", "01", "23"))
COI <- list("IMGDEGIS" = c("0", "1", "2", "3"),
            "IMGDEGIS_01_2_3" = c("01", "2", "3"),
            "IMGDEGIS_01_23" = c("01", "23"))

for (COInum in seq_len(length(COI))) {
    
    ## Select the category to analyze
    COIsel <- names(COI)[COInum]
    metasel <- inmetatable[,COIsel,drop=FALSE]
    metasel[,1] <- as.character(metasel[,1])
    # metasel[is.na(metasel[,1]),1] <- "NA"
    
    dir.create(paste0(eigengenebpplot_path, COIsel, "/"), showWarnings = FALSE, recursive = TRUE)
    
    ## Then run through each eigenodule
    bptablist <- wilcox_statoutlist <- spearman_statoutlist <- list()
    for (eigengenenum in seq_len(ncol(eigengenes))) {
        eigengene_sel <- eigengenes[,eigengenenum,drop=FALSE]
        eigengene_sellabel <- colnames(eigengene_sel)
        
        bptab <- merge(metasel, eigengene_sel, by = "row.names")[,c(2,1,3)]
        bptab[,2] <- eigengene_sellabel
        bptab[,1] <- as.character(bptab[,1])
        colnames(bptab) <- c("cohort", "eigengene", "value")
        # bpout <- boxplot_plotter(boxplottable = bptab, xsplit = "feature", 
        #                          labsparam = list(title = paste0(eigengene_sellabel, " vs ", COIsel), x = COIsel, y = eigengene_sellabel, 
        #                                           catorder = eigengene_sellabel, featorder = COI[COInum][[1]]), 
        #                          plotstats = "inter", testtypeparam = "wilcox.test"
        #                          # comparisonparam = combn(COI[COInum][[1]], 2, simplify=F)
        #                          )
        bpout <- boxplot_plotter(boxplottable = bptab, xsplit = "category", 
                                 labsparam = list(title = paste0(eigengene_sellabel, " vs ", COIsel), x = COIsel, y = eigengene_sellabel, 
                                                  catorder = eigengene_sellabel, featorder = COI[COInum][[1]]), 
                                 plotstats = "intra", testtypeparam = "wilcox.test"
                                 # comparisonparam = combn(COI[COInum][[1]], 2, simplify=F)
        )
        
        pdf(paste0(eigengenebpplot_path, COIsel, "/", colnames(eigengene_sel), "_boxplot.pdf"))
        print(bpout)
        junk <- dev.off()
        
        scattertab <- bptab[,c(1,3,2)]
        scattertab[,1] <- as.numeric(scattertab[,1])
        scattertab[,3] <- gsub("ME", "", scattertab[,3])
        scatterout <- scatter_plotter(indata = scattertab[,c(1,2)], colorvar = scattertab[,3,drop=FALSE], 
                                      labsparam = list(title = paste0(eigengene_sellabel, " vs ", COIsel), x = COIsel, y = eigengene_sellabel), 
                                      plotstats = TRUE)
        pdf(paste0(eigengenebpplot_path, COIsel, "/", eigengene_sellabel, "_scatter.pdf"))
        print(scatterout)
        junk <- dev.off()
        
        ## Save out the table
        bptablist[[eigengenenum]] <- bptab
        
        ## Grab the stats into a table
        # Create the combinatorial table for all combos of the groups
        combtab = combn(as.character(unique(bptab[!is.na(bptab[,1]),1])), 2, simplify=F)
        # Apply a functiona cross each combo using the bptab and pulling out each group
        wilcoxstatsout <- lapply(combtab, function(x){
            group1 <- bptab[bptab[,1] %in% x[1], 3]
            group2 <- bptab[bptab[,1] %in% x[2], 3]
            pval_out <- c(COIsel, x[1], x[2], eigengene_sellabel, wilcox.test(group1, group2)$p.value)
            cohen_out <- cohens_d(group1, group2)
            sign_out <- sign(mean(group2, na.rm = TRUE)-mean(group1, na.rm = TRUE))
            c(pval_out, cohen_out, sign_out)
        })
        wilcox_statoutlist[[eigengenenum]] <- do.call(rbind, wilcoxstatsout)
        
        # Repeat for the spearman correlation
        spearmanstatout <- cor.test(scattertab[,c(1)], scattertab[,c(2)], method = "spearman", exact = FALSE)
        spearman_statoutlist[[eigengenenum]] <- c(COIsel, eigengene_sellabel, spearmanstatout$p.value, spearmanstatout$estimate)
        
        names(bptablist)[eigengenenum] <- names(wilcox_statoutlist)[eigengenenum] <- names(spearman_statoutlist)[eigengenenum] <- eigengene_sellabel

    }
    wilcox_summarytable <- do.call(rbind, wilcox_statoutlist)
    colnames(wilcox_summarytable) <- c("Category", "Feature1", "Feature2", "Eigengene", "wilcox_pval", "cohens_d", "cohen_sign")
    spearman_summarytable <- do.call(rbind, spearman_statoutlist)
    colnames(spearman_summarytable) <- c("Category", "Eigengene", "spearman_pval", "spearman_rval")
    write.table(wilcox_summarytable, paste0(eigengenebpplot_path, COIsel, "/", "wilcox_summary_table.csv"), 
                sep = ",", col.names = TRUE, row.names = FALSE)
    write.table(spearman_summarytable, paste0(eigengenebpplot_path, COIsel, "/", "spearman_summary_table.csv"), 
                sep = ",", col.names = TRUE, row.names = FALSE)
    
}


wilcox_summarytable <- read.table(outfilepathwgcna, "eigengene_boxplots_noNAs/IMGDEGIS_01_2_3/wilcox_summary_table.csv", header = TRUE, sep = ",")

wilcox_plottable_temp <- wilcox_summarytable[apply(wilcox_summarytable[,c("Feature1", "Feature2")], 1, function(x) 
    paste0(x, collapse = "__")) %in% c("1__3","3__1"),]
wilcox_plottable <- data.frame(value = wilcox_plottable_temp[,"cohens_d"] * ifelse(wilcox_plottable_temp[,"Feature1"] == 3, -wilcox_plottable_temp[,"cohen_sign"], wilcox_plottable_temp[,"cohen_sign"]), row.names = wilcox_plottable_temp[,"Eigengene"])

heatmapcolorparam = colorRamp2(breaks = c(1, 0, -1), c("darkred", "white", "darkblue"))
rowmetatable <- data.frame(sigpval = wilcox_plottable_temp[,"wilcox_pval"] < 0.05, row.names = wilcox_plottable_temp[,"Eigengene"])
rowannotationlist <- annotationlist_builder(rowmetatable)
hmout <- create_heatmap(counttab = wilcox_plottable, scale_data = FALSE, 
                        rowmetatable = rowmetatable, rowannotationlist = rowannotationlist, heatmapcolorparam = heatmapcolorparam)
pdf(paste0(eigengenebpplot_path, "IMGDEGIS_01_2_3/", "summary_heatmap.pdf"))
draw(hmout[[1]])
junk <- dev.off()


## Ok what about WGCNA with custom modules as well:
WGCNA_metatable <- inmetatable[,c(grepl("IMGDEGIS", colnames(inmetatable)))]
WGCNA_counttable <- as.data.frame(t(normcounttable))
WGCNA_counttable <- WGCNA_counttable[rownames(WGCNA_metatable),]

custommodule_outfilepath <- paste0(outfilepathwgcna, "custom_modules/")
dir.create(custommodule_outfilepath, showWarnings = FALSE, recursive = TRUE)

speciesparam = "Homo sapiens"
# keggpathwaytab = as.data.frame(msigdbr(species = speciesparam, subcategory = c("CP:KEGG"))[,c("gs_name", "gene_symbol")])
# hallmarkpathwaytab = as.data.frame(msigdbr(species = speciesparam, category = c("H"))[,c("gs_name", "gene_symbol")])
gopathwaytab = as.data.frame(msigdbr(species = speciesparam, category = c("C5"))[,c("gs_name", "gene_symbol")])
pathwaytab = gopathwaytab

# Let's make some custom gene sets from GO terms:
make_supergs <- function(pathwaytab, termvector) {
    outlist <- list()
    for (termnum in seq_len(length(termvector))) {
        term <- termvector[termnum]
        outlist[[termnum]] <- cbind(gs_name = paste0(term, "_SUPERGS"), 
                         gene_symbol = unique(pathwaytab[pathwaytab[,1] %in% unique(pathwaytab[grepl(term, pathwaytab[,1]),1]),2]))
    }
    return(do.call(rbind.data.frame, outlist))
}
SUPERGS_termvector <- c("INFLAMMATORY", "THROMBOSIS", "FIBROSIS", "COAGULATION", "COMPLEMENT", "CHOLESTEROL")
superpathwaytab <- make_supergs(pathwaytab, termvector = SUPERGS_termvector)
POIlist <- paste0(SUPERGS_termvector, "_SUPERGS")

# Run over each superGS and run WGCNA
POIoutlist_R <- list()
POIoutlist_P <- list()
eigengeneoutlist <- list()
for (POInum in seq_len(length(POIlist))) {
    POIlabel <- POIlist[POInum]
    POI <- superpathwaytab[superpathwaytab[,1] == POIlabel,]
    
    #### genetocolors - has rownames genes, moduleLabels as number module, and moduleColors as color name
    genestocolors <- data.frame(moduleLabels = ifelse(colnames(WGCNA_counttable) %in% POI[,2], 1, 0), 
                                moduleColors = ifelse(colnames(WGCNA_counttable) %in% POI[,2], POIlabel, "grey"), 
                                row.names = colnames(WGCNA_counttable), stringsAsFactors = FALSE)
    moduleColors <- genestocolors[,"moduleColors"]
    # write.table(genestocolors, paste0(outfilepath, "genestocolors.csv"), sep = ",", col.names = NA, row.names = TRUE)
    
    # WGCNA_out3 <- WGNCA_module_feature_correlation(counttable, moduleColors, metatable, outfilepath)
    WGCNA_out3 <- WGNCA_module_feature_correlation(WGCNA_counttable, moduleColors, WGCNA_metatable)
    moduleTraitCor <- WGCNA_out3$moduleTraitCor
    moduleTraitPvalue <- WGCNA_out3$moduleTraitPvalue
    eigengenes = orderMEs(moduleEigengenes(WGCNA_counttable, moduleColors)$eigengenes)
    
    POIoutlist_R[[POInum]] <- moduleTraitCor[paste0("ME", POIlabel),]
    POIoutlist_P[[POInum]] <- moduleTraitPvalue[paste0("ME", POIlabel),]
    eigengeneoutlist[[POInum]] <- eigengenes[,paste0("ME", POIlabel),drop=FALSE]
    names(eigengeneoutlist)[POInum] <- names(POIoutlist_R)[POInum] <- names(POIoutlist_P)[POInum] <- POIlabel
}

moduleTraitCor <- POI_CorTab <- do.call(rbind, POIoutlist_R)
moduleTraitPvalue <- POI_PTab <- do.call(rbind, POIoutlist_P)
eigengene_table <- do.call(cbind, eigengeneoutlist)
colnames(eigengene_table) <- gsub("^ME", "", colnames(eigengene_table))
write.table(POI_CorTab, paste0(custommodule_outfilepath, "supergs_rval_tab.csv"), sep = ",", col.names = NA, row.names = TRUE)
write.table(POI_PTab, paste0(custommodule_outfilepath, "supergs_pval_tab.csv"), sep = ",", col.names = NA, row.names = TRUE)
write.table(eigengene_table, paste0(custommodule_outfilepath, "supergs_eigengene_tab.csv"), sep = ",", col.names = NA, row.names = TRUE)

# eigengenebpplot_path <- paste0(outfilepathwgcna, "eigengene_boxplots/")
eigengenebpplot_path <- paste0(custommodule_outfilepath, "eigengene_boxplots/")
dir.create(eigengenebpplot_path, showWarnings = FALSE, recursive = TRUE)

## I want to do the bptab for each comparison and each module... So need to add a forloop around this with our categories of interest...
# FIXME: Update path - # inmetafile <- "# PATH_UPDATED: output/run1c_rmoutliers2/rna_processing/metatable_in.csv"
# inmetatable <- read.table(inmetafile, sep = ",", header = TRUE, stringsAsFactors = FALSE, row.names = 1)

## I want to add some more categories here - namely combining some metrics
inmetatable[,"IMGDEGIS_01_2_3"] <- inmetatable[,"IMGDEGIS",drop=FALSE]
inmetatable[,"IMGDEGIS_01_2_3"] <- inmetatable[,"IMGDEGIS",drop=FALSE]
inmetatable[inmetatable[,"IMGDEGIS_01_2_3"] %in% c("0", "1"),"IMGDEGIS_01_2_3"] <- "01"
inmetatable[,"IMGDEGIS_01_23"] <- inmetatable[,"IMGDEGIS",drop=FALSE]
inmetatable[inmetatable[,"IMGDEGIS_01_23"] %in% c("0", "1", "2", "3"),"IMGDEGIS_01_23"] <- ifelse(
    inmetatable[inmetatable[,"IMGDEGIS_01_23"] %in% c("0", "1", "2", "3"),"IMGDEGIS_01_23"] %in% c("0", "1"), "01", "23")

# COI <- list("IMGDEGIS" = c(NA, "None", "Mild", "Moderate", "Severe"))
# COI <- list("IMGDEGIS" = c("NA", "0", "1", "2", "3"),
#             "IMGDEGIS_01_2_3" = c("NA", "01", "2", "3"),
#             "IMGDEGIS_01_23" = c("NA", "01", "23"))
COI <- list("IMGDEGIS" = c("0", "1", "2", "3"),
            "IMGDEGIS_01_2_3" = c("01", "2", "3"),
            "IMGDEGIS_01_23" = c("01", "23"))

for (COInum in seq_len(length(COI))) {
    
    ## Select the category to analyze
    COIsel <- names(COI)[COInum]
    metasel <- inmetatable[,COIsel,drop=FALSE]
    metasel[,1] <- as.character(metasel[,1])
    # metasel[is.na(metasel[,1]),1] <- "NA"
    
    dir.create(paste0(eigengenebpplot_path, COIsel, "/"), showWarnings = FALSE, recursive = TRUE)
    
    ## Then run through each eigenodule
    bptablist <- wilcox_statoutlist <- spearman_statoutlist <- list()
    for (eigengenenum in seq_len(ncol(eigengene_table))) {
        eigengene_sel <- eigengene_table[,eigengenenum,drop=FALSE]
        eigengene_sellabel <- colnames(eigengene_sel)
        
        bptab <- merge(metasel, eigengene_sel, by = "row.names")[,c(2,1,3)]
        bptab[,2] <- eigengene_sellabel
        bptab[,1] <- as.character(bptab[,1])
        colnames(bptab) <- c("cohort", "eigengene", "value")
        # bpout <- boxplot_plotter(boxplottable = bptab, xsplit = "feature", 
        #                          labsparam = list(title = paste0(eigengene_sellabel, " vs ", COIsel), x = COIsel, y = eigengene_sellabel, 
        #                                           catorder = eigengene_sellabel, featorder = COI[COInum][[1]]), 
        #                          plotstats = "inter", testtypeparam = "wilcox.test"
        #                          # comparisonparam = combn(COI[COInum][[1]], 2, simplify=F)
        #                          )
        bpout <- boxplot_plotter(boxplottable = bptab, xsplit = "category", 
                                 labsparam = list(title = paste0(eigengene_sellabel, " vs ", COIsel), x = COIsel, y = eigengene_sellabel, 
                                                  catorder = eigengene_sellabel, featorder = COI[COInum][[1]]), 
                                 plotstats = "intra", testtypeparam = "wilcox.test"
                                 # comparisonparam = combn(COI[COInum][[1]], 2, simplify=F)
        )
        
        pdf(paste0(eigengenebpplot_path, COIsel, "/", colnames(eigengene_sel), "_boxplot.pdf"))
        print(bpout)
        junk <- dev.off()
        
        scattertab <- bptab[,c(1,3,2)]
        scattertab[,1] <- as.numeric(scattertab[,1])
        scattertab[,3] <- gsub("ME", "", scattertab[,3])
        scatterout <- scatter_plotter(indata = scattertab[,c(1,2)], colorvar = scattertab[,3,drop=FALSE], 
                                      labsparam = list(title = paste0(eigengene_sellabel, " vs ", COIsel), x = COIsel, y = eigengene_sellabel), 
                                      plotstats = TRUE)
        pdf(paste0(eigengenebpplot_path, COIsel, "/", eigengene_sellabel, "_scatter.pdf"))
        print(scatterout)
        junk <- dev.off()
        
        ## Save out the table
        bptablist[[eigengenenum]] <- bptab
        
        ## Grab the stats into a table
        # Create the combinatorial table for all combos of the groups
        combtab = combn(as.character(unique(bptab[!is.na(bptab[,1]),1])), 2, simplify=F)
        # Apply a functiona cross each combo using the bptab and pulling out each group
        wilcoxstatsout <- lapply(combtab, function(x){
            group1 <- bptab[bptab[,1] %in% x[1], 3]
            group2 <- bptab[bptab[,1] %in% x[2], 3]
            out <- c(COIsel, x[1], x[2], eigengene_sellabel, wilcox.test(group1, group2)$p.value)
            out
        })
        wilcox_statoutlist[[eigengenenum]] <- do.call(rbind, wilcoxstatsout)
        
        # Repeat for the spearman correlation
        spearmanstatout <- cor.test(scattertab[,c(1)], scattertab[,c(2)], method = "spearman", exact = FALSE)
        spearman_statoutlist[[eigengenenum]] <- c(COIsel, eigengene_sellabel, spearmanstatout$p.value, spearmanstatout$estimate)
        
        names(bptablist)[eigengenenum] <- names(wilcox_statoutlist)[eigengenenum] <- names(spearman_statoutlist)[eigengenenum] <- eigengene_sellabel
        
    }
    wilcox_summarytable <- do.call(rbind, wilcox_statoutlist)
    colnames(wilcox_summarytable) <- c("Category", "Feature1", "Feature2", "Eigengene", "wilcox_pval")
    spearman_summarytable <- do.call(rbind, spearman_statoutlist)
    colnames(spearman_summarytable) <- c("Category", "Eigengene", "spearman_pval", "spearman_rval")
    write.table(wilcox_summarytable, paste0(eigengenebpplot_path, COIsel, "/", "wilcox_summary_table.csv"), 
                sep = ",", col.names = TRUE, row.names = FALSE)
    write.table(spearman_summarytable, paste0(eigengenebpplot_path, COIsel, "/", "spearman_summary_table.csv"), 
                sep = ",", col.names = TRUE, row.names = FALSE)
}













# ## I need to think about how to add a side annotation or the genesets that are a part of each module
# fullannotationtable <- rbind(WGCNA_out2_H)
# # fullannotationtable <- rbind(WGCNA_out2_H, WGCNA_out2_C5)
# # modulesel <- WGCNA_out2[grepl("PLATELET", WGCNA_out2[,"ID"]),c("module", "p.adjust", "ID")]
# modulesel <- na.omit(fullannotationtable[fullannotationtable[,"p.adjust"] < 0.05,c("module", "p.adjust", "ID")])
# 
# # add on the other colors with blanks
# if (sum(!unique(fullannotationtable[,"module"]) %in% modulesel[,"module"]) != 0){
#     blankmoduletab <- data.frame(module = unique(fullannotationtable[,"module"])[!unique(fullannotationtable[,"module"]) %in% modulesel[,"module"]],
#                                  p.adjust = NA, ID = modulesel[1,"ID"])
# } else {blankmoduletab <- NULL}
# 
# dotplottab <- rbind(modulesel, blankmoduletab)
# dotplottab[,3] <- gsub("HALLMARK_|GO_", "", dotplottab[,3])
# 
# ## Reorder to match the order of eigengenes
# colorOrder <- gsub("ME", "", colnames(eigengenes))
# dotplottab <- dotplottab[order(match(dotplottab[,1], colorOrder)),]
# 
# dotplottab[,1] <- factor(dotplottab[,1], levels = colorOrder)
# dotplottab <- dotplottab[order(dotplottab[,"module"], dotplottab[,"p.adjust"]),]
# dotplottab[,3] <- factor(dotplottab[,3], levels = unique(dotplottab[,3]))
# 
# 
# 
# labsparam <- list(title = "module dot plot", x = "modules", y = "GO terms", color = "module", size = "adjusted p value")
# colorvec <- as.character(unique(dotplottab[,1]))
# names(colorvec) <- colorvec
# 
# pout <- ggplot(dotplottab, aes(x=dotplottab[,1], y=dotplottab[,3], size=dotplottab[,2], color=as.character(dotplottab[,1])))
# # pout <- ggplot(modulesel, aes(x=modulesel[,1], y=modulesel[,3], size=modulesel[,2], color=modulesel[,2]))
# pout <- pout + geom_point(alpha = 0.8, shape = 16)
# pout <- pout + scale_size_continuous(name = "adjusted p value", range = c(8, 2), limits = (c(0,0.1)))
# pout <- pout + scale_color_manual(values = colorvec)
# pout <- pout + labs(title = labsparam$title, x = labsparam$x, y = labsparam$y, color = labsparam$color, shape = labsparam$shape, size = labsparam$size)
# pout <- pout + theme_bw()
# pout <- pout + theme(legend.position = "bottom", axis.text.x = element_text(angle = 90))
# 
# pdf(paste0(outfilepathwgcna, "module_enrichment_dotplot.pdf"), 8, 6)
# print(pout)
# junk <- dev.off()
# 
# 
# 
# ## Ordering gene modules by number of genes
# temporder <- names(table(moduleColors)[order(table(moduleColors))])
# 
# ## Subset the dotplot for only the modules that were sig dif
# subdotplottab <- na.omit(dotplottab[dotplottab[,"module"] %in% outstattab[outstattab[,"pval_PACE_v_THR"] < 0.05,"Color"],])
# subdotplottab[,1] <- factor(subdotplottab[,1], levels = colorOrder)
# subdotplottab <- subdotplottab[order(match(subdotplottab[,"module"], temporder), subdotplottab[,"p.adjust"]),]
# subdotplottab[,3] <- factor(subdotplottab[,3], levels = unique(subdotplottab[,3]))
# 
# labsparam <- list(title = "module dot plot", x = "modules", y = "GO terms", color = "module", size = "adjusted p value")
# colorvec <- as.character(unique(dotplottab[,1]))
# names(colorvec) <- colorvec
# 
# subpout <- ggplot(subdotplottab, aes(x=subdotplottab[,1], y=subdotplottab[,3], size=subdotplottab[,2], color=as.character(subdotplottab[,1])))
# # pout <- ggplot(modulesel, aes(x=modulesel[,1], y=modulesel[,3], size=modulesel[,2], color=modulesel[,2]))
# subpout <- subpout + geom_point(alpha = 0.8, shape = 16)
# subpout <- subpout + scale_size_continuous(name = "adjusted p value", range = c(8, 2), limits = (c(0,0.1)))
# subpout <- subpout + scale_color_manual(values = colorvec)
# subpout <- subpout + labs(title = labsparam$title, x = labsparam$x, y = labsparam$y, color = labsparam$color, shape = labsparam$shape, size = labsparam$size)
# subpout <- subpout + theme_bw()
# subpout <- subpout + theme(legend.position = "bottom", axis.text.x = element_text(angle = 90))
# 
# pdf(paste0(outfilepathwgcna, "subset_module_enrichment_dotplot.pdf"), 8, 6)
# print(subpout)
# junk <- dev.off()
# 
# 
# 
# ## The dotplots dont work really with that few genesets
# # So I need to just grab genesets and then plot with a barplot
# # fullannotationtable <- rbind(WGCNA_out2_H)
# fullannotationtable <- rbind(WGCNA_out2_C5)
# 
# brownPOI <- c("brown.GO_MYELOID_LEUKOCYTE_MEDIATED_IMMUNITY",
#               "brown.GO_SECRETION",
#               "brown.GO_SECRETORY_GRANULE",
#               "brown.GO_CELLULAR_RESPONSE_TO_OXYGEN_CONTAINING_COMPOUND",
#               "brown.GO_ACTIVATION_OF_INNATE_IMMUNE_RESPONSE",
#               "brown.GO_REGULATION_OF_CELL_DEATH",
#               "brown.GO_LYMPHOCYTE_ACTIVATION",
#               "brown.GO_REGULATION_OF_RESPONSE_TO_STRESS",
#               "brown.GO_T_CELL_ACTIVATION",
#               "brown.GO_COAGULATION")
# brownGSEAplottab <- cbind.data.frame(fullannotationtable[brownPOI,c("ID"),drop=FALSE], -log10(fullannotationtable[brownPOI,c("p.adjust"),drop=FALSE]))
# brownGSEAplottab[,1] <- gsub("_", " ", brownGSEAplottab[,1])
# rownames(brownGSEAplottab) <- fullannotationtable[brownPOI,c("ID")]
# 
# 
# brownpout <- ggplot(brownGSEAplottab, mapping = aes(x = str_wrap(brownGSEAplottab[,1], 20), y = brownGSEAplottab[,2]))
# brownpout <- brownpout + geom_bar(stat = "identity", fill = "brown")
# # brownpout <- brownpout + scale_fill_viridis(limits = c(0,max(0.1, max(gseaplotin[,3]))), direction = -1, option = "magma")
# brownpout <- brownpout + scale_x_discrete(limits=rev(str_wrap(brownGSEAplottab[,1], 20)))
# brownpout <- brownpout + coord_flip()
# brownpout <- brownpout + labs(x="Geneset", y = "NES", title = "Brown POI GSEA")
# pdf(paste0(outfilepathwgcna, "brown_gsea_barplot.pdf"))
# plot(brownpout)
# junk <- dev.off()
# 
# tanPOI <- c("tan.GO_B_CELL_ACTIVATION",
#             "tan.GO_LYMPHOCYTE_ACTIVATION",
#             "tan.GO_ANTIGEN_RECEPTOR_MEDIATED_SIGNALING_PATHWAY",
#             "tan.GO_CELL_SURFACE",
#             "tan.GO_REGULATION_OF_IMMUNE_RESPONSE",
#             "tan.GO_B_CELL_PROLIFERATION",
#             "tan.GO_HUMORAL_IMMUNE_RESPONSE",
#             "tan.GO_ACTIVATION_OF_IMMUNE_RESPONSE",
#             "tan.GO_IMMUNOGLOBULIN_COMPLEX",
#             "tan.GO_IMMUNOGLOBULIN_SECRETION")
# tanGSEAplottab <- cbind.data.frame(fullannotationtable[tanPOI,c("ID"),drop=FALSE], -log10(fullannotationtable[tanPOI,c("p.adjust"),drop=FALSE]))
# tanGSEAplottab[,1] <- gsub("_", " ", tanGSEAplottab[,1])
# rownames(tanGSEAplottab) <- fullannotationtable[tanPOI,c("ID")]
# 
# 
# 
# tanpout <- ggplot(tanGSEAplottab, mapping = aes(x = str_wrap(tanGSEAplottab[,1], 20), y = tanGSEAplottab[,2]))
# tanpout <- tanpout + geom_bar(stat = "identity", fill = "tan")
# # tanpout <- tanpout + scale_fill_viridis(limits = c(0,max(0.1, max(gseaplotin[,3]))), direction = -1, option = "magma")
# tanpout <- tanpout + scale_x_discrete(limits=rev(str_wrap(tanGSEAplottab[,1], 20)))
# tanpout <- tanpout + coord_flip()
# tanpout <- tanpout + labs(x="Geneset", y = "NES", title = "tan POI GSEA")
# pdf(paste0(outfilepathwgcna, "tan_gsea_barplot.pdf"))
# plot(tanpout)
# junk <- dev.off()
# 
# # fullannotationtable <- rbind(WGCNA_out2_H, WGCNA_out2_C5)
# # modulesel <- na.omit(fullannotationtable[fullannotationtable[,"p.adjust"] < 0.05,c("module", "p.adjust", "ID")])
# # 
# # # add on the other colors with blanks
# # if (sum(!unique(fullannotationtable[,"module"]) %in% modulesel[,"module"]) != 0){
# #     blankmoduletab <- data.frame(module = unique(fullannotationtable[,"module"])[!unique(fullannotationtable[,"module"]) %in% modulesel[,"module"]],
# #                                  p.adjust = NA, ID = modulesel[1,"ID"])
# # } else {blankmoduletab <- NULL}
# # 
# # dotplottab <- rbind(modulesel, blankmoduletab)
# # dotplottab[,3] <- gsub("HALLMARK_|GO_", "", dotplottab[,3])
# # 
# # ## Reorder to match the order of eigengenes
# # colorOrder <- gsub("ME", "", colnames(eigengenes))
# # dotplottab <- dotplottab[order(match(dotplottab[,1], colorOrder)),]
# # 
# # dotplottab[,1] <- factor(dotplottab[,1], levels = colorOrder)
# # dotplottab <- dotplottab[order(dotplottab[,"module"], dotplottab[,"p.adjust"]),]
# # dotplottab[,3] <- factor(dotplottab[,3], levels = unique(dotplottab[,3]))
# 
# ## Ordering gene modules by number of genes
# temporder <- names(table(moduleColors)[order(table(moduleColors))])
# 
# ## Subset the dotplot for only the modules that were sig dif
# subdotplottab <- na.omit(dotplottab[dotplottab[,"module"] %in% outstattab[outstattab[,"pval_PACE_v_THR"] < 0.05,"Color"],])
# subdotplottab[,1] <- factor(subdotplottab[,1], levels = colorOrder)
# subdotplottab <- subdotplottab[order(match(subdotplottab[,"module"], temporder), subdotplottab[,"p.adjust"]),]
# subdotplottab[,3] <- factor(subdotplottab[,3], levels = unique(subdotplottab[,3]))
# 
# 
# 
# 
# 
# 
# 
# # ## Can we collapse GO terms
# # # Step 1 - take GO terms of interest
# # POI <- WGCNA_out2_C5[WGCNA_out2_C5[,"module"] %in% "brown","ID"]
# # 
# # # Step 2 - convert goterms to goIDs
# # ingoterms <- POI
# # gotermconvertout <- goterm_to_goID(ingoterms)
# # goids_found <- gotermconvertout[[1]]
# # goids_closest <- gotermconvertout[[2]]
# # 
# # # Step 3 - take our goIDs, and then find ancestors
# # goid_conv_table <- rbind(goids_found[,c("GOID", "Input_GOterm")], goids_closest[,c("GOID", "Input_GOterm")])
# # ingoids <- goid_conv_table[,1]
# # ancestral_GOterm_out <- find_goterm_ancestors(ingoids)
# # ancestral_GOterm_counttable <- ancestral_GOterm_out[[1]]
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
