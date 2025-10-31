###########################################################################
#
#                       Integrative
#
###########################################################################
# Author: ISCHEMIA Study Team
# Date: Updated 2025-08-08
# Description: Integrative analysis for ISCHEMIA study
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
# outfilepathintegration = create_output_dir("analysis_output")
outfilepathintegration = create_output_dir("analysis_output")
dir.create(outfilepathintegration, recursive = TRUE, showWarnings = FALSE)

## Infiles
inmetafile <- file.path(DATA_DIR, "metatable_in.csv")
inmetatable <- read.table(inmetafile, sep = ",", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
# SOI <- rownames(na.omit(inmetatable[,"comp_ischemia__Sev_v_MildNone",drop=FALSE]))
SOI <- rownames(inmetatable)

## Grab the matched samples
incountfile <- file.path(DATA_DIR, "normcounttab.txt")
normcounttable <- read.table(incountfile, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

## So lets start simple, do the WGCNA info and nmf cluster info look interesting at all when put together????????
# wgcna_eigengenes_tablefile <- "# PATH_UPDATED: output/run2b_rmoutliers2_controlage/WGCNA_power14/wgcna_eigengenes.csv"
wgcna_eigengenes_tablefile <- "# PATH_UPDATED: output/run2b_rmoutliers2_controlage/WGCNA/WGCNA_power16_size40/wgcna_eigengenes.csv"
wgcna_eigengenes_table <- read.table(wgcna_eigengenes_tablefile, sep = ",", header = TRUE, row.names = 1)
eigengenes <- colnames(wgcna_eigengenes_table)

nmf_clustermembership_tablefile <- "# PATH_UPDATED: output/run1c_rmoutliers2/NMF/test_run2/nmf_clustermembership_table.csv"
# nmf_clustermembership_tablefile <- "# PATH_UPDATED: output/run2b_rmoutliers2_controlage/NMF_bp/try5_newmethylation_rmDMPsex1_rnamethonly/rank_4_nmf/nmf_clustermembership_table.csv"
nmf_clustermembership_table <- read.table(nmf_clustermembership_tablefile, sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
nmf_clustercores_tablefile <- "# PATH_UPDATED: output/run1c_rmoutliers2/NMF/test_run2/nmf_clustercores_table.csv"
# nmf_clustercores_tablefile <- "# PATH_UPDATED: output/run2b_rmoutliers2_controlage/NMF_bp/try5_newmethylation_rmDMPsex1_rnamethonly/rank_4_nmf/nmf_clustercores_table.csv"
nmf_clustercores_table <- read.table(nmf_clustercores_tablefile, sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

nmf_clustermembership_table[,"combined"] <- apply(nmf_clustermembership_table, 1, function(x) paste(na.omit(x)))
nmf_clustercores_table[,"combined"] <- nmf_clustermembership_table[,"combined"]
nmf_clustercores_table[rowSums(is.na(nmf_clustercores_table[,!grepl("combined", colnames(nmf_clustercores_table))])) == 4,"combined"] <- NA
nmfcluster_allsamples = nmf_clustermembership_table[,"combined",drop=FALSE]

## Full metafile read in - to grab some more info...
# Apparently there is no CBC? So lets use what we have at least that may be relevant:
inbioreptablefile <- "# PATH_UPDATED: data/biorep_10_27.csv"
inbioreptable <- read.table(inbioreptablefile, sep = ",", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("NA", "", NA))
extraCOI <- c("PATNUM",
              "HEMOGLOB", "PLATELET", "WBC",
              "IMGDEGIS",
              "CTNDV50"
)
addonmetatable <- inbioreptable[inbioreptable[,"PATNUM"] %in% colnames(normcounttable),colnames(inbioreptable) %in% extraCOI]
rownames(addonmetatable) <- addonmetatable[,"PATNUM"]
addonmetatable <- addonmetatable[,!grepl("PATNUM", colnames(addonmetatable))]

## Adding some custom categories of IMGDEGIS
addonmetatable[,"IMGDEGIS_01_2_3"] <- addonmetatable[,"IMGDEGIS"]
addonmetatable[addonmetatable[,"IMGDEGIS_01_2_3"] %in% c("None", "Mild"),"IMGDEGIS_01_2_3"] <- "NoneMild"
addonmetatable[,"IMGDEGIS_01_23"] <- addonmetatable[,"IMGDEGIS_01_2_3"]
addonmetatable[addonmetatable[,"IMGDEGIS_01_23"] %in% c("Moderate", "Severe"),"IMGDEGIS_01_23"] <- "ModerateSevere"

## What about having analyses SPECIFIC to each cluster, so what is happening to ischemia vs WGCNA modules IN Cluster1
addonmetatable[,"CLUSTER1_IMGDEGIS_01_2_3"] <- addonmetatable[,"IMGDEGIS_01_2_3"]
addonmetatable[!rownames(addonmetatable) %in% rownames(nmf_clustermembership_table[nmf_clustermembership_table[,"combined"] %in% "nmf_cluster_1",]),
               "CLUSTER1_IMGDEGIS_01_2_3"] <- NA
addonmetatable[,"CLUSTER2_IMGDEGIS_01_2_3"] <- addonmetatable[,"IMGDEGIS_01_2_3"]
addonmetatable[!rownames(addonmetatable) %in% rownames(nmf_clustermembership_table[nmf_clustermembership_table[,"combined"] %in% "nmf_cluster_2",]),
               "CLUSTER2_IMGDEGIS_01_2_3"] <- NA
addonmetatable[,"CLUSTER3_IMGDEGIS_01_2_3"] <- addonmetatable[,"IMGDEGIS_01_2_3"]
addonmetatable[!rownames(addonmetatable) %in% rownames(nmf_clustermembership_table[nmf_clustermembership_table[,"combined"] %in% "nmf_cluster_3",]),
               "CLUSTER3_IMGDEGIS_01_2_3"] <- NA
addonmetatable[,"CLUSTER4_IMGDEGIS_01_2_3"] <- addonmetatable[,"IMGDEGIS_01_2_3"]
addonmetatable[!rownames(addonmetatable) %in% rownames(nmf_clustermembership_table[nmf_clustermembership_table[,"combined"] %in% "nmf_cluster_4",]),
               "CLUSTER4_IMGDEGIS_01_2_3"] <- NA

## I think I also want to add on the ISCHEMIA severity data here and do a similar analysis
# And I think I want to treat it like the CBC data - as continuous, but then when we do the "raw" section", make sure it doesnt get split
# isch_addonmetatable <- inmetatable[,"IMGDEGIS",drop=FALSE]
# isch_addonmetatable <- within( inmetatable[,"IMGDEGIS",drop=FALSE], IMGDEGIS <- factor( IMGDEGIS, 
#                                                                                         levels = c("None", "Mild", "Moderate", "Severe"), 
#                                                                                         labels = c("0", ) ))


## For each eigengenecluster - we want to see if the clusters have different values
dir.create(paste0(outfilepathintegration, "wgcna_nmf_integration/"), showWarnings = FALSE, recursive = TRUE)

## The cleanest way to run this is over a list of cohort TABLES, and then work with each of those
COItablist <- list(
    nmfcluster_allsamples = nmf_clustermembership_table[,"combined",drop=FALSE],
    nmfcluster_coresamples = nmf_clustercores_table[,"combined",drop=FALSE],
    hemoglobin = addonmetatable[,"HEMOGLOB", drop=FALSE],
    platelet = addonmetatable[,"PLATELET", drop=FALSE],
    wbc = addonmetatable[,"WBC", drop=FALSE],
    IMGDEGIS = addonmetatable[,"IMGDEGIS", drop=FALSE],
    IMGDEGIS_01_2_3 = addonmetatable[,"IMGDEGIS_01_2_3", drop=FALSE],
    IMGDEGIS_01_23 = addonmetatable[,"IMGDEGIS_01_23", drop=FALSE],
    CLUSTER1_IMGDEGIS_01_2_3 = addonmetatable[,"CLUSTER1_IMGDEGIS_01_2_3",drop=FALSE],
    CLUSTER2_IMGDEGIS_01_2_3 = addonmetatable[,"CLUSTER2_IMGDEGIS_01_2_3",drop=FALSE],
    CLUSTER3_IMGDEGIS_01_2_3 = addonmetatable[,"CLUSTER3_IMGDEGIS_01_2_3",drop=FALSE],
    CLUSTER4_IMGDEGIS_01_2_3 = addonmetatable[,"CLUSTER4_IMGDEGIS_01_2_3",drop=FALSE],
    CTNDV50 = addonmetatable[,"CTNDV50",drop=FALSE]
)

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
        wgcna_sel <- wgcna_eigengenes_table[,eigen_sel,drop=FALSE]
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
spearmanpval_plottab <- spearmanpval_plottab[,!grepl("Eigengene", colnames(spearmanpval_plottab))]
keepeigennames <- rownames(spearmanpval_plottab)
spearmanpval_plottab <- t(apply(spearmanpval_plottab, 2, function(x) as.numeric(as.character(x))))
colnames(spearmanpval_plottab) <- keepeigennames

spearmanrval_plottab <- dcast(temptab[,c("Cohort", "Eigengene", "spearman_rval")], Eigengene ~ Cohort, value.var = "spearman_rval")
rownames(spearmanrval_plottab) <- spearmanrval_plottab[,"Eigengene"]
spearmanrval_plottab <- spearmanrval_plottab[,!grepl("Eigengene", colnames(spearmanrval_plottab))]
spearmanrval_plottab <- t(apply(spearmanrval_plottab, 2, function(x) as.numeric(as.character(x))))
colnames(spearmanrval_plottab) <- keepeigennames

spearmanhmout <- WGCNA_custom_heatmap(moduleTraitPvalue = spearmanpval_plottab, moduleTraitCor = spearmanrval_plottab)
pdf(paste0(outfilepathintegration, "wgcna_nmf_integration/", "spearman_values.pdf"), 10, 20, useDingbats = FALSE)
print(spearmanhmout)
junk <- dev.off()
















### Maybe do a survival analysis for EVERY EVENT with our clusters??? Maybe... Worth a shot?
# outfilepathsurvival = paste0(outfilepathintegration, "nmf_survival_analysis/run1_nmfcluster_IMGDEGIS/")
# outfilepathsurvival = paste0(outfilepathintegration, "nmf_survival_analysis/run2_nmfcluster_IMGDEGIS_nmfclustervrest/")
outfilepathsurvival = paste0(outfilepathintegration, "nmf_survival_analysis/run3_nmfcluster_IMGDEGIS_nmfclustervrest_wlabs/")
dir.create(outfilepathsurvival, recursive = TRUE, showWarnings = FALSE)

EOI <- gsub("^C_", "", colnames(inbioreptable)[grepl("^C_", colnames(inbioreptable))])
EOIcols <- apply(expand.grid(c("T_", "C_"), EOI, stringsAsFactors = FALSE), 1, paste, collapse = "")

## I am also going to compare each cluster against the rest of samples as well - that may also be interesting (and make the HR more understandable)
expand_nmfcomps <- make_comparison_columns(nmf_clustermembership_table[,"combined",drop=FALSE])

## Also want to output the labs data to see if its predictive for anything (hopefully not....)
addonmetatable_quartiled <- addonmetatable[,c("HEMOGLOB", "PLATELET", "WBC")]
for (metavarnum in seq_len(ncol(addonmetatable_quartiled))) {
    colsel <- addonmetatable_quartiled[,metavarnum,drop=FALSE]
    collable <- colnames(colsel)
    cattab <- cut2(t(colsel[,1]),g=4)
    levels(cattab)[match(levels(cattab)[1],levels(cattab))] <- paste0(collable, "_Q1min")
    levels(cattab)[match(levels(cattab)[2],levels(cattab))] <- paste0(collable, "_Q2")
    levels(cattab)[match(levels(cattab)[3],levels(cattab))] <- paste0(collable, "_Q3")
    levels(cattab)[match(levels(cattab)[4],levels(cattab))] <- paste0(collable, "_Q4max")
    addonmetatable_quartiled[,metavarnum] <- as.character(cattab)
}
expand_addonmetatable_quartiled <- make_comparison_columns(na.omit(addonmetatable_quartiled))




COItablist <- list(
    # nmfcluster_allsamples = nmf_clustermembership_table[,"combined",drop=FALSE],
    # nmfcluster_coresamples = nmf_clustercores_table[,"combined",drop=FALSE],
    
    # IMGDEGIS = addonmetatable[,"IMGDEGIS", drop=FALSE],
    # IMGDEGIS_01_2_3 = addonmetatable[,"IMGDEGIS_01_2_3", drop=FALSE],
    IMGDEGIS_01_23 = addonmetatable[,"IMGDEGIS_01_23", drop=FALSE],

    nmfcluster1vrest = data.frame(expand_nmfcomps[,"combined_nmf_cluster_1"]),
    nmfcluster2vrest = data.frame(expand_nmfcomps[,"combined_nmf_cluster_2"]),
    nmfcluster3vrest = data.frame(expand_nmfcomps[,"combined_nmf_cluster_3"]),
    nmfcluster4vrest = data.frame(expand_nmfcomps[,"combined_nmf_cluster_4"]),
    
    # hemoglobin = addonmetatable[,"HEMOGLOB", drop=FALSE],
    # platelet = addonmetatable[,"PLATELET", drop=FALSE],
    # wbc = addonmetatable[,"WBC", drop=FALSE],
    # CLUSTER1_IMGDEGIS_01_2_3 = addonmetatable[,"CLUSTER1_IMGDEGIS_01_2_3",drop=FALSE],
    # CLUSTER2_IMGDEGIS_01_2_3 = addonmetatable[,"CLUSTER2_IMGDEGIS_01_2_3",drop=FALSE],
    # CLUSTER3_IMGDEGIS_01_2_3 = addonmetatable[,"CLUSTER3_IMGDEGIS_01_2_3",drop=FALSE],
    # CLUSTER4_IMGDEGIS_01_2_3 = addonmetatable[,"CLUSTER4_IMGDEGIS_01_2_3",drop=FALSE],
    # CTNDV50 = addonmetatable[,"CTNDV50",drop=FALSE]
    
    HEMOGLOBQ1vrest = data.frame(expand_addonmetatable_quartiled[,"HEMOGLOB_HEMOGLOB_Q1min"]),
    HEMOGLOBQ2vrest = data.frame(expand_addonmetatable_quartiled[,"HEMOGLOB_HEMOGLOB_Q2"]),
    HEMOGLOBQ3vrest = data.frame(expand_addonmetatable_quartiled[,"HEMOGLOB_HEMOGLOB_Q3"]),
    HEMOGLOBQ4vrest = data.frame(expand_addonmetatable_quartiled[,"HEMOGLOB_HEMOGLOB_Q4max"]),
    PLATELETQ1vrest = data.frame(expand_addonmetatable_quartiled[,"PLATELET_PLATELET_Q1min"]),
    PLATELETQ2vrest = data.frame(expand_addonmetatable_quartiled[,"PLATELET_PLATELET_Q2"]),
    PLATELETQ3vrest = data.frame(expand_addonmetatable_quartiled[,"PLATELET_PLATELET_Q3"]),
    PLATELETQ4vrest = data.frame(expand_addonmetatable_quartiled[,"PLATELET_PLATELET_Q4max"]),
    WBCQ1vrest = data.frame(expand_addonmetatable_quartiled[,"WBC_WBC_Q1min"]),
    WBCQ2vrest = data.frame(expand_addonmetatable_quartiled[,"WBC_WBC_Q2"]),
    WBCQ3vrest = data.frame(expand_addonmetatable_quartiled[,"WBC_WBC_Q3"]),
    WBCQ4vrest = data.frame(expand_addonmetatable_quartiled[,"WBC_WBC_Q4max"])
)

## Ok - so we want to iterate over every GOI, and every event, and see what comes out
# survival_intable <- merge(inmetatable[grepl("PACE", inmetatable[,"PATNUM"]),c("PATNUM", EOIcols)], 
#                           GOIcat_counttable, by = "PATNUM")


fullpvaloutlist <- fullcumhazoutlist <- list()
for (COInum in seq_len(length(COItablist))) {
    ## Select gene
    COIsel <- COItablist[[COInum]]
    COIlabel <- names(COItablist)[COInum]
    
    survival_intable <- na.omit(merge(inbioreptable[inbioreptable[,"PATNUM"] %in% colnames(normcounttable),c("PATNUM", EOIcols)],
                              cbind(PATNUM = rownames(COIsel), COIsel), by = "PATNUM"))
    colnames(survival_intable)[ncol(survival_intable)] <- COIlabel
    
    survpvallist <- survHRlist <- list()
    for (eventnum in seq_len(length(EOI))) {
    # for (eventnum in seq_len(40)) {
        ## Select event
        eventsel <- EOI[eventnum]
        
        ## Create the subsurvivaltab
        survivaldata <- data.frame(na.omit(survival_intable[,c(paste0(c("T_", "C_"), eventsel), COIlabel)]))
        rownames(survivaldata) <- survival_intable[,"PATNUM"]
        
        ## Ok, for the competing events - we have to clean those out cause it screws with the analysis - just remove I guess?
        survivaldata <- survivaldata[!survivaldata[,grepl("^C_", colnames(survivaldata))] %in% 2,]
        
        # ## Apply a max time of 5 years - it cleans up the data a little - NOPE, dont need it
        # maxtimelimit <- 5*365
        # survivaldata <- survivaldata[survivaldata[,grepl("^T_", colnames(survivaldata))] < maxtimelimit,]
        
        ## Need to properly factorize our cohort data
        if (COIlabel %in% c("nmfcluster_allsamples", "nmfcluster_coresamples")) {survivaldata[,3] <- factor(survivaldata[,3], levels = c("nmf_cluster_1", "nmf_cluster_2", "nmf_cluster_3", "nmf_cluster_4"))}
        if (COIlabel %in% c("IMGDEGIS")) {survivaldata[,3] <- factor(survivaldata[,3], levels = c("None", "Mild", "Moderate", "Severe"))}
        if (COIlabel %in% c("IMGDEGIS_01_2_3")) {survivaldata[,3] <- factor(survivaldata[,3], levels = c("NoneMild", "Moderate", "Severe"))}
        if (COIlabel %in% c("IMGDEGIS_01_23")) {survivaldata[,3] <- factor(survivaldata[,3], levels = c("NoneMild", "ModerateSevere"))}
        if (COIlabel %in% c("nmfcluster1vrest", "nmfcluster2vrest", "nmfcluster3vrest", "nmfcluster4vrest")) {
            cohortlabels <- unique(survivaldata[,3])
            survivaldata[,3] <- factor(survivaldata[,3], levels = c(as.character(cohortlabels[grepl("not_", cohortlabels)]), 
                                                                    as.character(cohortlabels[!grepl("not_", cohortlabels)])))
        }
        # if (COIlabel %in% c("HEMOGLOBQ1vrest", "HEMOGLOBQ2vrest", "HEMOGLOBQ3vrest", "HEMOGLOBQ4vrest",
        #                     "PLATELETQ1vrest", "PLATELETQ2vrest", "PLATELETQ3vrest", "PLATELETQ4vrest",
        #                      WBCQ1vrest", "WBCQ2vrest", "WBCQ3vrest", "WBCQ4vrest")) {
        #     cohortlabels <- unique(survivaldata[,3])
        #     survivaldata[,3] <- factor(survivaldata[,3], levels = c(as.character(cohortlabels[grepl("not_", cohortlabels)]), 
        #                                                             as.character(cohortlabels[!grepl("not_", cohortlabels)])))
        # }
        
        ## Run the analysis
        survivalanalysisout <- create_survival_plot(survivaldata = survivaldata, timebreakparam = NULL, ylimitparam = c(0.25,1))
        # survivalanalysisout <- create_survival_plot(survivaldata = survivaldata, timebreakparam = NULL)
        outsurvtable <- survivalanalysisout$outsurvtable
        outsurvplot <- survivalanalysisout$outsurvplot
        outsurvpvalue <- survivalanalysisout$outsurvpvalue
        
        ## Save out the pvalue of the log-rank test
        survpvallist[[eventnum]] <- outsurvpvalue[,"pval"]
        names(survpvallist)[eventnum] <- eventsel
        
        # Ok, I think we need to return the HR as well for coxph to get some kind of effect size for this analysis
        outcoxphobject <- survivalanalysisout$outcoxphobject
        ## Put a hack in here to get the orientation correct. If theres only one comp, then coerce to a 1x3 dataframe:
        coxtempout <- summary(outcoxphobject)[["conf.int"]][,c("exp(coef)", "lower .95", "upper .95")]
        if(is.null(dim(coxtempout))){
            coxtempout <- t(data.frame(coxtempout))
            rownames(coxtempout) <- levels(survivaldata[,3])[2] ## I can do this because I know the first level is always the ref and I set this earlier above ^^
        } 
        hrvalue <- cbind(coxtempout, Event = eventsel, Cohort = COIlabel)
        
        ## Save out the coxph HR of the coxph test
        survHRlist[[eventnum]] <- hrvalue
        names(survHRlist)[eventnum] <- eventsel
        
        ## Write out the plot and table
        outsubdir <- paste0(outfilepathsurvival, COIlabel, "/", eventsel, "/")
        dir.create(outsubdir, showWarnings = FALSE, recursive = TRUE)
        
        write.table(outsurvtable, paste0(outsubdir, COIlabel, "_", eventsel, "_survival_stat_table.csv"),
                    sep = ",", col.names = TRUE, row.names = FALSE)
        pdf(paste0(outsubdir, COIlabel, "_", eventsel, "_survival_plot.pdf"), width = 8, height = 5, onefile=FALSE)
        print(outsurvplot)
        junk <- dev.off()
        
    }
    survpvaltab <- do.call(rbind, survpvallist)
    colnames(survpvaltab) <- paste0(COIlabel, "_pval")
    survhazardtab <- do.call(rbind, survHRlist)
    colnames(survhazardtab) <- paste0(COIlabel, "_", colnames(survhazardtab))
    
    fullpvaloutlist[[COInum]] <- survpvaltab
    names(fullpvaloutlist[COInum]) <- paste0(COIlabel)
    fullcumhazoutlist[[COInum]] <- survhazardtab
    names(fullcumhazoutlist)[COInum] <- paste0(COIlabel)
    
}
fullpvalouttable <- do.call(cbind, fullpvaloutlist)
write.table(fullpvalouttable, paste0(outfilepathsurvival, "full_survival_pval_table.csv"), sep = ",", col.names = NA, row.names = TRUE)
fullcumhazouttable <- do.call(rbind, fullcumhazoutlist)
colnames(fullcumhazouttable) <- c("HR", "lower_CI", "upper_CI", "Event", "Cohort")
write.table(fullcumhazouttable, paste0(outfilepathsurvival, "full_survival_cumhaz_table.csv"), sep = ",", col.names = NA, row.names = TRUE)


## Summary heatmap of the pvals
# fullstatouttablefile <- "# PATH_UPDATED: output/run1c_rmoutliers2/integrative_analyses/nmf_survival_analysis/full_survival_pval_table.csv"
# fullstatouttable <- read.table(fullstatouttablefile, sep = ",", header = TRUE, row.names = 1)
fullstatouttable <- fullpvalouttable

## Turn the event outcome into a heatmap
eventpvalcutoff <- 0.1
# eventpvalcutoff <- 0.99

# GOI heatmap with events
## Lets create a heatmap with the GOI, and then use a side annotation for if its indicative of events
# eventpvaltabfile <- "/Users/tosh/Desktop/Ruggles_Lab/projects/newman-pace-2020/output_hg19_newcontrols/run1_nooutlier1/survival_analysis/topbotquartile/full_survival_stat_table.csv"
# outeventpvaltab <- read.table(eventpvaltabfile, sep = ",", header = TRUE, row.names = 1)

## TOSH - selecting just nmf for now, the cell stuff is important, but not worth it atm
outeventpvaltab <- fullstatouttable[,c("IMGDEGIS_01_23_pval", "nmfcluster1vrest_pval", "nmfcluster2vrest_pval", "nmfcluster3vrest_pval", "nmfcluster4vrest_pval")]

## Create the initial maptab
# maptab <- data.frame(outeventpvaltab[rowSums(outeventpvaltab[,EOIvec] < eventpvalcutoff) > 0, EOIvec])
maptab <- data.frame(outeventpvaltab)

## MULTIPLE HYPOTHESIS CORRECTION - OMITTING FOR NOW BECAUSE WE HAVE WAY TOO MANY BAD HYPOTHESES AND ITS CLOGGING IT UP
# p.adjust(unlist(maptab), method = "fdr") ## Sanity check - it works
# “holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”

# correctedmat <- matrix(p.adjust(as.vector(as.matrix(maptab)), method='fdr'),ncol=ncol(maptab))
# rownames(correctedmat) <- rownames(maptab)
# colnames(correctedmat) <- colnames(maptab)
# maptab <- correctedmat

# Create a heatmap of just the event pvalues
EOIvec <- colnames(outeventpvaltab)
## Add in geneordering
# geneorder <- indeseqtable[rownames(outeventpvaltab)[rowSums(outeventpvaltab[,EOIvec] < eventpvalcutoff) > 0],]
# geneorder <- rownames(geneorder[order(geneorder[,"padj"], decreasing = FALSE),]) ## order by padj
# geneorder <- rownames(geneorder[order(geneorder[,"log2FoldChange"], decreasing = TRUE),]) ## order by log2fc
# geneorder <- rownames(correctedmat[order(indeseqtable[rownames(correctedmat),"log2FoldChange"], decreasing = TRUE),]) ## order by log2fc
## Order by avg pval across all events
# geneorder <- names(rowMeans(correctedmat[rownames(geneorder),])[order(rowMeans(correctedmat[rownames(geneorder),]), decreasing = FALSE)])
# geneorder <- names(rowMeans(correctedmat[order(rowMeans(correctedmat), decreasing = FALSE),]))
geneorder <- rownames(outeventpvaltab[order(outeventpvaltab[,1], decreasing = FALSE),])

# maptab <- data.frame(outeventpvaltab[rowSums(outeventpvaltab[,EOIvec] < eventpvalcutoff) > 0, EOIvec])


# maptab <- maptab[geneorder,][rowSums(maptab[geneorder,] < eventpvalcutoff) > 0,]

## Add in coloring here based on directionality:
hazardmaptabtemp <- dcast(data.frame(fullcumhazouttable), Event ~ Cohort, value.var = "HR")
rownames(hazardmaptabtemp) <- hazardmaptabtemp[,1]
# hazardmaptab <- hazardmaptabtemp[rownames(maptab),!grepl("Event", colnames(hazardmaptabtemp))]
hazardmaptabtemp <- hazardmaptabtemp[rownames(maptab),gsub("_pval", "", colnames(maptab))]
# hazardmaptab <- hazardmaptabtemp[rownames(maptab),gsub("_pval", "", colnames(maptab))]
hazardmaptab <- apply(hazardmaptabtemp, 2, as.numeric)
rownames(hazardmaptab) <- hazardmaptabtemp[,1]
hazardsigntab <- ifelse(hazardmaptab > 1, 1, -1)

maptab <- maptab * hazardsigntab

heatmapcolorparam = colorRamp2(c(-1, -eventpvalcutoff - 0.0001, -eventpvalcutoff, -0.0001, 0, 0.0001, eventpvalcutoff, eventpvalcutoff + 0.0001, 1), 
                               c("grey", "grey", "#bcb2fd", "#11007d", "white", "#7d0000", "#fd9999", "grey", "grey"))
heatmapcolorLEGEND = colorRamp2(c(-eventpvalcutoff, 0, eventpvalcutoff), c("blue", "white", "red"))



# heatmapcolorparam = colorRamp2(c(1, eventpvalcutoff + 0.0001, eventpvalcutoff, 0), c("white", "white", "#c9e0dc", "#00705E"))
# heatmapcolorLEGEND = colorRamp2(c(eventpvalcutoff, 0), c("#c9e0dc", "#00705E"))

## LEAVE IN ALL EVENTS:
# heatmapcolorparam = colorRamp2(c(1, 0), c("white", "#00705E"))
# heatmapcolorLEGEND = colorRamp2(c(1, 0), c("#white", "#00705E"))

# heatmapcolorparam = colorRamp2(c(1, 0.1 + 0.0001, 0.1, 0), c("white", "white", "#c9e0dc", "#00705E"))
# heatmapcolorLEGEND = colorRamp2(c(0.1, 0), c("#c9e0dc", "#00705E"))

ht1 = Heatmap(as.matrix(maptab), 
              col = heatmapcolorparam,    ## Define the color scale for the heatmap
              row_title = "Genes",                                       ## Name the rows
              column_title = "Events",                                  ## Name the columns
              border = TRUE,
              na_col = "white",
              rect_gp = gpar(col = "black", lwd = 0.5),
              
              cluster_columns = FALSE,                         ## Cluster the columns or leave as is
              cluster_rows = FALSE,                            ## Cluster the rows or leave as is
              clustering_method_rows = "ward.D2",
              clustering_method_columns = "ward.D2",
              #"ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", centroid"
              
              show_column_names = TRUE,                                  ## Show the Column Names
              column_names_gp = gpar(fontsize = 6),                      ## Change the size of the column names
              show_row_names = TRUE,                                    ## Show the row names
              row_names_side = "left",                                   ## Place the row names on the Left side of the heatmap
              row_names_gp = gpar(fontsize=6),
              
              show_row_dend = TRUE, #nrow(maptab) <=500,                                     ## Show the dendrogram on the rows
              show_column_dend = TRUE,                                   ## Show the dendrogram on the columns
              
              heatmap_legend_param = list(
                  col_fun = heatmapcolorLEGEND,
                  # at = seq(0, eventpvalcutoff, 0.01),
                  # at = seq(-eventpvalcutoff, eventpvalcutoff, 0.01),
                  # title = ifelse(samplesample==FALSE, "Zscore", "Spearman\nCorrelation"),
                  title = "pvalue",
                  legend_height = unit(2.5, "cm"),
                  title_gp = gpar(fontsize = 8, fontface = "bold")),
              height = unit(min((nrow(maptab)/2), 12),"cm"),
              width = unit(min(ncol(maptab), 18),"cm")
              
)
draw(ht1)


# hmoutfile <- paste0(outfilepathsurvival, "event_heatmap_005.pdf")
hmoutfile <- paste0(outfilepathsurvival, "event_heatmap_010.pdf")
# hmoutfile <- paste0(outfilepathsurvival, "event_heatmap_099.pdf")
# hmoutfile <- paste0(outfilepathsurvival, "event_heatmap_all_005.pdf")
# hmoutfile <- paste0(outfilepathsurvival, "event_heatmap_all_005_wcolor.pdf")
pdf(hmoutfile, 11, 9)
draw(ht1)
junk <- dev.off()


## Heatmap of eigengene values - with annotation of nmf cluster (and maybe annotation for nmf cluster driver for modules??)
# wgcna_eigengenes_table ## table we need

# nmfcluster_allsamples ## sample annotation
nmfeigenranktable_file <- "# PATH_UPDATED: output/run1c_rmoutliers2/integrative_analyses/wgcna_nmf_integration/boxplots_nmfcluster_allsamples/nmfcluster_allsamples_eigenrank_summary_table.csv"
eigenrank_summarytable <- read.table(nmfeigenranktable_file, sep = ",", header = TRUE, row.names = 1)

## If I want to only do drivers:
nmfeigenranktable <-dcast(eigenrank_summarytable, eigengene ~ Group, value.var = "Eigen_Score_Rank")
rownames(nmfeigenranktable) <- nmfeigenranktable[,"eigengene"]
nmfeigenranktable <- nmfeigenranktable[do.call(order, nmfeigenranktable[,2:ncol(nmfeigenranktable)]),2:ncol(nmfeigenranktable)]
nmfeigenranktable_drivers <- nmfeigenranktable == 4
w <- which(nmfeigenranktable_drivers == TRUE,arr.ind=TRUE)
nmfeigenranktable_drivers[w] <- colnames(nmfeigenranktable_drivers)[w[,"col"]]
nmfeigenranktable_drivers[nmfeigenranktable_drivers == "FALSE"] <- NA
nmfeigenranktable_drivers_collapse<- data.frame(driver = apply(nmfeigenranktable_drivers, 1, function(x) paste(na.omit(x))))

# Adding in the DESeq information to sort the genes by the log2fc
# plottab = t(wgcna_eigengenes_table)
plottab = t(wgcna_eigengenes_table[,colnames(wgcna_eigengenes_table)[!colnames(wgcna_eigengenes_table) %in% "MEgrey"]]) ## Need to remove the grey eigengene - not important
colmeta <- nmfcluster_allsamples[colnames(plottab),,drop=FALSE]
rowmeta <- nmfeigenranktable[rownames(plottab),]

## Custom sorting by cols
clusterOrder = c("nmf_cluster_1", "nmf_cluster_2","nmf_cluster_3", "nmf_cluster_4")
colmeta <- colmeta[ order(match(colmeta[,1], clusterOrder)),,drop=FALSE]
plottab <- plottab[,rownames(colmeta)]


rowannotationlist <- annotationlist_builder(rowmeta)
colannotationlist <- annotationlist_builder(colmeta)
# colannotationlist1 <- annotationlist_builder(rowmeta, customcolorlist = list(
#     log2FoldChange = colorRamp2(c(ceiling(max(rowmeta[,1])), 0, floor(min(rowmeta[,1]))), brewer.pal(3,"RdBu")),
#     pvalue = colorRamp2(seq(0,1,length = 5), brewer.pal(5,"Reds"))
# ))
# heatmapcolorparam <- colorRamp2(
#     c(lowvalue, midvalue, highvalue), c(lowcolor, midcolor, highcolor))
heatmapcolorparam <- colorRamp2(
    c(-5, 0, 5), c("blue", "white", "red"))
outhm1 <- create_heatmap(counttab = plottab, scale_data = TRUE,
                         colmetatable = colmeta, colannotationlist = colannotationlist,
                         rowmetatable = rowmeta, rowannotationlist = rowannotationlist,
                         colclusterparam = TRUE, rowclusterparam = TRUE, heatmapcolorparam = heatmapcolorparam)
outheatmapfile1 = paste0(outfilepathintegration, "wgcna_nmf_integration/", "eigengene_nmfcluster_heatmap.pdf")
pdf(outheatmapfile1, useDingbats = FALSE, height = 7, width = 13)
draw(outhm1$heatmap)
junk <- dev.off()


## What if we make this simpler - and just heatmap the avg eigengene values
nmfeigenvaluetable <-dcast(eigenrank_summarytable, eigengene ~ Group, value.var = "Eigengene")
rownames(nmfeigenvaluetable) <- nmfeigenvaluetable[,"eigengene"]
nmfeigenvaluetable <- nmfeigenvaluetable[do.call(order, nmfeigenvaluetable[,2:ncol(nmfeigenvaluetable)]),2:ncol(nmfeigenvaluetable)]
nmfeigenvaluetable <- nmfeigenvaluetable[!rownames(nmfeigenvaluetable) %in% "MEgrey",] ## MEgrey removal

rowmeta2 <- nmfeigenranktable_drivers_collapse[rownames(nmfeigenvaluetable),,drop=FALSE]
rowannotationlist2 <- annotationlist_builder(rowmeta2)
outhm2 <- create_heatmap(counttab = nmfeigenvaluetable, scale_data = TRUE,
                         colmetatable = NULL, colannotationlist = NULL,
                         rowmetatable = rowmeta2, rowannotationlist = rowannotationlist2,
                         colclusterparam = TRUE, rowclusterparam = TRUE, heatmapcolorparam = NULL)
outheatmapfile2 = paste0(outfilepathintegration, "wgcna_nmf_integration/", "eigengeneavg_nmfcluster_heatmap.pdf")
pdf(outheatmapfile2, useDingbats = FALSE, height = 7, width = 7)
draw(outhm2$heatmap)
junk <- dev.off()



## For each NMF Cluster - I want to run over our metadata and do an in vs out comparison for enrichment of whatever variable we are looking at
outfilepathtabulation <- paste0(outfilepathintegration, "nmf_cluster_tabulations/")
dir.create(outfilepathtabulation, showWarnings = FALSE, recursive = TRUE)

# Created a key, use the coltype for the character classes when reading in the biorep table
inbioreptable_keyfile <- "# PATH_UPDATED: data/biorep_10_27_question_key.txt"
inbioreptable_key <- read.table(inbioreptable_keyfile, sep = "\t", header = TRUE, stringsAsFactors = FALSE, na.strings = c("NA", "", NA))

inbioreptablefile <- "# PATH_UPDATED: data/biorep_10_27.csv"
inbioreptable <- read.table(inbioreptablefile, sep = ",", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE,
                            na.strings = c("NA", "", NA),
                            colClasses = inbioreptable_key[,"coltype"]
                            )

nmfcluster_allsamples = nmf_clustermembership_table[,"combined",drop=FALSE]

nmfcluster_metadata_analysis_intab <- merge(cbind(PATNUM = rownames(nmfcluster_allsamples), nmfclusterlabel = nmfcluster_allsamples[,1]),
                                            inbioreptable[rownames(nmfcluster_allsamples) %in% inbioreptable[,"PATNUM"],],
                                            by = "PATNUM")
# tabulateCOI <- colnames(nmfcluster_metadata_analysis_intab)
## List out all of the variables we actually want to test.
tabulateCOI <- c("PATNUM", "nmfclusterlabel", inbioreptable_key[!is.na(inbioreptable_key[,"summarystats"]),"variable"])
nmfcluster_metadata_analysis_intab <- nmfcluster_metadata_analysis_intab[,tabulateCOI]

## Create summarized stats for all clusters, not A vs not A, but chi sq across full pop
summarize_intable <- nmfcluster_metadata_analysis_intab
# summarize_intable[,"nmfclusterlabel"] <- ifelse(summarize_intable[,"nmfclusterlabel"] == clusterlabel_sel, clusterlabel_sel, paste0("not_", clusterlabel_sel))

summarystatfile <- paste0(outfilepathtabulation, "summarystats_ischmemia_by_nmfcluster.csv")
sumstattab_seq <- summarize_table(intable = summarize_intable, groupvar = "nmfclusterlabel", outfile = summarystatfile, calc_stats = TRUE)


simplified_outlist <- simplified_subcluster_outlist <- list()
for (clusternum in seq_len(length(unique(nmfcluster_allsamples[,1])))) {
    ## Grab cluster
    clusterlabel_sel <- unique(nmfcluster_allsamples[,1])[clusternum]
    
    ## create summarize intable
    summarize_intable <- nmfcluster_metadata_analysis_intab
    summarize_intable[,"nmfclusterlabel"] <- ifelse(summarize_intable[,"nmfclusterlabel"] == clusterlabel_sel, clusterlabel_sel, paste0("not_", clusterlabel_sel))
    
    summarystatfile <- paste0(outfilepathtabulation, "summarystats_ischmemia_by_", clusterlabel_sel, ".csv")
    sumstattab_seq <- summarize_table(intable = summarize_intable, groupvar = "nmfclusterlabel", outfile = summarystatfile, calc_stats = TRUE)
    
    ## I need to simplify my job a little by outputting sig pvals or something:
    # so lets run over the whole list, and grab all of the labels and p values, and then sort by pvalue?
    outtablist <- list()
    for (statnum in seq_len(length(names(sumstattab_seq[!names(sumstattab_seq) %in% "PATNUM"])))) {
    # for (statnum in seq_len(5)) {
        statlabel_sel <- names(sumstattab_seq[!names(sumstattab_seq) %in% "PATNUM"])[statnum]
        stattab_sel <- sumstattab_seq[[statlabel_sel]]
        
        # If numerical, then just grab the pval
        if (as.character(stattab_sel[1,1]) == "Min." & as.character(stattab_sel[2,1]) == "1st Qu.") {
            outtab <- cbind(statlabel_sel, nmfclusterlabel_cat = "numeric", stattab_sel[1,"statval",drop=FALSE])
        } else {
            outtab <- cbind(statlabel_sel, stattab_sel[,1,drop=FALSE], stattab_sel[,"statval",drop=FALSE])
        }
        
        # If categorical, than grab the pval and the cat label and save each of those
        outtablist[[statnum]] <- outtab
    }
    statssummarytable <- do.call(rbind, outtablist)
    write.table(statssummarytable, paste0(outfilepathtabulation, "simplified_stats_ischmemia_by_", clusterlabel_sel, ".csv"),
                sep = ",", col.names = TRUE, row.names = FALSE)
    simplified_outlist[[clusternum]] <- cbind(clusterlabel = clusterlabel_sel, statssummarytable)
    
    
    
    
    ## Then I also want to run intervention analysis for each cluster
    subcluster_summarize_intable <- summarize_intable[!grepl("not_", summarize_intable[,"nmfclusterlabel"]),]
    subcluster_summarystatfile <- paste0(outfilepathtabulation, "summarystats_", clusterlabel_sel, "_by_intervention", ".csv")
    subcluster_sumstattab_seq <- summarize_table(intable = subcluster_summarize_intable, 
                                                 groupvar = "TREATMNT", outfile = subcluster_summarystatfile, calc_stats = TRUE)
    
    ## I need to simplify my job a little by outputting sig pvals or something:
    # so lets run over the whole list, and grab all of the labels and p values, and then sort by pvalue?
    subcluster_outtablist <- list()
    for (statnum in seq_len(length(names(subcluster_sumstattab_seq[!names(subcluster_sumstattab_seq) %in% "PATNUM"])))) {
        # for (statnum in seq_len(5)) {
        statlabel_sel <- names(subcluster_sumstattab_seq[!names(subcluster_sumstattab_seq) %in% "PATNUM"])[statnum]
        stattab_sel <- subcluster_sumstattab_seq[[statlabel_sel]]
        
        # If numerical, then just grab the pval
        if (as.character(stattab_sel[1,1]) == "Min." & as.character(stattab_sel[2,1]) == "1st Qu.") {
            outtab <- cbind(statlabel_sel, TREATMNT_cat = "numeric", stattab_sel[1,"statval",drop=FALSE])
        } else {
            outtab <- cbind(statlabel_sel, stattab_sel[,1,drop=FALSE], stattab_sel[,"statval",drop=FALSE])
        }
        
        # If categorical, than grab the pval and the cat label and save each of those
        subcluster_outtablist[[statnum]] <- outtab
    }
    subcluster_statssummarytable <- do.call(rbind, subcluster_outtablist)
    write.table(subcluster_statssummarytable, paste0(outfilepathtabulation, "simplified_stats_by_", clusterlabel_sel, "_and_intervention", ".csv"),
                sep = ",", col.names = TRUE, row.names = FALSE)
    simplified_subcluster_outlist[[clusternum]] <- cbind(clusterlabel = clusterlabel_sel, subcluster_statssummarytable)
    
}
## Rerun for the intervention analysis too
full_simplified_stattable <- do.call(rbind, simplified_outlist)
full_simplified_stattable[,"combined_labels"] <- apply(full_simplified_stattable[,c("statlabel_sel", "nmfclusterlabel_cat")], 1, paste, collapse = "__")

simplified_stat_matrix <- dcast(full_simplified_stattable, combined_labels ~ clusterlabel, value.var = "statval")
rownames(simplified_stat_matrix) <- simplified_stat_matrix[,"combined_labels"]
simplified_stat_matrix <- simplified_stat_matrix[,!grepl("combined_labels", colnames(simplified_stat_matrix))]

write.table(simplified_stat_matrix,  paste0(outfilepathtabulation, "simplified_stats_ischmemia_by_cluster_matix", ".csv"),
            sep = ",", col.names = NA, row.names = TRUE)

## Rerun for the intervention analysis too
subcluster_simplified_stattable <- do.call(rbind, simplified_subcluster_outlist)
subcluster_simplified_stattable[,"combined_labels"] <- apply(subcluster_simplified_stattable[,c("statlabel_sel", "TREATMNT_cat")], 1, paste, collapse = "__")

subcluster_simplified_stat_matrix <- dcast(subcluster_simplified_stattable, combined_labels ~ clusterlabel, value.var = "statval")
rownames(subcluster_simplified_stat_matrix) <- subcluster_simplified_stat_matrix[,"combined_labels"]
subcluster_simplified_stat_matrix <- subcluster_simplified_stat_matrix[,!grepl("combined_labels", colnames(subcluster_simplified_stat_matrix))]

write.table(subcluster_simplified_stat_matrix,  paste0(outfilepathtabulation, "simplified_stats_ischmemia_by_cluster_intervention_matix", ".csv"),
            sep = ",", col.names = NA, row.names = TRUE)



## This was good for exploration - but I think I need odds ratios for these things too, so I have to go manual for each of these things - outputting p vals and odds ratios
dir.create(paste0(outfilepathtabulation, "COI_tabulations/"), showWarnings = FALSE, recursive = TRUE)
COI <- c("DEGRISCH", "IMGDEGIS", "CTMULT50", "CTANYD50", "DUKESCORE")
clusterlabels <- as.character(unique(nmfcluster_metadata_analysis_intab[,"nmfclusterlabel"]))
# nmfcluster_metadata_analysis_intab[,c("PATNUM", "nmfclusterlabel", COI)]

# Run over each category
statoutlist_category <- list()
for (COInum in seq_len(length(COI))) {
    COIlabel_sel <- COI[COInum]
    
    # For each category, run for each cluster
    statoutlist_cluster <- list()
    for (clusternum in seq_len(length(clusterlabels))) {
        cluster_sel <- clusterlabels[clusternum]
        datatable_sel <- nmfcluster_metadata_analysis_intab[,c("PATNUM", "nmfclusterlabel", COIlabel_sel)]
        datatable_sel[,"nmfclusterlabel"] <- ifelse(datatable_sel[,"nmfclusterlabel"] == cluster_sel, cluster_sel, paste0("not__", cluster_sel))
        
        featurelabels <- as.character(na.omit(unique(datatable_sel[,COIlabel_sel])))
        statoutlist_feature <- list()
        for (featurenum in seq_len(length(featurelabels))) {
            # Grab the info I need
            feature_sel <- featurelabels[featurenum]
            feature_datatable_sel <- datatable_sel
            feature_datatable_sel[,COIlabel_sel] <- ifelse(datatable_sel[,COIlabel_sel] == feature_sel, feature_sel, paste0("not__", feature_sel))
            
            # Create the fisher table for this analysis
            fishertab_sel <- dcast(data.frame(table(feature_datatable_sel[,2:3])), as.formula(paste0("nmfclusterlabel  ~ ", COIlabel_sel)), value.var = "Freq")
            rownames(fishertab_sel) <- fishertab_sel[,1]
            fishertab_sel <- fishertab_sel[c(cluster_sel, paste0("not__", cluster_sel)), c(feature_sel, paste0("not__", feature_sel))]
            
            # Run the fisher test, and then grab the results and save them out
            fisherout <- fisher.test(fishertab_sel)
            fisherout$p.value
            fisherout$estimate
            statoutlist_feature[[featurenum]] <- c(category = COIlabel_sel, feature = feature_sel, clusterlabel = cluster_sel, 
                                                   pvalue = fisherout$p.value, oddsratio = fisherout$estimate)
        }
        statoutlist_cluster[[clusternum]] <- do.call(rbind, statoutlist_feature)
    }
    statoutlist_category[[COInum]] <- do.call(rbind, statoutlist_cluster)
}
statouttable_full <- do.call(rbind, statoutlist_category)
write.table(statouttable_full, paste0(outfilepathtabulation, "COI_tabulations/", "COI_stat_table.csv"), col.names = TRUE, row.names = FALSE, sep = ",")




## What about comparing it to the biomarker data...
inbiomarkertable_file <- BIOMARKER_FILE
inbiomarkertable <- read.table(inbiomarkertable_file, sep = ",", header = TRUE)
biomarkerlabels <- colnames(inbiomarkertable)[grepl("_clean", colnames(inbiomarkertable))]

nmfcluster_biomarker_analysis_intab <- merge(cbind(PATNUM = rownames(nmfcluster_allsamples), nmfclusterlabel = nmfcluster_allsamples[,1]),
                                             inbiomarkertable[rownames(nmfcluster_allsamples) %in% inbioreptable[,"PATNUM"], grepl("_clean|PATNUM", colnames(inbiomarkertable))],
                                             by = "PATNUM")
# Run over each category
statoutlist_biomarker <- list()
for (biomarkernum in seq_len(length(biomarkerlabels))) {
    biomarkerlabel_sel <- biomarkerlabels[biomarkernum]
    
    # For each category, run for each cluster
    statoutlist_cluster <- list()
    for (clusternum in seq_len(length(clusterlabels))) {
        cluster_sel <- clusterlabels[clusternum]
        datatable_sel <- nmfcluster_biomarker_analysis_intab[,c("PATNUM", "nmfclusterlabel", biomarkerlabel_sel)]
        datatable_sel[,"nmfclusterlabel"] <- ifelse(datatable_sel[,"nmfclusterlabel"] == cluster_sel, cluster_sel, paste0("not__", cluster_sel))
        
        ## Now just a simple t.test between the two groups:
        pstatout <- t.test(datatable_sel[datatable_sel[,"nmfclusterlabel"] == cluster_sel, biomarkerlabel_sel],
                          datatable_sel[datatable_sel[,"nmfclusterlabel"] != cluster_sel, biomarkerlabel_sel])$p.value
        ratiostatout <- mean(datatable_sel[datatable_sel[,"nmfclusterlabel"] == cluster_sel, biomarkerlabel_sel], na.rm = TRUE) /
                        mean(datatable_sel[datatable_sel[,"nmfclusterlabel"] != cluster_sel, biomarkerlabel_sel], na.rm = TRUE)
        
        statoutlist_cluster[[clusternum]] <- c(category = biomarkerlabel_sel, clusterlabel = cluster_sel, 
                                                   pvalue = pstatout, oddsratio = ratiostatout)
    }
    statoutlist_biomarker[[biomarkernum]] <- do.call(rbind, statoutlist_cluster)
}
biomarker_statouttable_full <- do.call(rbind, statoutlist_biomarker)
write.table(biomarker_statouttable_full, paste0(outfilepathtabulation, "COI_tabulations/", "biomarker_stat_table.csv"), col.names = TRUE, row.names = FALSE, sep = ",")









## Plot out the stats for a quick viz
# Have to omit the time variables - they are always going to be really different
# plotCOI <- rownames(simplified_stat_matrix)[!grepl("^T_", rownames(simplified_stat_matrix))]
plotCOI <- rownames(simplified_stat_matrix)[grepl(paste(c("DEGRISCH", "IMGDEGIS", "CTMULT50", "CTANYD50", "DUKESCORE"), collapse = "|"), rownames(simplified_stat_matrix))]

# Create plottab
pvalcutoff <- 0.05
plottab <- simplified_stat_matrix[plotCOI,]
plottab <- apply(plottab, 2, function(x) as.numeric(as.character(x)))
rownames(plottab) <- rownames(simplified_stat_matrix[plotCOI,])
plottab[plottab > pvalcutoff] <- 1

# Make heatmapcolorparam
heatmapcolorparam <- colorRamp2(c(1, eventpvalcutoff + 0.0001, eventpvalcutoff, 0), c("white", "white", "#c9e0dc", "#00705E"))

create_heatmap(counttab = plottab, scale_data = FALSE, rowclusterparam = TRUE, colclusterparam = TRUE, heatmapcolorparam = heatmapcolorparam)








## Compare these imaging metrics: NUCPCTISC, SUMDIFFSC, SUMSTRSC, SUMRSTSC, ECHOISCSG, ECHOINFSG, CMRPERDF, CMRPERDI, ETTST1DP, ETTST2DP, ETTAGEPHR
dir.create(paste0(outfilepathintegration, "IMGMETRIC_analysis/"), showWarnings = FALSE, recursive = TRUE)

# Quick table of how many patients have what imagmod and how many severity
metasummarytable <- dcast(data.frame(table(inbioreptable[inbioreptable[,"PATNUM"] %in% colnames(normcounttable),c("IMAGMOD", "DEGRISCH")])), formula=IMAGMOD~DEGRISCH)
rownames(metasummarytable) <- metasummarytable[,"IMAGMOD"]
metasummarytable <- metasummarytable[,!grepl("IMAGMOD", colnames(metasummarytable))]
metasummarytable <- metasummarytable[,c("None", "Mild", "Moderate", "Severe", "Uninterpretable")]
metasummarytable["total",] <- colSums(metasummarytable)
metasummarytable[,"total"] <- rowSums(metasummarytable)
write.table(metasummarytable, paste0(outfilepathintegration, "IMGMETRIC_analysis/", "IMAGMOD_ISCHEMIA_summarytable.csv"), sep = ",", col.names = NA, row.names = TRUE)

# Ok - there are a lot of variables here, and each imaging modality has continuous metrics associated with it:
# NUCPCTISC	SUMDIFFSC	SUMSTRSC	SUMRSTSC - Nuclear (PET/SPECT)
# ECHOISCSG	ECHOINFSG - ECHO
# CMRPERDF	CMRPERDI - CMR
# ETTST1DP	ETTST2DP	ETTAGEPHR - ETT
## with raw gene expression, and also EIGENgene expression
## Grab the matched samples
incountfile <- file.path(DATA_DIR, "normcounttab.txt")
normcounttable <- read.table(incountfile, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

## eigengenes
# wgcna_eigengenes_tablefile <- "# PATH_UPDATED: output/run2b_rmoutliers2_controlage/WGCNA_power14/wgcna_eigengenes.csv"
wgcna_eigengenes_tablefile <- "# PATH_UPDATED: output/run2b_rmoutliers2_controlage/WGCNA/WGCNA_power16_size40/wgcna_eigengenes.csv"
wgcna_eigengenes_table <- read.table(wgcna_eigengenes_tablefile, sep = ",", header = TRUE, row.names = 1)
eigengenes <- colnames(wgcna_eigengenes_table)

## POI Eigengenes
# POI_eigengenes_tablefile <- "# PATH_UPDATED: output/figures/raw/run2_20210423/imaging/wgcna/POI_eigengene_tab.csv"
POI_eigengenes_tablefile <- "# PATH_UPDATED: output/figures/raw/run2_20210423/imaging/wgcna/POI_expanded_eigengene_tab.csv"
POI_eigengenes_table <- read.table(POI_eigengenes_tablefile, sep = ",", header = TRUE, row.names = 1)
POI_eigengenes <- colnames(POI_eigengenes_table)

## biorep table
inbioreptablefile <- "# PATH_UPDATED: data/biorep_10_27.csv"
inbioreptable <- read.table(inbioreptablefile, sep = ",", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("NA", "", NA))

# continuous imaging metrics
cont_img_metrics <- c("NUCPCTISC", "SUMDIFFSC", "SUMSTRSC", "SUMRSTSC", "ECHOISCSG", "ECHOINFSG", 
                      "CMRPERDF", "CMRPERDI", "ETTST1DP", "ETTST2DP", "ETTAGEPHR")

##
IMGMETRIC_cor_intab <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "PATNUM", all = TRUE, sort = FALSE), list(
    inbioreptable[,c("PATNUM", cont_img_metrics),drop=FALSE],
    cbind.data.frame(PATNUM = rownames(wgcna_eigengenes_table), wgcna_eigengenes_table),
    cbind.data.frame(PATNUM = rownames(POI_eigengenes_table), POI_eigengenes_table),
    cbind.data.frame(PATNUM = colnames(normcounttable), t(normcounttable))
))
rownames(IMGMETRIC_cor_intab) <- IMGMETRIC_cor_intab[,"PATNUM"]
IMGMETRIC_cor_intab <- IMGMETRIC_cor_intab[,!colnames(IMGMETRIC_cor_intab) %in% "PATNUM"]
# NUCPCT_cortab <- na.omit(NUCPCT_cortab)
IMGMETRIC_cor_intab <- IMGMETRIC_cor_intab[colnames(normcounttable),]


IMGMETRIC_cortab <- apply(as.matrix(IMGMETRIC_cor_intab[,!colnames(IMGMETRIC_cor_intab) %in% cont_img_metrics]), 2, function(x) {
    out <- corr.test(as.matrix(IMGMETRIC_cor_intab[,cont_img_metrics,drop=FALSE]), x, method = "spearman", use="pairwise.complete.obs")
    c(out$r[cont_img_metrics,1], out$p[cont_img_metrics,1])
})
## Ok, so a little hacky - but this double the rows, and for R val and then P val, so we have to add that to the names
rownames(IMGMETRIC_cortab) <- paste0(rownames(IMGMETRIC_cortab), c(rep("_Rval", nrow(IMGMETRIC_cortab)/2), rep("_Pval", nrow(IMGMETRIC_cortab)/2)))

# NUCPCT_cortabout <- merge(cbind.data.frame(gene = rownames(temp_cortab_out), temp_cortab_out), 
#                           cbind.data.frame(gene = rownames(deseq_imaging_table), deseq_imaging_table[,c("log2FoldChange", "padj")]), by = "gene")
# write.table(NUCPCT_cortabout, paste0(outfilepathintegration, "NUCPCTISC_analysis/", "NUCPCTISC_corr_and_DGE_outtab.csv"),
#             sep = ",", col.names = TRUE, row.names = FALSE)

# Run through this table and return all genes per metric - one for the p value and one for the r value... but has to be sig
IMGMETRIC_cortab_R <- IMGMETRIC_cortab[grepl("_Rval", rownames(IMGMETRIC_cortab)),]
IMGMETRIC_cortab_P <- IMGMETRIC_cortab[grepl("_Pval", rownames(IMGMETRIC_cortab)),]

## Write out these tables
write.table(t(IMGMETRIC_cortab_R), paste0(outfilepathintegration, "IMGMETRIC_analysis/", "IMGMETRIC_correlation_analyses/", "allmetric_geneandeigen_cortab_R.csv"), sep = ",", col.names = NA, row.names = TRUE)
write.table(t(IMGMETRIC_cortab_P), paste0(outfilepathintegration, "IMGMETRIC_analysis/", "IMGMETRIC_correlation_analyses/", "allmetric_geneandeigen_cortab_P.csv"), sep = ",", col.names = NA, row.names = TRUE)


## Write out a heatmap of the correlation values
IMGMETRIC_corR_hmtab <- t(IMGMETRIC_cortab_R)
IMGMETRIC_corP_hmtab <- t(IMGMETRIC_cortab_P)
colmetatable <- data.frame(IMAGMOD = c("NUC", "NUC", "NUC", "NUC", "ECHO", "ECHO", "CMR", "CMR", "ETT", "ETT", "ETT"),
                           row.names = rownames(IMGMETRIC_cortab_R))
colannotationlist <- annotationlist_builder(metatable = colmetatable)
rowmetatable <- data.frame(IMAGMOD = c(rep("Eigengene", 16), 
                                       rep("Pathway", sum(grepl("HALLMARK_|GO_|HP_", rownames(IMGMETRIC_corR_hmtab)))), 
                                       rep("Gene", (length(rownames(IMGMETRIC_corR_hmtab))-16-sum(grepl("HALLMARK_|GO_|HP_", rownames(IMGMETRIC_corR_hmtab)))))
                                       ),
                           row.names = colnames(IMGMETRIC_cortab_R))
rowannotationlist <- annotationlist_builder(metatable = rowmetatable)
hmout_Rval <- create_heatmap(IMGMETRIC_corR_hmtab, colmetatable = colmetatable, colannotationlist = colannotationlist, 
                             rowmetatable = rowmetatable, rowannotationlist = rowannotationlist,
                             scale_data = FALSE, rowclusterparam = TRUE)
pdf(paste0(outfilepathintegration, "IMGMETRIC_analysis/", "IMGMETRIC_correlation_Rval_summaryheatmap.pdf"))
draw(hmout_Rval[[1]])
junk <- dev.off()

heatmapcolorparam <- colorRamp2(breaks = c(1,0.05,0), c("white", "white", "darkred"))
hmout_Pval <- create_heatmap(IMGMETRIC_corP_hmtab, colmetatable = colmetatable, colannotationlist = colannotationlist, 
                             rowmetatable = rowmetatable, rowannotationlist = rowannotationlist,
                             scale_data = FALSE, rowclusterparam = TRUE, 
                             heatmapcolorparam = heatmapcolorparam)
pdf(paste0(outfilepathintegration, "IMGMETRIC_analysis/", "IMGMETRIC_correlation_Pval_summaryheatmap.pdf"))
draw(hmout_Pval[[1]])
junk <- dev.off()


IMGMETRIC_cortab_R_filtered <- IMGMETRIC_cortab_R
IMGMETRIC_cortab_R_filtered[IMGMETRIC_cortab_P > 0.05] <- NA

IMGMETRIC_cortab_R_filtered <- t(IMGMETRIC_cortab_R_filtered[,colSums(is.na(IMGMETRIC_cortab_R_filtered)) < nrow(IMGMETRIC_cortab_R_filtered)])

out1 <- apply(IMGMETRIC_cortab_R_filtered, 2, function(x) rownames(IMGMETRIC_cortab_R_filtered)[!is.na(x)])
pout <- create_upset_plot(upset_indata = out1, comboparam = NULL, transposeparam = FALSE)
## Ok, so there are 0 features that overlap across all...

dir.create(paste0(outfilepathintegration, "IMGMETRIC_analysis/", "IMGMETRIC_correlation_analyses/"), showWarnings = FALSE, recursive = TRUE)
## This is crude, but lets just do a hypergeo test for genes that go up and down for each variable, maybe we will get something...
for (imgmetricnum in seq_len(ncol(IMGMETRIC_cortab_R_filtered))) {
    
    imgmetricsel <- colnames(IMGMETRIC_cortab_R_filtered)[imgmetricnum]
    datasel <- na.omit(IMGMETRIC_cortab_R_filtered[,imgmetricsel,drop=FALSE])
    
    ## Create sub folder to pipe out to:
    dir.create(paste0(outfilepathintegration, "IMGMETRIC_analysis/", "IMGMETRIC_correlation_analyses/", imgmetricsel, "/"), showWarnings = FALSE, recursive = TRUE)
    
    UPgenes <- data.frame(rownames(datasel[datasel[,1] > 0,,drop=FALSE]))
    UPgenes <- UPgenes[UPgenes[,1] %in% rownames(normcounttable),,drop=FALSE]
    DOWNgenes <- data.frame(rownames(datasel[datasel[,1] < 0,,drop=FALSE]))
    DOWNgenes <- DOWNgenes[DOWNgenes[,1] %in% rownames(normcounttable),,drop=FALSE]
    
    hypergeo_genetest_out_GO_UP = hypergeo_genetest(UPgenes,
                                                 statcutoffparam = c("stattype" = "pvalue", "pstatcutoff" = 0.01, "log2fccutoff" = 0), ## Dummy list
                                                 genesetparam = c("C5"), speciesparam = "Homo sapiens")
    hypergeo_genetest_out_GO_DOWN = hypergeo_genetest(DOWNgenes,
                                                    statcutoffparam = c("stattype" = "pvalue", "pstatcutoff" = 0.01, "log2fccutoff" = 0), ## Dummy list
                                                    genesetparam = c("C5"), speciesparam = "Homo sapiens")
    
    write.table(datasel, paste0(outfilepathintegration, "IMGMETRIC_analysis/", "IMGMETRIC_correlation_analyses/", imgmetricsel, "/", imgmetricsel, "_corrvalues.csv"), 
                sep = ",", row.names = TRUE, col.names = NA)
    write.table(hypergeo_genetest_out_GO_UP$enricherUPout, 
                file = paste0(outfilepathintegration, "IMGMETRIC_analysis/", "IMGMETRIC_correlation_analyses/", imgmetricsel, "/", imgmetricsel, "_hypergeo_UP_GO.csv"), 
                sep = ",", row.names = TRUE, col.names = NA)
    write.table(hypergeo_genetest_out_GO_DOWN$enricherUPout, ## Purposeful hack
                file = paste0(outfilepathintegration, "IMGMETRIC_analysis/", "IMGMETRIC_correlation_analyses/", imgmetricsel, "/", imgmetricsel, "_hypergeo_DOWN_GO.csv"), 
                sep = ",", row.names = TRUE, col.names = NA)
    
    ## I also want to make a histo of the values - just for reference
    imgmetric_values <- na.omit(IMGMETRIC_cor_intab[,gsub("_Rval", "", imgmetricsel), drop = FALSE])
    pout <- plot_histogram(imgmetric_values, binparam = 20, 
                           labsparam = list(title = gsub("_Rval", "", imgmetricsel), x = "number of samples", y = paste0(gsub("_Rval", "", imgmetricsel), " value")))
    pdf(paste0(outfilepathintegration, "IMGMETRIC_analysis/", "IMGMETRIC_correlation_analyses/", imgmetricsel, "/", imgmetricsel, "_values_histogram.pdf"))
    print(pout)
    junk <- dev.off()
    
    ## And get a snapshot of the correlation values that come back:
    pout2 <- plot_histogram(data.frame(datasel), binparam = 20, 
                           labsparam = list(title = gsub("_Rval", "", imgmetricsel), x = "number of samples", y = paste0(gsub("_Rval", "", imgmetricsel), " value")))
    pdf(paste0(outfilepathintegration, "IMGMETRIC_analysis/", "IMGMETRIC_correlation_analyses/", imgmetricsel, "/", imgmetricsel, "_corvalues_histogram.pdf"))
    print(pout2)
    junk <- dev.off()
    
}




## Scatterplots of Sev vs MildNone avg for each modality...
dir.create(paste0(outfilepathintegration, "IMGMETRIC_analysis/", "IMAGMOD_comps/"), showWarnings = FALSE, recursive = TRUE)
IMAGMODparam <- unique(inbioreptable[,"IMAGMOD"])
plottab_outlist <- list()
for (IMAGMODnum in seq_len(length(IMAGMODparam))) {
    # Select modality and samples
    IMAGMODsel <- IMAGMODparam[IMAGMODnum]
    IMAGMODsel_cleaned <- gsub("\\(|\\)|\\/|\\ ", "", IMAGMODsel)
    samplesel <- intersect(inbioreptable[inbioreptable[,"IMAGMOD"] == IMAGMODsel,"PATNUM"], colnames(normcounttable))
    SEVEREsamples <- inbioreptable[inbioreptable[,"DEGRISCH"] %in% "Severe" & inbioreptable[,"PATNUM"] %in% samplesel, "PATNUM"]
    NONEMILDsamples <- inbioreptable[inbioreptable[,"DEGRISCH"] %in% c("None", "Mild") & inbioreptable[,"PATNUM"] %in% samplesel, "PATNUM"]
    
    # Select our count data, and split it by severe and none/mild
    SEVEREdatasel <- normcounttable[,SEVEREsamples,drop=FALSE]
    NONEMILDdatasel <- normcounttable[,NONEMILDsamples,drop=FALSE]
    
    plottab <- merge(data.frame(SEVERE = rowMeans(SEVEREdatasel)), data.frame(NONEMILD = rowMeans(NONEMILDdatasel)), by = "row.names")
    plottab[,"log2foldchange"] <- log2(plottab[,"SEVERE"] / plottab[,"NONEMILD"])
    plottab[,"wilcox_p"] <- apply(data.frame(seq_len(nrow(normcounttable))), 1, function(x) 
        wilcox.test(unlist(SEVEREdatasel[x,]), unlist(NONEMILDdatasel[x,]))$p.value)
    plottab[,"wilcox_adjp"] <- p.adjust(plottab[,"wilcox_p"], method = "fdr")
    plottab[,"significant"] <- ifelse(plottab[,"wilcox_adjp"] < 0.05, "blue", "darkgrey")
    # plottab[,c("SEVERE", "NONEMILD")] <- apply(plottab[,c("SEVERE", "NONEMILD")], 2, log2)
    
    scatout <- scatter_plotter(indata = plottab[,c(3,2)], colorvar = plottab[,"significant",drop=FALSE], 
                            labsparam = list(title = IMAGMODsel, x = "NONEMILD log2 expression", y = "SEVERE log2 expression"))
    scatout <- scatout + geom_abline(slope = 1, intercept = 0, linetype = 2)
    histout <- plot_histogram(data = plottab[,"log2foldchange",drop=FALSE], binparam = 30, labsparam = list(title = IMAGMODsel, x = "Log2FC", y = "Number of Genes"))
    
    pdf(paste0(outfilepathintegration, "IMGMETRIC_analysis/", "IMAGMOD_comps/", IMAGMODsel_cleaned, "_Sev_v_MildNone_scatter.pdf"))
    print(scatout)
    junk <- dev.off()
    
    pdf(paste0(outfilepathintegration, "IMGMETRIC_analysis/", "IMAGMOD_comps/", IMAGMODsel_cleaned, "_Sev_v_MildNone_log2fc_hist.pdf"))
    print(histout)
    junk <- dev.off()
    
    plottab_outlist[[IMAGMODnum]] <- plottab
    names(plottab_outlist)[IMAGMODnum] <- IMAGMODsel
    
}

plottab_outlist





## OK.. so no eigengenes came out :( so insteads, lets get a viz, because maybe I can bucket these things
# corrplotout <- ggpairs(data = IMGMETRIC_cor_intab[,colnames(IMGMETRIC_cor_intab) %in% c(eigengenes, cont_img_metrics)])
# tt1 <- temp_cortab_out[,colSums(temp_cortab_out[grepl("_Pval", rownames(temp_cortab_out)),] < 0.05) > 0]


# ## Ok, now what about just the combination of our NMF info with the ISCHEMIA severity info
# # Simple Calculations
# nmfcluster_allsamples = nmf_clustermembership_table[,"combined",drop=FALSE]
# nmfcluster_coresamples = nmf_clustercores_table[,"combined",drop=FALSE]
# inmetatable[,"IMGDEGIS",drop=FALSE]
# 
# nmf_isch_combtab = Reduce(function(dtf1, dtf2)
#     merge(dtf1, dtf2, by = "PATNUM", all.x = TRUE, sort = FALSE),
#     list(cbind(PATNUM = rownames(nmfcluster_allsamples), nmfcluster_allsamples),
#          cbind(PATNUM = rownames(nmfcluster_coresamples), nmfcluster_coresamples),
#          cbind(PATNUM = rownames(inmetatable), inmetatable[,"IMGDEGIS",drop=FALSE])
#     )
# )
# rownames(nmf_isch_combtab) = nmf_isch_combtab[,1]
# colnames(nmf_isch_combtab) <- c("PATNUM", "nmfclusters", "nmfcoreclusters", "IMGDEGIS")
# nmf_isch_combtab <- nmf_isch_combtab[,!grepl("PATNUM", colnames(nmf_isch_combtab))]
# 
# 
# nmf_isch_sumtab <- table(nmf_isch_combtab[,c("nmfclusters","IMGDEGIS")])
# nmf_isch_sumtab_perc <- nmf_isch_sumtab / rowSums(nmf_isch_sumtab)
# 
# nmfcore_isch_sumtab <- table(nmf_isch_combtab[,c("nmfcoreclusters","IMGDEGIS")])
# nmfcore_isch_sumtab_perc <- nmfcore_isch_sumtab / rowSums(nmfcore_isch_sumtab)
# 
# 
# # Simple Calculations for events
# nmfcluster_allsamples = nmf_clustermembership_table[,"combined",drop=FALSE]
# nmfcluster_coresamples = nmf_clustercores_table[,"combined",drop=FALSE]
# # inmetatable[,"C_PRIMARY",drop=FALSE]
# # inmetatable[,"T_PRIMARY",drop=FALSE]
# 
# ## Adding for less than a year
# ineventtable <- inmetatable[,"C_PRIMARY",drop=FALSE]
# # ineventtable[inmetatable[,"T_PRIMARY"] > 365 & inmetatable[,"C_PRIMARY"] != 0,] <- 2
# 
# nmf_event_combtab = Reduce(function(dtf1, dtf2)
#     merge(dtf1, dtf2, by = "PATNUM", all.x = TRUE, sort = FALSE),
#     list(cbind(PATNUM = rownames(nmfcluster_allsamples), nmfcluster_allsamples),
#          cbind(PATNUM = rownames(nmfcluster_coresamples), nmfcluster_coresamples),
#          # cbind(PATNUM = rownames(inmetatable), inmetatable[,"C_PRIMARY",drop=FALSE])
#          cbind(PATNUM = rownames(ineventtable), ineventtable[,"C_PRIMARY",drop=FALSE])
#     )
# )
# rownames(nmf_event_combtab) = nmf_event_combtab[,1]
# colnames(nmf_event_combtab) <- c("PATNUM", "nmfclusters", "nmfcoreclusters", "C_PRIMARY")
# nmf_event_combtab <- nmf_event_combtab[,!grepl("PATNUM", colnames(nmf_event_combtab))]
# 
# 
# nmf_event_sumtab <- table(nmf_event_combtab[,c("nmfclusters","C_PRIMARY")])
# nmf_event_sumtab_perc <- nmf_event_sumtab / rowSums(nmf_event_sumtab)
# 
# nmfcore_event_sumtab <- table(nmf_event_combtab[,c("nmfcoreclusters","C_PRIMARY")])
# nmfcore_event_sumtab_perc <- nmfcore_event_sumtab / rowSums(nmfcore_event_sumtab)



# 
# 
# 
# ## Can I create my own genesets - to specifically look for some processes that I suspect are there
# tt1 <- as.data.frame(msigdbr(species = "Homo sapiens", category = "C5")[,c("gs_name", "gene_symbol")])
# unique(tt1[grepl("INFLAMM", tt1[,1]),1])
# unique(tt1[grepl("THROMB", tt1[,1]),1])
# unique(tt1[grepl("COAG", tt1[,1]),1])
# unique(tt1[grepl("FIBR", tt1[,1]),1])







## Incorporate this into our data?
