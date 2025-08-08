################################
# Name: MacIntosh Cornwell
# Email: mcornwell1957@gmail.com
################################
## ISCHEMIA Cohort Characteristic Tabulations

## Load in Libraries
packagelist = c("Hmisc", "tools")
junk <- lapply(packagelist, function(xxx) suppressMessages(
  require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))

source("/Users/tosh/Desktop/Ruggles_Lab/code/mgc_plotting_functions.R")
source("/Users/tosh/Desktop/Ruggles_Lab/code/summarize_table_function.R")
source("/Users/tosh/Desktop/Ruggles_Lab/code/overlap_finder_function.R")

## Outpath
outfilepathstats = paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/ischemia2021/output/cohort_tabulations/")
dir.create(outfilepathstats, recursive = TRUE, showWarnings = FALSE)

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

dir.create(paste0(outfilepathstats, "sample_datatype_overlap_analysis/"), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(outfilepathstats, "sample_datatype_overlap_analysis/venndiagram.pdf"))
grid.draw(vennplot)
junk <- dev.off()
write.table(overlaptable, paste0(outfilepathstats, "sample_datatype_overlap_analysis/overlaptable.csv"), sep = ",", row.names = FALSE, col.names = TRUE)
write.table(overlapgrouptab, paste0(outfilepathstats, "sample_datatype_overlap_analysis/overlapgrouptab.csv"), sep = ",", row.names = FALSE, col.names = TRUE)
write.table(overlapsummary, paste0(outfilepathstats, "sample_datatype_overlap_analysis/overlapsummary.csv"), sep = ",", row.names = TRUE, col.names = TRUE)



## I think an upset plot here would be much more attractive
upset_indata <- overlap_inlist
upsetout <- create_upset_plot(upset_indata, comboparam = "all", transposeparam = TRUE)
upsetobject <- upsetout$upsetobject
upsetplot <- upsetout$upsetplot

pdf(paste0(outfilepathstats, "sample_datatype_overlap_analysis/upset_plot.pdf"))
draw(upsetplot)
junk <- dev.off()


## COI
COI <- c("SEX", "RACE", "BMI", "AGE_RAND", "SMOKSTAT", "HDLC", "LDLC", "RVBPDIA", "RVBPSYS", "HYPTENSE", "DIABETES", "IMGDEGIS", "CTNDV50", "CTMULT50")

## Create table to summarize by
summaryintable <- merge(overlapgrouptab, inbioreptable[,c("PATNUM", COI)], by.x = "features", by.y = "PATNUM", all.x = TRUE)
colnames(summaryintable)[1] <- "PATNUM"
rownames(summaryintable) <- summaryintable[,1]

## Now get the stats for each type on its own, and then the methyl/rna group
GOI <- list(RNA = summaryintable[grepl("RNAseq_PATNUMS", summaryintable[,"overlapgroup"]),"PATNUM"],
            Methyl = summaryintable[grepl("Methyl_PATNUMS", summaryintable[,"overlapgroup"]),"PATNUM"],
            Metab = summaryintable[grepl("Metabalomics_PATNUMS", summaryintable[,"overlapgroup"]),"PATNUM"],
            Biomarker = summaryintable[grepl("Biomarker_PATNUMS", summaryintable[,"overlapgroup"]),"PATNUM"],
            RNAMethyl = summaryintable[grepl("RNAseq_PATNUMS", summaryintable[,"overlapgroup"]) & 
                                       grepl("Methyl_PATNUMS", summaryintable[,"overlapgroup"]),"PATNUM"]
            )

summaryintablelist <- list()
for (GOInum in seq_len(length(GOI))) {
  GOIsel <- GOI[[GOInum]]
  GOIlabel <- names(GOI[GOInum])
  summarytable_sel <- summaryintable[GOIsel,]
  summaryintablelist[[GOInum]] <- cbind(summarytable_sel, cohort = GOIlabel)
  
  summarystatfile <- paste0(outfilepathstats, "summarystats_", GOIlabel, ".csv")
  sumstattab_seq <- summarize_table(intable = summarytable_sel, groupvar = NULL, outfile = summarystatfile, calc_stats = FALSE)
}
fullsummaryintable <- do.call(rbind, summaryintablelist)

summarystatfile <- paste0(outfilepathstats, "summarystats_across_cohorts.csv")
sumstattab_seq <- summarize_table(intable = fullsummaryintable, groupvar = "cohort", outfile = summarystatfile, calc_stats = TRUE)

