# ISCHEMIA Multiomics Analysis Scripts

Publication-ready analysis pipeline for ISCHEMIA multiomics study. This repository contains R scripts for data preprocessing, differential expression analysis, multiomic subtype identification, biomarker associations, drug target identification, and additive risk modeling. The scripts are organized to facilitate reproducibility and clarity, while keeping to the original functionality.

## Directory Structure

### Configuration & Utilities

- **`config.R`** - Central configuration (paths, parameters, random seeds)
- **`utils.R`** - Shared utility functions

### Core Analysis Scripts

Run these scripts in order for the main publication analyses:

1. **`preprocessing_differential_expression.R`** - RNA-seq QC, normalization, DESeq2, GSEA
2. **`subtypes_biomarker_associations.R`** - Biomarker associations with multiomics subtypes
3. **`subtypes_drug_interactions.R`** - Drug target identification using DGiDb
4. **`subtypes_multiomics_associations.R`** - Multi-omic feature associations (plaque, deconvolution, risk scores)
5. **`subtypes_additive_risk_models.R`** - Additive risk modeling

### Supporting Analysis Scripts

Additional analyses for specific components:

- **`cohort_tabulations.R`** - Cohort descriptive statistics
- **`dge_comparisons.R`** - Differential expression comparisons
- **`exploratory_analysis.R`** - Exploratory data analysis
- **`exploratory_analysis_functions.R`** - Helper functions for EDA
- **`integrative_analysis.R`** - Integrative analysis (version 1)
- **`integrative_analysis_v2.R`** - Integrative analysis (version 2)
- **`methylation_model_step1.R`** - Methylation preprocessing step 1
- **`methylation_processing.R`** - Methylation data processing
- **`multiomic_subtype_associations.R`** - Additional subtype associations
- **`nmf_multiomics.R`** - NMF-based multiomics clustering
- **`pace_validation.R`** - PACE cohort validation
- **`wgcna_analysis.R`** - WGCNA network analysis

## Contact

Matthew Muller - <mm12865@nyu.edu>
