# Available Resource and Workflow of Paper *Antigen presentation and tumor immunogenicity in cancer immunotherapy response prediction*

**There are two directories under this repository, which store two types of resource in the paper: data and code**.


## Data

* summary_data_new_20180906.csv - curated data for ORR linear model and prediction analysis.
* curated_TCGA_APM_TMB_Clinical.RData - curated TCGA pan-cancer data which clinical information, APM, TMB etc. available in RData format (which can be loaded by R). This is a main result of TCGA data process and used to analyze and plot figures showed in paper.
* df_combine_gsva_clinical.RData - TCGA pan-cancer data which combine clinical data and GSVA result from GSVA in RData format.

## Code

* download-TCGA-datasets-from-UCSCXena.R - for TCGA pancan clinical data and RNASeq data download.
* APM-pancan-explore.R - TCGA pan-cancer data preprocessing and EDA.
* APM-pancan-analysis.R - TCGA pan-cancer analysis (most of figures related to TCGA come from this).
* APM-enrichment.R - GSEA based APM score level using gene sets from MSigDB (hallmark, reactome, KEGG etc.).
* preprocess_data_with_both_TMB_and_GE.R - preprocess ICI datasets with both TMB and Gene Expression available.
* plotROC.R - ROC and Survival analysis for ICI datasets with both TMB and GE available.
* APM_TMB_and_ORR.R - pan-cancer ORR prediction analysis using pooled APM, TMB and TIGS.
* processTable.R - process ORR reference information to more concise summary table.

Other resource are available as supplementary files/tables of paper.

>NOTE: The code are not tidy, which means you cannot run automatically to reproduce all results (include figures and tables) of paper. If you are really interested in reproducing such a result, please run code one by one, read comments carefully and check input of code required. If you have any question, please open an issue or email to us.