# Available Resource and Workflow of Paper *Antigen presentation and tumor immunogenicity in cancer immunotherapy response prediction*


## Repository structure

* __code__ directory - tidy R functions
* __data__ directory - preprocessed data used for analysis and share
* __report__ directory - Rmarkdown analysis report and corresponding html version
* __manuscript__ directory - figures and related data used in manuscript
* __LICENSE__ file 
* This __README.md__ file

## How to use this repository

### For reader

For reader want to know more detail about methods, procedure and results in manuscript, you can download or clone this repository and open html file named **report_main.html** under report directory with your favorite browser.
 
### For reproducing results

TO DO.

## File description of data and code

### data

* summary_data_new_20180906.csv - curated data for ORR linear model and prediction analysis.
* curated_TCGA_APM_TMB_Clinical.RData - curated TCGA pan-cancer data which clinical information, APM, TMB etc. available in RData format (which can be loaded by R). This is a main result of TCGA data process and used to analyze and plot figures showed in paper.
* df_combine_gsva_clinical.RData - TCGA pan-cancer data which combine clinical data and GSVA result from GSVA in RData format.

### src （这里要修改和删除）

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

## Acknowledgement

We thank all the referees for providing valuable comments to promote the quality of this manuscript. We thank the authors and participating patients of immunotherapy publications for providing the data for this analysis. Thank TCGA project for making cancer genomics data available for analysis. Thank members of Liu lab for helpful discussion. This work was supported in part by the Shanghai Pujiang Program (16PJ1407400), The National Natural Science Foundation of China (31771373), and startup funding from ShanghaiTech University.


## Citation

If you use data, results or conclusion from this work, please cite:

*Antigen presentation and tumor immunogenicity in cancer immunotherapy response prediction, __Nature Communication__* (**Under Review**)


## LICENSE

Please note this work under __Apache License v2__ license, copyright belongs to Shixiang Wang and Xue-Song Liu, more detail please see [description of license](LICENSE).
