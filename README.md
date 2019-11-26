# Tumor Immunogenicity Score (TIGS)

## Contents

* [Overview](#overview)
* [Repo Contents](#repo-contents)
* [Instructions for Use](#instructions-for-use)
  * [Read analysis report](#read-analysis-report)
  * [Obtain data](#obtain-data)
  * [Reproduce analysis](#reproduce-analysis)
* [Test Environment](#test-environment)
* [Citation](#citation)
* [Acknowledgement](#acknowledgement)
* [LICENSE](#license)

## Overview

Purpose of this repository is to share analysis procedure, data and help readers or reviewers to know more detail of this work, reproduce or make use of results they are interested in.

## Repo Contents

* [code](./code): tidy R functions and R script for  GSVA score calculation with randomly selected APS/IIS genes.
* [data](./data): preprocessed data used for analysis and share
* [docs](./docs): website pages & data/figure files for showing analysis reports
* [report](./report): Rmarkdown analysis report and results of manuscript
  * [report/results](./report/results): important middle results and final results, most of them are in form of `.RData`, which can be easily loaded and operated by R. 
* __LICENSE__ file 

## Instructions for Use

### Read analysis report

For readers who want to know more details about methods, procedures or results in manuscript, you can read our analysis report, which is availble at <https://xsliulab.github.io/tumor-immunogenicity-score/>. 

### Obtain data

For readers who want to obtain raw/result data, please read corresponding preprocessing/analysis part in our [analysis report](https://xsliulab.github.io/tumor-immunogenicity-score/), locate data file, then download it with one of following ways:

* In Github, download file by clicking either `Download` button or `Raw` button at corresponding data page

* Use linux command `wget` or `curl`, fo example, you can download APM gene list by

  `wget https://github.com/XSLiuLab/tumor-immunogenicity-score/blob/master/data/APM.csv`

Or you can download whole respository with one of following ways:

* Clone this repository with `git clone https://github.com/XSLiuLab/tumor-immunogenicity-score.git`
* Download whole respository by clicking `Download` button at top right of url page <https://github.com/XSLiuLab/tumor-immunogenicity-score>

### Reproduce analysis

For readers who want to reproduce analysis shown in manuscript, please [install R](https://cran.r-project.org) in your computer, then install required R packages described at __Dependencies__ part of our [analysis report](https://xsliulab.github.io/tumor-immunogenicity-score/), followed by rendering `report_main.Rmd` file using [knitr package](https://github.com/yihui/knitr).

## Test Environment

* System: __MacOS__

* Software: __R v3.6.0__

* R packages:

  * [UCSCXenaTools](https://github.com/ShixiangWang/UCSCXenaTools) - download data from UCSC Xena
  * GEOquery - download data from NCBI GEO database
  * [tidyverse](https://www.tidyverse.org/) - operate data, plot
  * data.table - operate data
  * survival - built in R, used to do survival analysis 
  * metafor, metawho - meta-analysis
  * forestmodel - generate forestplot for meta-analysis model
  * forestplot - plot forestplot
  * survminer - plot survival fit
  * pROC - ROC analysis and visualization
  * [TCGAmutations](https://github.com/PoisonAlien/TCGAmutations) - download TCGA mutation data
  * [DT](https://cran.r-project.org/web/packages/DT/index.html) - show data table as a table in html
  * [GSVA](https://github.com/rcastelo/GSVA) - GSVA algorithm implementation
  * [ggstatsplot](https://github.com/IndrajeetPatil/ggstatsplot) - plot scatter with linear fit
  * [corrplot](https://cran.r-project.org/web/packages/corrplot/) - plot correlation 
  * knitr, rmdformats - used to compile this file
  * readxl - read xlsx data

These R packages are easily searched by internet and have no strict version requirements to reproduce the analyses.

## Citation

If you use data, results or conclusion from this work, please cite:

*Antigen presentation and tumor immunogenicity in cancer immunotherapy response prediction, __eLife__*. https://doi.org/10.7554/eLife.49020

## Acknowledgement

We thank all the reviewers for providing valuable comments to promote the quality of this manuscript. We thank the authors and participating patients of immunotherapy publications for providing the data for this analysis. Thank TCGA project for making cancer genomics data available for analysis. Thank members of Liu lab for helpful discussion. This work was supported in part by the Shanghai Pujiang Program (16PJ1407400), The National Natural Science Foundation of China (31771373), and startup funding from ShanghaiTech University.

## LICENSE

Code and documents of this work are made available for non commercial research purposes only under __Apache License v2__ license, more detail please see [description of license](LICENSE). However, the study has been applied for a national patent in China, notwithstanding any provision of the __Apache License v2__ License, the code/idea currently may not be used for commercial purposes without explicit written permission after contacting Xue-Song Liu <liuxs@shanghaitech.edu.cn> or Shixiang Wang <wangshx@shanghaitech.edu.cn>.

***

**Cancer Biology Group @ShanghaiTech**

**Research group led by Xue-Song Liu in ShanghaiTech. University**
