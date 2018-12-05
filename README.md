# Available Resource and Workflow of Paper *Antigen presentation and tumor immunogenicity in cancer immunotherapy response prediction*

Purpose of this repository is to share analysis procedure, data and help readers or reviewers to know more detail of this work, reproduce or make use of results they are interested in.

## Repository structure

* __code__ directory - tidy R functions
* __data__ directory - preprocessed data used for analysis and share
* __report__ directory - Rmarkdown analysis report and corresponding html version
  * `report/results` directory - important middle results and final results, most of them are in form of `.RData`, which can be easily loaded and operated by R. 
* __LICENSE__ file 

## Usage

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

For readers who want to reproduce analysis shown in manuscript, please [install R](https://cran.r-project.org) in your computer, then install required R packages described at __Dependencies__ part of our [analysis report](https://xsliulab.github.io/tumor-immunogenicity-score/), followed by rendering `report_main.Rmd` file using [`knitr` package](https://github.com/yihui/knitr).

## Citation

If you use data, results or conclusion from this work, please cite:

*Antigen presentation and tumor immunogenicity in cancer immunotherapy response prediction, __Nature Communication__* (**Under Review**)

## Acknowledgement

We thank all the reviewers for providing valuable comments to promote the quality of this manuscript. We thank the authors and participating patients of immunotherapy publications for providing the data for this analysis. Thank TCGA project for making cancer genomics data available for analysis. Thank members of Liu lab for helpful discussion. This work was supported in part by the Shanghai Pujiang Program (16PJ1407400), The National Natural Science Foundation of China (31771373), and startup funding from ShanghaiTech University.



## LICENSE

Please note this work under __Apache License v2__ license, copyright belongs to Shixiang Wang and Xue-Song Liu, more detail please see [description of license](LICENSE).
