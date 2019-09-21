# CIBERSORT result for TCGA patients
# were obtained from TCIA <https://tcia.at/home>

# library(tidyverse)
# df = data.table::fread("/Volumes/data/biodata/TCGA/TCGA-all-cellTypeFractionsAll-20180702.tsv")
# cibersort = df %>% 
#     as_tibble() %>% 
#     select(-quanTIseq_lsei_TIL10) %>% 
#     filter(!is.na(cibersort_LM22))
# 
# save(cibersort, file = "data/TCGA_cibersort.RData")

# This data only contains 6 cell types
# We also obtained cibersort result from TCGA publications https://gdc.cancer.gov/about-data/publications/panimmune
cibersort = data.table::fread("data/TCGA.Kallisto.fullIDs.cibersort.relative.tsv.gz")
# This data does not need clean, just read it into R

