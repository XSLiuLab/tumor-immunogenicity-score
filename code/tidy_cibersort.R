# CIBERSORT result for TCGA patients
# were obtained from TCIA <https://tcia.at/home>
library(tidyverse)
df = data.table::fread("/Volumes/data/biodata/TCGA/TCGA-all-cellTypeFractionsAll-20180702.tsv")
cibersort = df %>% 
    as_tibble() %>% 
    select(-quanTIseq_lsei_TIL10) %>% 
    filter(!is.na(cibersort_LM22))

save(cibersort, file = "data/TCGA_cibersort.RData")
