# There are more than 10 MHC II genes listed in Wikipedia
# https://en.wikipedia.org/wiki/MHC_class_II
#
# According to info from TCIA database,
# there are 7 classic genes
library(tidyverse)
source("code/functions.R")

MHC.all = c("HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB",
            "HLA-DPA1","HLA-DPB1", "HLA-DQA1", "HLA-DQA2",	"HLA-DQB1", "HLA-DQB2",
            "HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5")

MHC.classic = c("HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DRA", "HLA-DRB1")

load("report/results/merged_geneList.RData")
APM = merged_geneList %>% 
    dplyr::filter(Cell_type == "APM")

APM
load("report/results/TCGA_RNASeq_PanCancer.RData")

all(MHC.all %in% RNASeq_pancan$sample)  # not all genes in RNASeq data
all(MHC.classic %in% RNASeq_pancan$sample)

# Just using classic genes for GSVA

APM_MHC_II = APM %>% 
    filter(! Symbol %in% c("B2M", "HLA-A", "HLA-B", "HLA-C")) %>% # Remove MHC-I genes
    bind_rows(
        dplyr::tibble(
            Cell_type = "APM",
            Symbol = MHC.classic,
            Name = MHC.classic,
            inRNAseq = "Yes"
        )
    ) %>% 
    mutate(Cell_type = "APM_MHC_II")

res_MHCII.GSVA = applyGSVA(APM_MHC_II, group_col = "Cell_type", gene_col = "Symbol", ExprMatList = list(RNASeq_pancan), method = "gsva")

res_MHCII.GSVA = res_MHCII.GSVA[[1]]

save(APM_MHC_II, file = "data/APS_MHC_II_genes.RData")
save(res_MHCII.GSVA, file = "report/results/res_APS_MHC_II.GSVA.RData")
