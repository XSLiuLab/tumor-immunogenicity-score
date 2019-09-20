# FUN: Clean PHBR scores
#
# PHBR for MHC-I was obtained from Prof. Carter
# PHBR for MHC-II was obtained from PHBR MHC-II paper Table S2
#
# Reference:
# #1 MHC-I Genotype Restricts the Oncogenic Mutational Landscape
# #2 Evolutionary Pressure against MHC Class II Binding Cancer Mutations

library(data.table)
library(readxl)
library(tidyverse)

PHBR_I = fread("data/patient_affinities.cancer.PHBR.csv.gz")
PHBR_I[1:5, 1:5]
colnames(PHBR_I)[1] = "Patient_ID"


PHBR_II = read_xlsx("data/1-s2.0-S0092867418311097-mmc3.xlsx")
PHBR_II[1:5, 1:5]
colnames(PHBR_II)[1] = "Patient_ID"

# To compare PHBR with APS and TIGS, 
# we need just one score for each patient,
# thus we use PHBR median to represent presentation ability
# for each patient
PHBR_I.score = PHBR_I %>% 
    as_tibble() %>% 
    column_to_rownames("Patient_ID") %>% 
    apply(1, function(x) median(x, na.rm = TRUE))

PHBR_II.score = PHBR_II %>% 
    column_to_rownames("Patient_ID") %>% 
    apply(1, function(x) median(x, na.rm = TRUE))

head(PHBR_I.score)    
head(PHBR_I.score)

save(PHBR_I.score, file = "data/TCGA.PHBR_I.score.RData")
save(PHBR_II.score, file = "data/TCGA.PHBR_II.score.RData")

