#-----------------------------------------------------------
# preprocess data with both TMB and RNAseq/expression data
#----------------------------------------------------------
# There are nonsync mutation in common across three datasets
# We use it as mutation load

# get min and max APM value based on TCGA data
load("data/df_all.RData")
min_APM = min(df_all$APM)
max_APM = max(df_all$APM)
rm(df_all); gc()

library(GEOquery)
library(readxl)
library(tidyverse)

load("data/immune_cellType.RData")

##-----
## Ref1: Genomic and Transcriptomic Features of Response to Anti-PD-1 Therapy in Metastatic Melanoma. Cell 2016 Mar 24;165(1):35-44. PMID: 26997480

geo_dir = "data/GEOdata"

gse = "GSE78220"
GSE_78220 = getGEO(gse, GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = geo_dir)
GSE_78220 = GSE_78220$GSE78220_series_matrix.txt.gz

# download.file(url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE78220&format=file&file=GSE78220%5FPatientFPKM%2Exlsx", destfile = "data/GSE78220_FPKM.xlsx")

GSE78220_FPKM = read_xlsx("data/GSE78220_FPKM.xlsx")
GSE78220_FPKM = GSE78220_FPKM %>% 
    as.data.frame() %>% 
    mutate(mean_expr = rowMeans(.[, -1], na.rm = TRUE)) %>% 
    arrange(Gene, desc(mean_expr)) %>%
    distinct(Gene, .keep_all = TRUE) %>% 
    dplyr::select(-mean_expr) %>% 
    as.tibble()

# pData(GSE_78220) %>% View()

# applyGSVA(immune_cellType, group_col = "Cell_type", 
#           gene_col = "Symbol", ExprMatList = list(GSE78220_FPKM), method = "gsva", kcdf = "Poisson") -> gsva.gse78220
# 
# APM1_1 = (gsva.gse78220[[1]]$APM)

GSE78220_Norm = GSE78220_FPKM
GSE78220_Norm[,-1] = log2(GSE78220_Norm[,-1] + 1)
applyGSVA(immune_cellType, group_col = "Cell_type", 
          gene_col = "Symbol", ExprMatList = list(GSE78220_Norm), method = "gsva") -> gsva.gse78220_2

APM1_2 = gsva.gse78220_2[[1]]$APM
APM1_2 = (APM1_2 - min_APM) / (max_APM - min_APM) 

names(APM1_2) = colnames(GSE78220_Norm)[-1]
APM1_2
names(APM1_2) = sub(pattern = "(Pt.*)\\.baseline", "\\1", names(APM1_2))
# keep biggest value
APM1_2 = APM1_2[-21]
names(APM1_2)[c(14, 20)] = c("Pt16", "Pt27")

ge.GSE78220 = tibble(PatientID = names(APM1_2), APM=APM1_2)

# read clincal and tmb data
GSE_78220.TMB = read_csv("data/Cell2016.csv")
GSE_78220.TMB$nTMB = GSE_78220.TMB$TotalNonSyn / 38
#GSE_78220.TMB$nTMB = log(GSE_78220.TMB$TotalNonSyn+1)

Cell2016 = full_join(GSE_78220.TMB, ge.GSE78220, by = "PatientID")
Cell2016 = Cell2016[,c(1:5, 14,15, 6:13)]

Cell2016 %>% filter(!is.na(APM)) %>% 
    summarise(MedianAPM = median(APM, na.rm=TRUE), MedianTMB = median(nTMB, na.rm = TRUE), 
              ORR = sum(Response=="R")/n(), N = n())

# calculate TIGS score
#Cell2016$TIGS = log(Cell2016$nTMB + 1) * Cell2016$APM 
Cell2016$TIGS = Cell2016$nTMB * Cell2016$APM 

# save(Cell2016, file = "data/Cell2016.RData")
write_csv(Cell2016, "results/Cell2016_results.csv")

#---------------------------------------------------------------------------------------
# Science 2015: Genomic correlates of response to CTLA4 blockade in metastatic melanoma
#---------------------------------------------------------------------------------------
library(readxl)
science2015_TMB = read_xlsx("data/Science2015_TMB_List_AllPatients.xlsx")
science2015_cli1 = read_xlsx("data/Science2015_Clinical.xlsx", sheet = 1)
science2015_cli2 = read_xlsx("data/Science2015_Clinical.xlsx", sheet = 2)
science2015_rna = read_tsv("data/MEL-IPI-Share.rpkm.gct")
# science2015_rna_raw = read_tsv("data/MEL-IPI-Share.reads.gct")

#---- Preprocess
vc.nonSilent = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
                 "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del",
                 "In_Frame_Ins", "Missense_Mutation")

science2015_TMB %>% group_by(patient) %>% 
    summarise(TotalMutation=n(), 
              TotalNonSyn=sum(Variant_Classification %in% vc.nonSilent), 
             nTMB = TotalNonSyn / 38) -> science2015_TMB_summary

# science2015_TMB = filter(science2015_TMB, Variant_Classification %in% vc.nonSilent)
# 
# science2015_TMB_summary = science2015_TMB %>% group_by(patient) %>% 
#     summarize(TotalMutation=n(), nTMB=TotalMutation/38)
# summary(science2015_TMB_summary)

#---- remove duplicate symbols of rna data
science2015_rna %>% mutate(symbol = sub("(.*)\\..*$", "\\1", Description)) %>% 
    dplyr::select(-Name, -Description) %>% 
    dplyr::select(symbol, everything()) %>% 
    as.data.frame() %>% 
    mutate(mean_expr = rowMeans(.[, -1], na.rm = TRUE)) %>% 
    arrange(symbol, desc(mean_expr)) %>%
    distinct(symbol, .keep_all = TRUE) %>% 
    dplyr::select(-mean_expr) %>% 
    as.tibble() -> science2015_ge

length(intersect(science2015_ge$symbol, GSE78220_Norm$Gene))
# the gene number of science paper is bigger than other two paper
# maybe miRNA and some other RNA included
# keep them almost equal\
science2015_ge = science2015_ge %>% 
    filter(symbol %in% GSE78220_Norm$Gene)

science2015_ge[, -1] = log2(science2015_ge[, -1] + 1)

applyGSVA(immune_cellType, group_col = "Cell_type", 
          gene_col = "Symbol", ExprMatList = list(science2015_ge), method = "gsva") -> gsva.science2015

science2015_APM = gsva.science2015[[1]]$APM
names(science2015_APM) = colnames(science2015_ge)[-1]
science2015_APM = (science2015_APM - min_APM) / (max_APM - min_APM)

science2015_APM_df = tibble(patient = names(science2015_APM), APM = science2015_APM)
science2015_APM_df = science2015_APM_df %>% 
    mutate(patient = sub("^MEL.IPI_(.*)\\.Tumor.*$", "\\1", patient))

science2015_cli = bind_rows(
    dplyr::select(science2015_cli1, patient, age_start, RECIST, overall_survival, progression_free, primary, group, 
           histology, stage, gender, dead, progression, neos50),
    dplyr::select(science2015_cli2, patient, age_start, RECIST, overall_survival, progression_free, primary,
           group, histology, stage, gender, dead, progression) %>% filter(patient %in% c("Pat20", "Pat91"))
)
# combine three datasets
science2015 = full_join(full_join(science2015_cli, science2015_APM_df, by = "patient"), 
                         science2015_TMB_summary, by="patient")

# to avoid TIGS <0, TIGS = log(TMB +1) *APM for this dataset
science2015$TIGS = science2015$APM * log(science2015$nTMB +1)
write_csv(x = science2015, path = "results/science2015_results.csv")


#----------------------------------------------------------------
# Contribution of systemic and somatic factors to clinical response and resistance in urothelial cancer: an exploratory multi-omic analysis
#---------------------------------------------------------------
# Data Source: https://github.com/XSLiuLab/multi-omic-urothelial-anti-pdl1

library(tidyverse)

uro_cli = read_csv("~/biodata/Own/multi-omic-urothelial-anti-pdl1/data_clinical.csv")
uro_counts = read_csv("~/biodata/Own/multi-omic-urothelial-anti-pdl1/data_counts.csv")
uro_variants = read_csv("~/biodata/Own/multi-omic-urothelial-anti-pdl1/data_variants.csv",
                        col_types = "c??????")
uro_rna = read_csv("~/biodata/Own/multi-omic-urothelial-anti-pdl1/data_kallisto.csv")

uro_rna_wide = reshape2::dcast(uro_rna, gene_name ~ patient_id, value.var = "est_counts")

# try remove duplicated genes
uro_rna_wide %>% 
    as.data.frame() %>% 
    mutate(mean_expr = rowMeans(.[, -1], na.rm = TRUE)) %>% 
    arrange(gene_name, desc(mean_expr)) %>%
    distinct(gene_name, .keep_all = TRUE) %>% 
    dplyr::select(-mean_expr) %>% 
    as.tibble() -> uro_rna_ge

length(intersect(uro_rna_ge$gene_name, GSE78220_FPKM$Gene))
uro_rna_ge = dplyr::filter(uro_rna_ge, gene_name %in% GSE78220_FPKM$Gene)
uro_rna_ge_Norm = uro_rna_ge
uro_rna_ge_Norm[, -1] = log2(uro_rna_ge_Norm[,-1]+1)
applyGSVA(immune_cellType, group_col = "Cell_type", 
          gene_col = "Symbol", ExprMatList = list(uro_rna_ge_Norm), method = "gsva") -> gsva.uro
APM.uro = gsva.uro[[1]]$APM
APM.uro = (APM.uro - min_APM) / (max_APM - min_APM)

APM.uro = tibble(patient_id = colnames(uro_rna_ge_Norm)[-1], APM = APM.uro)
uro_cli_2 = uro_cli %>% 
    dplyr::select(patient_id, Age, Sex, is_benefit, is_benefit_os, os, `Alive Status`) %>% 
    mutate(event = ifelse(`Alive Status` == 'Y', 0, 1))

uro_tmb = uro_variants %>% 
    group_by(patient_id) %>% 
    summarise(TMB = n(), Nonsyn = sum(!is.na(gene_name)), nTMB = TMB/38) %>% #/38 
    mutate(patient_id = as.character(patient_id))

uro2017 = dplyr::full_join(uro_cli_2, dplyr::full_join(uro_tmb, APM.uro))
uro2017$TIGS = uro2017$APM * log(uro2017$nTMB)
#uro2017$TIGS2 = uro2017$APM * log(uro2017$nTMB)
 

write_csv(x = uro2017, path = "results/urothelial2017_results.csv")
