library(tidyverse)
library(pROC)

load("report/results/data_for_ICB_ROC.RData")
load("data/ICB_PHBR_I.RData")
hugo_gender = read_csv("data/Hugo_gender_info.csv")

tidy_df = dplyr::bind_rows(
    info_VanAllen2 %>% 
        dplyr::select(patient, gender, Response, APS1:APS100, 
                      APS, IIS, nTMB, TIGS, CD8, IFNG, IFNG.GS, ISG.RS, PDL1, TIDE, APS7_Sig, APS_MHC_II, 
                      CD8_comb_TMB, IFNG_comb_TMB, IFNG.GS_comb_TMB, ISG.RS_comb_TMB,    
                      PDL1_comb_TMB, TIDE_comb_TMB, APS7_Sig_comb_TMB, APS_MHC_II_comb_TMB) %>% 
        dplyr::mutate(
            Study = "Van2015",
            gender = dplyr::case_when(
                gender == "female" ~ "F",
                gender == "male" ~ "M",
                TRUE ~ NA_character_
            )
        ), 
    info_hugo %>% 
        dplyr::select(PatientID, Response, APS1:APS100, 
                      APS, IIS, nTMB, TIGS, CD8, IFNG, IFNG.GS, ISG.RS, PDL1, TIDE, APS7_Sig, APS_MHC_II, 
                      CD8_comb_TMB, IFNG_comb_TMB, IFNG.GS_comb_TMB, ISG.RS_comb_TMB,    
                      PDL1_comb_TMB, TIDE_comb_TMB, APS7_Sig_comb_TMB, APS_MHC_II_comb_TMB) %>% 
        dplyr::mutate(
            Study = "Hugo2016"
        ) %>% 
        dplyr::left_join(
            hugo_gender
        ) %>% 
        dplyr::rename(patient = PatientID, 
                      gender = Gender),
    info_Snyder %>% 
        dplyr::select(patient_id, is_benefit, Sex, APS1:APS100, 
                      APS, IIS, nTMB, TIGS, CD8, IFNG, IFNG.GS, ISG.RS, PDL1, TIDE, APS7_Sig, APS_MHC_II, 
                      CD8_comb_TMB, IFNG_comb_TMB, IFNG.GS_comb_TMB, ISG.RS_comb_TMB,    
                      PDL1_comb_TMB, TIDE_comb_TMB, APS7_Sig_comb_TMB, APS_MHC_II_comb_TMB) %>% 
        dplyr::mutate(
            Study = "Snyder2017",
            is_benefit = dplyr::case_when(
                is_benefit == "True" ~ "R",
                is_benefit == "False" ~ "NR",
                TRUE ~ NA_character_
            )
        ) %>% 
        dplyr::rename(patient = patient_id, 
                      gender = Sex, 
                      Response = is_benefit)
        
) %>% 
    dplyr::left_join(ICB_PHBR_I %>% dplyr::select(-study)) %>% 
    dplyr::slice(1:170)  # Remove 3 patients containing all NAs


# Compare PHBR_I score in two available datasets --------------------------

table(tidy_df$Study, !is.na(tidy_df$PHBR_I))
PHBRI = list()
PHBRI$Van2015 = roc(Response ~ PHBR_I, data = tidy_df, subset = (Study == "Van2015"))
PHBRI$Snyder2017 = roc(Response ~ PHBR_I, data = tidy_df, subset = (Study == "Snyder2017"))

sapply(PHBRI, auc)

p = ggroc(list(VanAllen2015=PHBRI$Van2015, Snyder2017=PHBRI$Snyder2017),
      legacy.axes = TRUE
) +
    labs(color = NULL) + 
    cowplot::theme_cowplot() + 
    theme(legend.position = c(0.6, 0.3)) + 
    scale_color_manual(labels = paste(c("VanAllen", "Snyder"), "AUC =", round(sapply(PHBRI, auc), 2)), values = c("blue", "red"))

ggsave(filename = "PHBR_ROCs.pdf", p)


# Do biomarkers have gender difference? -----------------------------------

genROC2 = function(info, markers, response="Response", direction="<", APSrandom=1:100) {
    ROC_LIST = list()
    ROC_LIST$Male = list()
    ROC_LIST$Female = list()
    # Calc APS random
    APSrandom_list = paste0("APS", APSrandom)
    info_m = info %>% dplyr::filter(gender == "M")
    info_f = info %>% dplyr::filter(gender == "F")
    for (i in c(markers,APSrandom_list)) {
        message("Male:", i)
        ROC_LIST$Male[[i]] = tryCatch(
            {
                roc(info_m[[response]], info_m[[i]], direction="<")
            }, error = function(e) {
                NA
            }
        )
        
    }
    for (i in c(markers,APSrandom_list)) {
        message("Female", i)
        ROC_LIST$Female[[i]] = tryCatch(
            {
                roc(info_f[[response]], info_f[[i]], direction="<")
            }, error = function(e) {
                NA
            }
        )
    }
    ROC_LIST
}

table(tidy_df$Study, tidy_df$gender)
# F  M
# Hugo2016   11 27
# Snyder2017  3 22
# Van2015    30 72

ROC_van2015 = genROC2(tidy_df %>% filter(Study == "Van2015"), colnames(tidy_df)[c(104:123, 125)])
ROC_hugo2016 = genROC2(tidy_df %>% filter(Study == "Hugo2016"), colnames(tidy_df)[c(104:123, 125)])

# summary AUC
van2015_auc_f = lapply(ROC_van2015$Female, function(x) tryCatch(as.numeric(auc(x)), error = function(e) NA))
van2015_auc_m = lapply(ROC_van2015$Male, function(x) tryCatch(as.numeric(auc(x)), error = function(e) NA))

hugo2016_auc_f = lapply(ROC_hugo2016$Female, function(x) tryCatch(as.numeric(auc(x)), error = function(e) NA))
hugo2016_auc_m = lapply(ROC_hugo2016$Male, function(x) tryCatch(as.numeric(auc(x)), error = function(e) NA))

auc_gender_df = list(van2015_auc_f, van2015_auc_m,hugo2016_auc_f, hugo2016_auc_m) %>% 
    dplyr::bind_rows() %>% 
    dplyr::mutate(Group = c("VanAllen2015_Female", "VanAllen2015_Male", "Hugo2016_Female", "Hugo2016_Male"), 
                  APSr = rowMeans(select(., paste0("APS", 1:100)))) %>% 
    tidyr::separate(Group, c("Study", "Gender")) %>% 
    dplyr::select(Study, Gender, dplyr::everything()) %>% 
    dplyr::select(-c(APS1:APS100)) %>% 
    dplyr::rename(TMB = nTMB)

auc_gender_df_long = auc_gender_df %>% 
    tidyr::gather(key = "Biomarker", value = "AUC", -Study, -Gender)

auc_gender_df_long

library(ggpubr)

p2 = ggdotchart(auc_gender_df_long, x = "Biomarker", y ="AUC",
           group = "Gender", color = "Gender", facet.by = "Study", rotate = TRUE,
           ggtheme = theme_bw()) + 
    rotate_x_text(angle = 60, hjust = 1, vjust = 1)
ggsave(filename = "ICB_gender_comparison.pdf", plot = p2)
