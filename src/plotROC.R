#------------------------------------
# Plot ROC for data analysis results
#------------------------------------
library(pROC)
library(tidyverse)
theme_set(theme_bw())

cell2016 = read_csv("results/Cell2016_results.csv")
cell2017 = read_csv("results/BMS038_results.csv")
science2015 = read_csv("results/science2015_results.csv") # to avoid TIGS <0, TIGS = log(TMB +1) *APM for this dataset



#------------------------------------------
# cell 2016
cell2016_roc1 = roc(cell2016$Response, cell2016$nTMB)
cell2016_roc2 = roc(cell2016$Response, cell2016$APM)
cell2016_roc3 = roc(cell2016$Response, cell2016$TIGS)

ggroc(list(TMB=cell2016_roc1, TIGS=cell2016_roc3), legacy.axes = TRUE) + 
    labs(color = "Predictor") + ggpubr::theme_classic2(base_size = 12, base_family = "Arial")

auc(cell2016_roc1)
auc(cell2016_roc2)
auc(cell2016_roc3)

#---------------------------
# Science 2015

science2015_2 = science2015 %>%
    mutate(Response = case_when(
        group == "response" ~ "R",
        group == "nonresponse" ~"NR",
        TRUE ~ NA_character_
    )) %>% filter(!is.na(Response)) # ,TIGS = log(nTMB) * APM 


science2015_roc1 = roc(science2015_2$Response, science2015_2$nTMB)
science2015_roc2 = roc(science2015_2$Response, science2015_2$APM)
science2015_roc3 = roc(science2015_2$Response, science2015_2$TIGS)

ggroc(list(TMB=science2015_roc1, TIGS=science2015_roc3), legacy.axes = TRUE) + labs(color = "Predictor") + 
    ggpubr::theme_classic2(base_size = 12, base_family = "Arial")

auc(science2015_roc1)
auc(science2015_roc2)
auc(science2015_roc3)

#---------------------------------------
# plos medicine 2017
urothelial2017 = read_csv("results/urothelial2017_results.csv")

urothelial2017 = filter(urothelial2017, !is.na(APM), !is.na(nTMB))

urothelial2017_roc1 = roc(urothelial2017$is_benefit, urothelial2017$nTMB)
urothelial2017_roc2 = roc(urothelial2017$is_benefit, urothelial2017$APM)
urothelial2017_roc3 = roc(urothelial2017$is_benefit, urothelial2017$TIGS)
# urothelial2017_roc4 = roc(urothelial2017$is_benefit, urothelial2017$TIGS2)


ggroc(list(TMB=urothelial2017_roc1, TIGS=urothelial2017_roc3), legacy.axes = TRUE) + 
    labs(color = "Predictor") + ggpubr::theme_classic2(base_size = 12, base_family = "Arial")

auc(urothelial2017_roc1)
auc(urothelial2017_roc2)
auc(urothelial2017_roc3)
auc(urothelial2017_roc4)

# cox
library(survival)
library(survminer)

#----------------------
# Cell 2016
cell2016_os = read_csv("data/cell2016_OS.csv")
cell2016 = full_join(cell2016, cell2016_os, by="PatientID")

coxph(Surv(OverallSurvival, VitalStatus == "Dead") ~ nTMB, data = cell2016)
coxph(Surv(OverallSurvival, VitalStatus == "Dead") ~ TIGS, data = cell2016)
coxph(Surv(OverallSurvival, VitalStatus == "Dead") ~ APM, data = cell2016)

cell2016 %>% mutate(
    TMB_Status = ifelse(nTMB>median(nTMB), "High", "Low"),
    TIGS_Status = ifelse(TIGS>median(TIGS, na.rm = TRUE), "High",
                         ifelse(is.na(TIGS), NA, "Low"))) -> cell2016

fit1 = survfit(Surv(OverallSurvival, VitalStatus == "Dead") ~ TMB_Status, data = cell2016)
ggsurvplot(fit1, data = cell2016, pval = TRUE, fun = "pct", 
           xlab = "Time (in days)")

fit2 = survfit(Surv(OverallSurvival, VitalStatus == "Dead") ~ TIGS_Status, data = cell2016)
ggsurvplot(fit2, data = cell2016, pval = TRUE, fun="pct", xlab = "Time (in days)")

# science 2015
science2015 %>% mutate(TIGS = log(nTMB) * APM,
                    TMB_Status = ifelse(nTMB>median(nTMB, na.rm = TRUE), "High",
                                        ifelse(is.na(nTMB), NA, "Low")),
                    TIGS_Status = ifelse(TIGS>median(TIGS, na.rm = TRUE), "High",
                                         ifelse(is.na(TIGS), NA, "Low")))  %>% 
    filter(group != "long-survival") -> science2015

fit1 = survfit(Surv(overall_survival, dead== 1) ~ TMB_Status, data = science2015)
ggsurvplot(fit1, data = science2015, pval = TRUE, fun="pct", xlab = "Time (in days)")

fit2 = survfit(Surv(overall_survival, dead== 1) ~ TIGS_Status, data = science2015)
ggsurvplot(fit2, data = science2015, pval = TRUE, fun="pct", xlab = "Time (in days)")

#--- plos medicine 2017
urothelial2017 %>% mutate(TMB_Status = ifelse(nTMB>median(nTMB, na.rm = TRUE), "High",
                                           ifelse(is.na(nTMB), NA, "Low")),
                       TIGS_Status = ifelse(TIGS>median(TIGS, na.rm = TRUE), "High",
                                            ifelse(is.na(TIGS), NA, "Low")))  -> urothelial2017

fit1 = survfit(Surv(os, event) ~ TMB_Status, data = urothelial2017)
ggsurvplot(fit1, data = urothelial2017, pval = TRUE, fun="pct", xlab = "Time (in days)")

fit2 = survfit(Surv(os, event) ~ TIGS_Status, data = urothelial2017)
ggsurvplot(fit2, data = urothelial2017, pval = TRUE, fun="pct", xlab = "Time (in days)")
