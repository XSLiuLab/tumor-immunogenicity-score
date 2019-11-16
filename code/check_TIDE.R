# Check the survival analysis for TIDE
#
# In .Rmd file
# We scale and TIDE and reverse the direction
# by TIDE = -TIDE to do roc (we are sure is right) 
# and survival analysis.
# To make sure the corresponding result 
# is right, here we directly use TIDE values to 
# take a check
library(tidyverse)
library(survival)
library(survminer)

load("report/results/Hugo2016_Info.RData")
load("report/results/VanAllen2015_Info.RData")
load("report/results/Snyder2017_Info.RData")


TIDE_hugo <- read_csv("report/results/TIDE_output_hugo.csv", col_types = cols())
TIDE_Snyder <- read_csv("report/results/TIDE_output_Snyder.csv", col_types = cols())
TIDE_VanAllen <- read_csv("report/results/TIDE_output_VanAllen.csv", col_types = cols())


hugo_os <- read_csv("data/cell2016_OS.csv", col_types = cols())
info_hugo <- full_join(info_hugo, hugo_os, by = "PatientID")


info_hugo2 = info_hugo %>% 
    dplyr::left_join(TIDE_hugo %>% dplyr::select(Patient, TIDE), by = c("PatientID" = "Patient")) %>% 
    dplyr::mutate(TIDE_status = ifelse(TIDE > median(TIDE, na.rm = TRUE), "High", "Low"))

info_VanAllen2 = info_VanAllen %>% 
    dplyr::left_join(TIDE_VanAllen %>% dplyr::select(Patient, TIDE), by = c("patient" = "Patient")) %>% 
    dplyr::mutate(TIDE_status = ifelse(TIDE > median(TIDE, na.rm = TRUE), "High", "Low"))

info_Snyder2 = info_Snyder %>% 
    dplyr::left_join(TIDE_Snyder %>% dplyr::select(Patient, TIDE), by = c("patient_id" = "Patient")) %>% 
    dplyr::mutate(TIDE_status = ifelse(TIDE > median(TIDE, na.rm = TRUE), "High", "Low"))


# Survival analysis
fit_hugo = survfit(Surv(OverallSurvival, VitalStatus=="Dead") ~ TIDE_status, data = info_hugo2)
ggsurvplot(fit_hugo,
           data = info_hugo2, pval = TRUE, fun = "pct",
           xlab = "Time (in days)"
)

fit_van = survfit(Surv(overall_survival, dead) ~ TIDE_status, data = info_VanAllen2)
ggsurvplot(fit_van,
           data = info_VanAllen2, pval = TRUE, fun = "pct",
           xlab = "Time (in days)"
)

fit_sny = survfit(Surv(os, event) ~ TIDE_status, data = info_Snyder2)
ggsurvplot(fit_sny,
           data = info_Snyder2, pval = TRUE, fun = "pct",
           xlab = "Time (in days)"
)
