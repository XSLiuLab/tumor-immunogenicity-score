# using metafor from CRAN
# and metawho from https://github.com/ShixiangWang/metawho
load("report/results/unicox.RData")
library(metawho)
library(tidyverse)

df_summary = df_summary %>% 
    mutate_at(vars(medianAPMn, medianTMBn, medianTIGS), 
              list(status = ~ifelse(. > median(.), "High", "Low")))
    
cox_APM = cox_APM %>% 
    left_join(df_summary %>% select(Project, medianAPMn_status, medianAPM)) %>% 
    rename(status = medianAPMn_status) %>% 
    arrange(medianAPM)

cox_TMB = cox_TMB %>% 
    left_join(df_summary %>% select(Project, medianTMBn_status, medianTMBn)) %>% 
    rename(status = medianTMBn_status) %>% 
    arrange(medianTMBn) %>% 
    mutate(medianTMB = exp(medianTMBn) - 1)

cox_TIGS = cox_TIGS %>% 
    left_join(df_summary %>% select(Project, medianTIGS_status, medianTIGS)) %>% 
    rename(status = medianTIGS_status) %>% 
    arrange(medianTIGS)

# cox_APM %>% 
#     rename(hr = Coef,
#            ci.lb = Lower,
#            ci.ub = Upper,
#            ni = N,
#            subgroup = status) %>% 
#     mutate(trial = Project,
#         entry = paste(trial, subgroup, sep = "-")) %>% 
#     deft_prepare() -> tt
# 
# model = rma(yi = yi, sei = sei, ni = ni,  data = tt %>% filter(subgroup == "Low"))
# forestmodel::forest_rma(model, 
#                         panels = metawho:::deft_panel(model, 
#                                                       headings = list(study = "Study", 
#                                                            n = "N", measure = "log Hazard Ratio", ci = NULL)),
#                         study_labels = tt$Project[tt$subgroup == "Low"], 
#                         limits = c(-4, 5))
# 
# model = rma(yi = yi, sei = sei, ni = ni,  data = tt %>% filter(subgroup == "High"))
# forestmodel::forest_rma(model, 
#                         panels = metawho:::deft_panel(model, 
#                                                       headings = list(study = "Study", 
#                                                                       n = "N", measure = "log Hazard Ratio", ci = NULL)),
#                         study_labels = tt$Project[tt$subgroup == "High"], 
#                         limits = c(-4, 5))

APS_df = cox_APM %>% 
    rename(hr = Coef,
           ci.lb = Lower,
           ci.ub = Upper,
           ni = N,
           subgroup = status) %>% 
    mutate(trial = Project,
           entry = paste(trial, subgroup, sep = "-")) %>% 
    deft_prepare() 
APS_df %>% 
    rma(yi = yi, sei = sei, ni = ni,  data = .) %>% 
    forestmodel::forest_rma(., 
                            panels = metawho:::deft_panel(., 
                                                          headings = list(study = "Project", 
                                                                          n = "N", measure = "log Hazard Ratio", ci = "HR (95% CI)")),
                            study_labels = APS_df$Project, 
                            limits = c(-4, 5))
    
TMB_df = cox_TMB %>% 
    rename(hr = Coef,
           ci.lb = Lower,
           ci.ub = Upper,
           ni = N,
           subgroup = status) %>% 
    mutate(trial = Project,
           entry = paste(trial, subgroup, sep = "-")) %>% 
    deft_prepare() 
TMB_df %>% 
    rma(yi = yi, sei = sei, ni = ni,  data = .) %>% 
    forestmodel::forest_rma(., 
                            panels = metawho:::deft_panel(., 
                                                          headings = list(study = "Project", 
                                                                          n = "N", measure = "log Hazard Ratio", ci = "HR (95% CI)")),
                            study_labels = TMB_df$Project, 
                            limits = c(-2, 3))

TIGS_df = cox_TIGS %>% 
    rename(hr = Coef,
           ci.lb = Lower,
           ci.ub = Upper,
           ni = N,
           subgroup = status) %>% 
    mutate(trial = Project,
           entry = paste(trial, subgroup, sep = "-")) %>% 
    deft_prepare() 
TIGS_df %>% 
    rma(yi = yi, sei = sei, ni = ni,  data = .) %>% 
    forestmodel::forest_rma(., 
                            panels = metawho:::deft_panel(., 
                                                          headings = list(study = "Project", 
                                                                          n = "N", measure = "log Hazard Ratio", ci = "HR (95% CI)")),
                            study_labels = TIGS_df$Project, 
                            limits = c(-4, 5))

# Plot correlation between median value and HR
plot_df = dplyr::bind_rows(
    cox_APM %>% select(Project, Coef, medianAPM, N) %>% 
        rename(Median = medianAPM) %>% mutate(type = "APS"),
    cox_TMB %>% select(Project, Coef, medianTMB, N) %>% 
        rename(Median = medianTMB) %>% mutate(type = "TMB"),
    cox_TIGS %>% select(Project, Coef, medianTIGS, N) %>% 
        rename(Median = medianTIGS) %>% mutate(type = "TIGS")
) %>% rename(HR = Coef) %>% 
    mutate(type = factor(type, levels = c("APS", "TMB", "TIGS")))

library(ggpubr)
#ggpubr::ggscatter(plot_df, x = "Median", y = "HR", facet.by = "type")
ggplot(plot_df %>% filter(HR < 10), aes(x = Median, y = HR)) +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_wrap(~type, scales = "free") +
    stat_cor(method = "pearson") + 
    geom_hline(yintercept = 1, linetype = 2) + 
    cowplot::theme_cowplot() + ylab("Hazard ratio") + xlab(label = NULL) -> p
ggsave("HR_vs_median_APS_TMB_TIGS.pdf", plot = p)
