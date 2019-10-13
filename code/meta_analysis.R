# using metafor from CRAN
# and metawho from https://github.com/ShixiangWang/metawho
load("report/results/unicox.RData")
library(metawho)
library(forestmodel)
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

# Set custom forest panels
custom_panel = function(model = NULL, factor_separate_line = FALSE,
                        headings = list(study = "Study", n = "N", measure = "HR", ci = NULL, p = "Pvalue"),
                        pvalue=NULL) {
    if (inherits(model, "rma")) {
        
        panels <- list(
            forest_panel(width = 0.01),
            forest_panel(
                width = 0.01, display = study, fontface = "bold", heading = headings$study,
                width_group = 1
            ),
            forest_panel(
                width = 0.18, display = stat, parse = TRUE,
                width_group = 1
            ),
            forest_panel(width = 0.03, display = n, hjust = 1, heading = headings$n),
            forest_panel(width = 0.03, item = "vline", hjust = 0.5),
            forest_panel(
                width = 0.45, item = "forest", hjust = 0.5, heading = headings$measure,
                linetype = "dashed", line_x = 0
            ),
            forest_panel(width = 0.03, item = "vline", hjust = 0.5),
            forest_panel(
                width = 0.20,
                display = sprintf("%0.2f (%0.2f, %0.2f)", exp(estimate), exp(conf.low), exp(conf.high)),
                heading = headings$ci,
                display_na = NA
            ),
            forest_panel(width = 0.01, item = "vline", hjust = 0.5),
            forest_panel(
                width = 0.1, display = round(pvalue, digits = 3),
                heading = headings$p,
                display_na = NA
            )
        )
    } else {
        stop("This function only support rma object.")
    }
    panels
}
    
APS_df = cox_APM %>% 
    rename(hr = Coef,
           ci.lb = Lower,
           ci.ub = Upper,
           ni = N,
           subgroup = status) %>% 
    mutate(trial = Project,
           entry = paste(trial, subgroup, sep = "-")) %>% 
    deft_prepare() 
model_APS = APS_df %>% 
    rma(yi = yi, sei = sei, ni = ni,  data = .)

model_APS %>% forestmodel::forest_rma(., 
                            panels = custom_panel(., headings = list(study = "Project", p = "P.value",
                                                                          n = "N", measure = "log Hazard Ratio",
                                                                          ci = "HR (95% CI)"), 
                                                  pvalue = c(p.adjust(APS_df$Pvalue, method = "fdr"), model_APS$pval)),
                            study_labels = APS_df$Project, 
                            limits = c(-4, 5)) -> p_aps
    
TMB_df = cox_TMB %>% 
    rename(hr = Coef,
           ci.lb = Lower,
           ci.ub = Upper,
           ni = N,
           subgroup = status) %>% 
    mutate(trial = Project,
           entry = paste(trial, subgroup, sep = "-")) %>% 
    deft_prepare() 
model_TMB = TMB_df %>% 
    rma(yi = yi, sei = sei, ni = ni,  data = .)
model_TMB %>% 
    forestmodel::forest_rma(., 
                            panels = custom_panel(., headings = list(study = "Project", p = "P.value",
                                                                     n = "N", measure = "log Hazard Ratio",
                                                                     ci = "HR (95% CI)"), 
                                                  pvalue = c(p.adjust(TMB_df$Pvalue, method = "fdr"), model_TMB$pval)),
                            study_labels = TMB_df$Project, 
                            limits = c(-2, 3)) -> p_tmb

TIGS_df = cox_TIGS %>% 
    rename(hr = Coef,
           ci.lb = Lower,
           ci.ub = Upper,
           ni = N,
           subgroup = status) %>% 
    mutate(trial = Project,
           entry = paste(trial, subgroup, sep = "-")) %>% 
    deft_prepare() 
model_TIGS = TIGS_df %>% 
    rma(yi = yi, sei = sei, ni = ni,  data = .)
model_TIGS %>% 
    forestmodel::forest_rma(., 
                            panels = custom_panel(., headings = list(study = "Project", p = "P.value",
                                                                     n = "N", measure = "log Hazard Ratio",
                                                                     ci = "HR (95% CI)"), 
                                                  pvalue = c(p.adjust(TIGS_df$Pvalue, method = "fdr"), model_TIGS$pval)),
                            study_labels = TIGS_df$Project, 
                            limits = c(-4, 7)) -> p_tigs

ggsave("Meta_APS.pdf", plot = p_aps, width = 7, height = 7)
ggsave("Meta_TMB.pdf", plot = p_tmb, width = 7, height = 7)
ggsave("Meta_TIGS.pdf", plot = p_tigs, width = 7, height = 7)

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
    cowplot::theme_cowplot() + ylab("Hazard ratio per unit increase") + xlab(label = "Tumor type APS/TMB/TIGS median") -> p
ggsave("HR_vs_median_APS_TMB_TIGS.pdf", plot = p, width = 8, height = 3)
