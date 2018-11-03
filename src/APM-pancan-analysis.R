#--------------------------------------------------------------------
# analysis for APM and TIGS in cancer type level
#--------------------------------------------------------------------
load(file = "data/curated_TCGA_APM_TMB_Clinical.RData")

df_all = df_TMB %>% 
    mutate(nAPM = (APM - min(APM, na.rm = TRUE))/ (max(APM, na.rm = TRUE) - min(APM, na.rm = TRUE)),
           nTMB = TMB_NonsynVariants / 38, 
           TIGS = log(nTMB+1) * nAPM) %>% 
    rename(Event = OS, Time = OS.time)

hist(df_all$nAPM, breaks = 100)
hist(df_all$TMB_NonsynVariants, breaks = 1000, 
     xlab = "TMB before log transformation", main = "TMB Distribution",
     las = 1)
hist(log(df_all$nTMB), breaks=120, xlab = "TMB after log transformation",
     main = "Normalized TMB Distribution", las = 1)
hist(df_all$TIGS, breaks = 100)


df_os = df_all %>% 
    filter(!is.na(Time), !is.na(Event))


df_os %>% filter(!is.na(nAPM)) %>% 
    group_by(Project) %>% summarise(N=n()) %>% arrange(N)

df_os %>% filter(!is.na(TIGS)) %>% 
    group_by(Project) %>% summarise(N=n()) %>% arrange(N)

# calculate APM cox model by project
model_APM = df_os %>% 
    filter(!is.na(nAPM)) %>% 
    group_by(Project) %>% 
    dplyr::do(coxfit = coxph(Surv(time = Time, event = Event) ~ nAPM, data = .)) %>% 
    summarise(Project = Project,
              Coef = summary(coxfit)$conf.int[1],
              Lower = summary(coxfit)$conf.int[3],
              Upper = summary(coxfit)$conf.int[4],
              Pvalue = summary(coxfit)$logtest[3])

# use nTMB or log(nTMB)
model_TMB = df_os %>% 
    filter(!is.na(nTMB)) %>% 
    group_by(Project) %>% 
    dplyr::do(coxfit = coxph(Surv(time = Time, event = Event) ~ log(nTMB), data = .)) %>% 
    summarise(Project = Project,
              Coef = summary(coxfit)$conf.int[1],
              Lower = summary(coxfit)$conf.int[3],
              Upper = summary(coxfit)$conf.int[4],
              Pvalue = summary(coxfit)$logtest[3])

model_TIGS = df_os %>% 
    filter(!is.na(TIGS)) %>% 
    group_by(Project) %>% 
    dplyr::do(coxfit = coxph(Surv(time = Time, event = Event) ~ TIGS, data = .)) %>% 
    summarise(Project = Project,
              Coef = summary(coxfit)$conf.int[1],
              Lower = summary(coxfit)$conf.int[3],
              Upper = summary(coxfit)$conf.int[4],
              Pvalue = summary(coxfit)$logtest[3])

N_APM = df_os %>% filter(!is.na(nAPM)) %>% 
    group_by(Project) %>% summarise(N=n()) 

N_TMB = df_os %>% filter(!is.na(nTMB)) %>% 
    group_by(Project) %>% summarise(N=n()) 

N_TIGS = df_os %>% filter(!is.na(TIGS)) %>% 
    group_by(Project) %>% summarise(N=n()) 


cox_APM = full_join(model_APM, N_APM)
cox_TMB = full_join(model_TMB, N_TMB)
cox_TIGS = full_join(model_TIGS, N_TIGS)

# test
test = df_os %>% 
    filter(!is.na(nAPM)) %>% 
    filter(Project == "DLBC")
coxph(Surv(Time, Event) ~ nAPM, data = test) -> fit

summary(fit) -> fit_sm
fit_sm$coefficients
fit_sm$conf.int

########### Forest plot
options(digits = 2)
forest_APM = rbind(c("Project", NA, NA, NA, "p.value", "No."),
                       cox_APM) %>% as.data.frame()
forest_APM$HR = c("HR", format(as.numeric(forest_APM$Coef[-1]), digits = 2))
forest_APM$Pvalue = c("p.value", format(as.numeric(forest_APM$Pvalue[-1]), digits = 2))

p_load(forestplot)

# forest_APM$Coef = log(as.numeric(forest_APM$Coef))
# forest_APM$Lower = log(as.numeric(forest_APM$Lower))
# forest_APM$Upper = log(as.numeric(forest_APM$Upper))

forestplot(fn.ci_norm = fpDrawCircleCI,
           forest_APM[,c("Project", "N", "HR", "Pvalue")],
               mean = c(NA, log(cox_APM$Coef)), lower = c(NA, log(cox_APM$Lower)), upper = c(NA, log(cox_APM$Upper)),
               is.summary = c(TRUE, rep(FALSE, 32)),
               clip = c(-4, 4), zero = 0,
           col=fpColors(box="royalblue",line="black", summary="royalblue", hrz_lines = "black"),
           vertices = TRUE,
           xticks = c(-4, -2, -1, 0, 1, 2, 4),
           hrzl_lines = list("2" = gpar(lty=1, col = "black")),
           boxsize = 0.5,
           # graph.pos = 3,
           xlab = "log Hazard Ratio")


#----- TMB
forest_TMB = rbind(c("Project", NA, NA, NA, "p.value", "No."),
                    cox_TMB) %>% as.data.frame()
forest_TMB$HR = c("HR", format(as.numeric(forest_TMB$Coef[-1]), digits = 2))
forest_TMB$Pvalue = c("p.value", format(as.numeric(forest_TMB$Pvalue[-1]), digits = 2))


forestplot(fn.ci_norm = fpDrawCircleCI,
           forest_TMB[,c("Project", "N", "HR", "Pvalue")],
           mean = c(NA, log(cox_TMB$Coef)), lower = c(NA, log(cox_TMB$Lower)), upper = c(NA, log(cox_TMB$Upper)),
           is.summary = c(TRUE, rep(FALSE, 32)),
           clip = c(-4, 4), zero = 0,
           col=fpColors(box="royalblue",line="black", summary="royalblue", hrz_lines = "black"),
           vertices = TRUE,
           xticks = c(-4, -2, -1, 0, 1, 2, 4),
           hrzl_lines = list("2" = gpar(lty=1, col = "black")),
           boxsize = 0.5,
           # graph.pos = 3,
           xlab = "log Hazard Ratio")


#---------- TIGS
forest_TIGS = rbind(c("Project", NA, NA, NA, "p.value", "No."),
                   cox_TIGS) %>% as.data.frame()
forest_TIGS$HR = c("HR", format(as.numeric(forest_TIGS$Coef[-1]), digits = 2))
forest_TIGS$Pvalue = c("p.value", format(as.numeric(forest_TIGS$Pvalue[-1]), digits = 2))


forestplot(fn.ci_norm = fpDrawCircleCI,
           forest_TIGS[,c("Project", "N", "HR", "Pvalue")],
           mean = c(NA, log(cox_TIGS$Coef)), lower = c(NA, log(cox_TIGS$Lower)), upper = c(NA, log(cox_TIGS$Upper)),
           is.summary = c(TRUE, rep(FALSE, 32)),
           clip = c(-4, 4), zero = 0,
           col=fpColors(box="royalblue",line="black", summary="royalblue", hrz_lines = "black"),
           vertices = TRUE,
           xticks = c(-4, -2, -1, 0, 1, 2, 4),
           hrzl_lines = list("2" = gpar(lty=1, col = "black")),
           boxsize = 0.5,
           # graph.pos = 3,
           xlab = "log Hazard Ratio")


#------------ TIGS, TMB, APM pancan, sort by value
df_summary = df_os %>% 
    group_by(Project) %>% 
    summarise(medianAPM = median(APM, na.rm = TRUE),
              medianTMB = median(TMB_NonsynVariants, na.rm = TRUE),
              medianTIGS = median(TIGS, na.rm = TRUE),
              medianAPMn = median(nAPM, na.rm = TRUE),
              medianTMBn = log(median(nTMB, na.rm = TRUE) +1 ))

# predict ORR
ORR = function(TIGS) 26.9 * TIGS - 6.6

df_summary = df_summary %>% 
    mutate(TIGS = medianAPMn * medianTMBn) %>% 
    mutate(ORR = ORR(TIGS))

library(scales)
library(ggpubr)
df_os %>% filter(!is.na(Project), !is.na(TIGS)) %>% 
    ggboxplot(x="Project", y="TIGS",  color="Project", add="jitter", xlab = "TCGA Projects", 
              ylab = "TIGS", add.params = list(size=0.6),
              legend = "none") + 
    rotate_x_text(angle = 45) + 
    geom_hline(yintercept = mean(df_os$TIGS, na.rm=TRUE), linetype=2) +
    scale_x_discrete(limits = arrange(df_summary, medianTIGS) %>% .$Project) -> p_tigs
    #stat_compare_means(method = "anova", label.y=4, label.x = 1.5) 

#ggsave("plots/pancan/TIGS_status_pan_cancer.pdf", plot = p_tigs, width = 12, height = 8)
ggsave("/Users/wsx/Documents/APM paper collection/TIGS_status_pan_cancer_correct.pdf", plot = p_tigs, width = 12, height = 8)

df_os %>% filter(!is.na(Project), !is.na(TMB_NonsynVariants)) %>% 
    ggboxplot(x="Project", y="TMB_NonsynVariants",  color="Project", add="jitter", xlab = "TCGA Projects", 
              ylab = "No. of Coding Somatic Nonsynonymous Mutation", add.params = list(size=0.6),
              legend = "none") + 
    rotate_x_text(angle = 45) + 
    geom_hline(yintercept = mean(df_os$TMB_NonsynVariants, na.rm=TRUE), linetype=2) + 
    scale_y_log10(breaks= 10^(-1:4), labels = trans_format("log10", math_format(10^.x))) +
    scale_x_discrete(limits = arrange(df_summary, medianTMB) %>% .$Project) -> p_tmb

ggsave("/Users/wsx/Documents/APM paper collection/TMB_status_pan_cancer_correct.pdf", plot = p_tmb, width = 12, height = 8)


#-- replot APM by rank
df_os %>% filter(!is.na(Project), !is.na(APM)) %>% 
    ggboxplot(x="Project", y="APM",  color="Project", add="jitter", xlab = "TCGA Projects", 
              ylab = "APM Score", add.params = list(size=0.6),
              legend = "none") + 
    rotate_x_text(angle = 45) + 
    geom_hline(yintercept = mean(df_os$APM, na.rm=TRUE), linetype=2) +
    scale_x_discrete(limits = arrange(df_summary, medianAPM) %>% .$Project) -> p_apm

ggsave("/Users/wsx/Documents/APM paper collection/APM_status_pan_cancer_correct.pdf", plot = p_apm, width = 12, height = 8)
# #------ sort Project by APM
# library(tidyverse)
# df_os %>% 
#     filter(!is.na(APM)) %>% 
#     group_by(Project) %>% 
#     summarize(MedianAPM = median(APM), N=n()) %>% 
#     arrange(MedianAPM)

#------------- New analysis
load(file = "results/curated_TCGA_APM_TMB_Clinical.RData")
library(tidyverse)
df_TMB = df_TMB %>% filter(!is.na(IIS), !is.na(TMB_NonsynVariants)) %>% 
    mutate(nTMB = TMB_NonsynVariants / 38, logTMB = log(nTMB + 1))

ggstatsplot::ggscatterstats(
    data = df_TMB, 
    x = nTMB, 
    y = IIS,
    xlab = "No. of Coding Somatic Nonsynonymous Mutation",
    ylab = "IIS Score",
    # title = "Correlation between APM and IIS score in pancancer",
    messages = FALSE, type = "spearman"
) 


ggstatsplot::ggscatterstats(
    data = df_TMB, 
    x = logTMB, 
    y = IIS,
    xlab = "No. of Coding Somatic Nonsynonymous Mutation (log)",
    ylab = "IIS Score",
    # title = "Correlation between APM and IIS score in pancancer",
    messages = FALSE, type = "spearman"
) 

plot_scatter = function(data, x, y, xlab = "Median APM", ylab = "Median IIS", conf.int=TRUE, method="spearman", label.x=-0.5, label.y=0.2, label="Project", ...){
    ggscatter(data, x=x, y=y,
              xlab = xlab, ylab = ylab,
              shape = 21, size = 3, color = "black",
              add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
              conf.int = conf.int,
              cor.coef = TRUE,
              cor.coeff.args = list(method = method, label.x = label.x, label.y=label.y, label.sep = "\n"),
              label=label, repel = TRUE, ...)
}


df_project2 = df_TMB %>%
    group_by(Project) %>% 
    summarise(TMB_median = median(nTMB, na.rm = TRUE), 
              IIS_median = median(IIS, na.rm = TRUE))

mean(df_TMB$nTMB)
plot_scatter(df_project2, x="TMB_median", y="IIS_median",
                         xlab = "Median TMB", ylab = "Median IIS", label.x = 0.2) + 
    geom_hline(yintercept = 0, linetype=2) + geom_vline(xintercept = 4.35, linetype = 2)

ggsave("plots/pancan/Correlation_APM_IIS_acrossTCGA.pdf", plot = p_part2_1,
       width = 5, height = 4)
