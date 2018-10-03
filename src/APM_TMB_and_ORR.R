#---------------------------------
# merge APM and TMB data with ORR
#--------------------------------
load("data/df_combine_gsva_clinical.RData")
rm(df.gsea); gc()

require(pacman)
p_load(tidyverse)

APM_tumor.tcga = df.gsva %>% 
    filter(!is.na(APM)) %>% 
    filter(sample_type == "Primary Tumor")

#-----------------------------------------
# Add APM values from mRNA chip
# this data came from result of addAPM.R
#----------------------------------------
load("results/Add_gsva_scc_sclc_merkel.RData")


buildDF = function(x, tumor_type = NULL){
    stopifnot(!is.null(tumor_type))
    APM = x$APM
    res = data.frame(Project = tumor_type, APM = APM, stringsAsFactors = FALSE)
    res
}

lapply(gsva.merkel, FUN = buildDF, tumor_type = "Merkel Cell Carcinoma") %>% purrr::reduce(rbind) -> df_merkel
lapply(gsva.scc, FUN = buildDF, tumor_type = "Cutaneous Squamous Cell Carcinoma") %>% purrr::reduce(rbind) -> df_scc
lapply(gsva.sclc, FUN = buildDF, tumor_type = "Small-Cell Lung Cancer") %>% purrr::reduce(rbind) -> df_sclc

df_list = list(APM_tumor.tcga, df_merkel, df_scc, df_sclc)
purrr::reduce(df_list, bind_rows) -> df_all

#save(df_all, file = "results/df_all.RData")

#---------
# Normalize APM to [0, 1] and summary APM data
df_all %>%
    mutate(N.APM = (APM - min(APM)) / (max(APM) - min(APM))) %>% 
    group_by(Project) %>% 
    summarise(MedianAPM = median(N.APM), Patients_APM = n()) -> sm_APM


#---- data above summary by hand and reload it !!!
rm(list = ls()); gc()

# sm_data = read_csv("data/summary_data_new_20180831.csv")
sm_data = read_csv("data/summary_data_new_20180906.csv")
# sm_data = sm_data %>% filter(!Cancer_Type %in% c("Cutaneous Squamous Cell Carcinoma", "Merkel Cell Carcinoma",
#                                                 "Small-Cell Lung Cancer"))

sm_data = sm_data %>% filter(!is.na(Pool_APM))

lm(Pool_ORR ~ Pool_APM, data = sm_data) %>% summary()
lm(Pool_ORR ~ log(Pool_TMB), data = sm_data)  %>% summary()
lm(Pool_ORR ~ log(Pool_TMB):Pool_APM, data = sm_data)  %>% summary() 
#lm(Pool_ORR ~ log(Pool_TMB)*Pool_APM, data = sm_data)  %>% summary() # this model consider all combination of TMB and APM
               

lm(Pool_ORR ~ Pool_APM, data = sm_data) -> fit1
lm(Pool_ORR ~ log(Pool_TMB), data = sm_data)  -> fit2
lm(Pool_ORR ~ log(Pool_TMB +1):Pool_APM, data = sm_data)  -> fit3

summary(fit1)
summary(fit2)
summary(fit3)

# Following X is a representation of TIGS
sm_data$X = sm_data$Pool_APM * log(sm_data$Pool_TMB + 1)

lm(Pool_ORR ~ X, data = sm_data) -> fit

#----------------------------------------------------------------------------------------
# New ORR can be predicted by this model
# NOTE: this model built on the pooled APM score and pooled TMB value of tumor type
#       if you wanna use this model, please exactly follow this
#       DO NOT use median TIGS value of tumor type calculate from individuals
nd = data.frame(X = c(0.7, 0.986))
predict(fit, newdata = nd, interval = "confidence")
#-----------------------------------------------------------------------------------------


#------- plot
require(ggrepel)
require(scales)

ggplot(sm_data, aes(x=Pool_APM, y=Pool_ORR)) + 
    geom_point(aes(color=Patients_APM, size=Patients_ORR)) + 
    geom_smooth(method="lm", se=T)  + 
    geom_text_repel(aes(label=Cancer_Type), size=3)  +
    labs(x="Median Normalized APM Score ", y="Objective Response Rate (%)",
         size="Objective Response Rate\n(no. of patients evaluated)", 
         color = "APM Score\n(no. of tumor analyzed)") +
    scale_size_continuous(breaks = c(50, 100, 500, 1000)) +
    scale_color_gradientn(colours = RColorBrewer::brewer.pal(5, name="OrRd")[-1],
                          breaks = c(50, 200, 500, 1000)) +
    theme_bw() 


ggplot(sm_data, aes(x=Pool_TMB, y=Pool_ORR)) + 
    geom_point(aes(color=Patients_TMB, size=Patients_ORR)) + 
    geom_smooth(method="lm", se=T)  + 
    geom_text_repel(aes(label=Cancer_Type), size=3)  +
    labs(x="Median No. of Coding Somatic Mutation per MB", y="Objective Response Rate (%)",
         size="Objective Response Rate\n(no. of patients evaluated)", 
         color = "Tumor Mutational Burden\n(no. of tumor analyzed)") +
    scale_x_continuous(trans = log_trans(),
                       breaks = c(2, 10, 20, 30, 40, 50),
                       labels = c(2, 10, 20, 30, 40, 50)) +
    scale_size_continuous(breaks = c(50, 100, 500, 1000)) +
    scale_color_gradientn(colours = RColorBrewer::brewer.pal(5, name="OrRd")[-1],
                          breaks = c(100, 1000, 5000, 10000)) +
    theme_bw()


ggplot(sm_data, aes(x=X, y=Pool_ORR)) + 
    geom_point(aes(size=Patients_ORR)) + 
    geom_smooth(method="lm", se=T)  + 
    geom_text_repel(aes(label=Cancer_Type), size=3)  +
    labs(x="Tumor Immunogenicity Score ", y="Objective Response Rate (%)",
         size="Objective Response Rate\n(no. of patients evaluated)") +
    scale_size_continuous(breaks = c(50, 100, 500, 1000)) +
    theme_bw() 


#------------------------------------------------------------------------------
# To explore the limitation and appropriate coeficient in ORR prediction model
#----------- fit a model, train coeficient 
df = sm_data
df$APS = df$Pool_APM
df$TMB = df$Pool_TMB * 38 # transform unit of per Mb to whole exome
df$ORR = df$Pool_ORR

df = select(df, APS, TMB, ORR)

# explore trends of model ORR ~ APS : log((TMB/x + 1))
calcTrends2 = function(df, a){
    sapply(a, function(x, data){
        summary(lm(ORR ~ APS : log((TMB/x + 1)), data = df) )$r.squared
    }, data = df)
}

# explore trends of model ORR ~ log((TMB/x + 1))
calcTrends3 = function(df, a){
    sapply(a, function(x, data){
        summary(lm(ORR ~ log((TMB/x + 1)), data = df) )$r.squared
    }, data = df)
}

# explore trends of model ORR ~ log((TMB/x)), this study the effect of 1 compared with ORR ~ APS : log((TMB/x + 1))
calcTrends4 = function(df, a){
    sapply(a, function(x, data){
        summary(lm(ORR ~ APS : log((TMB/x)), data = df) )$r.squared
    }, data = df)
}

a = seq(0.1,10000,0.1)

res2 = calcTrends2(df, a)

res3 = calcTrends3(df, a)

res4 = calcTrends4(df, a)

plot(a, res2)
plot(a, res3)
plot(a, res4)

res2_df = data.frame(a = a, Rsquare = res2)
res3_df = data.frame(a = a, Rsquare = res3)
res4_df = data.frame(a = a, Rsquare = res4)

max(res1_df$Rsquare)
res1_df$a[which(res1_df$Rsquare == max(res1_df$Rsquare))]
res2_df$a[which(res2_df$Rsquare == max(res2_df$Rsquare))]

fit2 = lm(ORR ~ APS : log(TMB/61.5 + 1), data = df)
fit3 = lm(ORR ~ APS : log(TMB/38 + 1), data = df)
fit4 = lm(ORR ~ APS : log(TMB/30 + 1), data = df)
fit5 = lm(ORR ~ log(TMB/38 + 1), data = df)
fit6 = lm(ORR ~ TMB, data = df)

# TMB/38 is a good enough model

#----- plot
opar = par(no.readonly = TRUE)
par(mfrow=c(1,2))
plot(res3_df, las = 1, xlab = "a", ylab = "R Square", type = "l")
plot(res3_df[1:500,], las = 1, xlab = "a", ylab = "R Square", type = "l")
par(opar)

opar = par(no.readonly = TRUE)
par(mfrow=c(1,2))
plot(res2_df, las = 1, xlab = "a", ylab = "R Square", type = "l")
plot(res2_df[1:1000,], las = 1, xlab = "a", ylab = "R Square", type = "l")
abline(v = 38, lty = 2, col = "red")
par(opar)


