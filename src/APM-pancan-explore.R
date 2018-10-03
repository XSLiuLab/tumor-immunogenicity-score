# APM score pan-cancer analysis
p_load(UCSCXenaTools, tidyverse)

TCGA_DIR = "G:/biodata/TCGA/UCSC_Xena/TCGA" # set data directory where TCGA clinical and pancan RNAseq data stored
dir(paste0(TCGA_DIR,"/RNAseq_Norm")) -> RNAseq_filelist
sub("TCGA\\.(.*)\\.sampleMap.*", "\\1", RNAseq_filelist) -> project_code
dir(paste0(TCGA_DIR,"/Clinical")) -> Clinical_filelist
# from <https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations>
TCGA_Study = read_tsv("
                      LAML	Acute Myeloid Leukemia
                      ACC	Adrenocortical carcinoma
                      BLCA	Bladder Urothelial Carcinoma
                      LGG	Brain Lower Grade Glioma
                      BRCA	Breast invasive carcinoma 
                      CESC	Cervical squamous cell carcinoma and endocervical adenocarcinoma
                      CHOL	Cholangiocarcinoma
                      LCML	Chronic Myelogenous Leukemia
                      COAD	Colon adenocarcinoma
                      CNTL	Controls
                      ESCA	Esophageal carcinoma
                      FPPP	FFPE Pilot Phase II
                      GBM	Glioblastoma multiforme
                      HNSC	Head and Neck squamous cell carcinoma
                      KICH	Kidney Chromophobe
                      KIRC	Kidney renal clear cell carcinoma
                      KIRP	Kidney renal papillary cell carcinoma
                      LIHC	Liver hepatocellular carcinoma
                      LUAD	Lung adenocarcinoma
                      LUSC	Lung squamous cell carcinoma
                      DLBC	Lymphoid Neoplasm Diffuse Large B-cell Lymphoma
                      MESO	Mesothelioma
                      MISC	Miscellaneous
                      OV	Ovarian serous cystadenocarcinoma
                      PAAD	Pancreatic adenocarcinoma
                      PCPG	Pheochromocytoma and Paraganglioma
                      PRAD	Prostate adenocarcinoma
                      READ	Rectum adenocarcinoma
                      SARC	Sarcoma
                      SKCM	Skin Cutaneous Melanoma
                      STAD	Stomach adenocarcinoma
                      TGCT	Testicular Germ Cell Tumors
                      THYM	Thymoma
                      THCA	Thyroid carcinoma
                      UCS	Uterine Carcinosarcoma
                      UCEC	Uterine Corpus Endometrial Carcinoma
                      UVM	Uveal Melanoma")
colnames(TCGA_Study) = c("StudyAbbreviation", "StudyName")

intersect(project_code, TCGA_Study$StudyAbbreviation) -> project_code
Clinical_filelist[rowSums(sapply(paste0("__",project_code,"_"), function(x) grepl(x, Clinical_filelist)))>0] -> Clinical_filelist2

setdiff(Clinical_filelist, Clinical_filelist2)
# Therefore, remove TCGA.FPPP.sampleMap__FPPP_clinicalMatrix.gz which has no RNAseq data
# remove other 3 datasets may merge more than one datasets, only keep individual study from TCGA

select_cols = c("sampleID", "OS", "OS.time", "OS.unit", "RFS", "RFS.time", "RFS.unit", 
                "age_at_initial_pathologic_diagnosis", "gender", "tobacco_smoking_history",
                "tobacco_smoking_history_indicator", "sample_type", "pathologic_M",
                "pathologic_N", "pathologic_T", "pathologic_stage")

# read clinical files
Cli_DIR = paste0(TCGA_DIR,"/Clinical")

#----------------------
# Process Clinical data
#----------------------
Clinical_List = XenaPrepare(paste0(Cli_DIR,"/",Clinical_filelist2))
names(Clinical_List)
sub("TCGA\\.(.*)\\.sampleMap.*", "\\1", names(Clinical_List)) -> project_code

TCGA_Clinical = tibble()
# for loop
for (i in 1:length(project_code)){
    clinical = names(Clinical_List)[i]
    project = project_code[i]
    df = Clinical_List[[clinical]]
    col_exist = select_cols %in% colnames(df)
    res = tibble()
    if(!all(col_exist)){
        res = df[, select_cols[col_exist]]
        res[, select_cols[!col_exist]] = NA
    }else{
        res = df[, select_cols]
    }
    res$Project = project
    res %>% select(Project, select_cols) -> res
    TCGA_Clinical = bind_rows(TCGA_Clinical, res)
}

rm(res, df, i, clinical, project, col_exist);gc()

str(TCGA_Clinical)
table(TCGA_Clinical$OS.unit)
table(TCGA_Clinical$RFS.unit)
length(unique(TCGA_Clinical$sampleID)) == length(TCGA_Clinical$sampleID)
table(TCGA_Clinical$gender)
table(TCGA_Clinical$tobacco_smoking_history)
table(TCGA_Clinical$tobacco_smoking_history_indicator)
table(TCGA_Clinical$sample_type)
table(TCGA_Clinical$pathologic_stage)
# all unit are same in days, remove the two column
TCGA_Clinical = TCGA_Clinical %>% select(-c(OS.unit, RFS.unit))

# tidy dataframe
TCGA_Clinical.tidy = TCGA_Clinical %>%
    rename(Age = age_at_initial_pathologic_diagnosis, Gender = gender, 
                         Smoking_history = tobacco_smoking_history,
                         Smoking_indicator = tobacco_smoking_history_indicator,
                         Tumor_Sample_Barcode = sampleID) %>% 
    filter(sample_type %in% c("Solid Tissue Normal", "Primary Tumor", "Metastatic", "Recurrent Tumor",
                              "Primary Blood Derived Cancer - Peripheral Blood")) %>% # Additional - New Primary, Additional Metastatic, FFPE Scrolls total 14 sample removed
    mutate(Gender = case_when(
        Gender == "FEMALE" ~ "Female",
        Gender == "MALE"   ~ "Male",
        TRUE               ~ NA_character_
    ), Tumor_stage = case_when(
        pathologic_stage == "Stage 0" ~ "0",
        pathologic_stage %in% c("Stage I", "Stage IA", "Stage IB") ~ "I",
        pathologic_stage %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC") ~ "II",
        pathologic_stage %in% c("Stage IIIA", "Stage IIIB", "Stage IIIC") ~ "III",
        pathologic_stage %in% c("Stage IV", "Stage IVA", "Stage IVB", "Stage IVC") ~ "IV",
        pathologic_stage == "Stage X" ~ "X",
        TRUE ~ NA_character_)) %>% 
    mutate(Gender = factor(Gender, levels = c("Male", "Female")), 
           Tumor_stage = factor(Tumor_stage, levels = c("0", "I", "II", "III", "IV", "X"))) 

if(!file.exists("results/TCGA_tidy_Clinical.RData")){
    save(TCGA_Clinical.tidy, file = "results/TCGA_tidy_Clinical.RData")
}    


#--------------------------
# Process RNASeq.norm data
#--------------------------
RNAseq_filelist[rowSums(sapply(paste0("TCGA.",project_code,".sampleMap"), function(x) grepl(x, RNAseq_filelist)))>0] -> RNAseq_filelist2

RNASeqDIR.norm = paste0(TCGA_DIR,"/RNAseq_Norm")
RNASeq_List.norm = XenaPrepare(paste0(RNASeqDIR.norm,"/",RNAseq_filelist2))
names(RNASeq_List.norm)
sapply(RNASeq_List.norm, function(x)nrow(x))

names(RNASeq_List.norm) = sub("TCGA\\.(.*)\\.sampleMap.*", "\\1", names(RNASeq_List.norm))

if(!file.exists("results/TCGA_RNASeq_Norm.RData")){
    save(RNASeq_List.norm, file = "results/TCGA_RNASeq_Norm.RData")
}

#--------------------------------------------------------------
# Calclutate Immune cell type Enrichment Score by GSVA package
#--------------------------------------------------------------
APM_genes = read_csv("G:/biodata/GeneList/APM.csv", skip = 1)
APM_genes = APM_genes %>% filter(! Gene_Name %in% c("HLA-E", "HLA-F", "HLA-G", "HLA-H"))

immune_cellType = read_csv("G:/biodata/GeneList/Immune_Cell_type_List.csv", skip=1)
immune_cellType = immune_cellType %>% filter(inRNAseq == "YES")

# add APM genes to immune celltype
immune_cellType = bind_rows(immune_cellType,
          tibble(Cell_type="APM", Symbol=APM_genes$Gene_Name,
                                 Name = APM_genes$Protein_Name, inRNAseq = "YES")) 

applyGSVA = function(group_df, group_col, gene_col, ExprMatList, 
                     method=c("ssgsea", "gsva", "zscore", "plage"),
                              kcdf=c("Gaussian", "Poisson")){
    stopifnot(inherits(group_df, "tbl_df") &
                  inherits(group_col, "character") &
                  inherits(gene_col, "character") &
                  inherits(ExprMatList, "list"))
    if(!require(GSVA)){
        stop("GSVA package need to be installed!")
    }
    
    method = match.arg(method)
    kcdf = match.arg(kcdf)
    
    require(dplyr)
    
    i = 1
    resList = list()
    groups = names(table(group_df[, group_col]))
    gset_list = lapply(groups, function(x){
        group_df[group_df[,group_col] == x, gene_col] %>% unlist %>% as.character()
    })
    
    names(gset_list) = groups

    for(expr_mat in ExprMatList){
        if(!inherits(expr_mat, "tbl_df")){
            stop("element of ExprMatList should be tibble!")
        }
        expr_mat = as.data.frame(expr_mat)
        rownames(expr_mat) = expr_mat[, 1]
        expr_mat = expr_mat[, -1] %>% as.matrix()
        
        res = gsva(expr=expr_mat, gset.idx.list=gset_list, method = method, kcdf = kcdf)
        res = as.data.frame(t(res))
        # modify----
        res = cbind(rownames(res), res)
        colnames(res)[1] = 'tsb'
        #-----
        resList[[i]] = res
        names(resList)[i] = names(ExprMatList)[i]
        i = i + 1
    }   
    return(resList)
}

applyGSVA(immune_cellType, group_col = "Cell_type", gene_col = "Symbol", ExprMatList = RNASeq_List.norm) -> res.ssGSEA
applyGSVA(immune_cellType, group_col = "Cell_type", gene_col = "Symbol", ExprMatList = RNASeq_List.norm, method = "gsva") -> res.GSVA

rm(RNASeq_List.norm); gc()
#------------------------------------
# Pan-cancer level GSVA calculation
#------------------------------------
dir(paste0(TCGA_DIR,"/RNAseq_Pancan")) -> RNAseq_filelist.Pancan
RNAseq_filelist.Pancan[rowSums(sapply(paste0("TCGA.",project_code,".sampleMap"), function(x) grepl(x, RNAseq_filelist.Pancan)))>0] -> RNAseq_filelist.Pancan2

RNASeqDIR.pancan = paste0(TCGA_DIR,"/RNAseq_Pancan")
RNASeq_List.pancan = XenaPrepare(paste0(RNASeqDIR.pancan,"/",RNAseq_filelist.Pancan2))
names(RNASeq_List.pancan)
sapply(RNASeq_List.pancan, function(x)nrow(x))

names(RNASeq_List.pancan) = sub("TCGA\\.(.*)\\.sampleMap.*", "\\1", names(RNASeq_List.pancan))

RNASeq_pancan = purrr::reduce(RNASeq_List.pancan, full_join)
if(!file.exists("results/TCGA_RNASeq_PanCancer.RData")){
    save(RNASeq_pancan, file = "results/TCGA_RNASeq_PanCancer.RData")
}


class(RNASeq_pancan)
rm(RNASeq_List.pancan);gc()

if(!file.exists("results/immune_cellType.RData")){
    save(immune_cellType, file="results/immune_cellType.RData")
}


applyGSVA(immune_cellType, group_col = "Cell_type", gene_col = "Symbol", ExprMatList = list(RNASeq_pancan)) -> res_pancan.ssGSEA


# rm(list=ls());gc()
# load("results/immune_cellType.RData")
# load("results/TCGA_RNASeq_PanCancer.RData")
# This step should run on a machine with big memory
# applyGSVA(immune_cellType, group_col = "Cell_type", gene_col = "Symbol", ExprMatList = list(RNASeq_pancan), method = "gsva") -> res_pancan.GSVA

save(TCGA_Clinical.tidy, res.GSVA, res.ssGSEA, res_pancan.ssGSEA, file = "Rscript/cacheData.RData")

rm(list=ls());gc()
#------------------------------------
# Calculate TIS and IIS
#------------------------
# load("Rscript/cacheData.RData")
# load("Rscript/res_pancan.GSVA.RData")

# save tumor-specific (analysis by tumor type)
# save(res.GSVA, res.ssGSEA, file = "results/GSVA_Result_by_TumorType.RData")
# save(res_pancan.GSVA, res_pancan.ssGSEA, file = "results/GSVA_panCancer_Result.RData")
#rm(res.GSVA, res.ssGSEA);gc()

res_pancan.GSVA = res_pancan.GSVA[[1]]
res_pancan.ssGSEA = res_pancan.ssGSEA[[1]]

plot(density(res_pancan.GSVA$APM))
plot(density(res_pancan.ssGSEA$APM))
# the two method have different distribution

calc_TisIIs = function(df){
    df %>% 
        mutate(TIS = (`CD8 T cells` + `T helper cells` + `T cells` + `Tcm cells` + `Tem cells` + `Th1 cells` + `Th2 cells` + `Th17 cells` + `Treg cells`) / 9,
               IIS = (`CD8 T cells` + `T helper cells` + `T cells` + `Tcm cells` + `Tem cells` + `Th1 cells` + `Th2 cells` + `Th17 cells` + `Treg cells` + aDC + `B cells` + `Cytotoxic cells` + DC + Eosinophils + iDC + Macrophages + `Mast cells` + Neutrophils + `NK CD56bright cells` + `NK CD56dim cells` + `NK cells` + pDC) / 22) -> df
}

# there are miner error in calculation at gsva above
# correct it
colnames(res_pancan.GSVA)[1] <- colnames(res_pancan.ssGSEA)[1] <- 'aDC'
res_pancan.GSVA = rownames_to_column(res_pancan.GSVA, var = "tsb")
res_pancan.ssGSEA = rownames_to_column(res_pancan.ssGSEA, var = "tsb")

gsva.pac = calc_TisIIs(res_pancan.GSVA)
gsea.pac = calc_TisIIs(res_pancan.ssGSEA)

df.gsva = full_join(TCGA_Clinical.tidy, gsva.pac, by=c("Tumor_Sample_Barcode"="tsb"))
df.gsea = full_join(TCGA_Clinical.tidy, gsea.pac, by=c("Tumor_Sample_Barcode"="tsb"))

save(df.gsva, df.gsea, file="results/df_combine_gsva_clinical.RData")

#----------------------------------------------
# exploration of APM scores at pan-cancer level
#----------------------------------------------
load(file="results/df_combine_gsva_clinical.RData")
p_load(survival, survminer, tidyverse, ggstatsplot, ggfortify)

df.gsea$Smoking_indicator = factor(df.gsea$Smoking_indicator,
                                   levels = c("Lifelong Non-smoker", "Current reformed smoker for < or = 15 years", "Current reformed smoker for > 15 years", "Current smoker"), labels = c("Non-smoker", "Reformed-smoker <=15 years", "Reformed-smoker >15 years", "Current-smoker"))

df.gsea.tumor = df.gsea %>% 
    filter(sample_type == "Primary Tumor", !Tumor_stage%in%c("0", "X")) %>%
    mutate(Tumor_stage = factor(Tumor_stage, levels = c("I", "II", "III", "IV")))

df.gsva$Smoking_indicator = factor(df.gsva$Smoking_indicator,
                                   levels = c("Lifelong Non-smoker", "Current reformed smoker for < or = 15 years", "Current reformed smoker for > 15 years", "Current smoker"), labels = c("Non-smoker", "Reformed-smoker <=15 years", "Reformed-smoker >15 years", "Current-smoker"))

df.gsva.tumor = df.gsva %>% 
    filter(sample_type == "Primary Tumor", !Tumor_stage%in%c("0", "X")) %>%
    mutate(Tumor_stage = factor(Tumor_stage, levels = c("I", "II", "III", "IV")))

# pheatmap
# correlation of GSVA score between APM and cell types
df.gsva.heat = df.gsva.tumor %>% 
    select(-c(Tumor_Sample_Barcode:Tumor_stage, `Antigen presenting machinery`)) %>% 
    filter(!is.na(APM)) %>% select(Project, APM, everything())
heat_mat = sapply(unique(df.gsva.heat$Project), function(x){
    mat = filter(df.gsva.heat, Project == x)
    mat = mat[, -1]
    cor_mat = cor(mat, method = "spearman")
    cor_mat[,1]
    })
heat_mat = heat_mat[-1, ]

heat_mat = heat_mat[!rownames(heat_mat) %in% c("TIS", "IIS"), ]


library(pheatmap)
breaksList = seq(-1, 1, by = 0.01)

clust = pheatmap(heat_mat, 
                 color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
                 breaks = breaksList)
annotation_col = data.frame(
    ProjectGroup = factor(paste0('ColCluster',cutree(clust$tree_col,3)))
)
annotation_row = data.frame(
    ImmuneCellGroup = factor(paste0('RowCluster',cutree(clust$tree_row,2)))
)
rownames(annotation_col) = colnames(heat_mat)
rownames(annotation_row) = rownames(heat_mat)
pheatmap(heat_mat, color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
         breaks = breaksList,
         annotation_row = annotation_row,
         annotation_names_row = FALSE, 
         fontsize_row = 8, fontsize_col = 6,
         cellheight = 10, cellwidth = 10, filename = "plots/pancan/correlation_of_GSVA_score_and_cell_types.pdf")



# export


# select method
# according to GSVA paper, gsva and gsea have similar result but gsva is better
#--------------------------------
df = dplyr::select(df.gsva.tumor, Project, APM, TIS, IIS, Gender, Age, Tumor_stage, OS.time, OS, Tumor_Sample_Barcode)
# df = dplyr::select(gsea.tumor, Project, IIS, APM, Gender, Age, Tumor_stage, OS.time, OS)
#--------------------------------------------------------------------------------


#-------------------------------------------
# Part1: APM status across TCGA project
#-------------------------------------------
# p_load(yyplot, grid)
## APM
df %>% filter(!is.na(Project)) %>% 
ggboxplot(x="Project", y="APM",  color="Project", add="jitter", xlab = "TCGA Projects", 
          ylab = "APM Score", add.params = list(size=0.6),
          legend = "none") + 
    rotate_x_text(angle = 45) + 
    geom_hline(yintercept = mean(df$APM, na.rm=TRUE), linetype=2) + 
    stat_compare_means(method = "anova", label.y=1, label.x = 1.5) -> p_part1

ggsave("plots/pancan/part1_APM_status_pan_cancer.pdf", plot = p_part1, width = 12, height = 8)

## TIS
df %>% filter(!is.na(Project)) %>% 
    ggboxplot(x="Project", y="TIS",  color="Project", add="jitter", xlab = "TCGA Projects", 
              ylab = "TIS Score", add.params = list(size=0.6),
              legend = "none") + 
    rotate_x_text(angle = 45) + 
    geom_hline(yintercept = mean(df$TIS, na.rm=TRUE), linetype=2) + 
    stat_compare_means(method = "anova", label.y=0.6, label.x = 1.5) -> p_part1_2

ggsave("plots/pancan/part1_TIS_status_pan_cancer.pdf", plot = p_part1_2, width = 12, height = 8)

## IIS
df %>% filter(!is.na(Project)) %>% 
    ggboxplot(x="Project", y="IIS",  color="Project", add="jitter", xlab = "TCGA Projects", 
              ylab = "IIS Score", add.params = list(size=0.6),
              legend = "none") + 
    rotate_x_text(angle = 45) + 
    geom_hline(yintercept = mean(df$TIS, na.rm=TRUE), linetype=2) + 
    stat_compare_means(method = "anova", label.y=0.6, label.x = 1.5) -> p_part1_3

ggsave("plots/pancan/part1_IIS_status_pan_cancer.pdf", plot = p_part1_3, width = 12, height = 8)

# correlation heatmap of APM and other cell type across tumor type

#----------------------------------------------------------------------
# Part2: relationship of APM and TIL status (multiple methods result)
#----------------------------------------------------------------------
ggstatsplot::ggscatterstats(
    data = df %>% filter(!is.na(APM)), 
    x = APM, 
    y = IIS,
    xlab = "APM Score",
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


df_project = df %>%
    filter(!is.na(APM)) %>% 
    group_by(Project) %>% 
    summarise(APM_median = median(APM), 
              IIS_median = median(IIS))

p_part2_1 = plot_scatter(df_project, x="APM_median", y="IIS_median",
             xlab = "Median APM", ylab = "Median IIS") + 
    geom_hline(yintercept = 0, linetype=2) + geom_vline(xintercept = 0, linetype = 2)

ggsave("plots/pancan/Correlation_APM_IIS_acrossTCGA.pdf", plot = p_part2_1,
       width = 5, height = 4)


# ggscatter(df_project, x="APM_median", y="IIS_median",
#           xlab = "Median APM", ylab = "Median IIS",
#           shape = 21, size = 3, color = "black",
#           add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
#           conf.int = TRUE,
#           cor.coef = TRUE,
#           cor.coeff.args = list(method = "spearman", label.x = -0.5, label.y=0.2, label.sep = "\n"),
#           label="Project", repel = TRUE)


# timer
##############
timer = read_tsv("/Volumes/data/biodata/TCGA/TCGA_sixcell_population_TIMER.txt")
timer$sample %>% length()
timer$sample %>% unique() %>% length()
timer$sample %>% substr(start=1, stop=12) %>% unique() %>% length()

timer %>% filter(as.numeric(substr(sample, 14, 15)) %in% seq(1,9)) %>% 
    mutate(sample = substr(sample, 1, 15)) -> timer2
timer2 %>% arrange(desc(T_cell.CD8)) %>% distinct(sample, .keep_all = TRUE) -> timer2


df_timer = dplyr::left_join(x = df, y = timer2, by = c("Tumor_Sample_Barcode"="sample"))

#---------
# APM
#-------
# df_project = df_timer %>%
#     filter(!is.na(APM)) %>% 
#     filter(T_cell.CD8 > 0.01) %>% 
#     group_by(Project) %>% 
#     summarise_at(vars(c(APM, B_cell:DC)), median, na.rm=TRUE) %>% filter(!is.na(T_cell.CD8))
# 
# plot_scatter(df_project, x="APM", y="B_cell",
#                          xlab = "Median APM", ylab = "Median B cell Fraction") + 
#     geom_vline(xintercept = 0, linetype = 2)
# plot_scatter(df_project, x="APM", y="T_cell.CD4",
#              xlab = "Median APM", ylab = "Median CD4 T cell Fraction") + 
#     geom_vline(xintercept = 0, linetype = 2)
# plot_scatter(df_project, x="APM", y="T_cell.CD8",
#              xlab = "Median APM", ylab = "Median CD8 T cell Fraction", label.y = 0.3) + 
#     geom_vline(xintercept = 0, linetype = 2)
# plot_scatter(df_project, x="APM", y="Neutrophil",
#              xlab = "Median APM", ylab = "Median Neutrophil cell Fraction", label.y = 0.15) + 
#     geom_vline(xintercept = 0, linetype = 2)
# plot_scatter(df_project, x="APM", y="Macrophage",
#              xlab = "Median APM", ylab = "Median Macrophage cell Fraction", label.y = 0.1) + 
#     geom_vline(xintercept = 0, linetype = 2)
# plot_scatter(df_project, x="APM", y="DC",
#              xlab = "Median APM", ylab = "Median  DC cell Fraction", label.y = 0.6) + 
#     geom_vline(xintercept = 0, linetype = 2)
# 
# 
# #-------
# # IIS
# #-------
# df_project = df_timer %>%
#     filter(!is.na(IIS)) %>% 
#     filter(T_cell.CD8 > 0.01) %>% 
#     group_by(Project) %>% 
#     summarise_at(vars(c(IIS, B_cell:DC)), median, na.rm=TRUE) %>% filter(!is.na(T_cell.CD8))
# 
# plot_scatter(df_project, x="IIS", y="B_cell",
#              xlab = "Median IIS", ylab = "Median B cell Fraction",
#              label.x = -0.25, label.y = 0.15) + 
#     geom_vline(xintercept = 0, linetype = 2)
# plot_scatter(df_project, x="IIS", y="T_cell.CD4",label.x = -0.25,
#              xlab = "Median IIS", ylab = "Median CD4 T cell Fraction") + 
#     geom_vline(xintercept = 0, linetype = 2)
# plot_scatter(df_project, x="IIS", y="T_cell.CD8",label.x = -0.25,
#              xlab = "Median IIS", ylab = "Median CD8 T cell Fraction", label.y = 0.3) + 
#     geom_vline(xintercept = 0, linetype = 2)
# plot_scatter(df_project, x="IIS", y="Neutrophil",label.x = -0.25,
#              xlab = "Median IIS", ylab = "Median Neutrophil cell Fraction", label.y = 0.15) + 
#     geom_vline(xintercept = 0, linetype = 2)
# plot_scatter(df_project, x="IIS", y="Macrophage",label.x = -0.25,
#              xlab = "Median IIS", ylab = "Median Macrophage cell Fraction", label.y = 0.1) + 
#     geom_vline(xintercept = 0, linetype = 2)
# plot_scatter(df_project, x="IIS", y="DC",label.x = -0.25,
#              xlab = "Median IIS", ylab = "Median  DC cell Fraction", label.y = 0.6) + 
#     geom_vline(xintercept = 0, linetype = 2)


#----------------------------
# Correlation heatmap with TIL
#-----------------------------
# pheatmap
# correlation of GSVA score between APM/IIS and TIL
p_load(corrplot)
mat_timer_APM = df_timer %>% 
    select(-c(TIS:Tumor_Sample_Barcode)) %>% 
    filter(!is.na(APM) & !is.na(T_cell.CD8)) %>% 
    filter(T_cell.CD8 != 0)
mat_timer_APM.heat = sapply(unique(mat_timer_APM$Project), function(x){
    mat = filter(mat_timer_APM, Project == x)
    mat = mat[, -1]
    cor_mat = cor(mat, method = "spearman")
    cor_mat[,1]
})

mat_timer_APM.heat = t(mat_timer_APM.heat[-1, -22])

col = colorRampPalette(c("blue", "white", "red"))(200)
corrplot(mat_timer_APM.heat, method = "color", tl.col="black", tl.srt = 45,
         col = col)
#-----------------------------
mat_timer_IIS = df_timer %>% 
    select(-c(APM, TIS, Gender:Tumor_Sample_Barcode)) %>% 
    filter(!is.na(IIS) & !is.na(T_cell.CD8)) %>% 
    filter(T_cell.CD8 != 0)
mat_timer_IIS.heat = sapply(unique(mat_timer_IIS$Project), function(x){
    mat = filter(mat_timer_IIS, Project == x)
    mat = mat[, -1]
    cor_mat = cor(mat, method = "spearman")
    cor_mat[,1]
})

mat_timer_IIS.heat = mat_timer_IIS.heat[-1, -22]

p.mat = sapply(unique(mat_timer_IIS$Project), function(x){
    mat = filter(mat_timer_IIS, Project == x)
    mat = mat[, -1]
    tryCatch(expr = {
        p_mat = cor.mtest(mat, method = "spearman", conf.level = .95)
        p_mat[[1]][,1]
    }, error = function(e){
        rep(NA, 7)
    })
    
})

p.mat = p.mat[-1, -22]
corrplot(mat_timer_IIS.heat, method = "color", tl.col="black", tl.srt = 45,
         col = col, p.mat = p.mat, sig.level = 0.05)

breaksList = seq(-1, 1, by = 0.01)
pheatmap(mat_timer_IIS.heat, 
         color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
         breaks = breaksList,
         fontsize_row = 8, fontsize_col = 10)

# tcia
#################
tcia = read_tsv("/Volumes/data/biodata/TCGA/TCGA-all-cellTypeFractionsAll.tsv", col_types = "cccdd")
tcia = unique(tcia)
tcia %>% dplyr::select(-quanTIseq_lsei_TIL10) %>% 
    spread(key = cell_type, value = cibersort_LM22) -> tcia_cibersort

tcia %>% dplyr::select(-cibersort_LM22) %>% 
    spread(key = cell_type, value = quanTIseq_lsei_TIL10) -> tcia_quanTIseq
df2 = df %>% mutate(Tumor_Sample_Barcode = substr(Tumor_Sample_Barcode, 1, 12)) %>% 
    arrange(APM) %>% distinct(Tumor_Sample_Barcode, .keep_all = TRUE)
df_cibersort = dplyr::left_join(x = df2, y = tcia_cibersort, by = c("Tumor_Sample_Barcode"="patientBarcode")) 
df_quanTIseq = dplyr::left_join(x = df2, y = tcia_quanTIseq, by = c("Tumor_Sample_Barcode"="patientBarcode"))

# cibersort
mat_cibersort_IIS = df_cibersort %>% 
    select(-c(APM, TIS, Gender:`CD4 T cells`, `Dendritic cells`, `Natural killer cells`, `Uncharacterized cells`)) %>% 
    filter(!is.na(IIS) & !is.na(`CD8 T cells`)) %>% 
    filter(`CD8 T cells` != 0)
mat_cibersort_IIS.heat = sapply(unique(mat_cibersort_IIS$Project), function(x){
    mat = filter(mat_cibersort_IIS, Project == x)
    mat = mat[, -1]
    cor_mat = cor(mat, method = "spearman")
    cor_mat[,1]
})

mat_cibersort_IIS.heat = mat_cibersort_IIS.heat[-1, ]

p.mat = sapply(unique(mat_cibersort_IIS$Project), function(x){
    mat = filter(mat_cibersort_IIS, Project == x)
    mat = mat[, -1]
    tryCatch(expr = {
        p_mat = cor.mtest(mat, method = "spearman", conf.level = .95)
        p_mat[[1]][,1]
    }, error = function(e){
        rep(NA, 7)
    })
    
})

p.mat = p.mat[-1, ]
corrplot(mat_cibersort_IIS.heat, method = "color", tl.col="black", tl.srt = 45,
         col = col, p.mat = p.mat, sig.level = 0.05)

# quanTIseq
mat_quanTIseq_IIS = df_quanTIseq %>% 
    select(-c(APM, TIS, Gender:disease)) %>% 
    filter(!is.na(IIS) & !is.na(`CD8 T cells`)) %>% 
    filter(`CD8 T cells` != 0)
mat_quanTIseq_IIS.heat = sapply(unique(mat_quanTIseq_IIS$Project), function(x){
    mat = filter(mat_quanTIseq_IIS, Project == x)
    mat = mat[, -1]
    cor_mat = cor(mat, method = "spearman")
    cor_mat[,1]
})

mat_quanTIseq_IIS.heat = mat_quanTIseq_IIS.heat[-1, ]

p.mat = sapply(unique(mat_quanTIseq_IIS$Project), function(x){
    mat = filter(mat_quanTIseq_IIS, Project == x)
    mat = mat[, -1]
    p_mat = cor.mtest(mat, method = "spearman", conf.level = .95)
    p_mat[[1]][,1]
    # tryCatch(expr = {
    #     p_mat = cor.mtest(mat, method = "spearman", conf.level = .95)
    #     p_mat[[1]][,1]
    # }, error = function(e){
    #     rep(NA, 7)
    # })
    
})

p.mat = p.mat[-1, ]
corrplot(mat_quanTIseq_IIS.heat, method = "color", tl.col="black", tl.srt = 45,
         col = col, p.mat = p.mat, sig.level = 0.05, na.label = "NA")

#-----------------------------------------
# Part 3: Correlation between APM and TMB
#------------------------------------------

#---------------------------------
# Using mutation from maftools
# devtools::install_github(repo = "PoisonAlien/TCGAmutations")
require(TCGAmutations)

study_list = names(table(df.gsva$Project))
cohorts = system.file('extdata', 'cohorts.txt', package = 'TCGAmutations')
cohorts = data.table::fread(input = cohorts)

# calculate TMB
lapply(study_list, function(study){
    require(maftools)
        TCGAmutations::tcga_load(study)
        maf = eval(as.symbol(tolower(paste0("TCGA_", study, "_mc3"))))
        maf.silent <- maf@maf.silent
        sample.silent <- maf.silent[,.N, .(Tumor_Sample_Barcode)]
        sample.nonsilent <- getSampleSummary(maf)
        res <- dplyr::full_join(sample.silent, sample.nonsilent, by="Tumor_Sample_Barcode")
        res = res %>% dplyr::mutate(TMB_Total=ifelse(!is.na(N), N+total, total), 
                              TMB_NonsynSNP=Missense_Mutation+Nonsense_Mutation,
                              TMB_NonsynVariants=total) %>% 
            dplyr::select(TMB_Total:TMB_NonsynVariants, Tumor_Sample_Barcode)   
        res
    
}) -> test.mut
names(test.mut) = study_list

# 32 study available
# test.mut[names(test.mut) != "MESO"] -> tcga_tmb
tcga_tmb = test.mut
TCGA_TMB = purrr::reduce(tcga_tmb, bind_rows)
rm(list = grep("tcga_*", ls(), value = TRUE))

save(TCGA_TMB, file = "results/TCGA_TMB.RData")

df_TMB = full_join(df2, TCGA_TMB, by="Tumor_Sample_Barcode")

cor.test(df_TMB$APM, df_TMB$TMB_NonsynVariants, method="spearman")

mat_TMB_APM = df_TMB %>% 
    select(-c(TIS:Tumor_Sample_Barcode, TMB_Total, TMB_NonsynSNP)) %>% 
    filter(!is.na(APM) & !is.na(TMB_NonsynVariants)) 

mat_TMB_APM.heat = sapply(unique(mat_TMB_APM$Project), function(x){
    mat = filter(mat_TMB_APM, Project == x)
    mat = mat[, -1]
    cor_mat = cor(mat, method = "spearman")
    cor_mat[,1]
})

mat_TMB_APM.heat = mat_TMB_APM.heat[-1,]

p.mat = sapply(unique(mat_TMB_APM$Project), function(x){
    mat = filter(mat_TMB_APM, Project == x)
    mat = mat[, -1]
    tryCatch(expr = {
        p_mat = cor.mtest(mat, method = "spearman", conf.level = .95)
        p_mat[[1]][,1]
    }, error = function(e){
        rep(NA, 2)
    })
    
})

p.mat = p.mat[-1, ]
# corrplot(matrix(mat_TMB_APM.heat), method = "color", tl.col="black", tl.srt = 45,
#          col = col, p.mat = p.mat, sig.level = 0.05)

N_project = mat_TMB_APM %>% group_by(Project) %>% summarise(N_samples = n())

TMB_cor.mat = inner_join(N_project, 
           cbind(mat_TMB_APM.heat, p.mat) %>% as.data.frame() %>% rownames_to_column(var = "Project"),
           by="Project")
knitr::kable(TMB_cor.mat)
TMB_cor.mat %>% rename(Corr=mat_TMB_APM.heat, pval=p.mat) %>% 
    arrange(desc(Corr))

#----------------------------
# TMB and IIS
#----------------------------
mat_TMB_IIS = df_TMB %>% 
    select(-c(APM, TIS, Gender:Tumor_Sample_Barcode, TMB_Total, TMB_NonsynSNP)) %>% 
    filter(!is.na(IIS) & !is.na(TMB_NonsynVariants)) 

mat_TMB_IIS.heat = sapply(unique(mat_TMB_APM$Project), function(x){
    mat = filter(mat_TMB_IIS, Project == x)
    mat = mat[, -1]
    cor_mat = cor(mat, method = "spearman")
    cor_mat[,1]
})

mat_TMB_IIS.heat = mat_TMB_IIS.heat[-1,]

p.mat = sapply(unique(mat_TMB_IIS$Project), function(x){
    mat = filter(mat_TMB_IIS, Project == x)
    mat = mat[, -1]
    tryCatch(expr = {
        p_mat = cor.mtest(mat, method = "spearman", conf.level = .95)
        p_mat[[1]][,1]
    }, error = function(e){
        rep(NA, 2)
    })
    
})

p.mat = p.mat[-1, ]
# corrplot(matrix(mat_TMB_APM.heat), method = "color", tl.col="black", tl.srt = 45,
#          col = col, p.mat = p.mat, sig.level = 0.05)

N_project = mat_TMB_IIS %>% group_by(Project) %>% summarise(N_samples = n())

IIS_cor.mat = inner_join(N_project, 
                         cbind(mat_TMB_IIS.heat, p.mat) %>% as.data.frame() %>% rownames_to_column(var = "Project"),
                         by="Project")
knitr::kable(IIS_cor.mat)
IIS_cor.mat %>% rename(Corr=mat_TMB_IIS.heat, pval=p.mat) %>% 
    arrange(desc(Corr))

#=============================================
#-----------------------------------------
# Part 4: Correlation between APM and CNV
#------------------------------------------
# read CNA data
CNA_level = read_csv("/Volumes/data/biodata/TCGA/TCGA_SCNA_Levels.csv", comment = "#")
TCGA_vars = full_join(df_TMB, CNA_level, by = c("Tumor_Sample_Barcode"="TumorSample"))

TCGA_vars %>% 
    filter(!is.na(SCNA.Level.normalized.by.size) & !is.na(APM)) %>% 
    filter(TumorType != "SKCM") %>% 
    group_by(TumorType) %>% 
    summarise(N_samples = n(),
              Corr_Chrom = cor(APM, Chrom.SCNA.Level, method = "spearman"),
              Corr_Chrom.p = cor.test(APM, Chrom.SCNA.Level, method = "spearman")$p.value,
              Corr_Arm = cor(APM, Arm.SCNA.Level, method = "spearman"),
              Corr_Arm.p = cor.test(APM, Arm.SCNA.Level, method = "spearman")$p.value,
              Corr_Focal = cor(APM, Focal.SCNA.Level, method = "spearman"),
              Corr_Focal.p = cor.test(APM, Focal.SCNA.Level, method = "spearman")$p.value,
              Corr_SCNA = cor(APM, SCNA.Level, method = "spearman"),
              Corr_SCNA.p = cor.test(APM, SCNA.Level, method = "spearman")$p.value,
              Corr_nSCNA = cor(APM, SCNA.Level.normalized.by.size, method = "spearman"),
              Corr_nSCNA.p = cor.test(APM, SCNA.Level.normalized.by.size, method = "spearman")$p.value) -> vars_mat
    
knitr::kable(vars_mat)

#-------------------------------
# Part 5: Sex difference APM score ???
#-------------------------------
df %>% filter(!is.na(Project) & !is.na(Gender)) %>% 
    ggboxplot(x="Project", y="APM",  color="Gender", add="jitter", xlab = "TCGA Projects", 
              ylab = "APM Score", add.params = list(size=0.6),
              legend = "top") + 
    rotate_x_text(angle = 45) + 
    geom_hline(yintercept = mean(df$APM, na.rm=TRUE), linetype=2) -> p_part5_1


df %>% filter(!is.na(Project) & !is.na(Gender)) %>% 
    ggboxplot(x="Gender", y="APM",  color="Gender", add="jitter", xlab = "", 
              ylab = "APM Score", add.params = list(size=0.6),
              facet.by = "Project", palette = c("blue", "red"),
              legend = "top") + 
    rotate_x_text(angle = 45) + 
    geom_hline(yintercept = 0, linetype=2) + 
    stat_compare_means(label = "p.format") -> p_part5_2


ggsave("plots/pancan/part5_APM_gender_across_TCGA.pdf", plot = p_part5_1, width = 12, height = 8)
ggsave("plots/pancan/part5_gender_compare_APM_across_TCGA.pdf", plot = p_part5_2, width = 8, height = 12)


#----------------------------------------------------------------
# Part 6 APM and patients prognosis
#-----------------------------------------------------------------

save(df_TMB, file = "results/curated_TCGA_APM_TMB_Clinical.RData")
# gsva.tumor = df.gsva.tumor
# gsva.tumor$Tumor_Sample_Barcode = substr(gsva.tumor$Tumor_Sample_Barcode, 1, 12)
# gsea.tumor = df.gsea.tumor
# gsea.tumor$Tumor_Sample_Barcode = substr(gsea.tumor$Tumor_Sample_Barcode, 1, 12)
# 
# gsva.tumor = full_join(gsva.tumor, TCGA_vars, by=c("Tumor_Sample_Barcode"))
# gsea.tumor = full_join(gsea.tumor, TCGA_vars, by = c("Tumor_Sample_Barcode"))
 
load(file = "results/curated_TCGA_APM_TMB_Clinical.RData")

# fit = coxph(Surv(OS.time, OS) ~ Gender + Age + Tumor_stage + TMB_NonsynVariants + SCNA.Level.normalized.by.size, data=na.omit(TCGA_vars[, c("IIS", "APM", "Gender", "Age", "Tumor_stage", "OS.time", "OS", "TMB_NonsynVariants", "SCNA.Level.normalized.by.size", "CellCycle.Signature.Score")]))
# fit = coxph(Surv(OS.time, OS) ~ SCNA.Level.normalized.by.size,
#             data=df)

##### re-compute above
# TCGA_vars %>% 
#     filter(!is.na(TMB_NonsynVariants) & !is.na(IIS)) %>% 
#     group_by(Project) %>% 
#     summarise(N=n(), Corr = cor(TMB_NonsynVariants, IIS)) -> cor.IIS.TMB
# 
# TCGA_vars %>% 
#     filter(!is.na(IIS)) %>% 
#     group_by(Project) %>% 
#     summarise(N=n(), Corr = cor(APM, IIS)) -> cor.IIS.APM

# df_nona = TCGA_vars %>% 
#     dplyr::select(Project:Tumor_Sample_Barcode, TMB_NonsynVariants) %>% 
#     na.omit()
df_nona = df_TMB %>% 
    dplyr::mutate(TMB = TMB_Total / 38) %>% 
    na.omit()


    # dplyr::select(-c(Project, TumorType, age, `gender(1=male)`, stage, histologicaltype, ER.receptor.status, TotNofMutations.in.exons, TP53Mutation)) %>%
    # na.omit()

fit = coxph(Surv(OS.time, OS) ~  APM + Gender + Age + Tumor_stage + TMB,
            data=df_nona)

# summary(fit)
ggforest(fit,main = "", fontsize = 1)
 
# fit = coxph(Surv(OS.time, OS) ~  APM + Gender + Age + Tumor_stage,
#             data=df.gsva.tumor %>% filter(!is.na(APM)))
# ggforest(fit)
# 
# 
# fit = coxph(Surv(OS.time, OS) ~  APM,
#             data=df.gsva.tumor %>% filter(!is.na(APM)))
# ggforest(fit)


df.gsva.tumor.APM = df.gsva.tumor %>% filter(!is.na(APM))
res.cut = surv_cutpoint(data=df.gsva.tumor.APM, time = "OS.time", event = "OS", variables = "APM")
res.cat = surv_categorize(res.cut)
fit = survfit(Surv(OS.time, OS) ~ APM, data = res.cat)
ggsurvplot(fit, data =res.cat, pval = TRUE, fun = "pct")

df.gsva.tumor.APM$APM_Status = ifelse(df.gsva.tumor.APM$APM > 0, "High", "Low")
# df.gsva.tumor.APM = subset(df.gsva.tumor.APM, OS.time <= 3650)
fit = survfit(Surv(OS.time, OS) ~ APM_Status, data = df.gsva.tumor.APM)
ggsurvplot(fit, data =df.gsva.tumor.APM,  pval = TRUE, fun = "pct", 
           xlab = "Time (in days)")
fit

## quartile 
df.gsva.tumor.APM = df.gsva.tumor.APM %>% 
    mutate(APM_quartile = case_when(
        APM <=  quantile(APM)[2] ~ "<= 25%",
        APM > quantile(APM)[2] & APM <= quantile(APM)[3] ~ "> 25% and <= 50%",
        APM > quantile(APM)[3] & APM <= quantile(APM)[4] ~ "> 50% and <= 75%",
        APM > quantile(APM)[4] ~ "> 75%",
       TRUE ~ NA_character_ 
    ),
    APM_quartile = factor(APM_quartile, levels = c("<= 25%", "> 25% and <= 50%", "> 50% and <= 75%", "> 75%")))
fit = survfit(Surv(OS.time, OS) ~ APM_quartile, data = df.gsva.tumor.APM %>% filter(OS.time <=3650))
ggsurvplot(fit, data =df.gsva.tumor.APM %>% filter(OS.time <=3650), risk.table = TRUE, pval = TRUE)

## quantile 10
df.gsva.tumor.APM = df.gsva.tumor.APM %>% 
    mutate(APM_quan = case_when(
        APM <=  quantile(APM, probs = seq(0,1,0.1))[2] ~ "<= 10%",
        APM > quantile(APM, probs = seq(0,1,0.1))[2] & APM <= quantile(APM, , probs = seq(0,1,0.1))[3] ~ "> 10% and <= 20%",
        APM > quantile(APM, probs = seq(0,1,0.1))[3] & APM <= quantile(APM, , probs = seq(0,1,0.1))[4] ~ "> 20% and <= 30%",
        APM > quantile(APM, probs = seq(0,1,0.1))[4] & APM <= quantile(APM, , probs = seq(0,1,0.1))[5] ~ "> 30% and <= 40%",
        APM > quantile(APM, probs = seq(0,1,0.1))[5] & APM <= quantile(APM, , probs = seq(0,1,0.1))[6] ~ "> 40% and <= 50%",
        APM > quantile(APM, probs = seq(0,1,0.1))[6] & APM <= quantile(APM, , probs = seq(0,1,0.1))[7] ~ "> 50% and <= 60%",
        APM > quantile(APM, probs = seq(0,1,0.1))[7] & APM <= quantile(APM, , probs = seq(0,1,0.1))[8] ~ "> 60% and <= 70%",
        APM > quantile(APM, probs = seq(0,1,0.1))[8] & APM <= quantile(APM, , probs = seq(0,1,0.1))[9] ~ "> 70% and <= 80%",
        APM > quantile(APM, probs = seq(0,1,0.1))[9] & APM <= quantile(APM, , probs = seq(0,1,0.1))[10] ~ "> 80% and <= 90%",
        APM > quantile(APM, probs = seq(0,1,0.1))[10] ~ "> 90%",
        TRUE ~ NA_character_ 
    ))
fit = survfit(Surv(OS.time, OS) ~ APM_quan, data = df.gsva.tumor.APM %>% filter(OS.time <=3650) )
ggsurvplot(fit, data =df.gsva.tumor.APM %>% filter(OS.time <=3650), risk.table = TRUE, pval = TRUE)





length(which(!is.na(df.gsea.tumor$Tumor_stage)))
fit = coxph(Surv(OS.time, OS) ~  IIS + Gender + Age + Tumor_stage,
            data=df.gsea.tumor %>% filter(!is.na(APM)))
ggforest(fit)



# plot forest plot by project
group_forest = function(data){
    projects = unique(data$Project)
    pdf(file = "plots/pancan/part6_forestPlot_by_projects.pdf", width = 7)
    for (i in 1:length(projects)){
        fit_data = dplyr::filter(data, Project == projects[i])
        fit = coxph(Surv(OS.time, OS) ~  APM + Gender + Age + Tumor_stage,
                    data=fit_data)
        p = ggforest(fit, main = paste0("Hazard ratio of ", projects[i]), noDigits = 3, data = fit_data)
        print(p, newpage = TRUE)
    }
    dev.off()
}
group_forest(df.gsea.tumor %>% filter(!is.na(APM) & !is.na(Gender) & !is.na(Age) & !is.na(Tumor_stage)))


 
require(MASS)
fit = lm(OS.time ~ IIS + APM + Gender + Age + Tumor_stage + TMB_Total + TMB_NonsynVariants + TMB_NonsynSNP + Immune.Signature.Score + CellCycle.Signature.Score + SCNA.Level + Chrom.SCNA.Level + Arm.SCNA.Level + Focal.SCNA.Level + Chrom.Arm.SCNA.Level + SCNA.Level.normalized.by.size, data = TCGA_vars)
stepAIC(fit)
fit = lm(OS.time ~ IIS + APM + Gender + Age + Tumor_stage, data = df.gsva.tumor)
stepAIC(fit)
plot(na.omit(df.gsva.tumor[, c("IIS", "APM", "Gender", "Age", "Tumor_stage", "OS.time")])$OS.time, fitted(fit))


