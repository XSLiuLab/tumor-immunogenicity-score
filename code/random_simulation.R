setwd("/Users/wsx/Documents/GitHub/tumor-immunogenicity-score/report")

load("results/merged_geneList.RData")
load("results/TCGA_RNASeq_PanCancer.RData")

library(tidyverse)

# keep primary tumor
keep_samples = colnames(RNASeq_pancan)[which(substr(colnames(RNASeq_pancan[-1]), 14, 15) == "01")]
RNASeq_pancan = RNASeq_pancan[, keep_samples]

all_genes = RNASeq_pancan$sample

# keep only genes used to calculate APM score and IIS
keep_types = c("CD8 T cells", "T helper cells", "T cells", "Tcm cells", "Tem cells",
               "Th1 cells", "Th2 cells", "Th17 cells", "Treg cells", "aDC", "B cells",
               "Cytotoxic cells", "DC", "Eosinophils", "iDC", "Macrophages", "Mast cells", 
               "Neutrophils", "NK CD56bright cells",  "NK CD56dim cells", "NK cells", "pDC", "APM")
merged_geneList = dplyr::filter(merged_geneList, Cell_type %in% keep_types)

# GSVA apply function
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
    
    if (length(ExprMatList) > 1) {
        for(expr_mat in ExprMatList){
            if(!inherits(expr_mat, "tbl_df")){
                stop("element of ExprMatList should be tibble!")
            }
            expr_mat = as.data.frame(expr_mat)
            rownames(expr_mat) = expr_mat[, 1]
            expr_mat = expr_mat[, -1] %>% as.matrix()
            
            res = gsva(expr=expr_mat, gset.idx.list=gset_list, method = method, kcdf = kcdf)
            res = as.data.frame(t(res))
            colnames(res)[1] = 'tsb'
            resList[[i]] = res
            names(resList)[i] = names(ExprMatList)[i]
            i = i + 1
        } 
    } else {
        for(expr_mat in ExprMatList){
            if(!inherits(expr_mat, "tbl_df")){
                stop("element of ExprMatList should be tibble!")
            }
            expr_mat = as.data.frame(expr_mat)
            rownames(expr_mat) = expr_mat[, 1]
            expr_mat = expr_mat[, -1] %>% as.matrix()
            
            res = gsva(expr=expr_mat, gset.idx.list=gset_list, method = method, kcdf = kcdf)
            res = as.data.frame(t(res))
            resList[[i]] = res
            names(resList)[i] = names(ExprMatList)[i]
            i = i + 1
        }      
    }
    return(resList)
}


normal_res = applyGSVA(merged_geneList, group_col = "Cell_type", gene_col = "Symbol", ExprMatList = list(RNASeq_pancan), method = "gsva") 

save(normal_res, "results/randomRes/normal_results.RData")

library(foreach)
library(doParallel)
registerDoParallel(cores = 4)

random_geneList = replicate(10,sample(all_genes,nrow(merged_geneList),replace = FALSE))
# use randomly selected genes substitute original genes

random_res = list()

foreach(k=1:ncol(random_geneList)) %dopar% {
    merged_geneList2 = merged_geneList
    merged_geneList2$Symbol = random_geneList[, k]
    
    res = applyGSVA(merged_geneList2, group_col = "Cell_type", gene_col = "Symbol", ExprMatList = list(RNASeq_pancan), method = "gsva")
    res = res[[1]]
    save(res, file = paste0("results/randomRes/random_res", k, ".RData"))
    random_res[[k]] = res
}

save(random_res, file = "results/randomRes/random_res.RData")
