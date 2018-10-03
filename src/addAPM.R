#----------------------------------------------------------------
# Author: Shixiang Wang <wangshx@shanghaitech.edu.cn>
#
# Purpose:
# Add APM data from GEO for some tumor types
# Include Merkel Cell Carcinoma and Small Cell Lung Cancer
#         and Cutaneous Squamous Carcinoma
#----------------------------------------------------------------
library(GEOquery)

if(!dir.exists("data/GEOdata")){
    dir.create("data/GEOdata")
}

geo_dir = "data/GEOdata"

GSE_39612 = getGEO("GSE39612", GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = geo_dir)
GSE_22396 = getGEO("GSE22396", GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = geo_dir)
GSE_36150 = getGEO("GSE36150", GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = geo_dir)
GSE_50451 = getGEO("GSE50451", GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = geo_dir)
GSE_99316 = getGEO("GSE99316", GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = geo_dir)

exprs(GSE_22396$GSE22396_series_matrix.txt.gz) %>% nrow()
pData(GSE_22396$GSE22396_series_matrix.txt.gz)
fData(GSE_22396$GSE22396_series_matrix.txt.gz) %>% nrow()

#---------------------------------------------------------
# process GSE_39612
gset1 = GSE_39612$GSE39612_series_matrix.txt.gz
View(pData(gset1))

gset1_mcc = gset1[, grep("MCC", gset1$title)]
gset1_scc = gset1[, grep("SCC", gset1$title)]

rm(gset1)
# process GSE_22396
gset2 = GSE_22396$GSE22396_series_matrix.txt.gz
View(pData(gset2))

# process GSE_36150
gset3 = GSE_36150$GSE36150_series_matrix.txt.gz
View(pData(gset3))

# process GSE_50451
gset4_1 = GSE_50451$`GSE50451-GPL570_series_matrix.txt.gz`
gset4_2 = GSE_50451$`GSE50451-GPL571_series_matrix.txt.gz`

View(pData(gset4_1))
View(pData(gset4_2))
#gset4_1 all cell_lines, not use it
#
gset4_mcc = gset4_2[, grep("MCC tumor", gset4_2$source_name_ch1)]
gset4_sclc = gset4_2[, grep("SCLC tumor", gset4_2$source_name_ch1)]

rm(gset4_1, gset4_2)
# process GSE_99316
gset5_1 = GSE_99316$`GSE99316-GPL10999_series_matrix.txt.gz`
gset5_2 = GSE_99316$`GSE99316-GPL570_series_matrix.txt.gz`
gset5_3 = GSE_99316$`GSE99316-GPL96_series_matrix.txt.gz`
gset5_4 = GSE_99316$`GSE99316-GPL97_series_matrix.txt.gz`

View(pData(gset5_1))
View(pData(gset5_2))
View(pData(gset5_3))
View(pData(gset5_4))

# only gset5_2 has data we need
gset5_sclc = gset5_2[, grep("SCLC", gset5_2$title)]

rm(gset5_1, gset5_2, gset5_3, gset5_4)

#-------------------------------------------------
# Now summary data
merkel_data = list(gset1_mcc, gset2, gset3, gset4_mcc)
scc_data = list(gset1_scc)
sclc_data = list(gset4_sclc, gset5_sclc)

#--------------------
# apply GSVA method
load("~/biodata/Own/immune_cellType.RData")

# transform exprset 2 tibble
genTibbleList = function(gsetList){
    stopifnot(is.list(gsetList), require("Biobase"), require("tidyverse"), class(gsetList[[1]]) == "ExpressionSet")
    
    res = list()
    i = 1
    for (gset in gsetList){
        eset = exprs(gset)
        # find gene symbol column
        fdata = fData(gset)
        symbol_col = grep("^gene.?symbol", colnames(fdata) , value = TRUE, ignore.case = TRUE)
        
        if(length(symbol_col) == 0) {
            message("Find nothing about gene symbol in fData, try search it...")
            symbol_col2 = grep("^gene_assignment", colnames(fdata) , value = TRUE, ignore.case = TRUE)
            message("find ", symbol_col2)
            
            message("processing...")
            strlist = strsplit(fdata[, symbol_col2], split = " // ")
            rowname = sapply(strlist, function(x) trimws(x[2]))
            rownames(eset) = rowname
            
           # stop("Something wrong with your fData of input List, please check it")
        }
        if(length(symbol_col) > 1) {
            warning("Multiple columns of fData matach gene symbol, only use the first one")
            symbol_col = symbol_col[1]
            rownames(eset) = fdata[, symbol_col]
        }else{
            rownames(eset) = fdata[, symbol_col]
        }
        
        
        
        # remove duplicate rows, keep the one with biggest mean value
        eset %>% 
            as.data.frame() %>% 
            rownames_to_column() %>% 
            mutate(mean_expr = rowMeans(.[, -1], na.rm = TRUE), 
                   rowname = sub("^(\\w+)\\..*", "\\1", rowname)) %>% 
            arrange(rowname, desc(mean_expr)) %>%
                distinct(rowname, .keep_all = TRUE) %>% 
            select(-mean_expr) %>% 
            as.tibble() -> res[[i]]
            
            
        i = i + 1    
    }
    return(res)
}

# apply GSVA method
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
        colnames(res)[1] = 'tsb'
        resList[[i]] = res
        names(resList)[i] = names(ExprMatList)[i]
        i = i + 1
    }   
    return(resList)
}


## scc
tibble.scc = genTibbleList(scc_data)
gsva.scc = applyGSVA(immune_cellType, group_col = "Cell_type", 
                     gene_col = "Symbol", ExprMatList = tibble.scc, method = "gsva")
## sclc
tibble.sclc = genTibbleList(sclc_data)
gsva.sclc = applyGSVA(immune_cellType, group_col = "Cell_type", 
                      gene_col = "Symbol", ExprMatList = tibble.sclc, method = "gsva")

## merkel1
tibble.merkel = genTibbleList(merkel_data)

gsva.merkel = applyGSVA(immune_cellType, group_col = "Cell_type",
                                         gene_col = "Symbol", ExprMatList = tibble.merkel, method = "gsva")

save(gsva.scc, gsva.sclc, gsva.merkel, file = "results/Add_gsva_scc_sclc_merkel.RData")
