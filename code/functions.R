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
    
    #library(parallel)
    
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

calc_TisIIs = function(df){
    df %>% 
        mutate(TIS = (`CD8 T cells` + `T helper cells` + `T cells` + `Tcm cells` + `Tem cells` + `Th1 cells` + `Th2 cells` + `Th17 cells` + `Treg cells`) / 9,
               IIS = (`CD8 T cells` + `T helper cells` + `T cells` + `Tcm cells` + `Tem cells` + `Th1 cells` + `Th2 cells` + `Th17 cells` + `Treg cells` + aDC + `B cells` + `Cytotoxic cells` + DC + Eosinophils + iDC + Macrophages + `Mast cells` + Neutrophils + `NK CD56bright cells` + `NK CD56dim cells` + `NK cells` + pDC) / 22) -> df
}

# transform exprset to tibble
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


findDEGs = function(info_df=NULL, expr_df=NULL, col_sample="Tumor_Sample_Barcode", col_group="APM", 
                    col_subset = "Project", method="limma", threshold = 0.25, 
                    save=FALSE, filename=NULL){
    
    # method=c("DESeq2", "limma", "edgeR", "voom")
    stopifnot(!is.null(info_df), !is.null(expr_df))
    stopifnot(threshold >=0 & threshold <=1)
    
    if(!require(data.table)){
        install.packages("data.table", dependencies = TRUE)
    }
    
    info_df =  setDT(info_df[, c(col_sample, col_group, col_subset)])
    colnames(info_df) = c("sample", "groupV", "subset")
    info_df = info_df[sample %in% colnames(expr_df)]
    info_df = info_df[!is.na(subset)]
    all_sets = unique(info_df[, subset])
    
    
    if('DEGroup' %in% colnames(info_df)){
        stop('DEGroup column exists, please rename and re-run.')
    }
    
    # set threshold
    th1 = threshold
    th2 = 1 - th1
    
    info_df[, DEGroup:=ifelse(groupV<quantile(groupV, th1), "Low",
                              ifelse(groupV>quantile(groupV, th2), "High", NA)), by=subset]
    info_df = info_df[!is.na(DEGroup)]
    sets = unique(info_df[, subset])
    may_del = setdiff(all_sets, sets)
    if(length(may_del)!=0) {
        message("Following groups has been filtered out because of threshold setting, you better check")
        print(may_del)
    }
    
    info_list = lapply(sets, function(x) info_df[subset == x])
    names(info_list) = sets
    
    col1 = colnames(expr_df)[1]
    
    
    options(digits = 4)
    
    doDEG = function(input, method=NULL){
        #-- prepare
        exprSet = expr_df[, c(col1, input$sample)]
        exprSet = as.data.frame(exprSet)
        exprSet = na.omit(exprSet)
        input = as.data.frame(input)
        rownames(exprSet) = exprSet[,1]
        exprSet = exprSet[, -1]
        
        
        sample_tb = table(input$DEGroup)
        
        #--- make sure have some samples
        if(!length(sample_tb) < 2  & all(sample_tb >= 5)){
            group_list = input$DEGroup
            
            #-- make sure packages are installed
            if(!all(require("DESeq2"), require("limma"), require("edgeR"))){
                source("https://bioconductor.org/biocLite.R")
                if(!require("DESeq2")) biocLite("DESeq2", dependencies = TRUE)
                if(!require("limma"))  biocLite("limma", dependencies = TRUE)
                if(!require("edgeR"))  biocLite("edgeR", dependencies = TRUE)
            }
            
            if("limma" %in% method){
                suppressMessages(library(limma))    
                design=model.matrix(~0+factor(group_list))
                colnames(design) = c('High', 'Low')
                cont.matrix=makeContrasts('High-Low',levels = design)
                
                fit=lmFit(exprSet,design)
                fit2=contrasts.fit(fit,cont.matrix)
                fit2=eBayes(fit2)
                topTable(fit2, number = Inf, adjust.method ='BH')
            }}
        
    }
    
    res = lapply(info_list, doDEG, method = method)
    names(res) = sets
    return(res)
}


goGSEA = function(DEG, prefix=NULL, pvalue=0.01, adj.pvalue=0.05, destdir="~/projects/APM/results/GSEA_results"){
    library(clusterProfiler)
    library(openxlsx)
    library(tidyverse)
    
    filterDEG = subset(DEG, subset = P.Value < pvalue & adj.P.Val < adj.pvalue)
    filterDEG$SYMBOL = rownames(filterDEG)
    filterDEG = filterDEG %>% 
        arrange(desc(logFC),adj.P.Val)
    geneList = filterDEG$logFC
    names(geneList) = filterDEG$SYMBOL
    
    res = list()
    if (base::exists("hallmark")) res$hallmark = GSEA(geneList, TERM2GENE=hallmark, verbose=FALSE)
    #if (base::exists("c1")) res$c1 = GSEA(geneList, TERM2GENE=c1, verbose=FALSE)
    if (base::exists("c2_kegg")) res$c2_kegg = GSEA(geneList, TERM2GENE=c2_kegg, verbose=FALSE)
    if (base::exists("c2_reactome")) res$c2_reactome = GSEA(geneList, TERM2GENE=c2_reactome, verbose=FALSE)
    #if (base::exists("c3")) res$c3 = GSEA(geneList, TERM2GENE=c3, verbose=FALSE)
    #if (base::exists("c4")) res$c4 = GSEA(geneList, TERM2GENE=c4, verbose=FALSE)
    if (base::exists("c5_mf")) res$c5_mf = GSEA(geneList, TERM2GENE=c5_mf, verbose=FALSE)
    if (base::exists("c5_bp")) res$c5_bp = GSEA(geneList, TERM2GENE=c5_bp, verbose=FALSE)
    if (base::exists("c6")) res$c6 = GSEA(geneList, TERM2GENE=c6, verbose=FALSE)
    if (base::exists("c7")) res$c7 = GSEA(geneList, TERM2GENE=c7, verbose=FALSE)
    
    getResultList = lapply(res, function(x) x@result)
    if(!dir.exists(destdir)) dir.create(destdir)
    outpath = file.path(destdir, prefix)
    write.xlsx(x = getResultList, file = paste0(outpath,".xlsx"))
    return(res)
}

summariseGSEA = function(pathway = NULL){
    purrr::reduce(Map(function(x, y, pathway){
        if(nrow(x[[pathway]]) == 0){
            res = NULL
        }else{
            res = x[[pathway]]
            res$Project = y
        } 
        res
    }, gsea_results, names(gsea_results), pathway), rbind)
}

plotPathway = function(df, save=FALSE, path=NULL, silent=FALSE){
    df_wide = reshape2::dcast(df, ID ~ Project, value.var="NES", fill = 0)
    
    breaksList = seq(-max(df_wide[,-1], na.rm = TRUE), max(df_wide[, -1], na.rm = TRUE), by = 0.01)
    pheatmap(df_wide %>% column_to_rownames(var = "ID"), 
             color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
             breaks = breaksList,
             #annotation_col = annotation_col,
             fontsize_row = 8, fontsize_col = 8, silent = silent,
             cellheight = 10, cellwidth = 10, filename = ifelse(save==FALSE, NA, path))
}

plotPathway2 = function(df, save=FALSE, path=NULL, silent=FALSE){
    df_wide = reshape2::dcast(df, ID ~ Project, value.var="NES", fill = 0)
    df_wide = df_wide[order(rowMeans(df_wide[,-1]), decreasing = TRUE),]
    rownames(df_wide) = NULL
    
    breaksList = seq(-max(df_wide[,-1], na.rm = TRUE), max(df_wide[, -1], na.rm = TRUE), by = 0.01)
    
    pheatmap(df_wide %>% column_to_rownames(var = "ID"), 
             color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
             breaks = breaksList, cluster_rows = FALSE,
             #annotation_col = annotation_col,
             fontsize_row = 8, fontsize_col = 8, silent = silent,
             cellheight = 10, cellwidth = 10, filename = ifelse(save==FALSE, NA, path))
}


buildDF = function(x, tumor_type = NULL){
    stopifnot(!is.null(tumor_type))
    APM = x$APM
    res = data.frame(Project = tumor_type, APM = APM, stringsAsFactors = FALSE)
    res
}

calcTrends = function(df, a){
    sapply(a, function(x, data){
        summary(lm(ORR ~ APS : log((TMB/x + 1)), data = df) )$r.squared
    }, data = df)
}
