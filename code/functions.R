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
        #colnames(res)[1] = 'tsb'
        resList[[i]] = res
        names(resList)[i] = names(ExprMatList)[i]
        i = i + 1
    }   
    return(resList)
}
