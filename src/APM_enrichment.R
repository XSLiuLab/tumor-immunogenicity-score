#----------------------------------------
# APM DEG and Pathway Enrichment Analysis
#----------------------------------------
library(tidyverse)

load("data/TCGA_RNASeq_PanCancer.RData")
load("data/df_combine_gsva_clinical.RData")

tcga_info = df.gsva
rm(df.gsea, df.gsva); gc()

tcga_info = tcga_info %>% 
    select(Project:OS.time, Age, Gender, sample_type, Tumor_stage, APM) %>% 
    filter(!is.na(APM))

projects = unique(tcga_info[, "Project"])
samples = colnames(RNASeq_pancan)[-1]

#-----------------------------------
# Build a workflow to calculate DEGs
#-----------------------------------
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
        
        #--- 保证有少量的病人
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
        }
            
            # if("DESeq2" %in% method){
            #     suppressMessages(library(DESeq2)）
            #                      
            #     colData <- data.frame(row.names=colnames(exprSet), group_list=group_list)
            #     dds <- DESeqDataSetFromMatrix(countData = exprSet,
            #                                   colData = colData,
            #                                   design = ~ group_list)
            #                      
            #     dds <- DESeq(dds)
            #     res <- results(dds, contrast=c("group_list","High","Low"))
            #     resOrdered <- res[order(res$padj),]
            #                      
            #     DEG =as.data.frame(resOrdered)
            #     DESeq2_DEG = na.omit(DEG)
            #                      
            #     nrDEG=DESeq2_DEG[,c(2,6)]
            #     colnames(nrDEG)=c('log2FoldChange','pvalue') 
            #                      
            # }
            # 
            # if("edgeR" %in% method){
            #     suppressMessages(library(edgeR)）
            #                      
            #     d <- DGEList(counts=exprSet,group=factor(group_list))
            #     keep <- rowSums(cpm(d)>1) >= 2
            #                      
            #     d <- d[keep, , keep.lib.sizes=FALSE]
            #     d$samples$lib.size <- colSums(d$counts)
            #     d <- calcNormFactors(d)
            #    
            #     dge=d
            #     design <- model.matrix(~0+factor(group_list))
            #     rownames(design)<-colnames(dge)
            #     colnames(design)<-levels(factor(group_list))
            #     dge=d
            #     dge <- estimateGLMCommonDisp(dge,design)
            #     dge <- estimateGLMTrendedDisp(dge, design)
            #     dge <- estimateGLMTagwiseDisp(dge, design)
            #                      
            #     fit <- glmFit(dge, design)
            #     # https://www.biostars.org/p/110861/
            #     lrt <- glmLRT(fit,  contrast=c(-1,1)) 
            #     nrDEG=topTags(lrt, n=nrow(dge))
            #     nrDEG=as.data.frame(nrDEG)
            #     head(nrDEG)
            #     edgeR_DEG =nrDEG 
            #     nrDEG=edgeR_DEG[,c(1,5)]
            #     colnames(nrDEG)=c('log2FoldChange','pvalue')
            #                      
            # }
            # 
            # if("voom" %in% method){
            #     suppressMessages(library(limma))
            #     
            #     design <- model.matrix(~0+factor(group_list))
            #     colnames(design)=levels(factor(group_list))
            #     rownames(design)=colnames(exprSet)
            #     
            #     dge <- DGEList(counts=exprSet)
            #     dge <- calcNormFactors(dge)
            #     logCPM <- cpm(dge, log=TRUE, prior.count=3)
            #     
            #     v <- voom(dge,design,plot=TRUE, normalize="quantile")
            #     fit <- lmFit(v, design)
            #     
            #     cont.matrix=makeContrasts(contrasts=c('High-Low'),levels = design)
            #     fit2=contrasts.fit(fit,cont.matrix)
            #     fit2=eBayes(fit2)
            #     
            #     tempOutput = topTable(fit2, coef='High-Low', n=Inf)
            #     DEG_limma_voom = na.omit(tempOutput)
            #     
            #     nrDEG=DEG_limma_voom[,c(1,4)]
            #     colnames(nrDEG)=c('log2FoldChange','pvalue')
            #     
            # }
            
        }
        
    }
    
    res = lapply(info_list, doDEG, method = method)
    names(res) = sets
    return(res)
}

DEG_pancan2 = findDEGs(info_df = tcga_info, expr_df = RNASeq_pancan)
save(DEG_pancan2, file = "results/DEG_pancan.RData")

#----------------------------------------------------
# Pathway enrichment analysis using clusterProfiler
#--------------------------------------------------
load(file="results/DEG_pancan.RData")
library(clusterProfiler)
library(tidyverse)
library(openxlsx)

#---------------

##- setting
pvalue = 0.01
adj.pvalue = 0.05

##- process

#------- Reading GSEA genesets files
hallmark = read.gmt("~/biodata/MsigDB/h.all.v6.2.symbols.gmt")
c1 = read.gmt("~/biodata/MsigDB/c1.all.v6.2.symbols.gmt")
c2_kegg = read.gmt("~/biodata/MsigDB/c2.cp.kegg.v6.2.symbols.gmt")
c2_reactome = read.gmt("~/biodata/MsigDB/c2.cp.reactome.v6.2.symbols.gmt")
c3 = read.gmt("~/biodata/MsigDB/c3.all.v6.2.symbols.gmt")
c4 = read.gmt("~/biodata/MsigDB/c4.all.v6.2.symbols.gmt")
c5_mf = read.gmt("~/biodata/MsigDB/c5.mf.v6.2.symbols.gmt")
c5_bp = read.gmt("~/biodata/MsigDB/c5.bp.v6.2.symbols.gmt")
c6 = read.gmt("~/biodata/MsigDB/c6.all.v6.2.symbols.gmt")
c7 = read.gmt("~/biodata/MsigDB/c7.all.v6.2.symbols.gmt")
#-=---------

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


# remove STAD, CHOL, KICH and DLBC
DEG_pancan2 = DEG_pancan2[setdiff(names(DEG_pancan2),c("STAD", "CHOL", "KICH", "DLBC"))]
GSEA_list = Map(goGSEA, DEG_pancan2, names(DEG_pancan2))
save(GSEA_list, file = "results/GSEA_results_rm_notEnriched.RData")

#------- GSEA example plot
# load("results/GSEA_results.RData")
load("results/GSEA_results_rm_notEnriched.RData")
library(clusterProfiler)

gseaplot(GSEA_list$SKCM$hallmark, geneSetID = "HALLMARK_INTERFERON_GAMMA_RESPONSE")

# gseaplot(egmt2, geneSetID = "INTEGRAL_TO_PLASMA_MEMBRANE")
# gseaplot(res$hallmark, geneSetID = "HALLMARK_INTERFERON_GAMMA_RESPONSE")

#------ get main result cloumn
gsea_results = lapply(GSEA_list, function(x) {
    #x@result %>% select()
    lapply(x, function(gsea) {
        gsea@result %>% dplyr::select(ID, setSize, enrichmentScore, NES, pvalue, qvalues)
    })
})

save(gsea_results, file = "results/GSEA_results_notEnriched_small.RData")
rm(GSEA_list); gc()

#----- plot
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

gsea_summary = list()
gsea_summary$hallmark = summariseGSEA(pathway = "hallmark")
gsea_summary$c2_kegg = summariseGSEA(pathway = "c2_kegg")
gsea_summary$c2_reactome = summariseGSEA(pathway = "c2_reactome")
gsea_summary$c5_mf = summariseGSEA(pathway = "c5_mf")
gsea_summary$c5_bp = summariseGSEA(pathway = "c5_bp")
gsea_summary$c6 = summariseGSEA(pathway = "c6")
gsea_summary$c7 = summariseGSEA(pathway = "c7")

save(gsea_summary, file ="results/GSEA_summary_notEnriched.RData")

## plot really
load("data/df_combine_gsva_clinical.RData")

df.gsva %>% 
    filter(!is.na(APM)) %>% 
    group_by(Project) %>% 
    summarize(MedianAPM = median(APM), N=n()) %>% 
    arrange(MedianAPM) -> sortAPM
sortAPM %>% filter(Project != "DLBC") -> sortAPM

library(scales)
gsea_summary$hallmark %>% 
    ggplot(mapping = aes(x = Project, y = ID)) +
    geom_tile(aes(fill = NES)) + 
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    scale_x_discrete(limits = sortAPM$Project) +  theme_classic() + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) 

# df_wide = reshape2::dcast(gsea_summary$hallmark, ID ~ Project, value.var="NES", fill = 0)

library(pheatmap)
# clust = pheatmap(heat_mat, 
#                  color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
#                  breaks = breaksList)

# annotation_col= data.frame(
#     MedianAPMScore = sortAPM$MedianAPM,
#     row.names =  sortAPM$Project
# )
# 
# rownames(annotation_row) = rownames(heat_mat)
# 
# pheatmap(df_wide %>% column_to_rownames(var = "ID"), 
#          color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
#          breaks = breaksList,
#          #annotation_col = annotation_col,
#          fontsize_row = 8, fontsize_col = 8,
#          cellheight = 10, cellwidth = 10,
#          silent = F
#          , filename = "results/test2_hallmark.pdf")

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

gg = plotPathway(gsea_summary$hallmark, silent = F, path = "results/GSEA_results_plot/Hallmark_Enriched.pdf", save = TRUE)
gg = plotPathway(gsea_summary$c2_kegg, silent = F, path = "results/GSEA_results_plot/KEGG_Enriched.pdf", save = TRUE)
gg = plotPathway(gsea_summary$c2_reactome, silent = F, path = "results/GSEA_results_plot/Reactome_Enriched.pdf", save = TRUE)
gg = plotPathway(gsea_summary$c5_mf, silent = F, path = "results/GSEA_results_plot/GO_MolecuFunction_Enriched.pdf", save = TRUE)
gg = plotPathway(gsea_summary$c5_bp, silent = F, path = "results/GSEA_results_plot/GO_BiologicalProcess_Enriched.pdf", save = TRUE)
gg = plotPathway(gsea_summary$c6, silent = F, path = "results/GSEA_results_plot/C6_oncogenicsignatures_Enriched.pdf", save = TRUE)
gg = plotPathway(gsea_summary$c7, silent = F, path = "results/GSEA_results_plot/C7_immunologicsignatures_Enriched.pdf", save = TRUE)

gg


#_-------- rank pheatmap by NES value
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

gg = plotPathway2(gsea_summary$c2_reactome, silent = F, path = "results/GSEA_results_plot/Reactome_Enriched2.pdf", save = TRUE)

#-----------------

# ggsave(filename = "results/test1_hallmark.pdf", width = 10, height = 10, plot = pp)

# gsea_hallmark = purrr::reduce(Map(function(x, y, pathway){
#     if(nrow(x[[pathway]]) == 0){
#         res = NULL
#     }else{
#         res = x[[pathway]]
#         res$Project = y
#     } 
#     res
# }, gsea_results, names(gsea_results), "hallmark"), rbind)

#_---------------------------------------------------

# #------------------------------
# # Pathway enrichment analysis
# # install.packages("pathfindR", dependencies = TRUE)
# library(pathfindR)
# 
# #------ Reactome
# doPathfindR = function(input, project){
#     input = data.frame(Gene.symbol=rownames(input), logFC=input$logFC, adj.P.Val=input$adj.P.Val)
#     Res = run_pathfindR(input, 
#                         output_dir = paste0("/home/wsx/projects/APM/pathfind/APM_enrich_Reactome_", project),
#                         gene_sets = "Reactome", n_processes = 8)
#     Res
# }
# 
# pathEnrich_Reactome = Map(doPathfindR, DEG_pancan2, names(DEG_pancan2))
# 
# save(pathEnrich_Reactome, file = "results/pathEnrich_Reactome.RData")
# 
# 
# Reactome_clustered <- choose_clusters(pathEnrich_Reactome$ACC)
# 
# 
# #----- KEGG
# # doPathfindR = function(input, project){
# #     input = data.frame(Gene.symbol=rownames(input), logFC=input$logFC, adj.P.Val=input$adj.P.Val)
# #     Res = run_pathfindR(input, 
# #                         output_dir = paste0("/home/wsx/projects/APM/pathfind/APM_enrich_KEGG_", project),
# #                         n_processes = 8)
# #     Res
# # }
# # 
# # pathEnrich_KEGG = Map(doPathfindR, DEG_pancan2, names(DEG_pancan2))
# 
