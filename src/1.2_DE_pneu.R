
load("../results/sample_count_df_symbol.RData")
################################################################################
# Normlization
library(limma)
library(edgeR)

# All gene count noramlization
S_raw <- as.matrix(dplyr::select(sample_count_df_symbol, -SYMBOL) )
rownames(S_raw) <- sample_count_df_symbol$SYMBOL
dge <- DGEList(counts=S_raw)

################################################################################
comp2groups <- list(nCoV_Heal=c("Healthy", "COVID") )



DEA_list <- list()

for (comp_var in names(comp2groups) ) {
    # comp_var <- names(comp2groups)[[1]]
    comp_i <- comp2groups[[comp_var]]
    print (comp_i)
    subsample_pheno <- dplyr::filter(sample_meta, Types %in% comp_i )
    subsample_pheno$Types <- factor(subsample_pheno$Types, levels=comp_i )
    Expdesign <- model.matrix(~subsample_pheno$Types)
    
    subsample_dge <- dge[,subsample_pheno$NewID]
    
    # DEA
    keep <- filterByExpr(subsample_dge, Expdesign)
    subsample_dge <- subsample_dge[keep,,keep.lib.sizes=FALSE]
    subsample_dge <- calcNormFactors(subsample_dge)
    
    v <- voom(subsample_dge, Expdesign, plot=FALSE, normalize="quantile")
    
    # boxplot(v$E[,nCoV_samp])
    # boxplot(v$E)
    
    Expfit1 <- lmFit(v, Expdesign)
    Expfit2 <- eBayes(Expfit1)
    limma_res <- topTable(Expfit2, coef=tail(colnames(Expdesign), 1), number=Inf) %>% 
        tibble::rownames_to_column() 
    colnames(limma_res)[1] <- comp_var
    
    
    # only protein_coding genes
    limma_DEA <- dplyr::filter(limma_res, adj.P.Val<=0.05, abs(logFC)>=2 )
    readr::write_csv(limma_DEA, paste0("../results/limma/", comp_var, "_limma_res.csv") )
    
    limma_UP <- dplyr::filter(limma_DEA, logFC>0)[[comp_var]]
    limma_DN <- dplyr::filter(limma_DEA, logFC<0)[[comp_var]]
    
    limma_DEG <- limma_DEA[[comp_var]]
    print(length(limma_DEG))
    limma_DE <- v$E[limma_DEG,]
    
    # limma_res include lncRNA, other objs exclude lncRNA
    DEA_list[["E"]][[comp_var]] <- v$E
    DEA_list[["limma_res"]][[comp_var]] <- limma_res
    DEA_list[["limma_DEA"]][[comp_var]] <- limma_DEA
    DEA_list[["limma_DEG"]][[comp_var]] <- limma_DEG
    DEA_list[["limma_UPDN"]][[paste0(comp_var, "_UP")]] <- limma_UP
    DEA_list[["limma_UPDN"]][[paste0(comp_var, "_DN")]] <- limma_DN
    DEA_list[["limma_DE"]][[comp_var]] <- limma_DE
    DEA_list[["subsample_pheno"]][[comp_var]] <- subsample_pheno
    
}

################################################################################

save.image("../results/1.0_DE_pneu.RData")
save(S_raw, DEA_list, sample_meta, file="../results/DEA_pneu_list.RData")

library(ggfortify)
# all gene
autoplot(prcomp(t(v$E), scale=F), 
         data=subsample_pheno, colour = "Types", 
         size = 12, label = T, label.size=4, label.colour="black")


# DEG
autoplot(prcomp(t(v$E[DEA_list[["limma_DEG"]]$nCoV_Heal,]), scale=F), 
         data=subsample_pheno, colour = "Types", 
         size = 12, label = T, label.size=4, label.colour="black")

# corrplot
M <- cor(v$E[DEA_list[["limma_DEG"]]$nCoV_Heal,], method ="p")
corrplot::corrplot(M, order = "hclust", addrect = 2)

# heatmap

nCOV_DE <- DEA_list[["limma_DE"]]$nCoV_Heal
gplots::heatmap.2(nCOV_DE, trace = "none", col="bluered", scale="row", labRow =NA, tracecol=NA,
                  srtCol=60, keysize=1, cexCol=1.3, cexRow=1.1, Colv=TRUE,
                  ylab=paste(nrow(nCOV_DE), "DE Genes"), margins = c(5, 2))






