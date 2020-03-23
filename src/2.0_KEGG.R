
load("../results/DEA_pneu_list.RData")

symbol2entrezID <- function(gene_symbols) {
    symbol_ENTREZID <- clusterProfiler::bitr(gene_symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    return(symbol_ENTREZID$ENTREZID)
}

c1 <- c("nCoV_Heal_UP", "nCoV_Heal_DN", "Vir_Heal_UP", "Vir_Heal_DN", 
        "Others_Heal_UP", "Others_Heal_DN")

################################################################################
# Pathway & GO analysis of DEG
library(clusterProfiler)
genesEntrezID_3g <- sapply(DEA_list[["limma_UPDN"]][c1], symbol2entrezID )
sapply(genesEntrezID_3g, length)

genesEntrezID_3g_KEGG <- compareCluster(genesEntrezID_3g, fun='enrichKEGG')
dotplot(genesEntrezID_3g_KEGG, showCategory=20)

library(ReactomePA)
genesEntrezID_3g_Reactome <- compareCluster(genesEntrezID_3g, fun='enrichPathway')
dotplot(genesEntrezID_3g_Reactome, showCategory=10)


genesEntrezID_3g_BP <- compareCluster(genesEntrezID_3g, fun='enrichGO', OrgDb='org.Hs.eg.db', ont="BP")
dotplot(genesEntrezID_3g_BP, showCategory=10)

save.image("../results/2.0_KEGG.RData")
