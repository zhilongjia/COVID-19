

CmapUp_raw <- cogena::gmt2list("../data/CmapUp100.gmt")


drug_DEA_list <- list()
# saquinavir@MCF7#5.2E-06M_6246
drug_DEA_list[["Saquinavir"]] <- CmapUp_raw[[4910]]

# ribavirin@PC3#1.64E-05M_7316
drug_DEA_list[["Ribavirin"]] <- CmapUp_raw[[5858]]


# names(CmapUp_raw)[grep("dinoprost@HL60#8.4e-06M_2446", names(CmapUp_raw))]
# names(CmapUp_raw)[grep("dinoprost@MCF7#8.4e-06M_5409", names(CmapUp_raw))]
# drug_DEA_list[["Upgene_dinoprost_2446"]] <- CmapUp_raw[[1661]]
# drug_DEA_list[["Upgene_dinoprost_5409"]] <- CmapUp_raw[[4177]]

devtools::load_all("../../cogena")
# load("../results/DEA_list.RData")
load("../results/DEA_pneu_list.RData")

load("../results/2.1_cogena_KEGG.RData")
DE <- geneExpInCluster(cogena_result, "dia", "3")
DE <- DE$clusterGeneExp
DE <- as.data.frame(DE)
rownames(DE[DE$cluster_id==2,])

drug_DEA_list[["COVID_C2"]] <- rownames(DE[DE$cluster_id==2,])


library(venn)

venn::venn(drug_DEA_list, ilab=TRUE, zcolor = "style", 
           box=FALSE, ilcs=1.5, sncs=1.4, cexil=10, cexsn=10)


drug_gene_venn_res <- gplots::venn(drug_DEA_list, show.plot=FALSE)

aa <- attr(drug_gene_venn_res,"intersections")
aa_df <- as.data.frame(sapply(aa, "length<-", max(lengths(aa))))

aa_tibble <- tidyr::gather(aa_df, key="set", value="SYMBOL", na.rm=TRUE)

source("./gene2KEGG.R")
aa_annotation <- gene2KEGG(aa_tibble$SYMBOL)

venn_annotation <- dplyr::left_join(aa_tibble, aa_annotation)

readr::write_tsv(venn_annotation, path="../results/drug_gene_venn_res.tsv", na="")

save.image("../results/3.1_drug_pathway.RData")


################################################################################
# check gene exp
tmp1 <- t(DEA_list$E$nCoV_Heal)

library(ggpubr)
tmp2 <- as.data.frame(sampleLabel[rownames(tmp1)])
colnames(tmp2) <- "group"
tmp <- cbind(tmp2, tmp1)


boxplot(AAK1 ~ group, tmp)

boxplot(INPP5B ~ group, tmp)

p <- ggboxplot(tmp, x = "group", y = "TMPRSS2",
                palette = "jco",
               add = "jitter")
#  Add p-value
p + stat_compare_means(method = "t.test")

