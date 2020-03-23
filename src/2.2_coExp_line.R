

################################################################################
# co-exp line plot
devtools::load_all("../../cogena")
load("../results/2.1_cogena_KEGG.RData")
DE <- geneExpInCluster(cogena_result, "dia", "3")
DE <- DE$clusterGeneExp
DE <- as.data.frame(DE)
DE$cluster_id <- as.factor(DE$cluster_id)
DE_short <- tidyr::gather(tibble::as_tibble(tibble::rownames_to_column(DE)), "sample", "expression", 3:(ncol(DE)+1))


DE_short$sample <- factor(DE_short$sample, levels=unique(DE_short$sample) )

ggplot(DE_short, aes_(x=~sample, y=~expression)) + #geom_line(aes_(group=~rowname, color=~cluster_id)) + 
  geom_smooth(aes_(group=~rowname, color=~cluster_id), method = "auto", se = FALSE)+
  facet_grid(rows = vars(cluster_id)) 

################################################################################
# limma res with cluster information
limma_res_raw <- DEA_list[["limma_DEA"]][["nCoV_Heal"]]

DE_tibble <- tibble::rownames_to_column(DE, var="SYMBOL")[,1:2]

source("./gene2KEGG.R")

gene_annotation <- gene2KEGG(DE_tibble$rowname)

DEG_annotation <- dplyr::left_join(DE_tibble, gene_annotation)

limma_res_CoExp <- dplyr::left_join(limma_res_raw, DE_tibble, by=c("nCoV_Heal"="rowname"))

readr::write_csv(DEG_annotation, path="../results/DEG_CoEXp.csv") 

readr::write_csv(limma_res_CoExp, path="../results/limma_res_CoEXp.csv") 

################################################################################
# Neutrophil_degranulation pathway
Neutrophil_degranulation_reactome_gene <- readr::read_delim("../data/Participating Molecules [R-HSA-6798695].tsv", delim="\t| ")
nCOV_DE <- DEA_list[["limma_DE"]]$nCoV_Heal

overlapped_gene <- intersect(Neutrophil_degranulation_reactome_gene$Name, rownames(nCOV_DE))
nCOV_DE <- nCOV_DE[overlapped_gene,]

gplots::heatmap.2(nCOV_DE, trace = "none", col="bluered", scale="row", labRow =NA, tracecol=NA,
                  srtCol=60, keysize=1, cexCol=1.3, cexRow=1.1, Colv=TRUE,
                  ylab=paste(nrow(nCOV_DE), "DE Genes"), margins = c(4, 2))

################################################################################

# clusterColor <- sample(rainbow(5))
# tmp <- map2col(as.numeric(as.factor(DE$cluster_id)), clusterColor)
