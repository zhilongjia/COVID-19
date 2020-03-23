

load("../results/DEA_pneu_list.RData")

nCoV_Heal_DE <- DEA_list[["limma_DE"]][["nCoV_Heal"]]
nCoV_Heal_pheno <- DEA_list[["subsample_pheno"]][["nCoV_Heal"]] 


sampleLabel <- nCoV_Heal_pheno$Types
names(sampleLabel) <- nCoV_Heal_pheno$NewID

################################################################################
devtools::load_all("../../cogena")
# library(cogena)

nClust <- 2:10
ncore <- 7
# clMethods <- c("hierarchical","kmeans","diana","fanny","som","sota","pam","clara","agnes")
clMethods <- c("hierarchical","kmeans","diana", "pam")

genecl_result <- coExp(nCoV_Heal_DE, nClust=nClust, 
                       clMethods=clMethods, 
                       metric="correlation", 
                       method="complete", 
                       ncore=ncore, 
                       verbose=TRUE)



# Distmat <- amap::Dist(nCoV_Heal_DE, method="correlation", nbproc=ncore)
# tmp <- diana(Distmat)
# plot(tmp)

###########################################
##########################################
annoGMT <- "c2.cp.kegg.v7.01.symbols.gmt.xz"

annofile <- system.file("extdata", annoGMT, package="cogena")


cogena_result <- clEnrich(genecl_result, annofile=annofile, sampleLabel=sampleLabel, ncore=ncore)

summary(cogena_result)
# library("scales") 
heatmapCluster(cogena_result, "diana", "3", maintitle="COVID-19", add2=F,  heatmapcol="bluered",
               cexCol=1.0, clusterColor=scales::hue_pal()(3)  )

heatmapPEI(cogena_result, "diana", "3", maintitle="COVID-19", add2=F,
           CutoffNumGeneset=20, orderMethod = "each")

heatmapPEI(cogena_result, "diana", "3", maintitle="COVID-19", add2=F,
           CutoffNumGeneset=20, geom="circle", orderMethod = "each")

save.image("../results/2.1_cogena_KEGG.RData")

##########################################
devtools::load_all("../../cogena")
annoGMT <- "c5.all.v7.0.symbols.gmt"

annofile <- system.file("extdata", annoGMT, package="cogena")
cogena_go_result <- clEnrich(genecl_result, annofile=annofile, sampleLabel=sampleLabel, ncore=ncore)

summary(cogena_go_result)

heatmapPEI(cogena_go_result, "diana", "3", maintitle="COVID-19", add2=F,
           CutoffNumGeneset=20, orderMethod = "each")

heatmapPEI(cogena_go_result, "diana", "3", maintitle="COVID-19", add2=TRUE,
           CutoffNumGeneset=20, orderMethod = "each", geom="circle")

save.image("../results/2.1_cogena_go.RData")

##########################################
devtools::load_all("../../cogena")
annoGMT <- "c2.cp.reactome.v7.01.symbols.gmt"

annofile <- system.file("extdata", annoGMT, package="cogena")
cogena_reactome_result <- clEnrich(genecl_result, annofile=annofile, sampleLabel=sampleLabel, ncore=ncore)

summary(cogena_reactome_result)

heatmapPEI(cogena_reactome_result, "diana", "3", maintitle="COVID-19", add2=F,
           CutoffNumGeneset=20, orderMethod = "each")

heatmapPEI(cogena_reactome_result, "diana", "3", maintitle="COVID-19", add2=TRUE,
                       CutoffNumGeneset=20, orderMethod = "each", geom="circle")

save.image("../results/2.1_cogena_reactome.RData")

#########################################################################
devtools::load_all("../../cogena")
annoGMT <- "CmapUp100.gmt.xz";
annofile <- system.file("extdata", annoGMT, package="cogena")
drugUP_result <- clEnrich(genecl_result, annofile=annofile, sampleLabel=sampleLabel, ncore=ncore)

heatmapPEI(drugUP_result, "diana", "3", maintitle="COVID-19", add2=F,
           CutoffNumGeneset=20, orderMethod = "2", printGS = TRUE)

heatmapPEI(drugUP_result, "diana", "3", maintitle="COVID-19", add2=FALSE,
           CutoffNumGeneset=20, orderMethod = "3", printGS = TRUE)

save.image("../results/2.1_cogena_CmapUP.RData")

#########################################################################
annoGMT <- "CmapDn100.gmt.xz"; 
annofile <- system.file("extdata", annoGMT, package="cogena")
drugDN_result <- clEnrich(genecl_result, annofile=annofile, sampleLabel=sampleLabel, ncore=ncore)

heatmapPEI(drugDN_result, "diana", "3", maintitle="COVID-19", add2=FALSE,
           CutoffNumGeneset=20, orderMethod = "1", printGS = TRUE)

save.image("../results/2.1_cogena_CmapDN.RData")


################################################################################
################################################################################
# drug annotation
load("../../cmap_annotation/results/CMap_MoA.RData")
# res_anno <- list()
res_PEI_C2 <- heatmapPEI(drugUP_result, "diana", "3", maintitle="COVID-19", add2=F,
                         CutoffNumGeneset=20, orderMethod = "2", printGS = TRUE)
drp_C2 <- dplyr::left_join(res_PEI_C2[,1], CMap_MoA)

res_PEI_C3 <- heatmapPEI(drugUP_result, "diana", "3", maintitle="COVID-19", add2=F,
                         CutoffNumGeneset=20, orderMethod = "3", printGS = TRUE)
drp_C3 <- dplyr::left_join(res_PEI_C3[,1], CMap_MoA)

res_PEI_C1 <- heatmapPEI(drugDN_result, "diana", "3", maintitle="COVID-19", add2=F,
           CutoffNumGeneset=20, orderMethod = "1", printGS = TRUE)
drp_C1 <- dplyr::left_join(res_PEI_C1[,1], CMap_MoA)


save.image("../results/2.1_cogena.RData")

