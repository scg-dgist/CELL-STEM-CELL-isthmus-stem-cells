##### Code for Stomach Single-cell stmn RNA-seq data
##### Author: Eunmin Lee (elee4472@dgist.ac.kr)
##### Last Update: 2019/07/04

library(Seurat)
library(scater)


data_path = "../data/STMN"
setwd(data_path)



PERPLEXITY =50
scesetFiltered = readRDS(file =paste0("scsetFilteredTSNEP",PERPLEXITY,".rds"))


load(file = "sc3set.RData")

seuset <- CreateSeuratObject(
  raw.data = counts(scesetFiltered)
)
seuset@data = assay(scesetFiltered,"logcounts")
seuset@scale.data = assay(scesetFiltered,"logcounts")

seuset@ident = sc3set@colData$sc3_9_clusters
names(seuset@ident) = colnames(seuset@data)


markers <- FindAllMarkers(object = seuset,
                          test.use = "wilcox")
# save(markers, file="markers_using_seurat.RData")


