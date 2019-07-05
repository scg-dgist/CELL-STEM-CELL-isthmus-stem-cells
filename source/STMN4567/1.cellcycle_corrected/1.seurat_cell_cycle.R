##### Code for Stomach Single-cell stmn RNA-seq data
##### Author: Eunmin Lee (elee4472@dgist.ac.kr)
##### Last Update: 2019/07/04


library(Seurat)
library(scater)
data_path = "../data/STMN4567"
setwd(data_path)

load("ensemblGenes2017-10-30.RData")
scesetFiltered = readRDS(file =paste0("scsetFiltered4567TSNEP50.rds"))

seuset <- CreateSeuratObject(
  raw.data = counts(scesetFiltered)
)
seuset@data = assay(scesetFiltered,"logcounts")
seuset@scale.data = assay(scesetFiltered,"logcounts")



# CELL CYCLE GENES 

G1.genes = c("Camk2a",
             "Camk2b",
             "Gpr132",
             "Itgb1",
             "Mtbp",
             "Myb",
             "Nfatc1",
             "Ppp2r3a",
             "Ppp3ca",
             "Skp2",
             "Taf10",
             "Slfn1")

S.genes = c("Dnajc2",
            "Mcm2",
            "Mcm3",
            "Mcm4",
            "Mki67",
            "Mre11a",
            "Msh2",
            "Pcna",
            "Rad17",
            "Rad51",
            "Sumo1")

#http://saweb2.sabiosciences.com/rt_pcr_product/HTML/PAMM-020A.html
G2M.genes = c(
  "Dnajc2",
  "Chek1",
  "Ppm1d",
  "Brca2",
  "Ccna1",
  "Ccnb1",
  "Cdc25a",
  "Cdc25b",
  "Cdk2",
  "Nek2",
  "Npm2",
  "Pes1",
  "Prm1",
  "Rad21",
  "Ran",
  "Shc1",
  "Smc1a",
  "Stag1",
  "Terf1",
  "Psmg2",
  "Wee1"
)

s.genes_mouse = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% S.genes]
g2m.genes_mouse = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% G2M.genes]




#Alternative CELL CYCLE GENE from Seurat tutorial

cc_gene_list = read.table("regev_lab_cell_cycle_genes.txt", header = FALSE)  # converting human cell cylce gene -> mouse
require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("hgnc_symbol","ensembl_gene_id"), filters = "hgnc_symbol", values = cc_gene_list$V1 , mart = human, attributesL = c("mgi_symbol", "ensembl_gene_id"), martL = mouse, uniqueRows=T)

s.genes_human <- cc_gene_list[1:43,"V1"]
g2m.genes_human <- cc_gene_list[44:97,"V1"]
s.genes_mouse <- genesV2[genesV2$HGNC.symbol %in% s.genes_human,]$Gene.stable.ID.1
g2m.genes_mouse <- genesV2[genesV2$HGNC.symbol %in% g2m.genes_human,]$Gene.stable.ID.1


# I've used first cell cycle genes.
seuset <- CellCycleScoring(object = seuset, s.genes = s.genes_mouse, g2m.genes = g2m.genes_mouse, set.ident = TRUE)

# using seuset@data <- normalized 
# return seuset@scale.data <- cell cycle corrected count table
seuset <- ScaleData(object = seuset, vars.to.regress = c("S.Score", "G2M.Score"), 
                    display.progress = FALSE)



# Save CC count table
CC_scaled_from_normlog2cd = seuset@scale.data
# save(CC_scaled_from_normlog2cd,file="CC_scaled_from_normlog2cd_4567_rebuttal.RData")


# finding HVG in CC 
seuset@data = seuset@scale.data
seuset <- FindVariableGenes(object = seuset, mean.function = ExpMean, dispersion.function = LogVMR, 
                            x.low.cutoff = 0.0125, x.high.cutoff = 2, y.cutoff = 2)
hvg_cc =  seuset@var.genes
# save(hvg_cc, file="hvg4567_cc_rebuttal.RData")
