##### Code for Stomach Single-cell pgc RNA-seq data
##### Author: Eunmin Lee (elee4472@dgist.ac.kr)
##### Last Update: 2019/07/04

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scater)

data_path = "../data/PGC"
setwd(data_path)

# "Differentiated Troy+ Chief Cells Act as Reserve Stem Cells to Generate 
#  All Lineages of the Stomach Epithelium" (Cell, 2013)
#Troy Tnfrsf19
#Mist1 Bhlha15

setwd(data_path)
load("ensemblGenes2017-10-30.RData")
scesetFiltered = readRDS(file = "scsetFilteredTSNEP50.rds")

# genelists
# "ENSMUSG00000024682" "ENSMUSG00000052271" "ENSMUSG00000005087" "ENSMUSG00000028664" "ENSMUSG00000023987"
# "ENSMUSG00000041961" "ENSMUSG00000030029" "ENSMUSG00000000142" "ENSMUSG00000009248" "ENSMUSG00000060548"
# "ENSMUSG00000034177" "ENSMUSG00000020140" "ENSMUSG00000031004" "ENSMUSG00000028832"
# "ENSMUSG00000026558", "ENSMUSG00000030048", "ENSMUSG00000039084", "ENSMUSG00000000567"
# "ENSMUSG00000030545", "ENSMUSG00000033453", "ENSMUSG00000033453", "ENSMUSG00000039530"
# "ENSMUSG00000041895"
ensemblGenes[ensemblGenes$external_gene_name == 'Tspan1','ensembl_gene_id']


TGgene = "ENSMUSG00000028699"

pdf(paste0("Figure_rTSNE_",ensemblGenes[TGgene,"external_gene_name"],"_perplexity_50.pdf"))
df = data.frame(x=attributes(scesetFiltered@reducedDims)$listData$TSNE[,1], 
                y=attributes(scesetFiltered@reducedDims)$listData$TSNE[,2],   
                expression=assay(scesetFiltered,"logcounts")[TGgene,])
ggplot(df,aes(x=x, y=y, colour=expression)) + 
  geom_point(size=2.5) + 
  scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
  ylab("Component 2") + 
  xlab("Component 1") + 
  theme_bw() +
  theme(text = element_text(size=20),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        axis.line=element_line(size=1),
        axis.ticks=element_line(size=1),
        legend.text=element_text(size=20), 
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.text.x = element_text(size=20)
  ) 
dev.off()














