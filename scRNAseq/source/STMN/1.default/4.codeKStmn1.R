##### Code for Stomach Single-cell stmn RNA-seq data
##### Author: Eunmin Lee (elee4472@dgist.ac.kr)
##### Last Update: 2019/07/04

library(scater)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)


data_path = "../../../data/STMN"
setwd(data_path)
PERPLEXITY =50
scesetFiltered = readRDS(file =paste0("scsetFilteredTSNEP",PERPLEXITY,".rds"))


load("ensemblGenes2017-10-30.RData")

ensemblGenes[ensemblGenes$external_gene_name == 'Stmn1','ensembl_gene_id']
# ENSMUSG00000028832
ensemblGenes[ensemblGenes$external_gene_name == 'Mki67','ensembl_gene_id']
# ENSMUSG00000031004

logcounts(scesetFiltered)


stmn1 = exprs(scesetFiltered)["ENSMUSG00000028832",]

allgenes = exprs(scesetFiltered)[rowSums(exprs(scesetFiltered)) > 0,]

corrGene = apply(allgenes, 1, function(x) cor(x,stmn1, method="spearman"))
corrGene["ENSMUSG00000031004"] # stmn1, mki67 0.70
stmn1_corrGene_0.3 = rownames(allgenes)[corrGene>=0.3]
# hist(corrGene)
# length(stmn1_corrGene)
allgenes_stmn1_corr_0.3 = exprs(scesetFiltered)[stmn1_corrGene_0.3,]
rownames(allgenes)

corr_df = data.frame(ensembl= rownames(allgenes), symbol=ensemblGenes[rownames(allgenes), "external_gene_name"], corr=corrGene )
corr_df = corr_df[order(-corr_df$corr),]
# write.csv(corr_df, file= "stmn1_corr_gene_list_spearman.csv", row.names=F )

pdf(paste("Figure_heatmap_Perplexity",PERPLEXITY,"_hclust_k4_stmn1_allgene_corr_spearman_0.3_2.pdf",sep=""), onefile=F)
set.seed(127) # seed changed??? from 123?
res = pheatmap(t(allgenes_stmn1_corr_0.3), 
               kmeans_k = 4,
               cluster_rows = FALSE,
               clustering_distance_cols = "correlation")

dev.off()

df = data.frame(x=attributes(scesetFiltered@reducedDims)$listData$TSNE[,1], 
                y=attributes(scesetFiltered@reducedDims)$listData$TSNE[,2],  
                expression=factor(res$kmeans[[1]]))
ggplot(df,aes(x=x, y=y, colour=expression)) + 
  geom_point(size=1.5) + 
  scale_colour_manual(values=c("#878787", "#d73027", "#4575b4", "#7fbc41")) +
  ylab("Component 2") + 
  xlab("Component 1") + 
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(text = element_text(size=20),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        axis.line=element_line(size=1),
        axis.ticks=element_line(size=1),
        legend.text=element_text(size=10), 
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.text.x = element_text(size=20)
  ) 
dev.off()
dim(scesetFiltered)
scesetFiltered_stmn1 = scesetFiltered[,!res$kmeans[[1]] == 1]


scesetFiltered_stmn1_minus = scesetFiltered[,res$kmeans[[1]] == 1]

#612 cells -> 321 cells
# saveRDS(scesetFiltered_stmn1,file =paste0("scsetFiltered_stmn1.rds"))



