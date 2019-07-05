##### Code for Stomach Single-cell stmn RNA-seq data
##### Author: Eunmin Lee (elee4472@dgist.ac.kr)
##### Last Update: 2019/07/04

library(ggplot2)
library(SC3)
library(Rtsne)
library(scater)

data_path = "../data/STMN4567"
setwd(data_path)


load(file="CC_scaled_from_normlog2cd_4567_rebuttal.RData")
load("hvg4567_cc_rebuttal.RData")
load("ensemblGenes2017-10-30.RData")
## tSNE 
PERPLEXITY = 50
set.seed(601) #601
Y = Rtsne(t(CC_scaled_from_normlog2cd[hvg_cc,]), max_iter=5000, verbose=TRUE, pca_scale=TRUE, pca_center=TRUE, perplexity=PERPLEXITY)

scesetFiltered = SingleCellExperiment(assays=list(counts =CC_scaled_from_normlog2cd, logcounts=CC_scaled_from_normlog2cd))
rowData(scesetFiltered)$feature_symbol = ensemblGenes[rownames(scesetFiltered), "external_gene_name"]


attributes(scesetFiltered@reducedDims)$listData$TSNE = Y$Y

setwd(data_path)
# save tSNE information into scesetFiltered
# saveRDS(scesetFiltered, paste("scsetFiltered_CC_4567TSNEP",PERPLEXITY,"_rebuttal.rds",sep=""))
scesetFiltered = readRDS(file = paste("scsetFiltered_CC_4567TSNEP",PERPLEXITY,"_rebuttal.rds",sep=""))

setwd(figure_path)
pdf(paste("Figure_rTSNE_Perplexity",PERPLEXITY,"_rebuttal.pdf",sep=""))
cellColor = c("#d73027", "#4575b4")
par(mar=c(5,5,1,1), cex.axis=1.5, cex.lab=1.5)
plot(attributes(scesetFiltered@reducedDims)$listData$TSNE[,1],
     attributes(scesetFiltered@reducedDims)$listData$TSNE[,2], pch=20, col="black",
     xlab="Component 1", ylab="Component 2", cex=1.5)
points(attributes(scesetFiltered@reducedDims)$listData$TSNE[colData(scesetFiltered)[,"condition"]=="23002", 1],
       attributes(scesetFiltered@reducedDims)$listData$TSNE[colData(scesetFiltered)[,"condition"]=="23002", 2], pch=20, col=cellColor[1], cex=1.5)
points(attributes(scesetFiltered@reducedDims)$listData$TSNE[colData(scesetFiltered)[,"condition"]=="23003", 1],
       attributes(scesetFiltered@reducedDims)$listData$TSNE[colData(scesetFiltered)[,"condition"]=="23003", 2], pch=20, col=cellColor[2], cex=1.5)
legend(-60, -50, c("23002", "23003"),
       pch=20,
       col=cellColor, bty="n", cex=1)
dev.off()


sc3_estimate_k(scesetFiltered[hvg_cc,])
# ks = use estimate_k 
sc3set <- sc3(scesetFiltered[hvg_cc,],ks = 6, gene_filter= F,k_estimator = TRUE, biology = TRUE)

est_k = 6

setwd(figure_path)
pdf("Figure_sc3_consensus_isthmus_rebuttal.pdf", onefile = FALSE)
sc3_plot_consensus(sc3set, k= est_k)
dev.off()


# save(sc3set, file= "sc3set4567_CC_rebuttal.RData")

i = est_k
cluster_sc3 = as.factor(sc3set$sc3_6_clusters) # sc3_k_clusters #if cluster # is changed, this also should be changed.

pdf(paste("Figure_TSNE_CC_cluster4567_Perplexity50_clustering_SC3_K",i,"_rebuttal.pdf", sep="")) 
df = data.frame(x=attributes(scesetFiltered@reducedDims)$listData$TSNE[,1], 
                y=attributes(scesetFiltered@reducedDims)$listData$TSNE[,2], 
                expression=cluster_sc3)
ggplot(df,aes(x=x, y=y, colour=expression)) + 
  geom_point(size=2.5) + 
  # scale_colour_manual(values=clust.col) +
  ylab("Component 2") + 
  xlab("Component 1") + 
  theme_bw() +
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
