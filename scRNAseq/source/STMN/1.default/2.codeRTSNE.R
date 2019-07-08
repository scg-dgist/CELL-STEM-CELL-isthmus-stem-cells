##### Code for Stomach Single-cell stmn RNA-seq data
##### Author: Eunmin Lee (elee4472@dgist.ac.kr)
##### Last Update: 2019/07/04

library(Rtsne)
library(scater)
library(dplyr)
library(SC3)


data_path = "../data/STMN"
setwd(data_path)
load("hvg.RData")
scesetFiltered = readRDS(file = "scsetFiltered_PCTDROPOUT999.rds")
load("ensemblGenes2017-10-30.RData")
## tSNE 
PERPLEXITY = 50
dim(exprs(scesetFiltered))
set.seed(601) #601
Y = Rtsne(t(exprs(scesetFiltered)[rownames(hvg),]), max_iter=5000, verbose=TRUE, pca_scale=TRUE, pca_center=TRUE, perplexity=PERPLEXITY)

attributes(scesetFiltered@reducedDims)$listData$TSNE = Y$Y

# saveRDS(scesetFiltered, paste("scsetFilteredTSNEP",PERPLEXITY,".rds",sep=""))
scesetFiltered = readRDS(file = paste("scsetFilteredTSNEP",PERPLEXITY,".rds",sep=""))

pdf(paste("Figure_rTSNE_Perplexity",PERPLEXITY,".pdf",sep=""))
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

rowData(scesetFiltered)$feature_symbol = ensemblGenes[rownames(scesetFiltered), "external_gene_name"]


sc3set <- sc3(scesetFiltered[rownames(hvg),], ks = 9,k_estimator = TRUE, biology = TRUE)
# sc3set <- sc3(sc3set, ks = 16,k_estimator = TRUE, biology = TRUE, n_cores = 4)
str(sc3set) # get estimated K

pdf("Figure_sc3_consensus_isthmus.pdf", onefile = FALSE)
sc3_plot_consensus(sc3set, k= 9)
dev.off()


# save(sc3set, file= "sc3set.RData")
load("sc3set.RData")

i = 9
cluster_sc3 = as.factor(sc3set$sc3_9_clusters)
setwd(figure_path)
pdf(paste("Figure_TSNE_Perplexity",PERPLEXITY,"_clustering_SC3_K",i,".pdf", sep="")) 
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

