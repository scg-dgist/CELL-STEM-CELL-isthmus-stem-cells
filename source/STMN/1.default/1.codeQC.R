##### Code for Stomach Single-cell stmn RNA-seq data
##### Author: Eunmin Lee (elee4472@dgist.ac.kr)
##### Last Update: 2019/07/04

library(scran)
library(scater)
library(biomaRt)
library(Rtsne)
library(ggplot2)
library(Rtsne)
library(SC3)
library(RColorBrewer)

data_path = "../data/STMN"
setwd(data_path)
##### Load data


countData = read.csv("stomach_isthmus_count_table.csv", row.names=1, header=TRUE, check.names=FALSE) # 384 cells

countData = countData[rowSums(countData) > 0, ]
cellInfo = data.frame(cell=colnames(countData), 
                      condition=sapply(colnames(countData), function(x) substr(x, 1, 5)))


cd = as.matrix(countData)
logcd = log2(cd + 1)

## Create SCSET
stomach_sce <- SingleCellExperiment(assays=list(counts=cd, logcounts=logcd), colData = cellInfo)

# exprs(stomach_sce) = log2(calculateCPM(stomach_sce, use.size.factors = FALSE) + 1)



rm(countData)
rm(cd)
rm(logcd)

# date = Sys.Date()
# ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", host="www.ensembl.org")
# ensemblGenes <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name',  'chromosome_name', 'gene_biotype'), mart=ensembl)
# rownames(ensemblGenes) <- ensemblGenes[,1]
# mtGenes = ensemblGenes[ensemblGenes[,3]=="MT",]
# save(ensemblGenes, file=paste("ensemblGenes",date,".RData",sep=""))
load("ensemblGenes2017-10-30.RData")

# remove control wells
control_ind = grepl("(^.*_1$|^.*_381$|^.*_382$|^.*_383$|^.*_384$)",colnames(stomach_sce))

stomach_sce = stomach_sce[,!control_ind]
colnames(stomach_sce)



mtGenes = ensemblGenes[ensemblGenes[,3]=="MT",]
is.mito = rownames(stomach_sce) %in% mtGenes[,1]
length(is.mito[is.mito == TRUE])
stomach_sce = calculateQCMetrics(stomach_sce,  feature_controls=list(Mt=is.mito))
stomach_sce = runPCA(stomach_sce, pca_data_input = "colData")  


# plotPCA(
#   stomach_sce,
#   size_by = "total_features",
#   # colour_by = "use",
#   pca_data_input = "pdata"
# )



cellColor = c("#d73027", "#4575b4")
pdf("Figure_QC_totalGeneCount_MTProp.pdf")
par(cex.axis=1.5, cex.lab=1.5)
plot(NULL, xaxt="n", xlab="Number of reads mapped to exons",
     ylab="Percentage of reads mapped to mitochondrial genes", xlim=c(10^2,1*10^7), ylim=c(0,100), log="x")
axis( 1, 10^(seq(2, 7, 2)), c(expression(10^2), expression(10^4), expression(10^6)) )
points(colData(stomach_sce)[colData(stomach_sce)[,"condition"]=="23002", "total_counts"],
       colData(stomach_sce)[colData(stomach_sce)[,"condition"]=="23002", "pct_counts_Mt"], pch=20, cex=1, col=cellColor[1])
points(colData(stomach_sce)[colData(stomach_sce)[,"condition"]=="23003", "total_counts"],
       colData(stomach_sce)[colData(stomach_sce)[,"condition"]=="23003", "pct_counts_Mt"], pch=20, cex=1, col=cellColor[2])
abline(h=10, v=20000, col="red", lty=2, lwd=2)
legend(1000000, 100,c("23002", "23003"),
       pch=20,
       col=cellColor, bty="n", cex=1.5)
dev.off()


pdf("Figure_QC_PCA.pdf")
cellColor = c("#d73027", "#4575b4")
PC1 <- format(attributes(attributes(stomach_sce@reducedDims)$listData$PCA)$percentVar[1]*100,digits=3)
PC2 <- format(attributes(attributes(stomach_sce@reducedDims)$listData$PCA)$percentVar[2]*100,digits=3)
par(mar=c(5,5,1,1), cex.axis=1.5, cex.lab=1.5)
plot(attributes(stomach_sce@reducedDims)$listData$PCA[,"PC1"], 
     attributes(stomach_sce@reducedDims)$listData$PCA[,"PC2"], pch=20, col="black",
     xlab=paste("PC1 ", "(", PC1, "%)", sep=""), ylab=paste("PC2 ", "(", PC2, "%)", sep=""), cex=0.5)
points(attributes(stomach_sce@reducedDims)$listData$PCA[colData(stomach_sce)[,"condition"]=="23002", "PC1"],
       attributes(stomach_sce@reducedDims)$listData$PCA[colData(stomach_sce)[,"condition"]=="23002", "PC2"], pch=20, col=cellColor[1], cex=0.5)
points(attributes(stomach_sce@reducedDims)$listData$PCA[colData(stomach_sce)[,"condition"]=="23003", "PC1"],
       attributes(stomach_sce@reducedDims)$listData$PCA[colData(stomach_sce)[,"condition"]=="23003", "PC2"], pch=20, col=cellColor[2], cex=0.5)
legend(26, 14, c("23002", "23003"),
       pch=20,
       col=cellColor, bty="n", cex=1)
dev.off()


pMT = 10

pdf("Figure_QC_mtProp_PCA.pdf")
cellColor = c("#fc710a", "#01665e")
PC1 <- format(attributes(attributes(stomach_sce@reducedDims)$listData$PCA)$percentVar[1]*100,digits=3)
PC2 <- format(attributes(attributes(stomach_sce@reducedDims)$listData$PCA)$percentVar[2]*100,digits=3)
par(mar=c(5,5,1,1), cex.axis=1.5, cex.lab=1.5)
plot(attributes(stomach_sce@reducedDims)$listData$PCA[,"PC1"], 
     attributes(stomach_sce@reducedDims)$listData$PCA[,"PC2"], pch=20, col="black",
     xlab=paste("PC1 ", "(", PC1, "%)", sep=""), ylab=paste("PC2 ", "(", PC2, "%)", sep=""), cex=0.5)
points(attributes(stomach_sce@reducedDims)$listData$PCA[colData(stomach_sce)[,"pct_counts_Mt"]>=pMT, "PC1"],
       attributes(stomach_sce@reducedDims)$listData$PCA[colData(stomach_sce)[,"pct_counts_Mt"]>=pMT, "PC2"], pch=20, col=cellColor[1], cex=0.5)
points(attributes(stomach_sce@reducedDims)$listData$PCA[colData(stomach_sce)[,"pct_counts_Mt"]<pMT, "PC1"],
       attributes(stomach_sce@reducedDims)$listData$PCA[colData(stomach_sce)[,"pct_counts_Mt"]<pMT, "PC2"], pch=20, col=cellColor[2], cex=0.5)
legend(27, 12, c(">10% of MT", "<10% of MT"),
       pch=20,
       col=cellColor, bty="n", cex=1)
dev.off()


# calculate pct_dropout
assay(stomach_sce, "is_exprs") <- calcIsExprs(stomach_sce, 
                                               lowerDetectionLimit = 5)

colData(stomach_sce)$pct_dropout <- (1 - apply(assay(stomach_sce, "is_exprs")  , 2, mean))* 100


pDROPOUT = 97
table(rowData(stomach_sce)$pct_dropout_counts < 99.9)

pdf("Figure_QC_dropout_PCA.pdf")
cellColor = c("#fc710a", "#01665e")
PC1 <- format(attributes(attributes(stomach_sce@reducedDims)$listData$PCA)$percentVar[1]*100,digits=3)
PC2 <- format(attributes(attributes(stomach_sce@reducedDims)$listData$PCA)$percentVar[2]*100,digits=3)
par(mar=c(5,5,1,1), cex.axis=1.5, cex.lab=1.5)
plot(attributes(stomach_sce@reducedDims)$listData$PCA[,"PC1"], 
     attributes(stomach_sce@reducedDims)$listData$PCA[,"PC2"], pch=20, col="black",
     xlab=paste("PC1 ", "(", PC1, "%)", sep=""), ylab=paste("PC2 ", "(", PC2, "%)", sep=""), cex=0.5)
points(attributes(stomach_sce@reducedDims)$listData$PCA[colData(stomach_sce)[,"pct_dropout"]>=pDROPOUT, "PC1"],
       attributes(stomach_sce@reducedDims)$listData$PCA[colData(stomach_sce)[,"pct_dropout"]>=pDROPOUT, "PC2"], pch=20, col=cellColor[1], cex=0.5)
points(attributes(stomach_sce@reducedDims)$listData$PCA[colData(stomach_sce)[,"pct_dropout"]<pDROPOUT, "PC1"],
       attributes(stomach_sce@reducedDims)$listData$PCA[colData(stomach_sce)[,"pct_dropout"]<pDROPOUT, "PC2"], pch=20, col=cellColor[2], cex=0.5)
legend(24, 12, c(">97% of dropout", "<97% of dropout"),
       pch=20,
       col=cellColor, bty="n", cex=1)
dev.off()

## After QC
scesetFiltered = stomach_sce[rowData(stomach_sce)$pct_dropout_counts<99.9, colData(stomach_sce)[,"pct_counts_Mt"]<pMT & colData(stomach_sce)[,"pct_dropout"]<pDROPOUT]

is.mito = rownames(scesetFiltered) %in% mtGenes[,1]
scesetFiltered = calculateQCMetrics(scesetFiltered,  feature_controls=list(Mt=is.mito))
scesetFiltered = runPCA(scesetFiltered, pca_data_input = "colData")  

pdf("Figure_POSTQC_PCA.pdf")
cellColor = c("#d73027", "#4575b4")
PC1 <- format(attributes(attributes(scesetFiltered@reducedDims)$listData$PCA)$percentVar[1]*100,digits=3)
PC2 <- format(attributes(attributes(scesetFiltered@reducedDims)$listData$PCA)$percentVar[2]*100,digits=3)
par(mar=c(5,5,1,1), cex.axis=1.5, cex.lab=1.5)
plot(attributes(scesetFiltered@reducedDims)$listData$PCA[,"PC1"], 
     attributes(scesetFiltered@reducedDims)$listData$PCA[,"PC2"], pch=20, col="black",
     xlab=paste("PC1 ", "(", PC1, "%)", sep=""), ylab=paste("PC2 ", "(", PC2, "%)", sep=""), cex=0.5)
points(attributes(scesetFiltered@reducedDims)$listData$PCA[colData(scesetFiltered)[,"condition"]=="23002", "PC1"],
       attributes(scesetFiltered@reducedDims)$listData$PCA[colData(scesetFiltered)[,"condition"]=="23002", "PC2"], pch=20, col=cellColor[1], cex=0.5)
points(attributes(scesetFiltered@reducedDims)$listData$PCA[colData(scesetFiltered)[,"condition"]=="23003", "PC1"],
       attributes(scesetFiltered@reducedDims)$listData$PCA[colData(scesetFiltered)[,"condition"]=="23003", "PC2"], pch=20, col=cellColor[2], cex=0.5)
legend(-19,16 , c("23002", "23003"),
       pch=20,
       col=cellColor, bty="n", cex=1)
dev.off()

cellColor = c("#d73027", "#4575b4")
pdf("Figure_POSTQC_totalGeneCount_MTProp.pdf")
par(cex.axis=1.5, cex.lab=1.5)
plot(NULL, xaxt="n", xlab="Number of reads mapped to exons",
     ylab="Percentage of reads mapped to mitochondrial genes", xlim=c(10^2,1*10^7), ylim=c(0,100), log="x")
axis( 1, 10^(seq(2, 7, 2)), c(expression(10^2), expression(10^4), expression(10^6)) )
points(colData(scesetFiltered)[colData(scesetFiltered)[,"condition"]=="23002", "total_counts"],
       colData(scesetFiltered)[colData(scesetFiltered)[,"condition"]=="23002", "pct_counts_Mt"], pch=20, cex=1, col=cellColor[1])
points(colData(scesetFiltered)[colData(scesetFiltered)[,"condition"]=="23003", "total_counts"],
       colData(scesetFiltered)[colData(scesetFiltered)[,"condition"]=="23003", "pct_counts_Mt"], pch=20, cex=1, col=cellColor[2])
abline(h=10, v=20000, col="red", lty=2, lwd=2)
legend(1000000, 100,c("23002", "23003"),
       pch=20,
       col=cellColor, bty="n", cex=1.5)
dev.off()

## Normalization with scran

clusters <- quickCluster(scesetFiltered)
scesetFiltered <- computeSumFactors(scesetFiltered, cluster=clusters)
scesetFiltered <- normalize(scesetFiltered)

pdf("Figure_SF_LibrarySize.pdf")
plot(sizeFactors(scesetFiltered), scesetFiltered$total_counts/1e3, log="xy",
     ylab="Library size (kilos)", xlab="Size factor", pch=20)
dev.off()

# saveRDS(scesetFiltered, "scsetFiltered_PCTDROPOUT999.rds")

## Detecting highly variable genes


var.fit <- trendVar(scesetFiltered, method="spline", parametric=TRUE, 
                    use.spikes=FALSE, span=0.2)#, design = batchDesign)
var.out <- decomposeVar(scesetFiltered, var.fit)
hvg <- var.out[which(var.out$FDR <= 0.05 & var.out$bio > .5),]

pdf("Figure_MeanExp_VarExp_BIO_0.5.pdf")
plot(var.out$mean, var.out$total, pch=16, cex=0.3, xlab="Mean log-expression",
     ylab="Variance of log-expression")
o <- order(var.out$mean)
lines(var.out$mean[o], var.out$tech[o], col="dodgerblue", lwd=2)
points(var.out$mean[var.out$FDR <= 0.05 & var.out$bio > .5], var.out$total[var.out$FDR <= 0.05 & var.out$bio > .5], pch=16, cex=0.3, col="red")
dev.off()

# save(hvg, file="hvg.RData")
