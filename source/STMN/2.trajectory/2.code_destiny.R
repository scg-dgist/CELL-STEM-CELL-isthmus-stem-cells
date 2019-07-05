##### Code for Stomach Single-cell stmn RNA-seq data
##### Author: Eunmin Lee (elee4472@dgist.ac.kr)
##### Last Update: 2019/07/04

library(scater)
library(destiny)
library(zoo)
library(reshape2)

data_path = "../data/STMN"
setwd(data_path)


PERPLEXITY =50
scesetFiltered = readRDS(file =paste0("scsetFilteredTSNEP",PERPLEXITY,".rds"))


load(file = "sc3set.RData")

newclusters = sc3set@colData$sc3_9_clusters

load("markers_using_seurat.RData")
deGenes_using_seurat = unique(markers[markers[,2]>log(1.5) & markers[,5] < 0.05,7])

deGenes_table = subset(markers, p_val_adj < 0.05 & avg_logFC > log(1.5) & gene %in% deGenes_using_seurat)
deGenes_table_df = data.frame(deGenes_table$cluster,
                              deGenes_table$gene,
                              ensemblGenes[deGenes_table$gene,"external_gene_name"],
                              deGenes_table$pct.1,
                              deGenes_table$pct.2,
                              deGenes_table$avg_logFC,
                              deGenes_table$p_val_adj)
colnames(deGenes_table_df) = c("Cluster","Ensembl gene id", "Gene Symbol","expressed pct.in","expressed pct.out", "AvgLogFC","Adj. p-value")
head(deGenes_table_df)
write.csv(deGenes_table_df, "Stmn single-cell data Marker gene List.csv", quote =F, row.names=F)

load("ensemblGenes2017-10-30.RData")


exprs.DE = exprs(scesetFiltered)[deGenes_using_seurat,]


dm = DiffusionMap(t(exprs.DE))
dpt = DPT(dm)
# save(dpt, file="dpt.RData")
plot(dpt, pch = 20)
plot(dpt, col_by = "branch")


levels(newclusters) = c("#F8766D",
                     "#D39200",
                     "#93AA00",
                     "#00BA38",
                     "#00C19F",
                     "#00B9E3",
                     "#619CFF",
                     "#DB72FB",
                     "#FF61C3")

pdf("FigureX_destiny_pseudotime.pdf", width = 9.6)
plot(dpt, root = 1, paths_to = c(2,3), pch = 20, col_path = c("grey","grey"), cex=1.5)
dev.off()


pdf("FigureX_destiny_cluster.pdf")
plot(dpt, root = 1, paths_to = c(2,3), pch = 20, col_path = c("grey","grey"), col = as.vector(newclusters), cex=1.5 )
dev.off()

pseudotime = dpt$DPT600
names(pseudotime) = colnames(scesetFiltered)

geneList_sym = c("Mki67",  "Stmn1" , "Muc6" ,  "Muc5ac")
geneList = rownames(ensemblGenes)[ensemblGenes$external_gene_name %in% geneList_sym]

############# path 6-> 1
ind_6_1 = colnames(scesetFiltered)[sc3set@colData$sc3_9_clusters %in% c(6,4,7,5,8,1)]
cds_6_1 = scesetFiltered[geneList,ind_6_1]

cd_zscore = t(scale(t(exprs(cds_6_1))))
pseudotime_6_1 = pseudotime[ind_6_1]

features.plot = geneList

control_cluster = factor(sc3set@colData$sc3_9_clusters[ind_6_1])

TSmat = matrix(data=NA, nrow=length(features.plot), 
               ncol=length(colnames(cd_zscore[,ind_6_1])), dimnames = list(features.plot, colnames(cd_zscore[,ind_6_1])))
rownames(TSmat) = ensemblGenes[features.plot,"external_gene_name"]
for(i in features.plot){
  y = cd_zscore[i,order(pseudotime_6_1)]
  Y_hat = rollapply(y,30, mean, align = 'center',fill = 'extend')
  TSmat[ensemblGenes[i,"external_gene_name"],] = Y_hat
}

df.TSmat = data.frame(t(TSmat), x = 1:length(colnames(cd_zscore[,ind_6_1])))
df.TSmat.long <- melt(df.TSmat, id=c("x"))

pdf("FigureX_path_6_1_destiny_pseudotime_plot_rollmean_selected_genes.pdf")
ggplot(data = df.TSmat.long, aes(x=x, y=value, colour=variable)) +
  geom_line() +
  ylab("auto. exprs") +
  xlab("Cells") +
  theme_bw() +
  theme(text = element_text(size=20),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        axis.line=element_line(size=1),
        axis.ticks=element_line(size=1),
        legend.position = "right",
        legend.text=element_text(size=10), 
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.text.x = element_text(size=20)
  ) + scale_color_manual(values=c("#619CFF","#00BA38","#F8766D","#D39200") )
dev.off()



############# path 6-> 2
ind_6_2 = colnames(scesetFiltered)[sc3set@colData$sc3_9_clusters %in% c(6,4,7,5,8,2)]
cds_6_2 = scesetFiltered[geneList,ind_6_2]

cd_zscore = t(scale(t(exprs(cds_6_2))))
pseudotime_6_2 = pseudotime[ind_6_2]

############# path 6-> 2
features.plot = geneList

control_cluster = factor(sc3set@colData$sc3_9_clusters[ind_6_2])

TSmat = matrix(data=NA, nrow=length(features.plot), 
               ncol=length(colnames(cd_zscore[,ind_6_2])), dimnames = list(features.plot, colnames(cd_zscore[,ind_6_2])))
rownames(TSmat) = ensemblGenes[features.plot,"external_gene_name"]
for(i in features.plot){
  y = cd_zscore[i,order(pseudotime_6_2)]
  Y_hat = rollapply(y,30, mean, align = 'center',fill = 'extend')
  TSmat[ensemblGenes[i,"external_gene_name"],] = Y_hat
}

df.TSmat = data.frame(t(TSmat), x = 1:length(colnames(cd_zscore[,ind_6_2])))
df.TSmat.long <- melt(df.TSmat, id=c("x"))
setwd(figure_path)
pdf("FigureX_path_6_2_destiny_pseudotime_plot_rollmean_selected_genes.pdf")
ggplot(data = df.TSmat.long, aes(x=x, y=value, color=variable)) +
  geom_line() +
  ylab("auto. exprs") +
  xlab("Cells") +
  theme_bw() +
  theme(text = element_text(size=20),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        axis.line=element_line(size=1),
        axis.ticks=element_line(size=1),
        legend.position = "right",
        legend.text=element_text(size=10), 
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.text.x = element_text(size=20)
  ) + scale_color_manual(values=c("#619CFF","#00BA38","#F8766D","#D39200") )
dev.off()


# annotation bar


publication_color = c("1"="#F8766D",
                      "2"="#D39200",
                      "3"="#93AA00",
                      "4"="#00BA38",
                      "5"="#00C19F",
                      "6"="#00B9E3",
                      "7"="#619CFF",
                      "8"="#DB72FB",
                      "9"="#FF61C3")
color_list = list(C = publication_color)
annotation = factor(sc3set@colData$sc3_9_clusters[sc3set@colData$sc3_9_clusters %in% c(6,4,7,5,8,1)])
annotdf = data.frame(row.names = ind_6_1[order(pseudotime_6_1)], C= annotation[order(pseudotime_6_1)])
annotdf
ha = ComplexHeatmap::HeatmapAnnotation(df = annotdf, col = color_list)
pdf("destiny_annotation_bar_6_1.pdf")
ComplexHeatmap::draw(ha,1:434)
dev.off()

publication_color = c("1"="#F8766D",
                      "2"="#D39200",
                      "3"="#93AA00",
                      "4"="#00BA38",
                      "5"="#00C19F",
                      "6"="#00B9E3",
                      "7"="#619CFF",
                      "8"="#DB72FB",
                      "9"="#FF61C3")
color_list = list(C = publication_color)
annotation = factor(sc3set@colData$sc3_9_clusters[sc3set@colData$sc3_9_clusters %in% c(6,4,7,5,8,2)])
annotdf = data.frame(row.names = ind_6_2[order(pseudotime_6_2)], C= annotation[order(pseudotime_6_2)])
annotdf
ha = ComplexHeatmap::HeatmapAnnotation(df = annotdf, col = color_list)
pdf("destiny_annotation_bar_6_2.pdf")
ComplexHeatmap::draw(ha,1:434)
dev.off()



