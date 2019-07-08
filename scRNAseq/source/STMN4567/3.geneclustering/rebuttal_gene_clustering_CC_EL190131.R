##### Code for Stomach Single-cell stmn RNA-seq data
##### Author: Eunmin Lee (elee4472@dgist.ac.kr)
##### Last Update: 2019/07/04


data_path = "../../../data/STMN4567"
setwd(data_path)

#library loading

# library(WGCNA) # dynamic cut X
# library(flashClust)          X
# library(slingshot) # slingshot
# library(dplyr) # slingshot
# library(plotly) # slingshot
library(ggplot2)
library(pheatmap)
library(Seurat)
library(amap)
library(dplyr)
library(RColorBrewer)
#data loading
scesetFiltered = readRDS(file = "scsetFiltered4567TSNEP50.rds")
load("CC_scaled_from_normlog2cd_4567_rebuttal.RData")

load("ensemblGenes2017-10-30.RData")
newclusters_45 = read.csv("newclusterid_cellname_181025.csv",header=TRUE, stringsAsFactors = F) # sc3 cluster info
newclusters = rep(0,ncol(scesetFiltered))
names(newclusters) = colnames(scesetFiltered)
newclusters[newclusters_45[,1]] = newclusters_45[,2]
newclusters = as.factor(newclusters)
newclusters_4567 = factor(newclusters[newclusters %in% c(1,2,3,4,5,6)])

seurat = CreateSeuratObject(CC_scaled_from_normlog2cd[,names(newclusters_4567)])
seurat@data = CC_scaled_from_normlog2cd[,names(newclusters_4567)]
seurat@scale.data = CC_scaled_from_normlog2cd[,names(newclusters_4567)]

#hirachial CC
deg_csv = read.csv("NL_PN_SL_markers.csv")
degene = as.vector(deg_csv[,2]) # EnsembleID
logcd = CC_scaled_from_normlog2cd[degene,names(newclusters_4567)]

logcd = logcd[rowMeans(logcd)>0,]
lognormcd_1 = rowMeans(logcd[,names(newclusters_4567)[newclusters_4567 == 1]])
lognormcd_2 = rowMeans(logcd[,names(newclusters_4567)[newclusters_4567 == 2]])
lognormcd_3 = rowMeans(logcd[,names(newclusters_4567)[newclusters_4567 == 3]])
lognormcd_4 = rowMeans(logcd[,names(newclusters_4567)[newclusters_4567 == 4]])
lognormcd_5 = rowMeans(logcd[,names(newclusters_4567)[newclusters_4567 == 5]])
lognormcd_6 = rowMeans(logcd[,names(newclusters_4567)[newclusters_4567 == 6]])
df = data.frame(c1=lognormcd_1,
                c2=lognormcd_2,
                c3=lognormcd_3,
                c4=lognormcd_4,
                c5=lognormcd_5,
                c6=lognormcd_6)
hc <- hclust(Dist(t(df),method = "spearman"))
pdf("cluster_hclust_spearman.pdf")
plot(hc)
dev.off()
# grouping 123, 6, 45
group1_ind = names(newclusters_4567)[newclusters_4567 %in% c(4,5)]
group2_ind = names(newclusters_4567)[newclusters_4567 %in% c(1,2,3)]
group3_ind = names(newclusters_4567)[newclusters_4567 %in% c(6)]
group_cluster = as.vector(newclusters_4567)
names(group_cluster) = names(newclusters_4567)
group_cluster[group1_ind] = "45"
group_cluster[group2_ind] = "123"
group_cluster[group3_ind] = "6"
seurat@ident = factor(group_cluster)

#
markers = FindAllMarkers(seurat)
degene =  subset(markers, p_val_adj < 0.05 & avg_logFC > log(1.5))$gene

degene = unique(subset(markers, p_val_adj < 0.05 & avg_logFC > log(1.5))$gene)

deGenes_table = subset(markers, p_val_adj < 0.05 & avg_logFC > log(1.5) & gene %in% degene)
deGenes_table_df = data.frame(deGenes_table$cluster,
                              deGenes_table$gene,
                              ensemblGenes[deGenes_table$gene,"external_gene_name"],
                              deGenes_table$pct.1,
                              deGenes_table$pct.2,
                              deGenes_table$avg_logFC,
                              deGenes_table$p_val_adj)
colnames(deGenes_table_df) = c("Cluster","Ensembl gene id", "Gene Symbol","expressed pct.in","expressed pct.out", "AvgLogFC","Adj. p-value")
head(deGenes_table_df)
write.csv(deGenes_table_df, "cell cycle corrected Stmn high single-cell data 3group Marker gene List.txt", quote =F, row.names=F)
write.csv(deGenes_table_df, "cell cycle corrected Stmn high single-cell data 3group Marker gene List.csv", quote =F, row.names=F)




load(file="umap_res.RData")
# write.csv(degene,file="backgroundgene_degene.txt", quote = F, row.names = F)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
colors = gg_color_hue(6)

umap_xy = data.frame(x=umap_res$layout[,2],
                     y=umap_res$layout[,3]
)


# color each Old cluster
colors = gg_color_hue(6)
df = data.frame(x=umap_xy$x, 
                y=umap_xy$y, 
                expression=newclusters_4567)

ggplot(df,aes(x=x, y=y, colour=expression)) +
  geom_point() + 
  scale_colour_manual(values = colors ) +
  ylab("Component 2") +
  xlab("Component 1") +
  # ggtitle(gene_names[target_genes[i],"external_gene_name"]) +
  theme_bw() +
  theme(    plot.title = element_text(hjust = 0.5,size = 40),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            panel.background=element_blank(),
            panel.border=element_blank(),
            plot.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(), 
            axis.line = element_blank(),
            axis.ticks=element_blank(),
            # legend.position="none",
            plot.margin=unit(c(0,3,0,0),"cm")
  )

# sling = slingshot(umap_res$layout,
#                   newclusters_4567, start.clus=4, end.clus = c("3","6"))
# save(sling, file="sling_20190201.RData")
# curve1_pseudotime_4_to_3 = slingPseudotime(sling)[,1]
# curve1_pseudotime_4_to_3 = curve1_pseudotime_4_to_3[!is.na(curve1_pseudotime_4_to_3)]
# save(curve1_pseudotime_4_to_3, file="curve1_pseudotime_4_to_3_20190201.RData")
# curve2_pseudotime_4_to_6 = slingPseudotime(sling)[,2]
# curve2_pseudotime_4_to_6 = curve2_pseudotime_4_to_6[!is.na(curve2_pseudotime_4_to_6)]
# save(curve2_pseudotime_4_to_6, file="curve2_pseudotime_4_to_6_20190201.RData")
load( file="curve1_pseudotime_4_to_3_20190201.RData")
load( file="curve2_pseudotime_4_to_6_20190201.RData")


#pseudotime
curve1_pseudotime_4_to_3 #path1
curve2_pseudotime_4_to_6 #path2

path1_cell_ind_sorted = names(sort(curve1_pseudotime_4_to_3))
path2_cell_ind_sorted = names(sort(curve2_pseudotime_4_to_6))

CC_scaled_gene_clustering = CC_scaled_from_normlog2cd[degene,] # CC scaled data using degene
# CC_scaled_gene_clustering = CC_scaled_from_normlog2cd # CC scaled data using all gene

CC_scaled_gene_clustering_path1 = CC_scaled_gene_clustering[,path1_cell_ind_sorted]
CC_scaled_gene_clustering_path2 = CC_scaled_gene_clustering[,path2_cell_ind_sorted]


#smoothing
#path1
df = data.frame(row.names= path1_cell_ind_sorted,
                pseudotime = sort(curve1_pseudotime_4_to_3)
)
cd = CC_scaled_gene_clustering_path1
target_pseudotime = df$pseudotime
x = target_pseudotime/max(target_pseudotime)
X = seq(0.0, 1.0, length.out = length(colnames(CC_scaled_gene_clustering_path1)))
TSmat = matrix(data=NA, nrow=length(degene),
               ncol=length(X), dimnames = list(degene, X))
for(i in degene){
  y = cd[i,rownames(df)]
  ls = loess(y ~ x, span=0.5)
  Y_hat = predict(ls, X)
  # TSmat[genename[i,"gene_symbol"],] = Y_hat  # use gene_symbol
  TSmat[i,] = Y_hat
}

plot(CC_scaled_gene_clustering_path1[1,])
plot(TSmat[1,])

gg_color_hue(6)

as.vector(newclusters_4567[rownames(df)])
annotation_col = data.frame(
  row.names = names(newclusters_4567[rownames(df)]),
  CC4567.cluster = newclusters_4567[rownames(df)])
gg_color_hue(6)
ann_colors = list("CC4567.cluster" = c("1"="#F8766D",
                              "2"="#B79F00",
                              "3"="#00BA38",
                              "4"="#00BFC4",
                              "5"="#619CFF",
                              "6"="#F564E3"))
TSmat_colnames = TSmat
colnames(TSmat_colnames) = rownames(annotation_col)
set.seed(1)
ph =pheatmap(TSmat_colnames,
             color = colorRampPalette(c("#4575B4","white","#D73027"))(100),
             cluster_cols = F, 
             clustering_distance_rows = "correlation",
             show_colnames = F,
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             kmeans_k = 5,
             filename = "CC 4567 path1_c4_to_c3_genecluster_heatmap.pdf")
dev.off()
df_path1_kmean = data.frame(kmeans = ph$kmeans$cluster, symbol = ensemblGenes[names(ph$kmeans$cluster),"external_gene_name"])
# write.csv(df_path1_kmean, file= "CC 4567 path1_kmean_cluster_info.csv")


#path2
df = data.frame(row.names= path2_cell_ind_sorted,
                pseudotime = sort(curve2_pseudotime_4_to_6)
)
cd = CC_scaled_gene_clustering_path2
target_pseudotime = df$pseudotime
x = target_pseudotime/max(target_pseudotime)
X = seq(0.0, 1.0, length.out = length(colnames(CC_scaled_gene_clustering_path2)))
TSmat = matrix(data=NA, nrow=length(degene),
               ncol=length(X), dimnames = list(degene, X))
for(i in degene){
  y = cd[i,rownames(df)]
  ls = loess(y ~ x, span=0.5)
  Y_hat = predict(ls, X)
  # TSmat[genename[i,"gene_symbol"],] = Y_hat  # use gene_symbol
  TSmat[i,] = Y_hat
}

plot(CC_scaled_gene_clustering_path2[1,])
plot(TSmat[1,])

gg_color_hue(6)

as.vector(newclusters_4567[rownames(df)])
annotation_col = data.frame(
  row.names = names(newclusters_4567[rownames(df)]),
  CC4567.cluster = newclusters_4567[rownames(df)])
gg_color_hue(6)
ann_colors = list("CC4567.cluster" = c("1"="#F8766D",
                                       "2"="#B79F00",
                                       "3"="#00BA38",
                                       "4"="#00BFC4",
                                       "5"="#619CFF",
                                       "6"="#F564E3"))
TSmat_colnames = TSmat
colnames(TSmat_colnames) = rownames(annotation_col)
set.seed(1)
ph =pheatmap(TSmat_colnames,
             color = colorRampPalette(c("#4575B4","white","#D73027"))(100),
             cluster_cols = F, 
             clustering_distance_rows = "correlation",
             show_colnames = F,
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             kmeans_k = 5,
             filename = "CC 4567 path2_c4_to_c6_genecluster_heatmap.pdf")
dev.off()
df_path2_kmean = data.frame(kmeans = ph$kmeans$cluster, symbol = ensemblGenes[names(ph$kmeans$cluster),"external_gene_name"])
# write.csv(df_path2_kmean, file= "CC 4567 path2_kmean_cluster_info.csv")

