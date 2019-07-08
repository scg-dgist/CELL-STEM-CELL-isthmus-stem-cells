##### Code for Stomach Single-cell stmn RNA-seq data
##### Author: Eunmin Lee (elee4472@dgist.ac.kr)
##### Last Update: 2019/07/04

library(scater)
library(reshape2)
library(ggplot2)
library(zoo)
library(slingshot)

data_path = "../../../data/STMN"
setwd(data_path)

PERPLEXITY =50
scesetFiltered = readRDS(file =paste0("scsetFilteredTSNEP",PERPLEXITY,".rds"))

load(file = "sc3set.RData")
load("ensemblGenes2017-10-30.RData")
newclusters = sc3set@colData$sc3_9_clusters

load("isthmus.lineages.RData")

geneList_sym = c("Mki67","Stmn1", "Muc6", "Muc5ac")
geneList = ensemblGenes$ensembl_gene_id[which(ensemblGenes$external_gene_name %in% geneList_sym)]

cd = assay(scesetFiltered,"logcounts")

cd_zscore = t(scale(t(cd)))

lambdas = slingPseudotime(lineages)

features.plot = geneList
# annotation bar

t=lambdas[,1]
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
annotation = factor(newclusters[order(t, na.last = NA)])
annotdf = data.frame(row.names = colnames(cd)[order(t, na.last = NA)], C= annotation)

ha = ComplexHeatmap::HeatmapAnnotation(df = annotdf, col = color_list)
pdf("annotation_path1.pdf")
ComplexHeatmap::draw(ha,1:dim(annotdf)[1])
dev.off()

# annotation bar

t=  lambdas[,3]
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
annotation = factor(newclusters[order(t, na.last = NA)])
annotdf = data.frame(row.names = colnames(cd)[order(t, na.last = NA)], C= annotation)

ha = ComplexHeatmap::HeatmapAnnotation(df = annotdf, col = color_list)
pdf("annotation_path3.pdf")
ComplexHeatmap::draw(ha,1:dim(annotdf)[1])
dev.off()

###### log rolling mean


t=  lambdas[,1]
sorted_cellname = colnames(cd)[order(t, na.last = NA)]
path1_cellname = sorted_cellname[150:190]
# save(path1_cellname, file="path1_cellname.RData")

features.plot = geneList

control_cluster = factor(sc3set@colData$sc3_9_clusters[!is.na(lambdas[,1])])

TSmat = matrix(data=NA, nrow=length(features.plot), 
               ncol=length(colnames(cd[,!is.na(lambdas[,1])])), dimnames = list(features.plot, colnames(cd_zscore[,!is.na(lambdas[,1])])))
rownames(TSmat) = ensemblGenes[features.plot,"external_gene_name"]
for(i in features.plot){
  y = cd[i,order(t, na.last = NA)]
  Y_hat = rollapply(y,30, mean, align = 'center',fill = 'extend')
  TSmat[ensemblGenes[i,"external_gene_name"],] = Y_hat
}

df.TSmat = data.frame(t(TSmat), x = 1:length(colnames(cd_zscore[,!is.na(lambdas[,1])])))
df.TSmat.long <- melt(df.TSmat, id=c("x"))

pdf("FigureX_path_totalpath1_psuedotime_plot_log_rollmean_selected_genes.pdf")
ggplot(data = df.TSmat.long, aes(x=x, y=value, color=variable)) +
  geom_line() +
  ylab("exprs") +
  xlab("Cell order") +
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
  )
dev.off()


t=  lambdas[,3]
sorted_cellname = colnames(cd)[order(t, na.last = NA)]
path2_cellname = sorted_cellname[150:200]
# save(path2_cellname, file="path2_cellname.RData")


features.plot = geneList

control_cluster = factor(sc3set@colData$sc3_9_clusters[!is.na(lambdas[,3])])

TSmat = matrix(data=NA, nrow=length(features.plot), 
               ncol=length(colnames(cd[,!is.na(lambdas[,3])])), dimnames = list(features.plot, colnames(cd_zscore[,!is.na(lambdas[,3])])))
rownames(TSmat) = ensemblGenes[features.plot,"external_gene_name"]
for(i in features.plot){
  y = cd[i,order(t, na.last = NA)]
  Y_hat = rollapply(y,30, mean, align = 'center',fill = 'extend')
  TSmat[ensemblGenes[i,"external_gene_name"],] = Y_hat
}

df.TSmat = data.frame(t(TSmat), x = 1:length(colnames(cd_zscore[,!is.na(lambdas[,3])])))
df.TSmat.long <- melt(df.TSmat, id=c("x"))

pdf("FigureX_path_totalpath3_psuedotime_plot_log_rollmean_selected_genes.pdf")
ggplot(data = df.TSmat.long, aes(x=x, y=value, color=variable)) +
  geom_line() +
  ylab("exprs") +
  xlab("Cell order") +
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
  )
dev.off()




#######

t=  lambdas[,1]

features.plot = geneList

control_cluster = factor(sc3set@colData$sc3_9_clusters[!is.na(t)])

TSmat = matrix(data=NA, nrow=length(features.plot), 
               ncol=length(colnames(cd_zscore[,!is.na(t)])), dimnames = list(features.plot, colnames(cd_zscore[,!is.na(t)])))
rownames(TSmat) = ensemblGenes[features.plot,"external_gene_name"]
for(i in features.plot){
  y = cd_zscore[i,order(t, na.last = NA)]
  Y_hat = rollapply(y,50, mean, align = 'center',fill = 'extend')
  TSmat[ensemblGenes[i,"external_gene_name"],] = Y_hat
}

df.TSmat = data.frame(t(TSmat), x = 1:length(colnames(cd_zscore[,!is.na(t)])))
df.TSmat.long <- melt(df.TSmat, id=c("x"))
setwd(figure_path)
pdf("FigureX_path_totalpath1_psuedotime_plot_rollmean_selected_genes.pdf")
ggplot(data = df.TSmat.long, aes(x=x, y=value, color=variable)) +
  geom_line() +
  ylab("exprs") +
  xlab("Cell order") +
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
  )
dev.off()

t=  lambdas[,3]

features.plot = geneList

control_cluster = factor(sc3set@colData$sc3_9_clusters[!is.na(t)])

TSmat = matrix(data=NA, nrow=length(features.plot), 
               ncol=length(colnames(cd_zscore[,!is.na(t)])), dimnames = list(features.plot, colnames(cd_zscore[,!is.na(t)])))
rownames(TSmat) = ensemblGenes[features.plot,"external_gene_name"]
for(i in features.plot){
  y = cd_zscore[i,order(t, na.last = NA)]
  Y_hat = rollapply(y,50, mean, align = 'center',fill = 'extend')
  TSmat[ensemblGenes[i,"external_gene_name"],] = Y_hat
}

df.TSmat = data.frame(t(TSmat), x = 1:length(colnames(cd_zscore[,!is.na(lambdas[,3])])))
df.TSmat.long <- melt(df.TSmat, id=c("x"))

pdf("FigureX_path_totalpath3_psuedotime_plot_rollmean_selected_genes.pdf")
ggplot(data = df.TSmat.long, aes(x=x, y=value, color=variable)) +
  geom_line() +
  ylab("exprs") +
  xlab("Cell order") +
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
  )
dev.off()
       

