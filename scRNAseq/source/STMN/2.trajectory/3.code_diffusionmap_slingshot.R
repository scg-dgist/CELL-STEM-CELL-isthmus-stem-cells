##### Code for Stomach Single-cell stmn RNA-seq data
##### Author: Eunmin Lee (elee4472@dgist.ac.kr)
##### Last Update: 2019/07/04

library(scater)
library(diffusionMap)
library(slingshot)
# library(M3Drop)
library(dplyr)

data_path = "../../../data/STMN"
setwd(data_path)

PERPLEXITY =50
scesetFiltered = readRDS(file =paste0("scsetFilteredTSNEP",PERPLEXITY,".rds"))


load(file = "sc3set.RData")

newclusters = sc3set@colData$sc3_9_clusters

# m3dgenes method
# m3dGenes <- as.character(
#   M3DropDifferentialExpression(assay(scesetFiltered,"logcounts")[rowSums(exprs(scesetFiltered)) > 0,])$Gene
# )
load("markers_using_seurat.RData")
deGenes_using_seurat = unique(markers[markers[,2]>log(1.5) & markers[,5] < 0.05,7])

isthmus.DE = assay(scesetFiltered,"logcounts")[deGenes_using_seurat,]

set.seed(42)
isthmus.DE.DM = diffuse(as.matrix(dist(t(isthmus.DE))))
load("isthmus.DE.DM.RData")
# save(isthmus.DE.DM, file="isthmus.DE.DM.RData")

df = data.frame(x=isthmus.DE.DM$X[, 1],
                y=isthmus.DE.DM$X[, 2],
                expression=newclusters)


centers = data.frame(df %>% dplyr::group_by(expression) %>% summarize(x = median(x = x), 
                                                                      y = median(x = y)))


ggplot() + 
  geom_point(data=df, size=1, aes(x=x, y=y, colour=expression)) +
  ylab("Component 2") + 
  xlab("Component 1") + 
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(
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
    legend.text=element_text(size=20), 
    legend.title=element_blank(),
    legend.key=element_blank()
  ) +
  geom_text(data = centers,
            mapping = aes(x = x, y = y, label = expression),
            size = 10, col = "black")


lineages <- getLineages(isthmus.DE.DM$X,
                        newclusters
                       , start.clus = '6', end.clus = c("1","2"))

# lineages <- getLineages(isthmus.DE.DM$X,
#                         newclusters
#                         , start.clus = '5')

lineages <- getCurves(lineages)
# save(lineages, file="isthmus.lineages.RData")
load("isthmus.lineages.RData")

publication_color = c("1"="#F8766D",
                      "2"="#D39200",
                      "3"="#93AA00",
                      "4"="#00BA38",
                      "5"="#00C19F",
                      "6"="#00B9E3",
                      "7"="#619CFF",
                      "8"="#DB72FB",
                      "9"="#FF61C3")
clust_col= newclusters
levels(clust_col)[levels(clust_col) == names(publication_color)] = publication_color

pdf("Figure_diffusionMap_slingshot_perplexity_50.pdf")
plot(isthmus.DE.DM$X[, 1],isthmus.DE.DM$X[, 2],col = as.vector(clust_col),pch=19)
lines(lineages@curves$curve1)
lines(lineages@curves$curve2)
lines(lineages@curves$curve4)

lines(lineages)
dev.off()

library(plotly)


df = data.frame(x=isthmus.DE.DM$X[ , 1],
                y=isthmus.DE.DM$X[, 2],
                z=isthmus.DE.DM$X[, 3],
                xl1=lineages@curves$curve3$s[,1], 
                yl1=lineages@curves$curve3$s[,2],
                zl1=lineages@curves$curve3$s[,3], 
                # xl2=lineages_3_2@curves$curve1$s[,1],
                # yl2=lineages_3_2@curves$curve1$s[,2],
                # zl2=lineages_3_2@curves$curve1$s[,3],
                cluster = factor(newclusters)
                )
plot_ly(df ) %>%
  add_trace(x=~x, y=~y, z= ~z, color = ~cluster, mode = "markers", marker = list( size = 4), colors = publication_color) %>%
  add_trace(x=~xl1, y=~yl1, z= ~zl1,marker = list( size = 4),color = I("red"))

