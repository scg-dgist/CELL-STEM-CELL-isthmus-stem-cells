##### Code for Stomach Single-cell stmn RNA-seq data
##### Author: Eunmin Lee (elee4472@dgist.ac.kr)
##### Last Update: 2019/07/04

library(umap)
library(zoo)
library(reshape2)
library(slingshot)
library(plotly)

data_path = "../../../data/STMN4567"
setwd(data_path)


load("CC_scaled_from_normlog2cd_4567_rebuttal.RData")

deg_csv = read.csv("NL_PN_SL_markers.csv")
degene = as.vector(deg_csv[,2])
load("ensemblGenes2017-10-30.RData")
newclusters_45 = read.csv("newclusterid_cellname_181025.csv",header=TRUE, stringsAsFactors = F)
newclusters = rep(0,ncol(CC_scaled_from_normlog2cd))
names(newclusters) = colnames(CC_scaled_from_normlog2cd)
newclusters[newclusters_45[,1]] = newclusters_45[,2]
newclusters = as.factor(newclusters)
newclusters_4567 = factor(newclusters[newclusters %in% c(1,2,3,4,5,6)])

custom.settings = umap.defaults
custom.settings$n_components =3
# set.seed(113)
# umap_res = umap(t(CC_scaled_from_normlog2cd[degene,names(newclusters_4567)]),
                # config=custom.settings)
# save(umap_res, file="umap_res.RData")
load(file="umap_res.RData")
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
colors = gg_color_hue(6)

umap_xy = data.frame(x=umap_res$layout[,2],
                     y=umap_res$layout[,3]
)

t(umap_res$layout)

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


sling = slingshot(umap_res$layout,
                  newclusters_4567, start.clus=4, end.clus = c("3","6"))
df = data.frame(x=umap_res$layout[,1],
                y=umap_res$layout[,2],
                z=umap_res$layout[,3],
                xl1=sling@curves$curve1$s[,1], 
                yl1=sling@curves$curve1$s[,2],
                zl1=sling@curves$curve1$s[,3], 
                xl2=sling@curves$curve2$s[,1],
                yl2=sling@curves$curve2$s[,2],
                zl2=sling@curves$curve2$s[,3],
                cluster = factor(newclusters_4567)
)

ax <- list(
  title = "",
  zeroline = FALSE,
  showline = FALSE,
  showticklabels = FALSE,
  showgrid = FALSE
)

publication_color = c("1"="#F8766D",
                      "2"="#B79F00",
                      "3"="#00BA38",
                      "4"="#00BFC4",
                      "5"="#619CFF",
                      "6"="#F564E3")
plot_ly(df) %>%
  add_trace(x=~x, y=~y, z= ~z, color = ~cluster, mode = "markers", marker = list( size = 4), colors = publication_color) %>%
  add_trace(x=~xl1, y=~yl1, z= ~zl1,marker = list( size = 4),color = I("red")) %>%
  add_trace(x=~xl2, y=~yl2, z= ~zl2,marker = list( size = 4),color = I("red")) %>%
  layout(scene = list(xaxis=ax,yaxis=ax,zaxis=ax,camera = list(eye = list(x = 1, y = 1, z =0.5))))

# 1 cluster
publication_color = c("1"="#F8766D",
                      "2"="grey",
                      "3"="grey",
                      "4"="grey",
                      "5"="grey",
                      "6"="grey")
plot_ly(df ) %>%
  add_trace(x=~x, y=~y, z= ~z, color = ~cluster, mode = "markers", marker = list( size = 4), colors = publication_color) %>%
  add_trace(x=~xl1, y=~yl1, z= ~zl1,marker = list( size = 4),color = I("red")) %>%
  add_trace(x=~xl2, y=~yl2, z= ~zl2,marker = list( size = 4),color = I("red"))%>%
  layout(scene = list(xaxis=ax,yaxis=ax,zaxis=ax,camera = list(eye = list(x = 1, y = 1, z =0.5))))
# 2 cluster
publication_color = c("1"="grey",
                      "2"="#B79F00",
                      "3"="grey",
                      "4"="grey",
                      "5"="grey",
                      "6"="grey")
plot_ly(df ) %>%
  add_trace(x=~x, y=~y, z= ~z, color = ~cluster, mode = "markers", marker = list( size = 4), colors = publication_color) %>%
  add_trace(x=~xl1, y=~yl1, z= ~zl1,marker = list( size = 4),color = I("red")) %>%
  add_trace(x=~xl2, y=~yl2, z= ~zl2,marker = list( size = 4),color = I("red")) %>%
  layout(scene = list(xaxis=ax,yaxis=ax,zaxis=ax,camera = list(eye = list(x = 1, y = 1, z =0.5))))

# 3 cluster
publication_color = c("1"="grey",
                      "2"="grey",
                      "3"="#00BA38",
                      "4"="grey",
                      "5"="grey",
                      "6"="grey")
plot_ly(df ) %>%
  add_trace(x=~x, y=~y, z= ~z, color = ~cluster, mode = "markers", marker = list( size = 4), colors = publication_color) %>%
  add_trace(x=~xl1, y=~yl1, z= ~zl1,marker = list( size = 4),color = I("red")) %>%
  add_trace(x=~xl2, y=~yl2, z= ~zl2,marker = list( size = 4),color = I("red"))%>%
  layout(scene = list(xaxis=ax,yaxis=ax,zaxis=ax,camera = list(eye = list(x = 1, y = 1, z =0.5))))

# 4 cluster
publication_color = c("1"="grey",
                      "2"="grey",
                      "3"="grey",
                      "4"="#00BFC4",
                      "5"="grey",
                      "6"="grey")
plot_ly(df ) %>%
  add_trace(x=~x, y=~y, z= ~z, color = ~cluster, mode = "markers", marker = list( size = 4), colors = publication_color) %>%
  add_trace(x=~xl1, y=~yl1, z= ~zl1,marker = list( size = 4),color = I("red")) %>%
  add_trace(x=~xl2, y=~yl2, z= ~zl2,marker = list( size = 4),color = I("red"))%>%
  layout(scene = list(xaxis=ax,yaxis=ax,zaxis=ax,camera = list(eye = list(x = 1, y = 1, z =0.5))))

# 5 cluster
publication_color = c("1"="grey",
                      "2"="grey",
                      "3"="grey",
                      "4"="grey",
                      "5"="#619CFF",
                      "6"="grey")
plot_ly(df ) %>%
  add_trace(x=~x, y=~y, z= ~z, color = ~cluster, mode = "markers", marker = list( size = 4), colors = publication_color) %>%
  add_trace(x=~xl1, y=~yl1, z= ~zl1,marker = list( size = 4),color = I("red")) %>%
  add_trace(x=~xl2, y=~yl2, z= ~zl2,marker = list( size = 4),color = I("red"))%>%
  layout(scene = list(xaxis=ax,yaxis=ax,zaxis=ax,camera = list(eye = list(x = 1, y = 1, z =0.5))))

# 6 cluster
publication_color = c("1"="grey",
                      "2"="grey",
                      "3"="grey",
                      "4"="grey",
                      "5"="grey",
                      "6"="#F564E3")
plot_ly(df ) %>%
  add_trace(x=~x, y=~y, z= ~z, color = ~cluster, mode = "markers", marker = list( size = 4), colors = publication_color) %>%
  add_trace(x=~xl1, y=~yl1, z= ~zl1,marker = list( size = 4),color = I("red")) %>%
  add_trace(x=~xl2, y=~yl2, z= ~zl2,marker = list( size = 4),color = I("red")) %>%
  layout(scene = list(xaxis=ax,yaxis=ax,zaxis=ax,camera = list(eye = list(x = 1, y = 1, z =0.5))))







slingPseudotime(sling)

curve1_pseudotime_4_to_3 = slingPseudotime(sling)[,1]
curve1_pseudotime_4_to_3 = curve1_pseudotime_4_to_3[!is.na(curve1_pseudotime_4_to_3)]
curve2_pseudotime_4_to_6 = slingPseudotime(sling)[,2]
curve2_pseudotime_4_to_6 = curve2_pseudotime_4_to_6[!is.na(curve2_pseudotime_4_to_6)]

geneList_sym = c("Mki67","Stmn1", "Muc6", "Muc5ac")
geneList = ensemblGenes$ensembl_gene_id[which(ensemblGenes$external_gene_name %in% geneList_sym)]


features.plot = geneList

cd =  CC_scaled_from_normlog2cd


#curve1_pseudotime_4_to_3
TSmat = matrix(data=NA, nrow=length(features.plot), 
               ncol=length(colnames(cd[,names(curve1_pseudotime_4_to_3)])), 
               dimnames = list(features.plot, colnames(cd[,names(curve1_pseudotime_4_to_3)])))
rownames(TSmat) = ensemblGenes[features.plot,"external_gene_name"]
for(i in features.plot){
  y = cd[i,names(curve1_pseudotime_4_to_3)[order(curve1_pseudotime_4_to_3)]]
  Y_hat = rollapply(y,20, mean, align = 'center',fill = 'extend')
  TSmat[ensemblGenes[i,"external_gene_name"],] = Y_hat
}

df.TSmat = data.frame(t(TSmat), x = 1:length(colnames(cd[,names(curve1_pseudotime_4_to_3)])))
df.TSmat.long <- melt(df.TSmat, id=c("x"))
pdf("rebuttal_UMAP_CC_4567_curve1_pseudotime_4_to_3_psuedotime_plot_rollmean_selected_genes.pdf")
ggplot(data = df.TSmat.long, aes(x=x, y=value, color=variable)) +
  geom_line() +
  ylab("Auto exprs") +
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

newclusters_4567[names(curve1_pseudotime_4_to_3)]

newclusters_45

# annotation bar

publication_color = c("1"="#F8766D",
                      "2"="#B79F00",
                      "3"="#00BA38",
                      "4"="#00BFC4",
                      "5"="#619CFF",
                      "6"="#F564E3")
color_list = list(C = publication_color)
annotation = factor(newclusters_4567[names(curve1_pseudotime_4_to_3)])
annotdf = data.frame(row.names = names(curve1_pseudotime_4_to_3)[order(curve1_pseudotime_4_to_3)], C= annotation[order(curve1_pseudotime_4_to_3)])
annotdf
ha = ComplexHeatmap::HeatmapAnnotation(df = annotdf, col = color_list)
pdf("annotation_curve1_pseudotime_4_to_3.pdf")
ComplexHeatmap::draw(ha,1:198)
dev.off

#curve2_pseudotime_4_to_6
TSmat = matrix(data=NA, nrow=length(features.plot), 
               ncol=length(colnames(cd[,names(curve2_pseudotime_4_to_6)])), 
               dimnames = list(features.plot, colnames(cd[,names(curve2_pseudotime_4_to_6)])))
rownames(TSmat) = ensemblGenes[features.plot,"external_gene_name"]
for(i in features.plot){
  y = cd[i,names(curve2_pseudotime_4_to_6)[order(curve2_pseudotime_4_to_6)]]
  Y_hat = rollapply(y,20, mean, align = 'center',fill = 'extend')
  TSmat[ensemblGenes[i,"external_gene_name"],] = Y_hat
}

df.TSmat = data.frame(t(TSmat), x = 1:length(colnames(cd[,names(curve2_pseudotime_4_to_6)])))
df.TSmat.long <- melt(df.TSmat, id=c("x"))
pdf("rebuttal_UMAP_CC_4567_curve2_pseudotime_4_to_6_psuedotime_plot_rollmean_selected_genes.pdf")
ggplot(data = df.TSmat.long, aes(x=x, y=value, color=variable)) +
  geom_line() +
  ylab("Auto exprs") +
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

newclusters_4567[names(curve2_pseudotime_4_to_6)]

newclusters_45

# annotation bar

publication_color = c("1"="#F8766D",
                      "2"="#B79F00",
                      "3"="#00BA38",
                      "4"="#00BFC4",
                      "5"="#619CFF",
                      "6"="#F564E3")
color_list = list(C = publication_color)
annotation = factor(newclusters_4567[names(curve2_pseudotime_4_to_6)])
annotdf = data.frame(row.names = names(curve2_pseudotime_4_to_6)[order(curve2_pseudotime_4_to_6)], C= annotation[order(curve2_pseudotime_4_to_6)])
annotdf
ha = ComplexHeatmap::HeatmapAnnotation(df = annotdf, col = color_list)
pdf("annotation_curve2_pseudotime_4_to_6.pdf")
ComplexHeatmap::draw(ha,1:180)
dev.off()



umap_xy = data.frame(x=umap_res$layout[,2],
                     y=umap_res$layout[,3]
)

t(umap_res$layout)

# color each Old cluster
colors = gg_color_hue(6)
df = data.frame(x=umap_xy$x, 
                y=umap_xy$y,
                xl1=sling@curves$curve1$s[,2],
                yl1=sling@curves$curve1$s[,3],
                xl2=sling@curves$curve2$s[,2],
                yl2=sling@curves$curve2$s[,3],
                expression=newclusters_4567)

ggplot(df) +
  geom_point(aes(x=x, y=y, colour=expression),size=3) + 
  geom_point(aes(x=xl1,y=yl1),size=0.5) +
  geom_point(aes(x=xl2,y=yl2),size=0.5) +
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

