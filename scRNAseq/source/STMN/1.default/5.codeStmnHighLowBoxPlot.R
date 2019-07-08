##### Code for Stomach Single-cell stmn RNA-seq data
##### Author: Eunmin Lee (elee4472@dgist.ac.kr)
##### Last Update: 2019/07/04


library(scater)
library(ggplot2)
data_path = "../../../data/STMN"
setwd(data_path)

PERPLEXITY =50
scesetFiltered = readRDS(file =paste0("scsetFilteredTSNEP",PERPLEXITY,".rds"))

scesetFiltered_stmn1_high = readRDS(file =paste0("scsetFiltered_stmn1.rds"))
scesetFiltered_stmn1_low = scesetFiltered[,!(colnames(scesetFiltered) %in% colnames(scesetFiltered_stmn1_high))]

# rowcount
genesList = rownames(scesetFiltered[rowMeans(assay(scesetFiltered)) > 10,])

genesList2 = rownames(scesetFiltered[rowMeans(2^(exprs(scesetFiltered))-1) > 10,])

genesList = genesList2


###################################



x = cor(exprs(scesetFiltered)[genesList,])
corr_all = 1-x[upper.tri(x)]

x = cor(exprs(scesetFiltered_stmn1_high)[genesList,])
corr_1 = 1-x[upper.tri(x)]

x = cor(exprs(scesetFiltered_stmn1_low)[genesList,])
corr_0 = 1-x[upper.tri(x)]


sample_name = c(rep("all", length(corr_all)),rep("stmn1_high", length(corr_1)), rep("stmn1_low", length(corr_0)))
corr_value = c(corr_all, corr_1, corr_0)
df = data.frame( sample = sample_name, value = corr_value)

write.csv(df, file = "1-corr_values.csv" )

setwd(figure_path)
pdf("Figure_stmn1_high_low_homogeneity_1-corr_boxplot_all.pdf")
ggplot(df, aes(x = sample, y = value)) + 
  geom_violin()  +
  geom_boxplot()  + 
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
  ) + ylab("1-Corr.") + xlab("")
dev.off()

