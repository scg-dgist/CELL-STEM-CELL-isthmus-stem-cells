##### Code for Stomach Single-cell stmn RNA-seq data
##### Author: Eunmin Lee (elee4472@dgist.ac.kr)
##### Last Update: 2019/07/04

library(scater)
library(ggplot2)


PERPLEXITY =50
data_path = "../data/STMN"
setwd(data_path)
scesetFiltered = readRDS(file =paste0("scsetFilteredTSNEP",PERPLEXITY,".rds"))


setwd(data_path)
load("ensemblGenes2017-10-30.RData")

ensemblGenes[ensemblGenes$external_gene_name == 'Stmn1','ensembl_gene_id']
# ENSMUSG00000028832
ensemblGenes[ensemblGenes$external_gene_name == 'Mki67','ensembl_gene_id']
# ENSMUSG00000031004
ensemblGenes[ensemblGenes$external_gene_name == 'Bhlha15','ensembl_gene_id'] # Mist1
# ENSMUSG00000052271
ensemblGenes[ensemblGenes$external_gene_name == 'Lgr5','ensembl_gene_id'] 
# ENSMUSG00000020140
ensemblGenes[ensemblGenes$external_gene_name == 'Sox2','ensembl_gene_id']
# ENSMUSG00000074637
ensemblGenes[ensemblGenes$external_gene_name == 'Runx1','ensembl_gene_id']
# ENSMUSG00000022952
ensemblGenes[ensemblGenes$external_gene_name == 'Lrig1','ensembl_gene_id']
# ENSMUSG00000030029 
ensemblGenes[ensemblGenes$external_gene_name == 'Tnfrsf19','ensembl_gene_id'] #Troy Tnfrsf19
# ENSMUSG00000060548
ensemblGenes[ensemblGenes$external_gene_name == 'Muc5ac','ensembl_gene_id'] 
# ENSMUSG00000037974
ensemblGenes[ensemblGenes$external_gene_name == 'Muc6','ensembl_gene_id'] 
# ENSMUSG00000048191

stmn1_1 = colnames(exprs(scesetFiltered))[exprs(scesetFiltered)["ENSMUSG00000028832",] > 0]
stmn1_0 = colnames(exprs(scesetFiltered))[!(exprs(scesetFiltered)["ENSMUSG00000028832",] > 0)]
mki67_1 = colnames(exprs(scesetFiltered))[exprs(scesetFiltered)["ENSMUSG00000031004",] > 0]
mki67_0 = colnames(exprs(scesetFiltered))[!(exprs(scesetFiltered)["ENSMUSG00000031004",] > 0)]
mist1_1 = colnames(exprs(scesetFiltered))[exprs(scesetFiltered)["ENSMUSG00000052271",] > 0]
mist1_0 = colnames(exprs(scesetFiltered))[!(exprs(scesetFiltered)["ENSMUSG00000052271",] > 0)]
lgr5_1 = colnames(exprs(scesetFiltered))[exprs(scesetFiltered)["ENSMUSG00000020140",] > 0]
lgr5_0 = colnames(exprs(scesetFiltered))[!exprs(scesetFiltered)["ENSMUSG00000020140",] > 0]
sox2_1 = colnames(exprs(scesetFiltered))[exprs(scesetFiltered)["ENSMUSG00000074637",] > 0]
sox2_0 = colnames(exprs(scesetFiltered))[!(exprs(scesetFiltered)["ENSMUSG00000074637",] > 0)]
runx1_1 = colnames(exprs(scesetFiltered))[exprs(scesetFiltered)["ENSMUSG00000022952",] > 0]
runx1_0 = colnames(exprs(scesetFiltered))[!(exprs(scesetFiltered)["ENSMUSG00000022952",] > 0)]
lrig1_1 = colnames(exprs(scesetFiltered))[exprs(scesetFiltered)["ENSMUSG00000030029",] > 0]
lrig1_0 = colnames(exprs(scesetFiltered))[!(exprs(scesetFiltered)["ENSMUSG00000030029",] > 0)]
troy_1 = colnames(exprs(scesetFiltered))[exprs(scesetFiltered)["ENSMUSG00000060548",] > 0]
troy_0 = colnames(exprs(scesetFiltered))[!(exprs(scesetFiltered)["ENSMUSG00000060548",] > 0)]
muc5ac_1 = colnames(exprs(scesetFiltered))[exprs(scesetFiltered)["ENSMUSG00000037974",] > 0]
muc5ac_0 = colnames(exprs(scesetFiltered))[!(exprs(scesetFiltered)["ENSMUSG00000037974",] > 0)]
muc6_1 = colnames(exprs(scesetFiltered))[exprs(scesetFiltered)["ENSMUSG00000048191",] > 0]
muc6_0 = colnames(exprs(scesetFiltered))[!(exprs(scesetFiltered)["ENSMUSG00000048191",] > 0)]

genesList = rownames(scesetFiltered[rowMeans(assay(scesetFiltered)) > 10,])


###############################################

#stmn1
x = cor(exprs(scesetFiltered)[genesList,])
corr_all = 1-x[upper.tri(x)]

x = cor(exprs(scesetFiltered)[genesList,stmn1_1])
corr_stmn1_1 = 1-x[upper.tri(x)]

x = cor(exprs(scesetFiltered)[genesList,stmn1_0])
corr_stmn1_0 = 1-x[upper.tri(x)]


sample_name = factor( c(rep("all", length(corr_all)),rep("stmn1_1", length(corr_stmn1_1)), rep("stmn1_0", length(corr_stmn1_0))) , levels = c("all","stmn1_1","stmn1_0"), ordered = TRUE)
corr_value = c(corr_all, corr_stmn1_1, corr_stmn1_0)
df = data.frame( sample = sample_name, value = corr_value)

ggplot(df, aes(x = sample, y = value)) + geom_boxplot()


#mki67
x = cor(exprs(scesetFiltered)[genesList,mki67_1])
corr_mki67_1 = 1-x[upper.tri(x)]

x = cor(exprs(scesetFiltered)[genesList,mki67_0])
corr_mki67_0 = 1-x[upper.tri(x)]


sample_name = factor( c(rep("all", length(corr_all)),rep("mki67_1", length(corr_mki67_1)), rep("mki67_0", length(corr_mki67_0))) , levels = c("all","mki67_1","mki67_0"), ordered = TRUE)
corr_value = c(corr_all, corr_mki67_1, corr_mki67_0)
df = data.frame( sample = sample_name, value = corr_value)

ggplot(df, aes(x = sample, y = value)) + geom_boxplot()


#mist1
x = cor(exprs(scesetFiltered)[genesList,mist1_1])
corr_mist1_1 = 1-x[upper.tri(x)]

x = cor(exprs(scesetFiltered)[genesList,mist1_0])
corr_mist1_0 = 1-x[upper.tri(x)]


sample_name = factor( c(rep("all", length(corr_all)),rep("mist1_1", length(corr_mist1_1)), rep("mist1_0", length(corr_mist1_0))) , levels = c("all","mist1_1","mist1_0"), ordered = TRUE)
corr_value = c(corr_all, corr_mist1_1, corr_mist1_0)
df = data.frame( sample = sample_name, value = corr_value)

ggplot(df, aes(x = sample, y = value)) + geom_boxplot()


#lgr5
x = cor(exprs(scesetFiltered)[genesList,lgr5_1])
corr_lgr5_1 = 1-x[upper.tri(x)]

x = cor(exprs(scesetFiltered)[genesList,lgr5_0])
corr_lgr5_0 = 1-x[upper.tri(x)]


sample_name = factor( c(rep("all", length(corr_all)),rep("lgr5_1", length(corr_lgr5_1)), rep("lgr5_0", length(corr_lgr5_0))) , levels = c("all","lgr5_1","lgr5_0"), ordered = TRUE)
corr_value = c(corr_all, corr_lgr5_1, corr_lgr5_0)
df = data.frame( sample = sample_name, value = corr_value)

ggplot(df, aes(x = sample, y = value)) + geom_boxplot()


#sox2
x = cor(exprs(scesetFiltered)[genesList,sox2_1])
corr_sox2_1 = 1-x[upper.tri(x)]

x = cor(exprs(scesetFiltered)[genesList,sox2_0])
corr_sox2_0 = 1-x[upper.tri(x)]


sample_name = factor( c(rep("all", length(corr_all)),rep("sox2_1", length(corr_sox2_1)), rep("sox2_0", length(corr_sox2_0))) , levels = c("all","sox2_1","sox2_0"), ordered = TRUE)
corr_value = c(corr_all, corr_sox2_1, corr_sox2_0)
df = data.frame( sample = sample_name, value = corr_value)

ggplot(df, aes(x = sample, y = value)) + geom_boxplot()


#runx1
x = cor(exprs(scesetFiltered)[genesList,runx1_1])
corr_runx1_1 = 1-x[upper.tri(x)]

x = cor(exprs(scesetFiltered)[genesList,runx1_0])
corr_runx1_0 = 1-x[upper.tri(x)]


sample_name = factor( c(rep("all", length(corr_all)),rep("runx1_1", length(corr_runx1_1)), rep("runx1_0", length(corr_runx1_0))) , levels = c("all","runx1_1","runx1_0"), ordered = TRUE)
corr_value = c(corr_all, corr_runx1_1, corr_runx1_0)
df = data.frame( sample = sample_name, value = corr_value)

ggplot(df, aes(x = sample, y = value)) + geom_boxplot()


#lrig1
x = cor(exprs(scesetFiltered)[genesList,lrig1_1])
corr_lrig1_1 = 1-x[upper.tri(x)]

x = cor(exprs(scesetFiltered)[genesList,lrig1_0])
corr_lrig1_0 = 1-x[upper.tri(x)]


sample_name = factor( c(rep("all", length(corr_all)),rep("lrig1_1", length(corr_lrig1_1)), rep("lrig1_0", length(corr_lrig1_0))) , levels = c("all","lrig1_1","lrig1_0"), ordered = TRUE)
corr_value = c(corr_all, corr_lrig1_1, corr_lrig1_0)
df = data.frame( sample = sample_name, value = corr_value)

ggplot(df, aes(x = sample, y = value)) + geom_boxplot()


#troy
x = cor(exprs(scesetFiltered)[genesList,troy_1])
corr_troy_1 = 1-x[upper.tri(x)]

x = cor(exprs(scesetFiltered)[genesList,troy_0])
corr_troy_0 = 1-x[upper.tri(x)]


sample_name = factor( c(rep("all", length(corr_all)),rep("troy_1", length(corr_troy_1)), rep("troy_0", length(corr_troy_0))) , levels = c("all","troy_1","troy_0"), ordered = TRUE)
corr_value = c(corr_all, corr_troy_1, corr_troy_0)
df = data.frame( sample = sample_name, value = corr_value)

ggplot(df, aes(x = sample, y = value)) + geom_boxplot()


#muc5ac
x = cor(exprs(scesetFiltered)[genesList,muc5ac_1])
corr_muc5ac_1 = 1-x[upper.tri(x)]

x = cor(exprs(scesetFiltered)[genesList,muc5ac_0])
corr_muc5ac_0 = 1-x[upper.tri(x)]


sample_name = factor( c(rep("all", length(corr_all)),rep("muc5ac_1", length(corr_muc5ac_1)), rep("muc5ac_0", length(corr_muc5ac_0))) , levels = c("all","muc5ac_1","muc5ac_0"), ordered = TRUE)
corr_value = c(corr_all, corr_muc5ac_1, corr_muc5ac_0)
df = data.frame( sample = sample_name, value = corr_value)

ggplot(df, aes(x = sample, y = value)) + geom_boxplot()


#muc6
x = cor(exprs(scesetFiltered)[genesList,muc6_1])
corr_muc6_1 = 1-x[upper.tri(x)]

x = cor(exprs(scesetFiltered)[genesList,muc6_0])
corr_muc6_0 = 1-x[upper.tri(x)]


sample_name = factor( c(rep("all", length(corr_all)),rep("muc6_1", length(corr_muc6_1)), rep("muc6_0", length(corr_muc6_0))) , levels = c("all","muc6_1","muc6_0"), ordered = TRUE)
corr_value = c(corr_all, corr_muc6_1, corr_muc6_0)
df = data.frame( sample = sample_name, value = corr_value)

ggplot(df, aes(x = sample, y = value)) + geom_boxplot()







sample_name = factor( c(
                        rep("stmn1_1", length(corr_stmn1_1)),
                        rep("mki67_1", length(corr_mki67_1)),
                        rep("mist1_1", length(corr_mist1_1)),
                        rep("lgr5_1", length(corr_lgr5_1)),
                        rep("sox2_1", length(corr_sox2_1)),
                        rep("runx1_1", length(corr_runx1_1)),
                        rep("lrig1_1", length(corr_lrig1_1)),
                        rep("troy_1", length(corr_troy_1)),
                        rep("muc5ac_1", length(corr_muc5ac_1)),
                        rep("muc6_1", length(corr_muc6_1))) , 
                      levels = c("all","stmn1_1","mki67_1","mist1_1","lgr5_1",
                                 "sox2_1","runx1_1","lrig1_1","troy_1","muc5ac_1", "muc6_1"),
                      ordered = TRUE)
corr_value = c(corr_stmn1_1,
               corr_mki67_1,
               corr_mist1_1,
               corr_lgr5_1,
               corr_sox2_1,
               corr_runx1_1,
               corr_lrig1_1,
               corr_troy_1,
               corr_muc5ac_1,
               corr_muc6_1)
df = data.frame( sample = sample_name, value = corr_value)
ggplot(df, aes(x = sample, y = value)) + geom_boxplot()

a = aov(value ~ sample, data =df)
summary(a)
