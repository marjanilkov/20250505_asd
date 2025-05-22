rm(list = ls())

library(ggplot2)
library(tidyverse)
library(WGCNA)
library(SuperExactTest)
library(lattice)
library(ggVennDiagram)

source("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/sinaiRepo/R-tomfunctions_static.R")
# We will load the marjanRfunctions.R set of functions directly from github from now on
devtools::source_url("https://raw.githubusercontent.com/marjanilkov/sinaiRepo/refs/heads/main/marjanRfunctions.R")

setwd("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/sample.heatmap/")
wkdir = getwd()

datExpr = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/data/UCLA-ASD/bulk_rnaSeq/limma.BA9/ucla.RNAseq.BA9.adj.RDS")
meta = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/clust/wina.BA9.UCLA.bulkRNAseq/gene.beta1/ucla.BA9.metadata.w.subtypes.RDS")

# Z-transform the data
zmatrix = apply(datExpr, 2, Ztransform, 4) # 1 means the apply function is done on rows and 2 would mean columns
rownames(zmatrix) = rownames(datExpr)
datExpr = zmatrix
dim(datExpr)
rm(zmatrix)
gc()

# create a common ID to be able to proterly kniw together these DFs
commonID = Reduce(intersect,list(rownames(meta), rownames(datExpr)))
meta = meta[commonID,]
datExpr = datExpr[commonID,]

pval_cutoff = 0.05 # The cutoff value for the pvalue we want to use
FC.cutoff = 1.1
# Read the gene data
h1row = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/sample.heatmap/gene_average/h1row.complete1.RDS")

gene_order = h1row
#plot(gene_order)
mydeepSplit       = FALSE# fine structure within module
myminModuleSize   = 3 # modules must have this minimum number of genes
myheightcutoff = 0.85
mcolcode2= cutTreeStatic(hiercluster=h1row, heightcutoff=myheightcutoff, minsize1=myminModuleSize)
#
colcode.reduced  = reassignModuleNames(mcolcode2, 
                                       minmodulesize=myminModuleSize, 
                                       anameallmodules=FALSE,
                                       auseblackwhite=FALSE,
                                       useNumberAsLabel=FALSE, 
                                       startlabel=1)
table(colcode.reduced)
plotDendrogramModuleLabels(mdendro=h1row, modulecolors=colcode.reduced, save2file=NULL, plotLabels=FALSE)

classification = as.data.frame(cbind(h1row$labels, as.character(colcode.reduced)))
colnames(classification) = c("id", "gene.module")
rownames(classification) = classification$id
gene_order = classification[order(classification$gene.module),]

i = "black"
j = "turquoise"

# let's try and use the absolute value of the datExpr for testing
#datExpr = abs(datExpr)

datExpr.mean = data.frame(matrix(ncol = 4, nrow = 0))
colnames(datExpr.mean) = c("gene.cluster", "subtype", "mean", "std")
for ( i in unique(gene_order$gene.module))
{
  tmp1.geneID = gene_order[gene_order$gene.module == i,]$id # extract the genes names in one gene cluster 
  comparing.df = list()
  for (j in unique(meta$cluster))
  {
    tmp2.sample = meta[meta$cluster == j,]$sample # extract the sample IDs for that particular SCZ subtype
    tmp3.df = reshape2::melt(datExpr[tmp2.sample, tmp1.geneID]) # extract the datExpr for the given subtype and gene cluster
    tmp3.df$gene.clust = i
    tmp3.df$sample.clust = j
    tmp.datExpr.mean = data.frame(i, j, mean(tmp3.df$value), sd(tmp3.df$value))
    datExpr.mean = rbind(datExpr.mean, tmp.datExpr.mean)
    # Find if there are significant differences between subtypes within each gene cluster
    comparing.df[[j]] = tmp3.df
  }
  
}
colnames(datExpr.mean) = c("gene.cluster", "subtype", "mean", "std")

ggplot(datExpr.mean, aes(x=gene.cluster, y=mean, fill=subtype)) + 
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=mean-std, ymax=mean+std), width=.2,position=position_dodge(.9))+ 
  scale_fill_manual("legend", values = c("blue" = "blue", "control" = "grey", "turquoise" = "turquoise"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20))
# ggsave(paste("datExpr.heatmap.barplot.png", sep = ""), 
#        width = 25, 
#        height = 20, 
#        units = "cm",
#        limitsize = F)

# Compare means
tmp.datExpr = t(datExpr)
tmp.datExpr = merge(tmp.datExpr, gene_order, by.x = "row.names", by.y = "id")
rownames(tmp.datExpr) = tmp.datExpr$Row.names
tmp.datExpr$Row.names = NULL

tmp.df = tmp.datExpr[tmp.datExpr$gene.module == "turquoise",]
tmp.df$gene.module = NULL

tmp.df = t(tmp.df)
tmp.meta = meta[, c("sample", "cluster"),]
tmp.df = merge(tmp.df, tmp.meta, by.x = "row.names", by.y = "sample" )
rownames(tmp.df) = tmp.df$Row.names
tmp.df$Row.names = NULL
library(tidyr)
# Reshape to long format using pivot_longer
tmp.df_long <- pivot_longer(
  tmp.df,
  cols = 1:(ncol(tmp.df)-1),
  names_to = "rnaseq",
  values_to = "expression"
)



library(ggplot2)




# Most basic error bar
#tmp.df$cluster = NULL
# ggplot(tmp.df) +
#   geom_bar( aes(x=cluster, y=expression), stat="identity", fill="skyblue", alpha=0.7) #+
  #geom_errorbar( aes(x=name, ymin=value-sd, ymax=value+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)
# library(ggplot2)
# library(dplyr)
library(hrbrthemes)
# data <- data.frame(
#   type = c( rep("variable 1", 1000), rep("variable 2", 1000) ),
#   value = c( rnorm(1000), rnorm(1000, mean=4) )
# )

# Represent it
tmp.df_long %>%
  ggplot( aes(x=expression, fill=cluster)) +
  geom_density( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("blue","grey", "turquoise")) +
  theme_ipsum() +
  labs(fill="")
# p
# 
# tmp.df_long$cluster = as.factor(tmp.df_long$cluster)
# compare_means(expression ~ cluster,  
#               data = tmp.df_long, 
#               p.adjust.method = "BH", method = "wilcox")
# library(dplyr)
# group_by(tmp.df_long, cluster) %>%
#   summarise(
#     count = n(),
#     mean = mean(expression, na.rm = TRUE),
#     sd = sd(expression, na.rm = TRUE)
#  )
library(ggpubr)
# Box plots
# ++++++++++++++++++++
# Plot weight by group and color by group
ggboxplot(tmp.df_long, 
          x = "cluster", 
          y = "expression", 
          bxp.errorbar = T, 
          #add = "jitter",
          color = "cluster", palette = c("grey", "blue", "turquoise"),
          order = c("control", "blue", "turquoise"),
          ylab = "Weight", xlab = "Treatment")
xxx = tmp.df_long[tmp.df_long$cluster == "blue" |tmp.df_long$cluster == "turquoise",]

t.test(expression ~ cluster, data = xxx, var.equal = F)
compare_means(expression ~ cluster, data = tmp.df_long, method = "t.test", paired = F)

font.size = 20
ggplot(tmp.df_long, aes(cluster, expression, colour = cluster)) + 
  scale_color_manual(values = c("blue","black", "turquoise"))+
  geom_violin(linewidth = 1) +
  geom_boxplot(width=0.1, linewidth = 1)+ # Add median and quartiles
  guides(fill = "none", color = "none", linetype = "none", shape = "none")+
  theme_classic() +
  xlab("SCZ subtype")+
  #ylab(aes.var)+
  theme( axis.text.x = element_text(size=font.size, colour = "black"),
         axis.text.y = element_text(size=font.size, colour = "black"),
         text = element_text(size=font.size))
  #geom_errorbar( aes(x=name, ymin=value-sd, ymax=value+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)
  
# Mean plots
# ++++++++++++++++++++
# Plot weight by group
# Add error bars: mean_se
# (other values include: mean_sd, mean_ci, median_iqr, ....)
# ggline(tmp.df_long, x = "cluster", y = "expression", 
#        add = c("mean_se", "jitter"), 
#        order = c("control", "blue", "turquoise"),
#        ylab = "Weight", xlab = "Treatment")

# Box plot
boxplot(expression ~ cluster, data = tmp.df_long,
        xlab = "Treatment", ylab = "Weight",
        frame = FALSE, col = c( "blue", "grey","turquoise"))
# plotmeans
library("gplots")
plotmeans(expression ~ cluster, data = tmp.df_long[c(1:50, 4500:4550, 9000:9050),], frame = FALSE,
          xlab = "Treatment", ylab = "Weight",
          main="Mean Plot with 95% CI") 
# Compute the analysis of variance
res.aov <- aov(expression ~ cluster, data = tmp.df_long[1:100,])
# Summary of the analysis
summary(res.aov)
plot(res.aov)




ggboxplot(tmp.df_long, x = "cluster", y = "expression", 
          color = "cluster", palette = c("grey", "turquoise", "blue"),
          ylab = "expression", xlab = "clusters")
res <- aov(expression ~ cluster,  
                   data = tmp.df_long)
res
summary(res)
#######################
# Testing, nothing of importance
s = 100
xx = sample(x = c("a", "b"), 
       prob = c(.5, .5),
       size = s, 
       replace = TRUE)
random_numbers <- runif(s)
df = data.frame(xx, random_numbers)
ggboxplot(df, 
          x = "xx", 
          y = "random_numbers", 
          bxp.errorbar = T, 
          palette = c("blue", "turquoise"),
          ylab = "Weight", xlab = "Treatment")

t.test(random_numbers ~ xx, data = df, var.equal = F)
compare_means(random_numbers ~ xx, data = df, method = "wilcox", paired = F)

x <- rnorm(12)
z.test(x,sigma.x=1)
# Two-sided one-sample z-test where the assumed value for
# sigma.x is one. The null hypothesis is that the population
# mean for 'x' is zero. The alternative hypothesis states
# that it is either greater or less than zero. A confidence
# interval for the population mean will be computed.

x <- c(7.8, 6.6, 6.5, 7.4, 7.3, 7., 6.4, 7.1, 6.7, 7.6, 6.8)
y <- c(4.5, 5.4, 6.1, 6.1, 5.4, 5., 4.1, 5.5)
x = tmp.df_long[tmp.df_long$cluster == "turquoise",]$expression
y = tmp.df_long[tmp.df_long$cluster == "blue",]$expression

z.test(x, sigma.x=1, y, sigma.y=1)
# Two-sided standard two-sample z-test where both sigma.x
# and sigma.y are both assumed to equal 0.5. The null hypothesis
# is that the population mean for 'x' less that for 'y' is 2.
# The alternative hypothesis is that this difference is not 2.
# A confidence interval for the true difference will be computed.

z.test(x, sigma.x=0.5, y, sigma.y=0.5, conf.level=0.90)
# Two-sided standard two-sample z-test where both sigma.x and
# sigma.y are both assumed to equal 0.5. The null hypothesis
# is that the population mean for 'x' less that for 'y' is zero.
# The alternative hypothesis is that this difference is not
# zero.  A 90% confidence interval for the true difference will
# be computed.
rm(x, y)


# library(gplots)
# show comparison with boxplot
data(state)
plotmeans(state.area ~ state.region)
# show some color and mean labels
plotmeans(state.area ~ state.region,
          mean.labels=TRUE, digits=-3,
          col="red", connect=FALSE)
# show how to specify which means should be connected
plotmeans(state.area ~ state.region, connect=list(1:2, 3:4),
          ccol="red", pch=7 )
# more complicated example showing how to show an interaction
data(esoph)
par(las=2, # use perpendicular axis labels
    mar=c(10.1,4.1,4.1,2.1), # create enough space for long x labels
    mgp=c(8,1,0) # move x axis legend down to avoid overlap
)
plotmeans(ncases/ncontrols ~ interaction(agegp , alcgp, sep =" "),
          connect=list(1:6,7:12,13:18,19:24),
          barwidth=2,
          col="dark green",
          data=esoph,
          xlab="Age Group and Alcohol Consumption",
          ylab="# Cases / # Controls",
          ylim = c(-.9,1.4),
          main=c("Fraction of Cases for by Age and Alcohol Consumption",
                 "Ile-et-Vilaine Esophageal Cancer Study")
)
abline(v=c(6.5, 12.5, 18.5), lty=2)