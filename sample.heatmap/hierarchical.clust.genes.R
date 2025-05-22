
rm(list=ls())	#destroys all the objects for a fresh start (better than restarting R)

# ---------------------------------- Parameters to be changed -----------------------------------------

# 1. specify the current working directory
wd = "C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/sample.heatmap/"
setwd(wd)
source("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew3/20220829_data/codes_for_data/R-tomfunctions_static.R")

clustering = "gene"

print(clustering)
outputDir  =paste(clustering,"_average", sep = "")
dir.create(paste(wd,"",outputDir, sep = ""), showWarnings = F)
setwd(paste(wd,"", outputDir, sep = ""))

cormethod="pearson";

# 6. choose a beta value (or power for power func) to derive a scalefree network 
#    if you set this variable as -1, the program will automatically choose a value to derive scale
#     free network. If corrlpower = 1 is chosen then scale free property is not sought after
corrlpower = 1 #

library(WGCNA)
library(lattice)
collect_garbage()

#*-------------------------------------------------------------------------------------
#* STEP 0: read in gene information, expression data and consolidate data
#*
#*
#*-------------------------------------------------------------------------------------
# we have differing numbers of AD samples in the protein, methylation and RNAseq
# data, so we will do an intersection of samples with AD that appear 
# simultaneously in all data sets. We use the standardized (Z transformed) data
datExpr = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/data/UCLA-ASD/bulk_rnaSeq/limma.BA9/ucla.RNAseq.BA9.adj.RDS")
meta = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/clust/wina.BA9.UCLA.bulkRNAseq/gene.beta1/ucla.BA9.metadata.w.subtypes.RDS")

# keep only the SCZ cases
meta = meta[meta$primaryDiagnosis == "Autism Spectrum Disorder",]
# create a common ID to be able to proterly kniw together these DFs
commonID = Reduce(intersect,list(rownames(meta), rownames(datExpr)))
meta = meta[commonID,]
datExpr = datExpr[commonID,]
# Z-transform the data
zmatrix = apply(datExpr, 2, Ztransform, 4) # 1 means the apply function is done on rows and 2 would mean columns
rownames(zmatrix) = rownames(datExpr)
datExpr = zmatrix
dim(datExpr)
rm(zmatrix)
gc()
# Choose randomly only 1% of genes for testing purtposes and to speed up things
#datExpr = datExpr[, sample(colnames(datExpr), ncol(datExpr)*0.01)]

#*-------------------------------------------------------------------------------------
#* compute correlation coefficient 
#*-------------------------------------------------------------------------------------
if (cormethod=="spearman") {
    rankExpr = apply(datExpr, 2, rank)
    corhelp <- cor(rankExpr, use = "pairwise.complete.obs", method="pearson") 
  } else{
    corhelp <- cor(datExpr, use = "pairwise.complete.obs", method=cormethod) 
  }

corhelp <- ifelse(is.na(corhelp), 0, corhelp)
corhelp = (corhelp +1)/2

dichotCor  <- abs(corhelp)^corrlpower
diag(dichotCor)<- 0

dist1 <- 1- abs(dichotCor)

#saveRDS(dist1, "dist1.RDS")

h1row <- hclust(as.dist(dist1), method = "complete")
saveRDS(h1row, "h1row.complete1.RDS")
collect_garbage()

# ----------------- output Hierarchical Clustering image -----------
openImgDev("imgHierClust.single.png",iwidth = 2048, iheight = 800, ipointsize = 8)
par(mfrow=c(1,1))
plot(h1row, labels=F, xlab="",ylab="",main="",sub="")
dev.off()

plot(h1row, labels=F, xlab="",ylab="",main="",sub="")

#myheightcutoff  ~ maximum dendrogam height allowed in modules
#mydeepSplit     ~ FALSE #TRUE will lead to iterative deep cut (resulting to many small modules), FALSE for normal dynamic cut
#myminModuleSize ~ minimal module size
mydeepSplit       = FALSE# fine structure within module
myminModuleSize   = 3 # modules must have this minimum number of genes
myheightcutoff = 0.85
mcolcode2= cutTreeStatic(hiercluster=h1row, heightcutoff=myheightcutoff, minsize1=myminModuleSize)

colcode.reduced  = reassignModuleNames(mcolcode2, minmodulesize=myminModuleSize, anameallmodules=FALSE, 
                    auseblackwhite=FALSE, useNumberAsLabel=FALSE, startlabel=1)
table(colcode.reduced)
plotDendrogramModuleLabels(mdendro=h1row, modulecolors=colcode.reduced, save2file=NULL, plotLabels=FALSE)


################################################################################
# Relabeling
# the colors of my clusters are different from the colors of Ryan's clusters
# We will try and match them here so as to have an easier way of comparing them 
# downstream
NewMatrix = as.data.frame(cbind(mcolcode2, as.character(colcode.reduced)))
colnames(NewMatrix) = c("x", "module")
NewMatrix$x = NULL
saveRDS(NewMatrix, "complete.gene.clustering.RDS")
# OldMatrix = read.delim("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew2/20220816_data/Ryans_results_subtype_definitions_12_2019/meta.MSBB.subtypes.BM36.121919.tsv", sep = "\t")
# outdf = merge(NewMatrix, OldMatrix, by.x = "row.names", by.y = "SynapseId", all.x = TRUE)
# rownames(outdf) = outdf$Row.names
# outdf = outdf[,c("Row.names", "module", "ADsubtype", "ADsubtypeclass")]
# colnames(outdf) = c("SampleID", "MarjanModule", "RyanModule", "RyanModuleClass")
# 
# # this line is creating a problem when comparing the module name to NA so we are
# # going to transform all NAs to 0
# outdf[is.na(outdf)] <- 0
# #overlapTable(outdf$MarjanModule, outdf$RyanModule)
# 
# # we need to do some manual reclustering of the protein clusters 
# # 
# # if (clustering == "prot")
# # {
# #   outdf$MarjanModule = ifelse(outdf$MarjanModule =="yellow" | outdf$MarjanModule =="green", "blue", outdf$MarjanModule)
# # }
# 
# # try and relabel Marjan's clusters to see what could happen then
# outdf$MarjanModuleRelabeled = matchLabels( outdf$MarjanModule, outdf$RyanModule)
# 
# # Marjan's module's relabeled to best match Ryan's modules and then checked how
# # they overlap
# modulesOverlap = overlapTable(outdf$MarjanModuleRelabeled, outdf$RyanModule)
# 
# overlapTable(outdf$MarjanModuleRelabeled, outdf$RyanModule)
# 
# ##########################################
# # projecting Marjan's 5 clusters onto Ryan's 3 subtypes typical, intermediate, atypical
# overlapTable(outdf$MarjanModuleRelabeled, outdf$RyanModuleClass)
# 
# moduleSubtypeOverlap = overlapTable(outdf$MarjanModuleRelabeled, outdf$RyanModuleClass)
# 
# # #regrouping Marjan's 5 clusters into 3 subtypes based on dendrogram
# # outdf$MarjanModuleClass = ifelse(outdf$MarjanModuleRelabeled=="turquoise" , "typical", outdf$MarjanModuleRelabeled)
# # outdf$MarjanModuleClass = ifelse(outdf$MarjanModuleClass=="blue" | outdf$MarjanModuleRelabeled=="yellow" , "atypical", outdf$MarjanModuleClass)
# # outdf$MarjanModuleClass = ifelse(outdf$MarjanModuleClass=="magenta" | outdf$MarjanModuleClass=="brown", "intermediate", outdf$MarjanModuleClass)
# 
# if(clustering == "rna"){
# outdf$MarjanModuleClass = ifelse(outdf$MarjanModuleRelabeled=="turquoise" , "typical", outdf$MarjanModuleRelabeled)
# outdf$MarjanModuleClass = ifelse(outdf$MarjanModuleClass=="blue", "atypical", outdf$MarjanModuleClass)
# outdf$MarjanModuleClass = ifelse(outdf$MarjanModuleClass=="yellow" | outdf$MarjanModuleClass=="brown", "intermediate", outdf$MarjanModuleClass)
# }
# 
# if (clustering =="prot")
# {
#   outdf$MarjanModuleClass = ifelse(outdf$MarjanModuleRelabeled=="turquoise" , "typical", outdf$MarjanModuleRelabeled)
#   outdf$MarjanModuleClass = ifelse(outdf$MarjanModuleClass=="yellow", "atypical", outdf$MarjanModuleClass)
#   outdf$MarjanModuleClass = ifelse(outdf$MarjanModuleClass=="green", "atypical", outdf$MarjanModuleClass)
# }
# 
# overlapTable( outdf$MarjanModuleClass, outdf$RyanModuleClass)
# 
# subtypeSubtypeOverlap = overlapTable(outdf$MarjanModuleClass, outdf$RyanModuleClass)
# 
# 
# 
# plotDendrogramModuleLabels(mdendro=h1row, modulecolors=outdf$MarjanModuleRelabeled, save2file=NULL, plotLabels=FALSE)
# 
# mcex =  0.05 + 0.5/log10(ncol(dist1))
# openImgDev("imgCorHeatMap.png")
# heatmapColorRG = ( rgcolors.func(50) )
# heatmapColor = heat.colors(50)
# par(mfrow=c(1,1))
# diag(dist1) = 0.5
# heatmap(dist1,
#         Rowv=as.dendrogram(h1row),
#         Colv="Rowv", 
#         scale="none",
#         revC=T,
#         ColSideColors=outdf$MarjanModuleRelabeled,
#         RowSideColors=outdf$MarjanModuleRelabeled,
#         cexRow = mcex, 
#         cexCol = mcex,
#         col=heatmapColor)
# dev.off()
# 
# save(list = ls(), file = paste(clustering,"_clustering_average.Rdata", sep = ""))
# #}
# 
