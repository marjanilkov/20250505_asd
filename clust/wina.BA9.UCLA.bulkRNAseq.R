# ## R version 3.5.3-3.6.3
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.10")
# install.packages(c('RcppEigen', 'Rcpp'))
# install.packages(c("http://cran.nexr.com/src/contrib/Rclusterpp_0.2.3.tar.gz"), repos=NULL, type="source")
# BiocManager::install(c("impute", "pcaMethods"))
# install.packages(c('dynamicTreeCut', 'flashClust', 'foreach'))
# install.packages("/sc/arion/projects/zhangb03a/neffr01/transition/subtyping_AD/MSBB_RNA_CDRadj/tools/WINA-ryan-0.1.5.2.tar.gz", type="source",repos=NULL)

rm(list = ls())
### input params ###
library("WINA")
library("jsonlite")
library(WGCNA)
setwd("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/clust/")
wkdir = getwd()
source("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew3/20220829_data/codes_for_data/R-tomfunctions_static.R")

num_cores = 7
beta = 1
RsquaredCut = 0.8
pcheicutoff = 0.20
minheight = 0.
myminModuleSize = 1
clustering = "gene"
out_dir = paste(wkdir,"/wina.BA9.UCLA.bulkRNAseq/",sep="")
dir.create(out_dir, showWarnings = F)
setwd(out_dir)

out_dir1 = paste(out_dir,clustering,".beta",beta,"/",sep="")
dir.create(out_dir1, showWarnings = F)

setwd(out_dir1)
#*-------------------------------------------------------------------------------------
#* STEP 0: read in gene information, expression data and consolidate data
#*
#*
#*-------------------------------------------------------------------------------------
# we have differing numbers of AD samples in the protein, methylation and RNAseq
# data, so we will do an intersection of samples with AD that appear 
# simultaneously in all data sets. We use the standardized (Z transformed) data
datExpr = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/data/UCLA-ASD/limma.BA9/ucla.RNAseq.BA9.adj.RDS")
#meta = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250415_scz/data/Metadata/METADATA.RDS")
meta1 = read.delim2("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/data/UCLA-ASD/Metadata/individual_UCLA-ASD_metadata.csv", sep = ",")
meta2 = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/data/UCLA-ASD/limma.BA9/sampleID.to.individualID.RDS")
meta2 = meta2[ grepl("ba9", meta2$ba, ignore.case = T),]

meta = merge(meta1, meta2, by.x = "individualID", by.y = "individualID")
# replace all _, - with .
meta$sample = gsub("-", ".", meta$sample)
meta$sample = gsub("_", ".", meta$sample)
# remove NAs in meta$sample and diagnosis
meta = meta[!is.na(meta$sample),]
meta = meta[!is.na(meta$primaryDiagnosis),]

rownames(meta) = meta$sample
rm(meta1, meta2)

# meta$PMI has 5 missing values. I will replace them with the mean.
# Apparently the PMI column is character
meta$PMI = as.numeric(meta$PMI)
meta$ageDeath = as.numeric(meta$ageDeath)

meta$PMI = ifelse(is.na(meta$PMI), mean(meta$PMI, na.rm = T), meta$PMI)

# Z-transform the data
zmatrix = apply(datExpr, 2, Ztransform, 4) # 1 means the apply function is done on rows and 2 would mean columns
rownames(zmatrix) = rownames(datExpr)
datExpr = zmatrix
dim(datExpr)
rm(zmatrix)
gc()


# Use only the SCZ patients
meta_ASD = meta[meta$primaryDiagnosis == "Autism Spectrum Disorder",]
datExpr = datExpr[rownames(datExpr) %in% meta_ASD$sample,]
dim(datExpr)

print("oneplus=TRUE")
print("running WINA. here we GO...")
result_true = wina(datExpr,
              headCol=NA,
              outputDir=out_dir1,
              fname='WINA',
              beta=beta,
              RsquaredCut=RsquaredCut,
              linkage="average",
              pcheicutoff=pcheicutoff,
              cormethod=c("pearson"),
              myminModuleSize=myminModuleSize,
              myheightcutoff=minheight,
              imagetype="png",
              gene.missing.rate.cutoff=0.5,
              sd.cutoff=NULL,
              sample.missing.rate.cutoff=0.5,
              impute=c('mean'),
              compute.connectivity.statistics=TRUE,
              plot.heatmap=TRUE,
              heatmap.downsample.size=5000,
              heatmap.enhancefold=6,
              heatmap.useRaster=T,
              ncores=num_cores,
              oneplus=TRUE)



# plotDendrogramModuleLabels(mdendro=h1row, modulecolors=colcode.reduced, save2file=NULL, plotLabels=T)
# 
# print("saving Rdata...")
# save.image(file=paste(out_dir1,"WINA_run.Rdata",sep=""))
# print("saved.")


#Use Static TreeCut  instead of dynamic
load(paste(out_dir1,"WINA_h1row_dendrogram.Rdata", sep = ""))
load(paste(out_dir1, "WINA_colcode-reduced_dendrogram.Rdata", sep = ""))
plot(h1row)
mydeepSplit       = FALSE# fine structure within module
myminModuleSize   = 3 # modules must have this minimum number of genes
myheightcutoff = 0.501
mcolcode2= cutTreeStatic(hiercluster=h1row, heightcutoff=myheightcutoff, minsize1=myminModuleSize)
#
colcode.reduced  = reassignModuleNames(mcolcode2, 
                                       minmodulesize=myminModuleSize, 
                                       anameallmodules=FALSE,
                                       auseblackwhite=FALSE,
                                       useNumberAsLabel=FALSE, 
                                       startlabel=1)
table(colcode.reduced)

png(file="dendro.png",width=600, height=350, units="mm", res=300, pointsize = 14)
plotDendrogramModuleLabels(mdendro=h1row, modulecolors=colcode.reduced, save2file=NULL, plotLabels=FALSE)
dev.off()

# we need to have the datExpr again to find which samples come from control and which from AD
# geneDF = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2024/20240305_scz_asd_ad/cm/scz/Release4/RNAseq.CommonMind.Msbb.Pitt.Penn.adj.std.RDS")
# meta = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2024/20240305_scz_asd_ad/cm/scz/Release4/Metadata/METADATA.RDS")

# if (clustering=="gene"){datExpr = as.data.frame(t(geneDF))}
# rm(geneDF)
# gc()

classification = as.data.frame(cbind(h1row$labels, as.character(colcode.reduced)))
colnames(classification) = c("id", "cluster")
rownames(classification) = classification$id
# Now do visual inspection of the heatmap and dendrogram and do any merging if needed
#classification$cluster = ifelse(classification$cluster=="grey", "blue", classification$cluster)
#classification$cluster = ifelse(classification$cluster=="brown", "blue", classification$cluster)
tmp1 = datExpr[,1:2]
classification = merge(classification, tmp1, by.x="id", by.y="row.names", all.y = T)
# remove the two unnecessary columns
classification = classification[, 1:2]
rownames(classification) = classification$id

# Add the controls and the MCI
#metaMSSM = meta[meta$Cohort == "MSSM-Penn-Pitt",]
classification = merge(meta, classification, by.x = "sample", by.y = "id", all.x = T)
classification$cluster = ifelse(classification$primaryDiagnosis == "control", "control", classification$cluster)
table(classification$cluster)

# check if the clusters overlap with Ryan's original clusters and relabel them accordingly
# overlapTable(classification$ADsubtype, classification$cluster)
# classification$cluster = ifelse(classification$cluster=="blue", "C1", classification$cluster)
# classification$cluster = ifelse(classification$cluster=="turquoise", "A", classification$cluster)
# classification = classification[,c("SynapseId", "cluster")]
# classification = classification[complete.cases(classification$cluster),]
# rownames(classification)=classification$SynapseId
colnames(classification)[ncol(classification)] = "mrna.subtype.b1"
saveRDS(classification, "ucla.BA9.metadata.w.subtypes.RDS")
