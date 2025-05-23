rm(list = ls())

options(stringsAsFactors=F)
library(DESeq2)
library(limma)
library(edgeR)
library(qvalue)
library(SuperExactTest)
set.seed(12345)
library(GOtest) ##from minghui
library(msigdb) ##from minghui
library(VennDiagram)
library(ggVennDiagram)
library(ggplot2)
library(scales)
# We will load the marjanRfunctions.R set of functions directly from github from now on
devtools::source_url("https://raw.githubusercontent.com/marjanilkov/sinaiRepo/refs/heads/main/marjanRfunctions.R")
################################################################################
# Variables set by the user
setwd("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/clust/wina.BA9.UCLA.bulkRNAseq/gene.beta1/")

p.val_cutoff =0.05
FC_cutoff = 1.5
tissue = "DLPFC"
################################################################################
# Read in the data
datExpr <- readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/data/UCLA-ASD/bulk_rnaSeq/limma.BA9/ucla.RNAseq.BA9.adj.RDS")
classification = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/clust/wina.BA9.UCLA.bulkRNAseq/gene.beta1/ucla.BA9.metadata.w.subtypes.RDS")

# have the name of the column with clusters called "cluster"
colnames(classification)[ncol(classification)] = "cluster"
modSize = ncol(datExpr)
datExpr = datExpr[complete.cases(datExpr),]

# create a common ID to be able to properly knit together these DFs
commonID = Reduce(intersect,list(rownames(classification), rownames(datExpr)))
# put everything together 
classification = classification[rownames(classification) %in% commonID,]
datExpr = datExpr[rownames(datExpr) %in% commonID,]

################################################################################
# DEA clusters
clusters = unique(classification$cluster)
design = model.matrix(~0+cluster, data=classification) 
colnames(design) = gsub("cluster", "", colnames(design))
fit = lmFit(as.data.frame(t(datExpr)), design) #fit all groups together!

# We will automate the makeContrasts function. Have the control cluster always 
# at the end of the vector
tmpClust = colnames(design)
tmpClust = tmpClust[!tmpClust %in% c("control")]
tmpClust = c(tmpClust, "control")

# make a combo of all clusters
tmpDesign = t(combn(tmpClust, 2))

# autoamtically create a text for the makeContrasts function
contrastForm = NULL
topTableForm = NULL
degForm = NULL

i=1

for (i in 1:nrow(tmpDesign))
{
  print(i)
  contrastForm = paste(contrastForm, tmpDesign[i,1],".",tmpDesign[i,2], " = ", tmpDesign[i,1]," - ",tmpDesign[i,2],",",sep = "")
  topTableForm = paste(topTableForm, tmpDesign[i,1],".",tmpDesign[i,2], " = topTable(efit, coef = '", tmpDesign[i,1],".",tmpDesign[i,2], "', number = Inf );",sep = "")
  degForm = paste(degForm, tmpDesign[i,1],".",tmpDesign[i,2] , "=", tmpDesign[i,1],".",tmpDesign[i,2] , ",", sep = "")
}

contrastForm = paste("contr <- makeContrasts(",contrastForm,"levels = colnames(coef(fit)))", sep = "")
eval(parse(text = contrastForm))

vfit <- contrasts.fit(fit, contr)
efit <- eBayes(vfit)

summary(decideTests(efit, p.value = p.val_cutoff, lfc = log2(FC_cutoff)))
# autoamtically execute created formulas
eval(parse(text = topTableForm))

degForm = paste("deg = list(", degForm, ")", sep = "")
# remove the last comma
library(stringi)
degForm = stri_replace_last(degForm, fixed = ',', '')
eval(parse(text = degForm))

# To reduce the names of deg elements we will reduce the name TURQUOISE to TURQ
# and control to ctrl
names(deg) = gsub("turquoise", "turq", names(deg))
names(deg) = gsub("control", "ctrl", names(deg))
################################################################################
# We will add also the DEGs from the overall cases VS controls
classification$case.ctrl = ifelse(classification$cluster != "control", "case", classification$cluster)
clusters = unique(classification$case.ctrl)
design = model.matrix(~0+case.ctrl, data=classification) 
colnames(design) = gsub("case.ctrl", "", colnames(design))
fit = lmFit(as.data.frame(t(datExpr)), design) #fit all groups together!

# We will automate the makeContrasts function. Have the control cluster always 
# at the end of the vector
tmpClust = colnames(design)
tmpClust = tmpClust[!tmpClust %in% c("control")]
tmpClust = c(tmpClust, "control")

# make a combo of all clusters
tmpDesign = t(combn(tmpClust, 2))

# autoamtically create a text for the makeContrasts function
contrastForm = NULL
topTableForm = NULL
degForm = NULL

i=1

for (i in 1:nrow(tmpDesign))
{
  print(i)
  contrastForm = paste(contrastForm, tmpDesign[i,1],".",tmpDesign[i,2], " = ", tmpDesign[i,1]," - ",tmpDesign[i,2],",",sep = "")
  topTableForm = paste(topTableForm, tmpDesign[i,1],".",tmpDesign[i,2], " = topTable(efit, coef = '", tmpDesign[i,1],".",tmpDesign[i,2], "', number = Inf );",sep = "")
  degForm = paste(degForm, tmpDesign[i,1],".",tmpDesign[i,2] , "=", tmpDesign[i,1],".",tmpDesign[i,2] , ",", sep = "")
}

contrastForm = paste("contr <- makeContrasts(",contrastForm,"levels = colnames(coef(fit)))", sep = "")
eval(parse(text = contrastForm))

vfit <- contrasts.fit(fit, contr)
efit <- eBayes(vfit)

summary(decideTests(efit, p.value = p.val_cutoff, lfc = log2(FC_cutoff)))
# autoamtically execute created formulas
eval(parse(text = topTableForm))

degForm = paste("deg.case.ctrl = list(", degForm, ")", sep = "")
# remove the last comma
library(stringi)
degForm = stri_replace_last(degForm, fixed = ',', '')
eval(parse(text = degForm))

# To reduce the names of deg elements we will reduce the name TURQUOISE to TURQ
# and control to ctrl
names(deg.case.ctrl) = gsub("control", "ctrl", names(deg.case.ctrl))

deg = c(deg, deg.case.ctrl)
#saveRDS(deg, paste(tissue,".DEG.RDS", sep = ""))
degList = DEGmachine(deg[2:3], FC.cutoff = FC_cutoff)
degList.case.ctrl = DEGmachine(deg.case.ctrl, FC.cutoff = FC_cutoff)
degList.case.ctrl = degList.case.ctrl[1:2]

degList = c(degList.case.ctrl, degList)
################################################################################
# F U N C T I O N A L    A N A L Y S I S 
################################################################################

library(org.Hs.eg.db)
library(GOtest) ##from minghui
library(msigdb) ##from minghui
GOsets = c('c5.go.bp','c5.go.cc', 'c5.go.mf')
gosets_genes = msigdb.genesets(sets=GOsets, type='symbols', species='human',return.data.frame=T)
universe = curated.genesets(c('HGNC_universe'))$Gene

# The names of genes now are in ENSEMBL, we will transform them to symbols
mapping <- readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/data/UCLA-ASD/bulk_rnaSeq/UCLA_ASD_bulk_rnaSeq.RDS")
mapping = mapping[,1:2]

enrichList = list()

tmp0 = degList

i = 1

for ( i in 1:length(tmp0))
{
  tmp1 = tmp0[[i]]
  print(i)
  # The gene names are with version at the end. Remove that (everything after the fullstop)
  tmp1 = gsub("\\..*","",tmp1)
  if (length(tmp1)>0)
  {
    # make it into a data frame for the same of the function working properly
    tmp1 = data.frame(tmp1, "group")
    tmp1 = merge(mapping, tmp1, by.x = "Geneid", by.y = "tmp1")
    
    tmp1 = tmp1[,c("Symbol", "X.group.")]
    colnames(tmp1) = c("gene","group")
    
    result_weight1 = GOtest(x = tmp1,
                            go = gosets_genes,
                            query.population = universe,
                            background = 'query',
                            #name.x = gsub(" & ", ".", names(tmp)), # This does not like special characters
                            name.x = names(tmp0)[i], # This does not like special characters
                            name.go = 'GOsets',
                            method = 'hypergeometric',
                            ncores = (detectCores(all.tests = FALSE, logical = TRUE) -1) )
    
    if(sum(result_weight1$P.adj<p.val_cutoff)>0)
    {
      tmp2 = result_weight1[result_weight1$Pvalue<p.val_cutoff,]
      enrichList[[paste(names(tmp0)[i], sep = "")]] = tmp2
    }
  }
}

saveRDS(enrichList, paste(tissue, ".FC",FC_cutoff, ".GOterms.RDS", sep = ""))
