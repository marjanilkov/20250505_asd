rm(list = ls())

library(ggplot2)
library(tidyverse)
library(WGCNA)
library(SuperExactTest)
library(lattice)
library(ggVennDiagram)

source("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/sinaiRepo/R-tomfunctions_static.R")
source("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/sinaiRepo/marjanRfunctions.R")

setwd("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/clust/wina.BA9.UCLA.bulkRNAseq/gene.beta1/")
wkdir = getwd()

pval_cutoff = 0.05 # The cutoff value for the pvalue we want to use
FC.cutoff = 1.3
# Read the differential expression data
deg = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/clust/wina.BA9.UCLA.bulkRNAseq/gene.beta1/DLPFCUCLA.b1.DEP.RDS")
degList = DEGmachine(deg, FC.cutoff = FC.cutoff, p.val.cutoff = pval_cutoff) 
n.background = nrow(deg[[1]])
# Extract the significant overlaps
cex.n = 1
deg.list.for.sTest = degList[1:4]
sTest = supertest(deg.list.for.sTest, degree = 2, n = n.background)
sTest = SuperTestOnlySignificant(sTest = sTest)

openImgDev(paste("gene.wina.beta1.supertest.p.",pval_cutoff,"FC",FC.cutoff ,".png", sep = ""),
           iwidth = 2048, iheight = 800, ipointsize = 40)
plot.msets(sTest, degree = 2, 
           sort.by = "size",
           show.overlap.size = T, 
           Layout = "landscape",
           minMinusLog10PValue = 1.3, legend.text.cex =cex.n, overlap.size.cex = cex.n, color.scale.cex = cex.n )
dev.off()

# Make a Venn diagram of the unique genes
###############################################################################  
set.Names = "y" # whether to include the set names or not. options are y/n
################################################################################

degs = deg.list.for.sTest

if (set.Names == "y"){NN = names(degs)} else{NN = ""}

ggVennDiagram(degs,  
              category.names = NN,
              label_alpha = 0,
              label = "count", 
              label_size = 15,
              set_size = 30) + 
  theme(text = element_text(size = 30),
        legend.title = element_text(size=50),
        legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(2, 'cm'),
        legend.position = "none",
        legend.text = element_text(colour="black", size=40),
        plot.title = element_text(size = 60, face = "bold"))+
  scale_fill_gradient(low = "#F4FAFE", high = "tan")+
  scale_x_continuous(expand = expansion(mult = .23))

ggsave(paste("Venn.scz.sub.",ifelse(set.Names=="y", "YESnames", "NOnames" ),".png", sep = ""), width = 50, height = 40, units = "cm",limitsize = F)

