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

setwd("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250415_scz/sample.heatmap/")
wkdir = getwd()

pval_cutoff = 0.05 # The cutoff value for the pvalue we want to use
FC.cutoff = 1.1
# Read the gene data
h1row = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250415_scz/sample.heatmap/gene_average/h1row.complete.RDS")

gene_order = h1row
#plot(gene_order)
mydeepSplit       = FALSE# fine structure within module
myminModuleSize   = 3 # modules must have this minimum number of genes
myheightcutoff = 0.75
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
#gene_order$id = NULL
deg.list = list()
i = "blue"
for ( i in unique(gene_order$gene.module))
{
  tmp1 = gene_order[gene_order$gene.module == i,]
  # remove the version tag
  tmp1$id = gsub("\\..*","",tmp1$id)
  
  deg.list[[i]] = tmp1$id  
}
# deg.list1 = DEGmachine(deg.list, FC.cutoff = 1.01)
# 
# deg.list1 =list()
# deg.list1[["blue.up"]] = rownames(deg.list$blue[deg.list$blue$adj.P.Val < 0.05 & deg.list$blue$logFC>0,])
# deg.list1[["blue.dn"]] = rownames(deg.list$blue[deg.list$blue$adj.P.Val < 0.05 & deg.list$blue$logFC<0,])
# deg.list1[["turquoise.up"]] = rownames(deg.list$turquoise[deg.list$turquoise$adj.P.Val < 0.05 & deg.list$turquoise$logFC>0,])
# deg.list1[["turquoise.dn"]] = rownames(deg.list$turquoise[deg.list$turquoise$adj.P.Val < 0.05 & deg.list$turquoise$logFC<0,])


# # Remove the version tag on he genes
# deg.list1[["blue.up"]] = gsub("\\..*","",deg.list1[["blue.up"]])
# deg.list1[["blue.dn"]] = gsub("\\..*","",deg.list1[["blue.dn"]])
# deg.list1[["turquoise.up"]] = gsub("\\..*","",deg.list1[["turquoise.up"]])
# deg.list1[["turquoise.dn"]] = gsub("\\..*","",deg.list1[["turquoise.dn"]])

# Extract the significant overlaps
# cex.n = 1
# sTest = supertest(deg.list1, degree = 2, n = n.background)
# openImgDev(paste("gene.wina.beta1.supertest.p.",pval_cutoff,"FC",FC.cutoff ,".png", sep = ""),
#            iwidth = 2048, iheight = 800, ipointsize = 40)
# plot.msets(sTest, degree = 2, 
#            sort.by = "size",
#            show.overlap.size = T, 
#            Layout = "landscape",
#            minMinusLog10PValue = 1.3, legend.text.cex =cex.n, overlap.size.cex = cex.n, color.scale.cex = cex.n )
# dev.off()
# The only significant intersections are blue.dn+turquoise.up and blue.up+turquoise.dn
# blue.up.turquoise.dn = intersect(deg.list1$blue.up, deg.list1$turquoise.dn)
# blue.dn.turquoise.up = intersect(deg.list1$blue.dn, deg.list1$turquoise.up)
# deg.list2 = list(blue.up.turquoise.dn=blue.up.turquoise.dn, blue.dn.turquoise.up=blue.dn.turquoise.up)
# deg.list1 = c(deg.list1, deg.list2)

# Make a Venn diagram of the unique genes
###############################################################################  
# ggVennDiagram(deg.list1[1:4],  
#               category.names = "",#c("turq.up", "turq.down", "blue.up", "blue.down"),
#               label_alpha = 0,
#               label = "count", 
#               label_size = 25,
#               set_size = 25) + 
#   scale_fill_gradient(low = "#F4FAFE", high = "#F4FAFE")+
#   # remove legend
#   theme(legend.position = "none")+
#   scale_x_continuous(expand = expansion(mult = .2))
# ggsave(paste("DEGs_Venn.",pval_cutoff,"FC",FC.cutoff,".nonames.png", sep = ""), width = 50, height = 30, units = "cm")


################################################################################
# F U N C T I O N A L    A N A L Y S I S
################################################################################
library(org.Hs.eg.db)
library(GOtest) ##from minghui
library(msigdb) ##from minghui
GOsets = c('c5.go.bp','c5.go.cc', 'c5.go.mf')
gosets_genes = msigdb.genesets(sets=GOsets, type='symbols', species='human',return.data.frame=T)
universe = curated.genesets(c('HGNC_universe'))$Gene
# transform the ENSEMBL gene IDs to symbols for the functional analysis
# library('biomaRt')
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# G_list <- getBM(attributes= c("ensembl_gene_id_version", "hgnc_symbol"), mart= mart)
G_list = read.delim2("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2024/20240305_scz/data/knownGenes.geneid.tsv", header = F)
colnames(G_list) = c("ENSEMBL.ID", "gene.symbol", "chr", "start", "end", "sign", "description")

# Choose a list of gene set names and do functional analysis on them
tmp0 = deg.list
#yp.miR195n = intersect(megaList$list.mir.pos$yellow.pos, megaList$list.mir.neg$mir.neg)
#yn.miR195p = intersect(megaList$list.mir.neg$yellow.neg, megaList$list.mir.pos$mir.pos)

#gene.set = list(yp.miR195n = yp.miR195n, yn.miR195p = yn.miR195p)

i = "blue"
enrichList = list()
for ( i in names(tmp0))
{
  tmp1 = tmp0[[i]]
  print(i)

  if (length(tmp1)>0)
  {
    # make it into a data frame for the same of the function working properly
    tmp1 = data.frame(tmp1, "group")
    # # change ENSEMBL to gene symbol
    #annots <- select(org.Hs.eg.db, keys=tmp1$tmp1, columns="SYMBOL", keytype="ENSEMBL")
    tmp1 = merge(tmp1, G_list, by.x = "tmp1", by.y = "ENSEMBL.ID")
    # tmp1 = merge(tmp1, annots, by.x = "tmp1", by.y = "ENSEMBL")

    tmp1 = tmp1[,c("gene.symbol", "X.group.")]
    colnames(tmp1) = c("gene","group")

    result_weight1 = GOtest(x = tmp1,
                            go = gosets_genes,
                            query.population = universe,
                            background = 'query',
                            #name.x = gsub(" & ", ".", names(tmp)), # This does not like special characters
                            name.x = gsub("-", "_", i), # This does not like special characters
                            name.go = 'GOsets',
                            method = 'hypergeometric',
                            ncores = (detectCores(all.tests = FALSE, logical = TRUE) -1) )

    if(sum(result_weight1$P.adj<pval_cutoff)>0)
    {
      tmp2 = result_weight1[result_weight1$P.adj<pval_cutoff,]
      enrichList[[paste(i,".p.", pval_cutoff,".FC.", FC.cutoff, sep = "")]] = tmp2
    }
  }
}

saveRDS(enrichList, "gene.complete.GOterms.sample.heatmap.RDS")

# go = gene.wina.beta.1.GOtermsFC1p2
# 
# # Get the top 10 enriched GO terms for each module
# top1 <- lapply(go, function(x) x[order(x$P.adj),]$GOsets[1])
# top1 <- as.data.frame(do.call(rbind, top1))
# 
# openxlsx::write.xlsx(top1, "commonMind.scz.GOterms.top1.xlsx", rowNames = T)
