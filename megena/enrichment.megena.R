rm(list = ls())

library(GeneOverlap)
library(MEGENA)
library(ggplot2)

source("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew3/20220829_data/codes_for_data/R-tomfunctions_static.R")
source("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/sinaiRepo/marjanRfunctions.R")

datExpr = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/data/UCLA-ASD/bulk_rnaSeq/limma.BA9/ucla.RNAseq.BA9.adj.RDS")

wd = ("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/megena/")
setwd(wd)
FC_cutoff = 1.3

########################################################################################
# For SCZvsCTRL DEGs
# # If we use FC=0
# degs <- readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2024/20240305_scz/wina.clust/MPP/initial_cluster_MPP/gene.beta1/SCZvsCTRL.MPP.DEGs.RDS")
# degs1 =list()
# degs1[["scz.up"]] = rownames(degs$scz[degs$scz$adj.P.Val < 0.05 & degs$scz$logFC>0,])
# degs1[["scz.dn"]] = rownames(degs$scz[degs$scz$adj.P.Val < 0.05 & degs$scz$logFC<0,])
# # Remove the version tag on he genes
# degs1[["scz.up"]] = gsub("\\..*","",degs1[["scz.up"]])
# degs1[["scz.dn"]] = gsub("\\..*","",degs1[["scz.dn"]])
########################################################################################

########################################################################################
# FOr SCZsubtypevsCTRL DEGs
degs <- readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/clust/wina.BA9.UCLA.bulkRNAseq/gene.beta1/DLPFC.DEG.RDS")

# remove the DEGs without ctrl contrast
degs = degs[2:4]
# If we are using FC=0
# degs1 =list()
# degs1[["blue.up"]] = rownames(degs$blue[degs$blue$adj.P.Val < 0.05 & degs$blue$logFC>0,])
# degs1[["blue.dn"]] = rownames(degs$blue[degs$blue$adj.P.Val < 0.05 & degs$blue$logFC<0,])
# degs1[["turquoise.up"]] = rownames(degs$turquoise[degs$turquoise$adj.P.Val < 0.05 & degs$turquoise$logFC>0,])
# degs1[["turquoise.dn"]] = rownames(degs$turquoise[degs$turquoise$adj.P.Val < 0.05 & degs$turquoise$logFC<0,])
degs1 = DEGmachine(degs, FC.cutoff = FC_cutoff)
# # Remove the version tag on he genes
#degs1 = version.remove(degs1)
########################################################################################
# We will do this for several scenarios
file.list = list.files(path = "C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/megena/from.minerva/", pattern = "\\.RData$")
i = file.list[1]             

#for (i in file.list)
#{
  print(i)
  load(paste("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/megena/from.minerva/", i, sep = ""))
  g = igraph::upgrade_graph(g)
  tmp0 = summary.output$modules
  
  #names(tmp0) = gsub("c1_", "M", names(tmp0))
  #tmp0 = version.remove(tmp0)
  n = ncol(datExpr)
  j="blue.ctrl.dn"
  ii = "c1_3"
  overlap.table = as.data.frame(matrix(1, nrow = 1, ncol = 3))
  colnames(overlap.table) = c("module.name", "deg.set", "pval")
  # Find the p.values for the overlaps between the SCZ or SCZsubtype signatures and MEGENA modules
  
  for (ii in names(tmp0))
  {
    for (j in names(degs1)[1:4])# dont use the overlaps and uninque DEGs
    {
      module.genes = tmp0[[ii]]
      go.obj <- newGeneOverlap(module.genes, degs1[[j]], genome.size=n)
      go.obj <- testGeneOverlap(go.obj)
      overlap.table = rbind(overlap.table, c(ii, j, as.numeric(go.obj@pval)))
    }
  }

  overlap.table = overlap.table[!(overlap.table$module.name == "1"),]
  overlap.table$pval = as.numeric(overlap.table$pval)
  overlap.table = overlap.table[order(overlap.table$pval),]
  # We will equate all p.val>0.05 to 1
  overlap.table$pval = ifelse(overlap.table$pval > 0.05, 1, overlap.table$pval)
  # get some coloring (with log transform option)
  # let's do it for each directionality separately
  deg.set = "turq.ctrl.dn"
  overlap.list = list()
  for (deg.set in names(degs1)[1:4])
    {
    print(deg.set)
    mdf = summary.output$module.table
   
    # Find the overlap subtable of only the DEG set in question in this loop
    overlap.table1 = overlap.table[overlap.table$deg.set == deg.set,]
    # Now we will have an mdf table that contains the overlap pvalue on top of 
    mdf = merge(mdf, overlap.table1, by.x = "module.id", by.y = "module.name")
    mdf.x = mdf
    mdf.x$module.id = gsub("c1_", "M", mdf.x$module.id)
    mdf.x$module.parent = gsub("c1_", "M", mdf.x$module.parent)
    
    sbobj = draw_sunburst_wt_fill(module.df = mdf.x,
                                  feat.col = "pval",
                                  log.transform = TRUE,
                                  fill.type = "continuous",
                                  fill.scale = scale_fill_gradient2(low = "white",
                                                                    mid = "white",
                                                                    high = "red",
                                                                    midpoint = -log10(0.05),
                                                                    na.value = "white"),
                                  id.col = "module.id",
                                  parent.col = "module.parent")
    ppi=300
    png(paste(i,".overlap.with.", deg.set,".png", sep = ""),res=ppi,width=10*ppi,height=8.5*ppi)
    sbobj
    dev.off()
    # remove the insignificant pvalues 
    mdf = mdf[mdf$pval < 0.05,]
    
    # print mdf1s first 10 columns
    mdf = mdf[order(mdf$pval),c("module.id", "pval")]
    mdf = mdf[1:10,]
    mdf = mdf[complete.cases(mdf),]
    
    overlap.list[[paste( i, deg.set, sep = ".")]] = mdf
    
    # List of genes for each significantly overlapping module with signatures from mdf
    gene.list = MEGENA.output$module.output$modules
    gene.list = gene.list[mdf$module.id]
    # Remove the version on the gene name
    gene.list = version.remove(gene.list)
    ################################################################################
    #    F U N C T I O N A L    A N A L Y S I S
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
    tmp0 = gene.list

    iiii = "c1_16"
    enrichList = list()
    pval_cutoff = 0.05
    for ( iiii in names(tmp0))
    {tryCatch({
      tmp1 = tmp0[[iiii]]
      print(iiii)

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
                            name.x = gsub("-", "_", iiii), # This does not like special characters
                            name.go = 'GOsets',
                            method = 'hypergeometric',
                            ncores = (detectCores(all.tests = FALSE, logical = TRUE) -1) )

        if(sum(result_weight1$P.adj<pval_cutoff)>0)
        {
          tmp2 = result_weight1[result_weight1$P.adj<pval_cutoff,]
          enrichList[[iiii]] = tmp2
        }
      }
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
    
  saveRDS(enrichList, paste("GO.",i,".overlap.with.", deg.set,".RDS", sep = ""),)
  
  # A part that can be run separately. It is here just because it is a natural flow
  # print mdf1s first 10 columns
  mdf2 = mdf[order(mdf$pval),c("module.id", "pval")]
  mdf2$go = ""
  #replace c1_ with M
  go = readRDS(paste("GO.",i,".overlap.with.", deg.set,".RDS", sep = ""))
  go.module = "c1_46"
  
  for (go.module in names(go))
  {tryCatch({
    tmp1 = go[[go.module]][1,"GOsets"]
    mdf2$go[match(go.module, mdf2$module.id )] = tmp1
  })
  }
  
  mdf2$module.id = gsub("c1_", "M", mdf2$module.id)
  
  write.csv(mdf2, paste("overlap", i, deg.set,"csv", sep = "."), row.names = F)
  }

  
  