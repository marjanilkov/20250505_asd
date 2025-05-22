rm(list=ls())
library(pheatmap)
options(expressions = 5e5)
library(e1071)
library(data.table)
library(viridis)
library("RColorBrewer")
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(as.matrix(xs), probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

### data loading ###

print("loading data...")
source("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew3/20220829_data/codes_for_data/R-tomfunctions_static.R")

setwd("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/sample.heatmap/")
datExpr = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/data/UCLA-ASD/bulk_rnaSeq/limma.BA9/ucla.RNAseq.BA9.adj.RDS")
# Z-transform the data
zmatrix = apply(datExpr, 2, Ztransform, 4) # 1 means the apply function is done on rows and 2 would mean columns
rownames(zmatrix) = rownames(datExpr)
datExpr = zmatrix
dim(datExpr)
rm(zmatrix)
gc()
datExpr = datExpr[complete.cases(datExpr),]
gnxp = as.data.frame(t(datExpr))
rm(datExpr)

#meta_br = read.delim2(file="msbb.meta.BM_36.tsv",header=T,sep="\t")
# classification = read.delim2(file="wina_subtypes_meta_with_new_2018_defs.tsv",header=T,sep="\t")
# classification$sample <- NULL
# classification = merge(meta_br,classification,by="BB")
# classification$cluster = classification$wina_cluster
#classification = classification[complete.cases(classification),]
classification = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/clust/wina.BA9.UCLA.bulkRNAseq/gene.beta1/ucla.BA9.metadata.w.subtypes.RDS")
colnames(classification)[ncol(classification)] = "cluster"

wina_modules = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd//sample.heatmap/gene_average/complete.gene.clustering.RDS") 
wina_modules$gene_name = rownames(wina_modules)
wina_modules = wina_modules[,c("gene_name", "module")]
colnames(wina_modules) = c("gene_name", "module")

# modules_genes_order = read.delim2("/sc/arion/projects/zhangb03a/neffr01/AD_subtyping_megena/MSBB_updated/MEGENA_network_genes_MSBB_BM36_cdr_gt_1/modules_genes_order.tsv",header=T)

megena_modules = wina_modules

clust.method = "complete"
if (clust.method == "complete"){h1row = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/sample.heatmap/gene_average/h1row.complete1.RDS")}

print("loaded data")

### preparing gene_exp_file ###
gene_exp = gnxp
#rownames(gene_exp) = make.unique(as.character(gene_exp$Symbol))
gene_exp_t = transpose(gene_exp)
rownames(gene_exp_t) = colnames(gene_exp)
colnames(gene_exp_t) = rownames(gene_exp)
#gene_exp_t = gene_exp_t[-1,]
gene_exp_t = data.frame(data.matrix(gene_exp_t))

### speedy settings ###
if (!exists("gene_exp_t_all")){
  gene_exp_t_all = gene_exp_t
  ###
  ##
  # USAGE OF MATCH::
  #
  # Y = Y[,match(X, Y)]
  # Y = Y[,match(X_df_in_order_you_want, Y_df_you_want_to_order)]
  #
  ##
  ###
  all_genes_order_master = match(h1row$labels[h1row$order],colnames(gene_exp_t_all)) ##THIS IS FROM WINA
  all_genes_order_master = all_genes_order_master[!is.na(all_genes_order_master)]
  gene_exp_t_all = gene_exp_t_all[,all_genes_order_master]
}
gene_exp_t = gene_exp_t_all#[,1:2000]
###

### preparing color bars ###
mat_col_samples = data.frame(
  cluster=classification[match(rownames(gene_exp_t),classification$sample),"cluster"],
  sample=classification[match(rownames(gene_exp_t),classification$sample),"sample"]
) #samples are rows
mat_col_genes = data.frame(
  gene_module=megena_modules[match(colnames(gene_exp_t),megena_modules$gene_name),"module"],
  gene = megena_modules[match(colnames(gene_exp_t),megena_modules$gene_name),"gene_name"]
) #genes are columns
mat_col_genes = mat_col_genes[complete.cases(mat_col_genes),]
mat_col_samples = mat_col_samples[complete.cases(mat_col_samples),]

rownames(mat_col_genes) = mat_col_genes$gene
rownames(mat_col_samples) = mat_col_samples$sample

gene_exp_t = gene_exp_t[mat_col_samples[order(mat_col_samples$cluster),]$sample,]

#mat_colors_samples = heat.colors(length(unique(mat_col_samples$cluster)))

### clustering within gene modules ###
#mat_col_genes_2 = mat_col_genes
#mat_col_genes_2$gene = colnames(gene_exp_t)
#all_genes_order = c()
#for(genemod in unique(mat_col_genes_2$gene_module)){
#  clust_cols = gene_exp_t[,mat_col_genes_2[mat_col_genes_2$gene_module==genemod,"gene"]]
#  order_genes = hclust(dist(as.matrix(transpose(clust_cols))))$order
#  all_genes_order = c(all_genes_order, colnames(clust_cols[order_genes,]))
#}
#gene_exp_t = gene_exp_t[,all_genes_order]
#mat_col_genes_2 = data.frame(gene_module=mat_col_genes[order(mat_col_genes$gene_module),"gene_module"]) #genes are columns

### let's order the samples by hclust ###
all_samples_order = c()
#for(clust in unique(mat_col_samples$cluster)){
for(clust in c("control", "turquoise", "blue")){
  
  selectsamp = match(mat_col_samples[mat_col_samples$cluster==clust,"sample"], rownames(gene_exp_t))
  clust_rows = gene_exp_t[selectsamp,]
  clust_rows = na.omit(clust_rows)
  #clust_rows = clust_rows[-1,]
  order_samples = hclust(dist(as.matrix(clust_rows)))$order
  all_samples_order = c(all_samples_order, rownames(clust_rows[order_samples,]))
}

gene_exp_t = gene_exp_t[all_samples_order,]

mat_col_samples = data.frame(
  cluster=classification[match(rownames(gene_exp_t),classification$sample),"cluster"],
  sample=classification[match(rownames(gene_exp_t),classification$sample),"sample"]
) 
rownames(mat_col_samples) = mat_col_samples$sample

### prepare to output
#mat_breaks <- quantile_breaks(gene_exp_t, n = 20)
gene_exp_t_2 = transpose(gene_exp_t) #for plotting
rownames(gene_exp_t_2) = colnames(gene_exp_t)
colnames(gene_exp_t_2) = rownames(gene_exp_t)
breakslist = c(-5,-2.32635,-1.28155,-0.841,-0.5244,-0.2535,0,0.2535,0.5244,0.841,1.28155,2.32635,5)
x = unique(mat_col_genes$gene_module)
y = unique(mat_col_genes$gene_module)
test = NULL
test$gene_module = unlist(mapply(function(x,y) { as.character(y) }, as.character(x), as.character(y), SIMPLIFY = FALSE,USE.NAMES = TRUE),use.names = TRUE)
x = unique(mat_col_samples$cluster)
y = brewer.pal(n=length(x),name="Set1")
test$cluster = unlist(mapply(function(x,y) { as.character(y) }, as.character(x), as.character(y), SIMPLIFY = FALSE,USE.NAMES = TRUE),use.names = TRUE)

### PLOTTING NOW
mat_col_genes$gene = NULL
mat_col_samples$sample=NULL
ppi = 300
png("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250415_scz/sample.heatmap/pheatmap_asd.png",res=ppi,width=8*ppi,height=8*ppi)
pheatmap(
  mat               = gene_exp_t_2,
  color             = colorRamps::blue2yellow(length(breakslist)),
  border_color      = NA,
  show_colnames     = FALSE,
  show_rownames     = FALSE,
  breaks            = breakslist,
  annotation_row    = mat_col_genes,
  annotation_col    = mat_col_samples,
  gaps_col          = match(unique(mat_col_samples$cluster),mat_col_samples$cluster)[2:length(unique(mat_col_samples$cluster))]-1,
  annotation_colors = test,
  cluster_cols      = FALSE,
  cluster_rows      = FALSE,
  drop_levels       = TRUE,
  legend            = TRUE,
  scale             = "row",
  fontsize          = 8,
  main              = "Gene expression of SCZ and controls in MSBB-Penn-Pitt by WINA cluster"
)
dev.off()

### DONE
