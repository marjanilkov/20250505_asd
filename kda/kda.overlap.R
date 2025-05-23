rm(list = ls())

library(SuperExactTest)

wd = "C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2024/20240305_scz_asd_ad/cm/scz/kda/"
setwd(wd)
source("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/R-tomfunctions_static.R")
p.val.thresh = 0.05
FC.thresh = 1.5

mpp.b = read.delim2("MEGENA-KeyDrivers-undirected-CDRadj-dn.MPP.cases.only/WINA_subnets_AD_L3_KD_blue_keydriver.xls")
mpp.t = read.delim2("MEGENA-KeyDrivers-undirected-CDRadj-dn.MPP.cases.only/WINA_subnets_AD_L3_KD_turquoise_keydriver.xls")
nihm = read.delim2("MEGENA-KeyDrivers-undirected-CDRadj-dn-NIHM.cases.only/WINA_subnets_AD_L3_KD_SCZ_keydriver.xls")

mpp.b = mpp.b$keydrivers
mpp.t = mpp.t$keydrivers
nihm = nihm$keydrivers

kda.list = list(MPP.blue=mpp.b, MPP.turquoise=mpp.t, NIHM=nihm)


# Extract the significant overlaps
sTest = supertest(kda.list, degree = 2, n = 19111)
#openImgDev(paste("kda.supertest.",p.val.thresh,"FC",FC.thresh ,".png", sep = ""),iwidth = 2048, iheight = 800, ipointsize = 20)
plot.msets(sTest, degree = 2, 
           sort.by = "size",
           show.overlap.size = T, 
           Layout = "landscape",
           minMinusLog10PValue = 1.3,
           legend.text.cex = 2,
           color.scale.cex = 2,
           cex.axis = 2, 
           cex=2, margin = c(0.5, 8, 1.5, 3))
#dev.off()

# The intersection between NIHM and MPP.blue
x = as.data.frame(intersect(kda.list$MPP.blue, kda.list$NIHM))
colnames(x) = "x"

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(attributes= c("ensembl_gene_id_version", "hgnc_symbol", "uniprotswissprot"), mart= mart)

y = merge(x, G_list, by.x = "x", by.y = "ensembl_gene_id_version")
y = y[!duplicated(y$x),]
colnames(y)[1] = "ensembl"
write.csv(y, "NIHM.MPP.blue.shared.kda.csv")
