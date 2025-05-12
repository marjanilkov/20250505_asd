# Transform the data from counts to CPM and normalize it for technical variables

# Here we will use specifically the ASD +countrols counts data from BA9 from UCLA cohort

rm(list = ls())

setwd("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/data/UCLA-ASD/")

datCounts = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/data/UCLA-ASD/bulk_rnaSeq/ucla.counts.BA9.ASD.and.ctrl.RDS")
rownames(datCounts) = datCounts$Geneid
datCounts = datCounts[7:ncol(datCounts)]

sampleInfo = read.delim2("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/data/UCLA-ASD/Metadata/rnaSeq_UCLA-ASD_metadata.csv", sep = ",")
# Choose only the Brodmann area 9
sampleInfo = sampleInfo[grep("ba9", sampleInfo$libraryID, ignore.case = T), ]

# Sampleinfo and datCounts countain the same identifyers, except that in many 
# places ".", "-", and "_" are used interchangeably. We will make everything 
# separated by "." and all letters will be capital
sampleInfo$libraryID = toupper(sampleInfo$libraryID)
sampleInfo$libraryID = gsub("_", ".", sampleInfo$libraryID)
sampleInfo$libraryID = gsub("-", ".", sampleInfo$libraryID)
rownames(sampleInfo) = sampleInfo$libraryID

colnames(datCounts) = toupper(colnames(datCounts))
colnames(datCounts) = gsub("_", ".", colnames(datCounts))
colnames(datCounts) = gsub("-", ".", colnames(datCounts))

# The sampleInfo libraryID is the one used as column names in datCounts
datCounts = as.data.frame(t(datCounts))

# just a check
datCounts1 = merge( sampleInfo, datCounts, by.x = "libraryID", by.y = "row.names")
datCounts = as.data.frame(t(datCounts))

# Do the analysis now
source("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/sinaiRepo/CountsQC.R")
design.model =  ~0+RIN+sequencingBatch
sampleInfo$group = sampleInfo$sequencingBatch

countsQC(Counts = datCounts, sampleInfo = sampleInfo, design.model = design.model)
