# Transform the data from counts to CPM and normalize it for technical variables

# Here we will use specifically the ASD +countrols counts data from BA9 from UCLA cohort

rm(list = ls())

setwd("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/data/UCLA-ASD/bulk_rnaSeq/")

datCounts = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/data/UCLA-ASD/bulk_rnaSeq/ucla.counts.BA9.ASD.and.ctrl.RDS")
rownames(datCounts) = datCounts$Geneid
datCounts = datCounts[7:ncol(datCounts)]

# Using the CountsQC.R script we found that one sample, AN01093.BA9.2012.150 had 
#only 345539 counts in total so we removed this before anything else could continue.
datCounts$AN01093_BA9_2012.150 = NULL

individualInfo = read.delim2("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/data/UCLA-ASD/Metadata/individual_UCLA-ASD_metadata.csv", sep = ",")
sampleID.to.individualID <- readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/data/UCLA-ASD/bulk_rnaSeq/sampleID.to.individualID.RDS")
# Choose only the Brodmann area 9
sampleID.to.individualID = sampleID.to.individualID[sampleID.to.individualID$ba == "BA9",]
sampleInfo = merge(sampleID.to.individualID, individualInfo, by.x = "individualID", by.y = "individualID")
rm(individualInfo, sampleID.to.individualID)
rownames(sampleInfo) = sampleInfo$sample
sampleInfo$group = as.factor(sampleInfo$primaryDiagnosis)

# Sampleinfo and datCounts countain the same identifiers, except that in many 
# places ".", "-", and "_" are used interchangeably. We will make everything 
# separated by "." and all letters will be capital
sampleInfo$sample = toupper(sampleInfo$sample)
sampleInfo$sample = gsub("_", ".", sampleInfo$sample)
sampleInfo$sample = gsub("-", ".", sampleInfo$sample)
rownames(sampleInfo) = sampleInfo$sample

colnames(datCounts) = toupper(colnames(datCounts))
colnames(datCounts) = gsub("_", ".", colnames(datCounts))
colnames(datCounts) = gsub("-", ".", colnames(datCounts))

# Remove NAs
sampleInfo = sampleInfo[!is.na(sampleInfo$group),]
# Knit everything properly
# create a common ID to be able to properly knit together these data frames
commonID = Reduce(intersect,list(rownames(sampleInfo), colnames(datCounts)))
# put everything together in a list
sampleInfo = sampleInfo[commonID,]
datCounts = datCounts[,commonID]

# Do the analysis now
source("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/sinaiRepo/CountsQC.R")
design.model = ~0+group


countsQC(Counts = datCounts, sampleInfo = sampleInfo, design.model = design.model)
