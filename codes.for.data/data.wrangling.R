# A script to explore the PsychEncode ASD data from UCLA and Yale
rm(list = ls())

setwd("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/data")

#save in RDS format for faster loading
# ucla.counts = read.delim2("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/data/UCLA-ASD/bulk_rnaSeq/UCLA_ASD_bulk_rnaSeq.tsv")
# saveRDS(ucla.counts, "C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/data/UCLA-ASD/bulk_rnaSeq/UCLA_ASD_bulk_rnaSeq.RDS")

ucla.counts = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/data/UCLA-ASD/bulk_rnaSeq/UCLA_ASD_bulk_rnaSeq.RDS")
meta.biospec = read.delim2("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/data/UCLA-ASD/Metadata/biospecimen_UCLA-ASD_metadata.csv", sep = ",")
meta.individual = read.delim2("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/data/UCLA-ASD/Metadata/individual_UCLA-ASD_metadata.csv", sep = ",")
meta.rnaSeq = read.delim2("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2025/20250505_asd/data/UCLA-ASD/Metadata/rnaSeq_UCLA-ASD_metadata.csv", sep = ",")

# The individual IDs are not uniform between meta.rnseq and the others, so we will make them so
meta.rnaSeq$specimenID = gsub("Sample_", "", meta.rnaSeq$specimenID)
# Since all the data in these columns is NA, we will remove them
meta.rnaSeq <- meta.rnaSeq[,colSums(is.na(meta.rnaSeq))<nrow(meta.rnaSeq)]
meta.individual <- meta.individual[,colSums(is.na(meta.individual))<nrow(meta.individual)]
meta.biospec <- meta.biospec[,colSums(is.na(meta.biospec))<nrow(meta.biospec)]

table(meta.rnaSeq$sequencingBatch)
# replace the dash divider with a point for uniformity with other places
meta.rnaSeq$sequencingBatch = gsub("-", ".", meta.rnaSeq$sequencingBatch)

# Biospec
table(meta.biospec$BrodmannArea)

# From the column names(samples ID) in ucla.counts we need to create a table and 
# find out how many samples from what brain regions etc
df.colnames = as.data.frame(colnames(ucla.counts[7:ncol(ucla.counts)]))
colnames(df.colnames) = "sample"
# I have found that the sample UMB5176.3.1.2.5_2016.9099 misses the string "BA"
# in front of the Brodmann area 3.1.2.5, so I will add it by hand since this is 
# the only case in the whole file with this error 
df.colnames$sample[767] = "UMB5176.BA3.1.2.5_2016.9099"

df.colnames$sample.BA = sub("_[^_]+$", "", df.colnames$sample)
df.colnames$sequencingBatch = gsub("^(.*)_", "", df.colnames$sample)
library(stringr)

# Take the BRodman area and put it in a separate column
separate.BA = function(a)
{
  if (grepl("BA", a)){b = str_remove(a, '.*(?=BA)')}
  if (grepl("CBL", a)){b = str_remove(a, '.*(?=CBL)')}
  return(b)
}
i = 767
for (i in 1:nrow(df.colnames))
{
  print(i)
  df.colnames$ba[i] = separate.BA(df.colnames$sample.BA[i])  
}

separate.sample = function(a)
{
  if (grepl("BA", a)){b = str_remove(a, '(?=_BA).*')}
  if (grepl("CBL", a)){b = str_remove(a, '(?=_CBL).*')}
  if (grepl("BA", a)){b = str_remove(a, '(?=.BA).*')}
  if (grepl("CBL", a)){b = str_remove(a, '(?=.CBL).*')}
  return(b)
}

i = 9
for (i in 1:nrow(df.colnames))
{
  print(i)
  df.colnames$individualID[i] = separate.sample(df.colnames$sample.BA[i])  
}

df.colnames = df.colnames[,c("sample", "individualID", "ba", "sequencingBatch")]
#saveRDS(df.colnames, "sampleID.to.individualID.RDS")
# # We need to figure out how many samples we have from each region
# table(df.colnames$ba)
# 
# # create a common ID to be able to proterly kniw together these DFs
# commonID = Reduce(intersect,list(df.colnames$individualID, meta.individual$individualID))
# 
# # We now have two tables, meta.individual binds each individualID to disease status
# # and df.colnames binds individualID to Broadmann area. We will merge these two
# df.asd.ba = merge(df.colnames, meta.individual, by = "individualID")
# 
# # Remove all the samples that don't have information on ASD disease status
# df.asd.ba = df.asd.ba[!is.na(df.asd.ba$primaryDiagnosis),]
# 
# tmp1 = df.asd.ba[df.asd.ba$ba == "BA9",]
# table(tmp1$primaryDiagnosis)
# 
# # We would like to have a df with only individualID and diagnosis regardless of 
# # brain region. This is another check to see if meta.individual is correctly formatted
# df.diagnosis = df.asd.ba[,c("individualID", "primaryDiagnosis")]
# df.diagnosis = df.diagnosis[!duplicated(df.diagnosis$individualID),]
# table(df.diagnosis$primaryDiagnosis)
# 
# # We now have the IDs of all individuals for which we have diagnosis and bulk RNAseq data
# df.colnames = df.colnames[df.colnames$individualID %in% commonID,]
# 
# # Save the UCLA BA9 counts
# cols <- grep("BA9", meta.individual$specimenID, value = TRUE)
# meta.biospec1 = meta.biospec[meta.biospec$specimenID %in% cols,]



# We will use only the samples from BA9=DLPFC
ucla.BA9.samples = grepl("BA9", colnames(ucla.counts))
sum(ucla.BA9.samples)
df.colnames = df.colnames[grepl("BA9", df.colnames$sample),]

ucla.counts.tmp1 = ucla.counts[, ucla.BA9.samples]
ucla.counts.tmp2 = ucla.counts[, colnames(ucla.counts)[1:6]]
ucla.counts = cbind(ucla.counts.tmp2, ucla.counts.tmp1)
rm(ucla.counts.tmp1, ucla.counts.tmp2)
saveRDS(ucla.counts, "ucla.counts.BA9.ASD.and.ctrl.RDS")
# We will now choose only the samples with ASD
meta.individual = meta.individual[meta.individual$primaryDiagnosis == "Autism Spectrum Disorder",]
meta.individual = merge(meta.individual, df.colnames, by = "individualID")

ucla.counts1 = ucla.counts[,meta.individual$sample]
ucla.counts2 = ucla.counts[,1:6]
ucla.counts = cbind(ucla.counts2, ucla.counts1)

#saveRDS(ucla.counts, "ucla.counts.BA9.ASD.RDS")
