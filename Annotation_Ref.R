# Code written by Micheal; revised by Xiaoya for Micheal's reference. 2022-04-23.


#Microarray agilent
#sureprint_g3_ge_8x60k rat from old and young control and Hearing loss
setwd("D:/microarray-firstCopy")

library(limma)
library(affy)
library(ggplot2)
library(EnhancedVolcano)
library(biomaRt)
require(biomaRt)

## Get annotation file.
rat <- useEnsembl(biomart = "genes",  dataset = "rnorvegicus_gene_ensembl")
tables <- listAttributes(rat)
tables[grep('agilent', tables[,1]),]
listAttributes(rat)
mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('rnorvegicus_gene_ensembl', mart)

# Used this 'agilent_sureprint_g3_ge_8x60k'
annotLookup <- getBM(
  mart = mart,
  attributes = c(
    'agilent_sureprint_g3_ge_8x60k',
    'wikigene_description',
    'ensembl_gene_id',
    'entrezgene_id',
    'gene_biotype',
    'external_gene_name',
    'rgd_symbol'))
colnames(annotLookup)[1] <- 'AgilentID'
annotLookup <- annotLookup[which(annotLookup$AgilentID != ""),]; dim(annotLookup) #66622     7



## Read 40 raw files in.
targets <- readTargets("targets.txt")

# convert the data to an EListRaw object: data object for single channel data
# specify green.only = TRUE for Agilent
# retain information about background via gIsWellAboveBG
project <- read.maimages(
  targets,
  source = 'agilent.median',
  green.only = TRUE,
  other.columns = c('gIsWellAboveBG', "gBGMeanSignal"))
colnames(project) <- gsub('raw\\/', '', colnames(project))   


## Get ID_Gene_Data matrix.
ID_Data <- cbind(as.data.frame(project$genes$ProbeName), project$other$gBGMeanSignal);dim(ID_Data); colnames(ID_Data)[1] <- "ProbeName"  # To get the numeric data which come from the file Project_Targets.
#  62976    41
ID_Gene_Data <- merge(x = annotLookup[,c("AgilentID", "external_gene_name", "rgd_symbol")], y = ID_Data, by.x = "AgilentID", by.y = "ProbeName", all = FALSE);dim(ID_Gene_Data)
# 64489    43

# ID_Gene_Data[which(ID_Gene_Data$AgilentID == "A_42_P453935"),1:4]
# length(which(ID_Gene_Data$AgilentID == "A_42_P453935")) # 10