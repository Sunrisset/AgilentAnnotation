# Code written by Micheal. 2022-04-22.


#Microarray agilent
#sureprint_g3_ge_8x60k rat from old and young control and Hearing loss
#

library(limma)
library(affy)
library(ggplot2)
library(EnhancedVolcano)

## Make annotation file
library(biomaRt)

setwd("G:/microarray-firstCopy")
require(biomaRt)

rat <- useEnsembl(biomart = "genes",  dataset = "rnorvegicus_gene_ensembl")
tables <- listAttributes(rat)
tables[grep('agilent', tables[,1]),]
# Used this 'agilent_sureprint_g3_ge_8x60k'
listAttributes(rat)
mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('rnorvegicus_gene_ensembl', mart)

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

# write.table(
#   annotLookup,
#   paste0('AgilentAnnotationfile', gsub("-", "_", as.character(Sys.Date())), '.tsv'),
#   sep = '\t',
#   row.names = FALSE,
#   quote = FALSE)

# ##Preprocessing data
# setwd("C:/Users/Physics/Desktop/Microarray-Data/microarray-firstCopy/")
targets <- readTargets("targets.txt")
# annotfile <- read.table('./AgilentAnnotationfile.tsv', sep = "\t")

# targetsinfo <- readTargets("targets.txt")
#targetsinfo <- readTargets(targetsfile, sep = '\t')

# convert the data to an EListRaw object: data object for single channel data
# specify green.only = TRUE for Agilent
# retain information about background via gIsWellAboveBG
##
project <- read.maimages(
  targets,
  source = 'agilent.median',
  green.only = TRUE,
  other.columns = c('gIsWellAboveBG', "gBGMeanSignal"))
colnames(project) <- gsub('raw\\/', '', colnames(project))
colnames(annotLookup)[1] <- 'AgilentID'

# overlap0 <- annotLookup[which(annotLookup$AgilentID %in% project$genes$ProbeName),]

Overlap_ID_Gene <- annotLookup[which(project$genes$ProbeName %in% annotLookup$AgilentID),]  # To get the (db) IDs that match our IDs. The gene name information also stored in Overlap_ID, which is what we need ultimatly.

ID_Data <- cbind(as.data.frame(project$genes$ProbeName), project$other$gBGMeanSignal); colnames(ID_Data)[1] <- "ProbeName"  # To get the numeric data which come from file Project_Targets.
Overlap_ID_data <- ID_Data[which(project$genes$ProbeName %in% annotLookup$AgilentID),]  # To get the numeric data of the matched IDs.

library(dplyr)


Overlap_ID_Gene_Data <- full_join(x = Overlap_ID_Gene,  y = Overlap_ID_data, by = "")


b <- merge(x = Overlap_ID_Gene$AgilentID, y = Overlap_ID_data[,1:2])





# overlap1 <- annotLookup[match(project$genes$ProbeName, annotLookup$AgilentID),]

# table(project$genes$ProbeName == annotLookup$AgilentID) # check that annots are aligned



# project$genes$AgilentID <- annotLookup$AgilentID
# project$genes$wikigene_description <- annotLookup$wikigene_description
# project$genes$ensembl_gene_id <- annotLookup$ensembl_gene_id
# project$genes$entrezgene_id <- annotLookup$entrezgene_id
# project$genes$gene_biotype <- annotLookup$gene_biotype
# project$genes$external_gene_name <- annotLookup$external_gene_name


# #perform background correction on the fluorescent intensities
# project.bgcorrect <- backgroundCorrect(project, method = 'normexp')
# # normalize the data with the 'quantile' method
# project.bgcorrect.norm <- normalizeBetweenArrays(project.bgcorrect, method = 'quantile')


# filter out control probes, those with no symbol, and those that fail:

#Control <- project.bgcorrect.norm$genes$ControlType==1L
#NoSymbol <- is.na(project.bgcorrect.norm$genes$external_gene_name)
#IsExpr <- rowSums(project.bgcorrect.norm$other$gIsWellAboveBG > 0) >= 1
#project.bgcorrect.norm.filt <- project.bgcorrect.norm[!Control & !NoSymbol & IsExpr, ]

# dim(project.bgcorrect.norm)
# #dim(project.bgcorrect.norm.filt)
# # remove annotation columns we no longer need
# 


project.bgcorrect.norm$genes <- project.bgcorrect.norm$genes[,c(
  'wikigene_description','ensembl_gene_id','entrezgene_id','gene_biotype',)] #'external_gene_name', 'ProbeName'
head(project.bgcorrect.norm$genes)

project.bgcorrect.norm.mean <- avereps(project.bgcorrect.norm,
                                       ID = project.bgcorrect.norm$genes$ProbeName)
dim(project.bgcorrect.norm.mean)
