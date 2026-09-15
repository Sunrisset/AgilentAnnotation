# Code written by Xiaoya for Micheal's reference. 2022-04-23.


# ====== Agilent annotation for raw txt files. ======

setwd("D:/microarray-firstCopy")

library(biomaRt)
library(dplyr)

# detremine the database we will use.
listMarts()

# find the name of the rat ensembl genes dataset under ensembl mart.
mart <- useMart('ENSEMBL_MART_ENSEMBL')
dl <- as.data.frame(listDatasets(mart=mart))
dl %>% filter(grepl("rat",tolower(description))) # mRatBN7.2   rnorvegicus_gene_ensembl.

# AgilentID from rnorvegicus_gene_ensembl.
mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset(dataset = 'rnorvegicus_gene_ensembl', mart)
la <- listAttributes(mart=mart);head(la);View(la)
searchAttributes(mart=mart,pattern="agilent")

# meta data from rnorvegicus_gene_ensembl.
mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset(dataset = 'rnorvegicus_gene_ensembl', mart)
rdata <- getBM(mart=mart,
               attributes=c('agilent_sureprint_g3_ge_8x60k',
                            'wikigene_description',
                            "ensembl_gene_id",
                            "entrezgene_id",
                            "gene_biotype",
                            "external_gene_name",
                            'rgd_symbol'),
               uniqueRows=TRUE, useCache=FALSE)
rdata_noEmptyID <- rdata[which(rdata$agilent_sureprint_g3_ge_8x60k != ""),]; colnames(rdata_noEmptyID)[1] <- "AgilentID"; dim(rdata_noEmptyID) #66622     7


# create list object G.
library(limma)
targets <- readTargets("targets.txt")
G <- read.maimages(targets, source = 'agilent.median', green.only = TRUE, other.columns = c('gIsWellAboveBG', "gBGMeanSignal"))


# Get ID_Gene_Data matrix.
ID_Data <- cbind(as.data.frame(G$genes$ProbeName), G$other$gBGMeanSignal);dim(ID_Data); colnames(ID_Data)[1] <- "ProbeName"
ID_Gene_Data <- merge(x = rdata_noEmptyID[,c("AgilentID", "external_gene_name", "rgd_symbol")], y = ID_Data, by.x = "AgilentID", by.y = "ProbeName", all = FALSE);dim(ID_Gene_Data)
# 64489    43
write.table(x = ID_Gene_Data, file = "ID_Gene_Data.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

# The end. :)



