###################################################################################
#	Bárbara Bitarello
#	Mast modified: 24.11.2016
#	This  is based on the old scripts, but hopefully with some improvements.
#
####################################################################################

library(parallel)
library(SOAR)
library(data.table)
library(dplyr)

setDT(read.table('Background.txt', header=T, sep="\t"))-> ALL.GENES
setDT(read.table('LWK_paralogs.txt', header=T, sep="\t"))-> LWK.cand.paralogs
setDT(read.table('OR_LWK_paralogs2.txt', header=T, sep="\t"))-> OR.paralogs

nrow(ALL.GENES) #69108

length(unique(ALL.GENES$Associated_Gene_Name)) # 13447
na.omit(ALL.GENES)-> ALL.GENES.1
#factor(ALL.GENES[,1])-> ALL.GENES[,1] #13439 genes IDs

length(unique(ALL.GENES.1$Associated_Gene_Name)) #9224  #these are the genes left after na.omit
nrow(ALL.GENES.1) #64871
remove(ALL.GENES)
gc()


nrow(LWK.cand.paralogs) #8832
na.omit(LWK.cand.paralogs)-> LWK.cand.paralogs1
remove(LWK.cand.paralogs); gc()
nrow(LWK.cand.paralogs1) #8528
nrow(OR.paralogs) #235
na.omit(OR.paralogs)-> OR.paralogs1
remove(OR.paralogs); gc()
nrow(OR.paralogs1) #234

gcinfo(FALSE)
gc()


split(ALL.GENES.1, ALL.GENES.1$Ensembl_Gene_ID)-> big.list
split(LWK.cand.paralogs, YRI.cand.paralogs$Ensembl_Gene_ID)-> LWK.cand.paralogs.split
split(OR.paralogs, OR.paralogs$Ensembl_Gene_ID)-> OR.paralogs.split

