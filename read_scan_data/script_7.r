############################
#	Bárbara Bitarello
#
#	Last modified: 28.10.2016
###########################


library(parallel)
library(SOAR)
library(ggplot2)
library(vioplot)
library(dplyr)
library(reshape)
library(reshape2)
Sys.setenv(R_LOCAL_CACHE="estsession")
library(data.table)

pops<-c("AWS","LWK","YRI","CEU", "FIN","GBR","TSI", "CHB","CHS" ,"JPT","MXL", "CLM","PUR")



my.function.improved<-function(B, E, df=LWK.chr, chr=6){
as.numeric(gsub("chr", "", chr))-> chr
rbind(filter(df[[chr]], Chr==chr,End.Win > B,End.Win < E), filter(df[[chr]], Chr==chr,Beg.Win > B,Beg.Win < E), filter(df[[chr]], Chr==chr,Beg.Win<B, End.Win>E), filter(df[[chr]], Chr==chr,Beg.Win>B ,End.Win<E))->res
setDT(res)
setkey(res, Win.ID)
#df[rownames(res[!duplicated(res),]),]-> res2
unique(res)-> res2
return(res2)
}

####
####

names(bed2)<-c('chr', 'beg', 'end','name', 'subtype', 'gene.id', 'type')
setDT(bed2)

mclapply(1:22, function(x) subset(bed2, chr==paste0('chr',x)))-> coding.per.chr.list

list.SCAN[[2]]-> LWK
setDT(LWK)

mclapply(1:22, function(x) setDT(filter(LWK, Chr==x)))-> LWK.chr

#test for chr 22:

system.time(mclapply(1:nrow(coding.per.chr.list[[22]]), function(x) my.function.improved(B=coding.per.chr.list[[22]]$beg[x], E=coding.per.chr.list[[22]]$end[x], chr=paste0("chr", 22), df=LWK.chr))-> test)  #1577.720
#
test-> test.22

system.time(mclapply(1:nrow(coding.per.chr.list[[21]]), function(x) my.function.improved(B=coding.per.chr.list[[21]]$beg[x], E=coding.per.chr.list[[21]]$end[x], chr=paste0("chr", 21), df=LWK.chr))-> test.21)  #462.091


system.time(mclapply(1:nrow(coding.per.chr.list[[21]]), function(x) my.function.improved(B=coding.per.chr.list[[21]]$beg[x], E=coding.per.chr.list[[21]]$end[x], chr=paste0("chr", 21), df=LWK.chr))-> test.21)  #462.091


#


