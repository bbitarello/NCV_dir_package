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
library(qqman)

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

Store(test.21, test.22)

system.time(mclapply(1:nrow(coding.per.chr.list[[21]]), function(x) my.function.improved(B=coding.per.chr.list[[21]]$beg[x], E=coding.per.chr.list[[21]]$end[x], chr=paste0("chr", 21), df=LWK.chr))-> test.21)  #462.091

testB<-vector('list', 20)

for ( i in 1:20){

system.time(mclapply(1:nrow(coding.per.chr.list[[i]]), function(x) my.function.improved(B=coding.per.chr.list[[i]]$beg[x], E=coding.per.chr.list[[i]]$end[x], chr=paste0("chr", i), df=LWK.chr))-> testB[[i]])

gc()
cat ('chromosome ',i, ' done\n')

} #chrs 17 and 20 failed. need to redo

system.time(mclapply(1:nrow(coding.per.chr.list[[17]]), function(x) my.function.improved(B=coding.per.chr.list[[17]]$beg[x], E=coding.per.chr.list[[17]]$end[x], chr=paste0("chr",17), df=LWK.chr))-> testB[[17]])


system.time(mclapply(1:nrow(coding.per.chr.list[[20]]), function(x) my.function.improved(B=coding.per.chr.list[[20]]$beg[x], E=coding.per.chr.list[[20]]$end[x], chr=paste0("chr", 20), df=LWK.chr))-> testB[[20]])






#manhattan plots g



mclapply(list.SCAN, function(x) select(x, Chr:End.Win, Dist.NCD.f0.5, Z.f0.5.P.val))->tes.manhattan.f0.5

mclapply(tes.manhattan.f0.5, function(x) cbind(x, BP=(x$Beg.Win+x$End.Win)/2))->tes.manhattan.2.f0.5


for (i in 1:7){
colnames(tes.manhattan.2.f0.5[[i]])[c(1,4,5)]<-c('CHR', 'SNP', 'P')}


mclapply(tes.manhattan.2.f0.5, function(x) setDT(arrange(x, SNP)))->tes.manhattan.f0.5


top829f0.5<-mclapply(tes.manhattan.f0.5,function(x) head(x,829))

mclapply(top829f0.5, function(x) arrange(x, CHR, Beg.Win))->sort.top829f0.5 



pdf('bedfiles/top829.my.man.test.pdf')
manhattan(tes.manhattan.f0.5[[3]], 
highlight=as.character(top829f0.5[[3]]$BP),suggestiveline=-log10(0.0001),genomewideline=-log10(0.0005001925))
dev.off()

