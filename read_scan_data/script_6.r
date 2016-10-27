###############################################
#
#	Barbara D Bitarello
#
#	Last modified: 20.10.2016
#
#	Make bedfiles for Joao
#################################################


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

BED.PATH<-'/mnt/sequencedb/PopGen/barbara/NCV_dir_package/read_scan_data/bedfiles/'

#generate a background bed file (all scanned genes)

###MAKE UNION OF TF FOR EACH POP ##########
mclapply(1:7, function(x) setDT(rbind(CANDf0.5[[x]], CANDf0.4[[x]], CANDf0.3[[x]])[-(which(duplicated(rbind(CANDf0.5[[x]], CANDf0.4[[x]], CANDf0.3[[x]])))),]))-> Union.CANDf0.5_0.4_0.3

mclapply(1:7, function(x) setDT(rbind(top829f0.5[[x]], top829f0.4[[x]], top829f0.3[[x]])[-(which(duplicated(rbind(top829f0.5[[x]], top829f0.4[[x]], top829f0.3[[x]])))),]))-> Union.top0.5_0.4_0.3

Store(Union.CANDf0.5_0.4_0.3)
Store(Union.top0.5_0.4_0.3)

Objects()

#######################
#a few sanity checks

####
#test
for(i in 1:7){
nrow(list.SCAN[[i]])-> n1;nrow(Union.CANDf0.5_0.4_0.3[[i]])-> n2;nrow(Union.top0.5_0.4_0.3[[i]])-> n3

cbind(rbind(select(list.SCAN[[i]], Nr.IS, Nr.FDs, Nr.SNPs, PtoD), select(Union.CANDf0.5_0.4_0.3[[i]], Nr.IS, Nr.FDs, Nr.SNPs, PtoD), select(Union.top0.5_0.4_0.3[[i]], Nr.IS, Nr.FDs, Nr.SNPs, PtoD)), Type=c(rep('genomic', n1), rep('Significant', n2), rep('Outliers', n3)))-> temp

tmpname<-paste0('Nr.IS.',pops[[i]], '.pdf')
tmpname1<-paste0('Nr.FDs.',pops[[i]], '.pdf')
tmpname2<-paste0('Nr.SNPs.',pops[[i]], '.pdf')
tmpname3<-paste0('PtoD.',pops[[i]], '.pdf')

ggplot(temp) + geom_density(aes(x = Nr.IS, colour = Type))
ggsave(paste0('figures/',tmpname))

ggplot(temp) + geom_density(aes(x = Nr.FDs, colour = Type))
ggsave(paste0('figures/',tmpname1))

ggplot(temp) + geom_density(aes(x = Nr.SNPs, colour = Type))
ggsave(paste0('figures/',tmpname2))

ggplot(temp) + geom_density(aes(x = PtoD, colour = Type))
ggsave(paste0('figures/',tmpname3))

cat('Finished', pops[i], '\n')
}


#

#plots SFS

source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/SFS_script.r')


#tf=0.5

unlist(mclapply(1:nrow(top829f0.5[[2]]), function(x) SFS.function(CHR=top829f0.5[[2]][x,]$Chr, BEG=top829f0.5[[2]][x,]$Beg.Win, END=top829f0.5[[2]][x,]$End.Win, POP=2)))-> LWK.top.f0.5

unlist(mclapply(1:nrow(top829f0.5[[3]]), function(x) SFS.function(CHR=top829f0.5[[3]][x,]$Chr, BEG=top829f0.5[[3]][x,]$Beg.Win, END=top829f0.5[[3]][x,]$End.Win, POP=3)))-> YRI.top.f0.5

unlist(mclapply(1:nrow(top829f0.5[[6]]), function(x) SFS.function(CHR=top829f0.5[[6]][x,]$Chr, BEG=top829f0.5[[6]][x,]$Beg.Win, END=top829f0.5[[6]][x,]$End.Win, POP=6)))-> GBR.top.f0.5

unlist(mclapply(1:nrow(top829f0.5[[7]]), function(x) SFS.function(CHR=top829f0.5[[7]][x,]$Chr, BEG=top829f0.5[[7]][x,]$Beg.Win, END=top829f0.5[[7]][x,]$End.Win, POP=7)))-> TSI.top.f0.5

#

unlist(mclapply(1:nrow(CANDf0.5[[2]]), function(x) SFS.function(CHR=CANDf0.5[[2]][x,]$Chr, BEG=CANDf0.5[[2]][x,]$Beg.Win, END=CANDf0.5[[2]][x,]$End.Win, POP=2)))-> LWK.cand.f0.5

unlist(mclapply(1:nrow(CANDf0.5[[3]]), function(x) SFS.function(CHR=CANDf0.5[[3]][x,]$Chr, BEG=CANDf0.5[[3]][x,]$Beg.Win, END=CANDf0.5[[3]][x,]$End.Win, POP=3)))-> YRI.cand.f0.5

unlist(mclapply(1:nrow(CANDf0.5[[6]]), function(x) SFS.function(CHR=CANDf0.5[[6]][x,]$Chr, BEG=CANDf0.5[[6]][x,]$Beg.Win, END=CANDf0.5[[6]][x,]$End.Win, POP=6)))-> GBR.cand.f0.5

unlist(mclapply(1:nrow(CANDf0.5[[7]]), function(x) SFS.function(CHR=CANDf0.5[[7]][x,]$Chr, BEG=CANDf0.5[[7]][x,]$Beg.Win, END=CANDf0.5[[7]][x,]$End.Win, POP=7)))-> TSI.cand.f0.5

#tf=0.3


unlist(mclapply(1:nrow(top829f0.3[[2]]), function(x) SFS.function(CHR=top829f0.3[[2]][x,]$Chr, BEG=top829f0.3[[2]][x,]$Beg.Win, END=top829f0.3[[2]][x,]$End.Win, POP=2)))-> LWK.top.f0.3

unlist(mclapply(1:nrow(top829f0.3[[3]]), function(x) SFS.function(CHR=top829f0.3[[3]][x,]$Chr, BEG=top829f0.3[[3]][x,]$Beg.Win, END=top829f0.3[[3]][x,]$End.Win, POP=3)))-> YRI.top.f0.3

unlist(mclapply(1:nrow(top829f0.3[[6]]), function(x) SFS.function(CHR=top829f0.3[[6]][x,]$Chr, BEG=top829f0.3[[6]][x,]$Beg.Win, END=top829f0.3[[6]][x,]$End.Win, POP=6)))-> GBR.top.f0.3

unlist(mclapply(1:nrow(top829f0.3[[7]]), function(x) SFS.function(CHR=top829f0.3[[7]][x,]$Chr, BEG=top829f0.3[[7]][x,]$Beg.Win, END=top829f0.3[[7]][x,]$End.Win, POP=7)))-> TSI.top.f0.3

#

unlist(mclapply(1:nrow(CANDf0.3[[2]]), function(x) SFS.function(CHR=CANDf0.3[[2]][x,]$Chr, BEG=CANDf0.3[[2]][x,]$Beg.Win, END=CANDf0.3[[2]][x,]$End.Win, POP=2)))-> LWK.cand.f0.3

unlist(mclapply(1:nrow(CANDf0.3[[3]]), function(x) SFS.function(CHR=CANDf0.3[[3]][x,]$Chr, BEG=CANDf0.3[[3]][x,]$Beg.Win, END=CANDf0.3[[3]][x,]$End.Win, POP=3)))-> YRI.cand.f0.3

unlist(mclapply(1:nrow(CANDf0.3[[6]]), function(x) SFS.function(CHR=CANDf0.3[[6]][x,]$Chr, BEG=CANDf0.3[[6]][x,]$Beg.Win, END=CANDf0.3[[6]][x,]$End.Win, POP=6)))-> GBR.cand.f0.3

unlist(mclapply(1:nrow(CANDf0.3[[7]]), function(x) SFS.function(CHR=CANDf0.3[[7]][x,]$Chr, BEG=CANDf0.3[[7]][x,]$Beg.Win, END=CANDf0.3[[7]][x,]$End.Win, POP=7)))-> TSI.cand.f0.3

#tf=0.4



unlist(mclapply(1:nrow(top829f0.4[[2]]), function(x) SFS.function(CHR=top829f0.4[[2]][x,]$Chr, BEG=top829f0.4[[2]][x,]$Beg.Win, END=top829f0.4[[2]][x,]$End.Win, POP=2)))-> LWK.top.f0.4

unlist(mclapply(1:nrow(top829f0.4[[3]]), function(x) SFS.function(CHR=top829f0.4[[3]][x,]$Chr, BEG=top829f0.4[[3]][x,]$Beg.Win, END=top829f0.4[[3]][x,]$End.Win, POP=3)))-> YRI.top.f0.4

unlist(mclapply(1:nrow(top829f0.4[[6]]), function(x) SFS.function(CHR=top829f0.4[[6]][x,]$Chr, BEG=top829f0.4[[6]][x,]$Beg.Win, END=top829f0.4[[6]][x,]$End.Win, POP=6)))-> GBR.top.f0.4

unlist(mclapply(1:nrow(top829f0.4[[7]]), function(x) SFS.function(CHR=top829f0.4[[7]][x,]$Chr, BEG=top829f0.4[[7]][x,]$Beg.Win, END=top829f0.4[[7]][x,]$End.Win, POP=7)))-> TSI.top.f0.4

#

unlist(mclapply(1:nrow(CANDf0.4[[2]]), function(x) SFS.function(CHR=CANDf0.4[[2]][x,]$Chr, BEG=CANDf0.4[[2]][x,]$Beg.Win, END=CANDf0.4[[2]][x,]$End.Win, POP=2)))-> LWK.cand.f0.4

unlist(mclapply(1:nrow(CANDf0.4[[3]]), function(x) SFS.function(CHR=CANDf0.4[[3]][x,]$Chr, BEG=CANDf0.4[[3]][x,]$Beg.Win, END=CANDf0.4[[3]][x,]$End.Win, POP=3)))-> YRI.cand.f0.4

unlist(mclapply(1:nrow(CANDf0.4[[6]]), function(x) SFS.function(CHR=CANDf0.4[[6]][x,]$Chr, BEG=CANDf0.4[[6]][x,]$Beg.Win, END=CANDf0.4[[6]][x,]$End.Win, POP=6)))-> GBR.cand.f0.4

unlist(mclapply(1:nrow(CANDf0.4[[7]]), function(x) SFS.function(CHR=CANDf0.4[[7]][x,]$Chr, BEG=CANDf0.4[[7]][x,]$Beg.Win, END=CANDf0.4[[7]][x,]$End.Win, POP=7)))-> TSI.cand.f0.4



#now neutral
filter(list.SCAN[[2]], Chr %in% c(21,22))-> chrs22_21
setDT(chrs22_21)

system.time(unlist(mclapply(1:nrow(chrs22_21), function(x) try(SFS.function(CHR=chrs22_21$Chr[x], BEG=chrs22_21$Beg.Win[x], END=chrs22_21$End.Win[x], POP=2))))-> genomicSFS) # 2004.006

genomicSFS[-(which(genomicSFS=="Error in read.table(out, header = F) : no lines available in input\n"))]-> test

as.numeric(test)-> genomic.SFS


Store(genomic.SFS);Store(LWK.top.f0.5, LWK.top.f0.3, LWK.top.f0.4)
Store(YRI.top.f0.5, YRI.top.f0.3, YRI.top.f0.4, GBR.top.f0.5, GBR.top.f0.4, GBR.top.f0.3, TSI.top.f0.5, TSI.top.f0.4, TSI.top.f0.3)
Store(YRI.cand.f0.5, YRI.cand.f0.3, YRI.cand.f0.4, GBR.cand.f0.5, GBR.cand.f0.4, GBR.cand.f0.3, TSI.cand.f0.5, TSI.cand.f0.4, TSI.cand.f0.3, LWK.cand.f0.5, LWK.cand.f0.4, LWK.cand.f0.3)
 

#TO DO: SFS FOR THE OTHER # POPS.


filter(list.SCAN[[3]], Chr %in% c(21,22))-> chrs22_21.YRI

setDT(chrs22_21.YRI)

system.time(unlist(mclapply(1:nrow(chrs22_21.YRI), function(x) try(SFS.function(CHR=chrs22_21.YRI$Chr[x], BEG=chrs22_21.YRI$Beg.Win[x], END=chrs22_21.YRI$End.Win[x], POP=3))))-> genomicSFS.YRI) # 

genomicSFS.YRI[-(which(genomicSFS.YRI=="Error in read.table(out, header = F) : no lines available in input\n"))]-> test

as.numeric(test)-> genomic.SFS.YRI



filter(list.SCAN[[6]], Chr %in% c(21,22))-> chrs22_21.GBR

setDT(chrs22_21.GBR)

system.time(unlist(mclapply(1:nrow(chrs22_21.GBR), function(x) try(SFS.function(CHR=chrs22_21.GBR$Chr[x], BEG=chrs22_21.GBR$Beg.Win[x], END=chrs22_21.GBR$End.Win[x], POP=6))))-> genomicSFS.GBR) # 

genomicSFS.GBR[-(which(genomicSFS.GBR=="Error in read.table(out, header = F) : no lines available in input\n"))]-> test

as.numeric(test)-> genomic.SFS.GBR



filter(list.SCAN[[7]], Chr %in% c(21,22))-> chrs22_21.TSI

setDT(chrs22_21.TSI)

system.time(unlist(mclapply(1:nrow(chrs22_21.TSI), function(x) try(SFS.function(CHR=chrs22_21.TSI$Chr[x], BEG=chrs22_21.TSI$Beg.Win[x], END=chrs22_21.TSI$End.Win[x], POP=7))))-> genomicSFS.TSI) # 2060.752

genomicSFS.TSI[-(which(genomicSFS.TSI=="Error in read.table(out, header = F) : no lines available in input\n"))]-> test

as.numeric(test)-> genomic.SFS.TSI


#actual plots
library(lattice)

pdf('figures/SFS.LWK.pdf')
par(mfrow=c(1,7))
histogram(genomicSFS[genomicSFS != 0 & genomicSFS !=100], col='lightgray', main='Neutral', xlab='DAF', breaks=seq(from=1,to=99,by=2)) #LWK
histogram(LWK.cand.f0.5[LWK.cand.f0.5 != 0 & LWK.cand.f0.5 !=100], col='cornflowerblue', main='Signficant tf=0.5', xlab='DAF',  breaks=seq(from=1,to=99,by=2))
histogram(LWK.top.f0.5[LWK.top.f0.5 != 0 & LWK.top.f0.5 !=100], col='cornflowerblue', main='Outlier tf=0.5', xlab='DAF', breaks=seq(from=1,to=99,by=2))

histogram(LWK.cand.f0.4[LWK.cand.f0.4 != 0 & LWK.cand.f0.4 !=100], col='sienna1', main='Signficant tf=0.4', xlab='DAF',  breaks=seq(from=1,to=99,by=2))
histogram(LWK.top.f0.4[LWK.top.f0.4 != 0 & LWK.top.f0.4 !=100], col='sienna1', main='Outlier tf=0.4', xlab='DAF', breaks=seq(from=1,to=99,by=2))

histogram(LWK.cand.f0.3[LWK.cand.f0.3 != 0 & LWK.cand.f0.3 !=100], col='violetred1', main='Signficant tf=0.3', xlab='DAF',  breaks=seq(from=1,to=99,by=2))
histogram(LWK.top.f0.3[LWK.top.f0.3 != 0 & LWK.top.f0.3 !=100], col='violetred1', main='Outlier tf=0.3', xlab='DAF', breaks=seq(from=1,to=99,by=2))
dev.off()

#
pdf('figures/SFS.YRI.pdf')
par(mfrow=c(1,7))
histogram(genomicSFS[genomicSFS != 0 & genomicSFS !=100], col='lightgray', main='Neutral', xlab='DAF', breaks=seq(from=1,to=99,by=2))  #LWK
histogram(YRI.cand.f0.5[YRI.cand.f0.5 != 0 & YRI.cand.f0.5 !=100], col='cornflowerblue', main='Signficant tf=0.5', xlab='DAF',  breaks=seq(from=1,to=99,by=2))
histogram(YRI.top.f0.5[YRI.top.f0.5 != 0 & YRI.top.f0.5 !=100], col='cornflowerblue', main='Outlier tf=0.5', xlab='DAF', breaks=seq(from=1,to=99,by=2))

histogram(YRI.cand.f0.4[YRI.cand.f0.4 != 0 & YRI.cand.f0.4 !=100], col='sienna1', main='Signficant tf=0.4', xlab='DAF',  breaks=seq(from=1,to=99,by=2))
histogram(YRI.top.f0.4[YRI.top.f0.4 != 0 & YRI.top.f0.4 !=100], col='sienna1', main='Outlier tf=0.4', xlab='DAF', breaks=seq(from=1,to=99,by=2))

histogram(YRI.cand.f0.3[YRI.cand.f0.3 != 0 & YRI.cand.f0.3 !=100], col='violetred1', main='Signficant tf=0.3', xlab='DAF',  breaks=seq(from=1,to=99,by=2))
histogram(YRI.top.f0.3[YRI.top.f0.3 != 0 & YRI.top.f0.3 !=100], col='violetred1', main='Outlier tf=0.3', xlab='DAF', breaks=seq(from=1,to=99,by=2))

dev.off()

#

pdf('figures/SFS.GBR.pdf')
par(mfrow=c(1,7))
histogram(genomicSFS[genomicSFS != 0 & genomicSFS !=100], col='lightgray', main='Neutral', xlab='DAF', breaks=seq(from=1,to=99,by=2)) #LWK
histogram(GBR.cand.f0.5[GBR.cand.f0.5 != 0 & GBR.cand.f0.5 !=100], col='cornflowerblue', main='Signficant tf=0.5', xlab='DAF',  breaks=seq(from=1,to=99,by=2))
histogram(GBR.top.f0.5[GBR.top.f0.5 != 0 & GBR.top.f0.5 !=100], col='cornflowerblue', main='Outlier tf=0.5', xlab='DAF', breaks=seq(from=1,to=99,by=2))

histogram(GBR.cand.f0.4[GBR.cand.f0.4 != 0 & GBR.cand.f0.4 !=100], col='sienna1', main='Signficant tf=0.4', xlab='DAF',  breaks=seq(from=1,to=99,by=2))
histogram(GBR.top.f0.4[GBR.top.f0.4 != 0 & GBR.top.f0.4 !=100], col='sienna1', main='Outlier tf=0.4', xlab='DAF', breaks=seq(from=1,to=99,by=2))

histogram(GBR.cand.f0.3[GBR.cand.f0.3 != 0 & GBR.cand.f0.3 !=100], col='violetred1', main='Signficant tf=0.3', xlab='DAF',  breaks=seq(from=1,to=99,by=2))
histogram(GBR.top.f0.3[GBR.top.f0.3 != 0 & GBR.top.f0.3 !=100], col='violetred1', main='Outlier tf=0.3', xlab='DAF', breaks=seq(from=1,to=99,by=2))

dev.off()



pdf('figures/SFS.f0.5.TSI.pdf')
par(mfrow=c(1,3))
histogram(genomicSFS[genomicSFS != 0 & genomicSFS !=100], col='lightgray', main='Neutral', xlab='DAF', breaks=seq(from=1,to=99,by=2)) #LWK
histogram(TSI.cand.f0.5[TSI.cand.f0.5 != 0 & TSI.cand.f0.5 !=100], col='cornflowerblue', main='Signficant tf=0.5', xlab='DAF',  breaks=seq(from=1,to=99,by=2))
histogram(TSI.top.f0.5[TSI.top.f0.5 != 0 & TSI.top.f0.5 !=100], col='cornflowerblue', main='Outlier tf=0.5', xlab='DAF', breaks=seq(from=1,to=99,by=2))


histogram(TSI.cand.f0.4[TSI.cand.f0.4 != 0 & TSI.cand.f0.4 !=100], col='sienna1', main='Signficant tf=0.4', xlab='DAF',  breaks=seq(from=1,to=99,by=2))
histogram(TSI.top.f0.4[TSI.top.f0.4 != 0 & TSI.top.f0.4 !=100], col='sienna1', main='Outlier tf=0.4', xlab='DAF', breaks=seq(from=1,to=99,by=2))

histogram(TSI.cand.f0.3[TSI.cand.f0.3 != 0 & TSI.cand.f0.3 !=100], col='violetred1', main='Signficant tf=0.3', xlab='DAF',  breaks=seq(from=1,to=99,by=2))
histogram(TSI.top.f0.3[TSI.top.f0.3 != 0 & TSI.top.f0.3 !=100], col='violetred1', main='Outlier tf=0.3', xlab='DAF', breaks=seq(from=1,to=99,by=2))


dev.off()


Store(genomicSFS, genomicSFS.GBR, genomicSFS.TSI, genomic.SFS.YRI)

##############
#write bed files

BED.PATH<-'/mnt/sequencedb/PopGen/barbara/NCV_dir_package/read_scan_data/bedfiles/'


write.table(select(Union.CANDf0.5_0.4_0.3[[3]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,"Union.CANDf0.5_0.4_0.3_YRI.bed"), quote=F, sep="\t", col.names=F, row.names=F) 

write.table(select(Union.CANDf0.5_0.4_0.3[[2]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,"Union.CANDf0.5_0.4_0.3_LWK.bed"), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(Union.CANDf0.5_0.4_0.3[[6]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,"Union.CANDf0.5_0.4_0.3_GBR.bed"), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(Union.CANDf0.5_0.4_0.3[[7]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,"Union.CANDf0.5_0.4_0.3_TSI.bed"), quote=F, sep="\t", col.names=F, row.names=F)




write.table(select(Union.top0.5_0.4_0.3[[3]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,"Union.top816.0.5_0.4_0.3_YRI.bed"), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(Union.top0.5_0.4_0.3[[2]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,"Union.top816.0.5_0.4_0.3_LWK.bed"), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(Union.top0.5_0.4_0.3[[6]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,"Union.top816.0.5_0.4_0.3_GBR.bed"), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(Union.top0.5_0.4_0.3[[7]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1),  file=paste0(BED.PATH,"Union.top816.0.5_0.4_0.3_TSI.bed"), quote=F, sep="\t", col.names=F, row.names=F)



write.table(select(top829f0.5[[2]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'top829f0.5.', pops[2], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(top829f0.5[[3]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'top829f0.5.', pops[3], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(top829f0.5[[6]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'top829f0.5.', pops[6], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(top829f0.5[[7]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'top829f0.5.', pops[7], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)


write.table(select(top829f0.4[[2]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'top829f0.4.', pops[2], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(top829f0.4[[3]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'top829f0.4.', pops[3], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(top829f0.4[[6]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'top829f0.4.', pops[6], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(top829f0.4[[7]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'top829f0.4.', pops[7], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)


write.table(select(top829f0.3[[2]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'top829f0.3.', pops[2], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(top829f0.3[[3]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'top829f0.3.', pops[3], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(top829f0.3[[6]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'top829f0.3.', pops[6], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(top829f0.3[[7]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'top829f0.3.', pops[7], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

#

write.table(select(CANDf0.5[[2]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'CANDf0.5.', pops[2], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(CANDf0.5[[3]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'CANDf0.5.', pops[3], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(CANDf0.5[[6]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'CANDf0.5.', pops[6], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(CANDf0.5[[7]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'CANDf0.5.', pops[7], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)


write.table(select(CANDf0.4[[2]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'CANDf0.4.', pops[2], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(CANDf0.4[[3]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'CANDf0.4.', pops[3], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(CANDf0.4[[6]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'CANDf0.4.', pops[6], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(CANDf0.4[[7]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'CANDf0.4.', pops[7], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)


write.table(select(CANDf0.3[[2]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'CANDf0.3.', pops[2], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(CANDf0.3[[3]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'CANDf0.3.', pops[3], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(CANDf0.3[[6]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'CANDf0.3.', pops[6], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(CANDf0.3[[7]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'CANDf0.3.', pops[7], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)








write.table(select(list.SCAN[[2]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'background_windows.bed'), quote=F, sep="\t", col.names=F, row.names=F)


####Scanned windows annotation ####


#sort ensembl bed file

read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/final_encode.bed')-> bed2

gsub("chr", "", bed2$V1)-> bed2$V1

arrange(bed2, V1, V2, V3)-> bed2

paste0("chr", bed2$V1)-> bed2$V1

#background windows

gsub("chr", "", list.SCAN[[2]]$Chr)-> list.SCAN[[2]]$Chr

arrange(list.SCAN[[2]], Chr, Beg.Win, End.Win)-> bed1

paste0("chr", bed1$Chr)-> bed1$Chr


setDT(bed1)




#this session below needs refinement...i mention proportion fo windows overlapping something, but these 'windows' are actually merged windows from the merge and then intersect sessions...need to find a way to make this clearer.


bedTools.merge(bed1=bed1)-> merge.scanned.windows

setDT(merge.scanned.windows)

with(merge.scanned.windows, cbind(merge.scanned.windows, Win.ID=paste0(V1,"|", V2, "|", V3)))-> merge.scanned.windows

length(unique(sort(merge.scanned.windows$Win.ID)))  #14,195 merged windows in the scan

system.time(bedTools.2in(bed1=merge.scanned.windows, bed2=bed2)-> intersect.scanned.windows) #64.045

setDT(intersect.scanned.windows)

#stopped here 27.10
length(unique(sort(intersect.scanned.windows$V8))) #48,254 number of 'coding elements' scanned

length(unique(sort(filter(intersect.scanned.windows, V11=="protein_coding")$V8))) #18,633 genes scanned

length(unique(intersect.scanned.windows$V4)) #number of merged windows overlapping coding elements 1,670, i,em 11670/14195=82% of the merged windows

length(unique(sort(filter(intersect.scanned.windows, V11=="protein_coding")$V4)))  #8514 number of merged windows scanning a prot. coding gene, i.e, 8514/14195=60%

length(unique(filter(intersect.scanned.windows, V11=="protein_coding")$V4))/length(unique(merge.scanned.windows$Win.ID)) # 60%  of (merged) background windows overlap protein_coding genes

mclapply(c(2,3,6,7), function(x) length(unique(filter(intersect.top829f0.5[[x]], V10=="protein_coding")$Win.ID))/length(unique(merge.top829f0.5[[x]]$Win.ID))) #for the top windows, around 70% overap prot coding genes

mclapply(c(2,3,6,7), function(x) length(unique(filter(intersect.CANDf0.5[[x]], V10=="protein_coding")$Win.ID))/length(unique(intersect.CANDf0.5[[x]]$Win.ID))) #and around 75% of the candidte windows
