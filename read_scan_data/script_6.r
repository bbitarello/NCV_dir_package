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

system.time(unlist(mclapply(1:nrow(chrs22_21.TSI), function(x) try(SFS.function(CHR=chrs22_21.TSI$Chr[x], BEG=chrs22_21.TSI$Beg.Win[x], END=chrs22_21.TSI$End.Win[x], POP=7))))-> genomicSFS.TSI) # 

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



#add Min.NCD.tf collun (which tf yields lowest NCD value)



test.col<-vector('list', 7)
for(i in 1:7){
mclapply(1:nrow(Union.CANDf0.5_0.4_0.3[[i]]),  function(y) colnames(select(Union.CANDf0.5_0.4_0.3[[i]][y,], Z.f0.5.P.val:Z.f0.3.P.val))[which(sapply(1:3, function(x) select(Union.CANDf0.5_0.4_0.3[[i]][y,], Z.f0.5.P.val:Z.f0.3.P.val)[,x] ==min(select(Union.CANDf0.5_0.4_0.3[[i]][y,], Z.f0.5.P.val:Z.f0.3.P.val))))])-> test.col[[i]]
which(unlist(lapply(test.col[[i]], function(x) length(x)==2)))-> repl2
for(j in 1: length(repl2)){
paste0(test.col[[i]][[repl2[j]]][1], "|", test.col[[i]][[repl2[j]]][2])-> test.col[[i]][[repl2[j]]]}
which(unlist(lapply(test.col[[i]], function(x) length(x)==3)))-> repl3
if(length(repl3)>=1){
for(j in 1: length(repl3)){
paste0(test.col[[i]][[repl3[j]]][1], "|", test.col[[i]][[repl3[j]]][2])-> test.col[[i]][[repl3[j]]]}}
which(unlist(lapply(test.col[[i]], function(x) length(x)==4)))-> repl4
if(length(repl4)>=1){
for(j in 1: length(repl4)){
paste0(test.col[[i]][[repl4[j]]][1], "|", test.col[[i]][[repl4[j]]][2])-> test.col[[i]][[repl4[j]]]}}
which(unlist(lapply(test.col[[i]], function(x) length(x)==5)))-> repl5
if(length(repl5)>=1){
for(j in 1: length(repl5)){
paste0(test.col[[i]][[repl5[j]]][1], "|", test.col[[i]][[repl5[j]]][2])-> test.col[[i]][[repl5[j]]]}}
cat('Done with', pops[i],'\n')
}

test.col2<-vector('list', 7)
for(i in 1:7){
mclapply(1:nrow(Union.top0.5_0.4_0.3[[i]]),  function(y) colnames(select(Union.top0.5_0.4_0.3[[i]][y,], Z.f0.5.P.val:Z.f0.3.P.val))[which(sapply(1:3, function(x) select(Union.top0.5_0.4_0.3[[i]][y,], Z.f0.5.P.val:Z.f0.3.P.val)[,x] ==min(select(Union.top0.5_0.4_0.3[[i]][y,], Z.f0.5.P.val:Z.f0.3.P.val))))])-> test.col2[[i]]
which(unlist(lapply(test.col2[[i]], function(x) length(x)==2)))-> repl2
for(j in 1: length(repl2)){
paste0(test.col2[[i]][[repl2[j]]][1], "|", test.col2[[i]][[repl2[j]]][2])-> test.col2[[i]][[repl2[j]]]}
which(unlist(lapply(test.col2[[i]], function(x) length(x)==3)))-> repl3
if(length(repl3)>=1){
for(j in 1: length(repl3)){
paste0(test.col2[[i]][[repl3[j]]][1], "|", test.col2[[i]][[repl3[j]]][2])-> test.col2[[i]][[repl3[j]]]}}
which(unlist(lapply(test.col2[[i]], function(x) length(x)==4)))-> repl4
if(length(repl4)>=1){
for(j in 1: length(repl4)){
paste0(test.col2[[i]][[repl4[j]]][1], "|", test.col2[[i]][[repl4[j]]][2])-> test.col2[[i]][[repl4[j]]]}}
which(unlist(lapply(test.col2[[i]], function(x) length(x)==5)))-> repl5
if(length(repl5)>=1){
for(j in 1: length(repl5)){
paste0(test.col[[i]][[repl2[j]]][1], "|", test.col[[i]][[repl2[j]]][2])-> test.col[[i]][[repl2[j]]]}
which(unlist(lapply(test.col[[i]], function(x) length(x)==3)))-> repl3
if(length(repl3)>=1){
for(j in 1: length(repl3)){
paste0(test.col[[i]][[repl3[j]]][1], "|", test.col[[i]][[repl3[j]]][2])-> test.col[[i]][[repl3[j]]]}}
which(unlist(lapply(test.col[[i]], function(x) length(x)==4)))-> repl4
if(length(repl4)>=1){
for(j in 1: length(repl4)){
paste0(test.col[[i]][[repl4[j]]][1], "|", test.col[[i]][[repl4[j]]][2])-> test.col[[i]][[repl4[j]]]}}
which(unlist(lapply(test.col[[i]], function(x) length(x)==5)))-> repl5
if(length(repl5)>=1){
for(j in 1: length(repl5)){
paste0(test.col[[i]][[repl5[j]]][1], "|", test.col[[i]][[repl5[j]]][2])-> test.col[[i]][[repl5[j]]]}}
}

test.col2<-vector('list', 7)
for(i in 1:7){
mclapply(1:nrow(Union.top0.5_0.4_0.3[[i]]),  function(y) colnames(select(Union.top0.5_0.4_0.3[[i]][y,], Z.f0.5.P.val:Z.f0.3.P.val))[which(sapply(1:3, function(x) select(Union.top0.5_0.4_0.3[[i]][y,], Z.f0.5.P.val:Z.f0.3.P.val)[,x] ==min(select(Union.top0.5_0.4_0.3[[i]][y,], Z.f0.5.P.val:Z.f0.3.P.val))))])-> test.col2[[i]]
which(unlist(lapply(test.col2[[i]], function(x) length(x)==2)))-> repl2
for(j in 1: length(repl2)){
paste0(test.col2[[i]][[repl2[j]]][1], "|", test.col2[[i]][[repl2[j]]][2])-> test.col2[[i]][[repl2[j]]]}
which(unlist(lapply(test.col2[[i]], function(x) length(x)==3)))-> repl3
if(length(repl3)>=1){
for(j in 1: length(repl3)){
paste0(test.col2[[i]][[repl3[j]]][1], "|", test.col2[[i]][[repl3[j]]][2])-> test.col2[[i]][[repl3[j]]]}}
which(unlist(lapply(test.col2[[i]], function(x) length(x)==4)))-> repl4
if(length(repl4)>=1){
for(j in 1: length(repl4)){
paste0(test.col2[[i]][[repl4[j]]][1], "|", test.col2[[i]][[repl4[j]]][2])-> test.col2[[i]][[repl4[j]]]}}
which(unlist(lapply(test.col2[[i]], function(x) length(x)==5)))-> repl5
if(length(repl5)>=1){
for(j in 1: length(repl5)){
paste0(test.col2[[i]][[repl4[j]]][1], "|", test.col2[[i]][[repl4[j]]][2])-> test.col2[[i]][[repl4[j]]]}}
which(unlist(lapply(test.col2[[i]], function(x) length(x)==5)))-> repl5
if(length(repl5)>=1){
for(j in 1: length(repl5)){
paste0(test.col2[[i]][[repl5[j]]][1], "|", test.col2[[i]][[repl5[j]]][2])-> test.col2[[i]][[repl5[j]]]}}
}

Store(test.col, test.col2)

mclapply(1:7, function(x) gsub(".P.val","", gsub("Z.f", "", unlist(test.col[[x]]))))-> extra.col

mclapply(1:7, function(x) gsub(".P.val","", gsub("Z.f", "", unlist(test.col2[[x]]))))-> extra.col.top


rbind(table(as.numeric(extra.col[[1]])), table(as.numeric(extra.col[[2]])), table(as.numeric(extra.col[[3]])), table(as.numeric(extra.col[[4]])), table(as.numeric(extra.col[[5]])), table(as.numeric(extra.col[[6]])), table(as.numeric(extra.col[[7]])))

for (i in 1:7){

cbind(Union.CANDf0.5_0.4_0.3[[i]], Min.ZPval.Feq=extra.col[[i]])-> Union.CANDf0.5_0.4_0.3[[i]]
cbind(Union.top0.5_0.4_0.3[[i]],  Min.ZPval.Feq=extra.col.top[[i]])-> Union.top0.5_0.4_0.3[[i]]}


Store(Union.CANDf0.5_0.4_0.3)
Store(Union.top0.5_0.4_0.4)


