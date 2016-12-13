##############################################
#	handling extreme windows
#	Barbara Bitarello
#	Creation: 12.10.2016
#	Last modified: 19.10.2016
#
##############################################

library(parallel)
library(SOAR)
library(ggplot2)
library(vioplot)
library(dplyr)
library(reshape)
library(reshape2)
Sys.setenv(R_LOCAL_CACHE="estsession")
library(data.table)

#this below I took from candidates_script_v2.r #stopped here 10.10.2016

pops<-c("AWS","LWK","YRI","CEU", "FIN","GBR","TSI", "CHB","CHS" ,"JPT","MXL", "CLM","PUR")

############################################## Part I #############################################
#This part is where I added info into the object list.SCAN. At the end of this part I stores this object, so 
#it can be skipped (go to 'Part II')
#RANK candidate windows

#system.time(mclapply(list.SCAN, function(x) cbind(x, Dist.NCV.f0.5=rep(NA, dim(x)[1]),Dist.NCV.f0.4=rep(NA, dim(x)[1]),Dist.NCV.f0.3=rep(NA, dim(x)[1]),Dist.NCV.f0.2=rep(NA, dim(x)[1])))-> list.SCAN.2)

Objects()

bin.list2<-vector('list', 7)

for (i in 1:3){
c(mclapply(bin.vec1, function(x) (which(list.SCAN[[i]]$Nr.IS==x))), list(which(list.SCAN[[i]]$Nr.IS>=bin.vec2)))->bin.list2[[i]]}

for (j in 4:7){
c(mclapply(bin.vec1.eu, function(x) (which(list.SCAN[[j]]$Nr.IS==x))), list(which(list.SCAN[[j]]$Nr.IS>=bin.vec2.eu)))->bin.list2[[j]]}

test.res<-vector('list', length(bin.list2[[1]]))

#calculate Z-scores
for (j in 1:3){ #AFRICA only
bin.list2[[j]][[length(bin.list2[[1]])]]->II  #first the last bin, which collapses all the remaining ones.
mean(l.bin.vec2[[1]]$ncvFD_f0.5)->mean5.II;sd(l.bin.vec2[[1]]$ncvFD_f0.5)-> sd5.II
mean(l.bin.vec2[[1]]$ncvFD_f0.4)->mean4.II;sd(l.bin.vec2[[1]]$ncvFD_f0.4)-> sd4.II;mean(l.bin.vec2[[1]]$ncvFD_f0.3)->mean3.II;sd(l.bin.vec2[[1]]$ncvFD_f0.3)-> sd3.II
mean(l.bin.vec2[[1]]$ncvFD_f0.2)->mean2.II;sd(l.bin.vec2[[1]]$ncvFD_f0.2)-> sd2.II;mean(l.bin.vec2[[1]]$ncvFD_f0.1)->mean1.II;sd(l.bin.vec2[[1]]$ncvFD_f0.1)-> sd1.II
gc()
if(length(II)>0){
((list.SCAN[[j]][II,]$NCDf5-mean5.II)/sd5.II)->list.SCAN[[j]]$Dist.NCD.f0.5[II]
((list.SCAN[[j]][II,]$NCDf4-mean4.II)/sd4.II)->list.SCAN[[j]]$Dist.NCD.f0.4[II];((list.SCAN[[j]][II,]$NCDf3-mean3.II)/sd3.II)->list.SCAN[[j]]$Dist.NCD.f0.3[II]
((list.SCAN[[j]][II,]$NCDf2-mean2.II)/sd2.II)->list.SCAN[[j]]$Dist.NCD.f0.2[II];((list.SCAN[[j]][II,]$NCDf1-mean1.II)/sd1.II)->list.SCAN[[j]]$Dist.NCD.f0.1[II]} #there was an error in the dist f0.1 before...
gc()

for (i in 1: (length(bin.list2[[1]])-1)){  #for all the other bins, except the last one
I<-bin.list2[[j]][[i]]
mean(l.bin.vec1[[i]]$ncvFD_f0.5)-> mean5.bin
mean(l.bin.vec1[[i]]$ncvFD_f0.4)-> mean4.bin
mean(l.bin.vec1[[i]]$ncvFD_f0.3)-> mean3.bin
mean(l.bin.vec1[[i]]$ncvFD_f0.2)-> mean2.bin
mean(l.bin.vec1[[i]]$ncvFD_f0.1)-> mean1.bin
sd(l.bin.vec1[[i]]$ncvFD_f0.5)-> sd5.bin
sd(l.bin.vec1[[i]]$ncvFD_f0.4)-> sd4.bin
sd(l.bin.vec1[[i]]$ncvFD_f0.3)-> sd3.bin
sd(l.bin.vec1[[i]]$ncvFD_f0.2)-> sd2.bin
sd(l.bin.vec1[[i]]$ncvFD_f0.1)-> sd1.bin
if(length(I)>0){
((list.SCAN[[j]][I,]$NCDf5-mean5.bin)/sd5.bin)->list.SCAN[[j]]$Dist.NCD.f0.5[I]
((list.SCAN[[j]][I,]$NCDf4-mean4.bin)/sd4.bin)->list.SCAN[[j]]$Dist.NCD.f0.4[I]
((list.SCAN[[j]][I,]$NCDf2-mean2.bin)/sd2.bin)->list.SCAN[[j]]$Dist.NCD.f0.2[I]
((list.SCAN[[j]][I,]$NCDf3-mean3.bin)/sd3.bin)->list.SCAN[[j]]$Dist.NCD.f0.3[I]
((list.SCAN[[j]][I,]$NCDf1-mean1.bin)/sd1.bin)->list.SCAN[[j]]$Dist.NCD.f0.1[I]
} 
cat('bin', i, 'done\n');
} 
cat('pop', j, 'done\n')}



#it works. now add this as as a " distance" collumn to the candidate windows data sets.
##############
#Now Europe#
##############

system.time(for (j in 4:7){ #Europe only
bin.list2[[j]][[length(bin.list2[[4]])]]->II
mean(l.bin.vec2.eu[[1]]$ncvFD_f0.5)->mean5.II;sd(l.bin.vec2.eu[[1]]$ncvFD_f0.5)-> sd5.II
mean(l.bin.vec2.eu[[1]]$ncvFD_f0.4)->mean4.II;sd(l.bin.vec2.eu[[1]]$ncvFD_f0.4)-> sd4.II
mean(l.bin.vec2.eu[[1]]$ncvFD_f0.3)->mean3.II;sd(l.bin.vec2.eu[[1]]$ncvFD_f0.3)-> sd3.II
mean(l.bin.vec2.eu[[1]]$ncvFD_f0.2)->mean2.II;sd(l.bin.vec2.eu[[1]]$ncvFD_f0.2)-> sd2.II
mean(l.bin.vec2.eu[[1]]$ncvFD_f0.1)->mean1.II;sd(l.bin.vec2.eu[[1]]$ncvFD_f0.1)-> sd1.II
gc()
if(length(II)>0){
((list.SCAN[[j]][II,]$NCDf5-mean5.II)/sd5.II)->list.SCAN[[j]]$Dist.NCD.f0.5[II]
((list.SCAN[[j]][II,]$NCDf4-mean4.II)/sd4.II)->list.SCAN[[j]]$Dist.NCD.f0.4[II]
((list.SCAN[[j]][II,]$NCDf3-mean3.II)/sd3.II)->list.SCAN[[j]]$Dist.NCD.f0.3[II]
((list.SCAN[[j]][II,]$NCDf2-mean2.II)/sd2.II)->list.SCAN[[j]]$Dist.NCD.f0.2[II]
((list.SCAN[[j]][II,]$NCDf1-mean1.II)/sd1.II)->list.SCAN[[j]]$Dist.NCD.f0.1[II]}
gc()
for (i in 1: (length(bin.list2[[4]])-1)){ #for all the bins except the last one (European pops)
I<-bin.list2[[j]][[i]]
mean(l.bin.vec1.eu[[i]]$ncvFD_f0.5)-> mean5.bin
mean(l.bin.vec1.eu[[i]]$ncvFD_f0.4)-> mean4.bin;mean(l.bin.vec1.eu[[i]]$ncvFD_f0.3)-> mean3.bin
mean(l.bin.vec1.eu[[i]]$ncvFD_f0.2)-> mean2.bin;mean(l.bin.vec1.eu[[i]]$ncvFD_f0.1)-> mean1.bin

sd(l.bin.vec1.eu[[i]]$ncvFD_f0.5)-> sd5.bin;sd(l.bin.vec1.eu[[i]]$ncvFD_f0.4)-> sd4.bin
sd(l.bin.vec1.eu[[i]]$ncvFD_f0.3)-> sd3.bin
sd(l.bin.vec1.eu[[i]]$ncvFD_f0.2)-> sd2.bin;sd(l.bin.vec1.eu[[i]]$ncvFD_f0.1)-> sd1.bin

gc()
if(length(I)>0){
((list.SCAN[[j]][I,]$NCDf5-mean5.bin)/sd5.bin)->list.SCAN[[j]]$Dist.NCD.f0.5[I]
((list.SCAN[[j]][I,]$NCDf4-mean4.bin)/sd4.bin)->list.SCAN[[j]]$Dist.NCD.f0.4[I]
((list.SCAN[[j]][I,]$NCDf3-mean3.bin)/sd3.bin)->list.SCAN[[j]]$Dist.NCD.f0.3[I]
((list.SCAN[[j]][I,]$NCDf2-mean2.bin)/sd2.bin)->list.SCAN[[j]]$Dist.NCD.f0.2[I]
((list.SCAN[[j]][I,]$NCDf1-mean1.bin)/sd1.bin)->list.SCAN[[j]]$Dist.NCD.f0.1[I]
}
cat('bin', i, 'done\n');
}
cat('pop', j, 'done\n');
})

#

save(list.SCAN, file='list.SCAN.2.RData')

Store(list.SCAN)
#IDEA for the future: save workspace and copy to darwin so I can use Debora's 1000G annotation.
# Sometimes it is necessary to store some stuff, otherwise the session crashes.
#####################################################
Objects()
#p-values based on Z scores.

mclapply(list.SCAN, function(x) cbind(arrange(x, Dist.NCD.f0.5),Z.f0.5.P.val=seq(1:nrow(x))/nrow(x)))-> tmp5
mclapply(1:7, function(x) cbind(arrange(tmp5[[x]],Dist.NCD.f0.4), Z.f0.4.P.val=seq(1:nrow(tmp5[[x]]))/nrow(tmp5[[x]])))->tmp4
remove(tmp5);gc()


mclapply(1:7, function(x) cbind(arrange(tmp4[[x]],Dist.NCD.f0.3), Z.f0.3.P.val=seq(1:nrow(tmp4[[x]]))/nrow(tmp4[[x]])))->tmp3
remove(tmp4);gc()
mclapply(1:7, function(x) cbind(arrange(tmp3[[x]],Dist.NCD.f0.2), Z.f0.2.P.val=seq(1:nrow(tmp3[[x]]))/nrow(tmp3[[x]])))->tmp2
remove(tmp3); gc()
mclapply(1:7, function(x) cbind(arrange(tmp2[[x]],Dist.NCD.f0.1), Z.f0.1.P.val=seq(1:nrow(tmp2[[x]]))/nrow(tmp2[[x]])))->list.SCAN.3
remove(tmp2);gc()

#check list.SCAN.2

list.SCAN.3-> list.SCAN
remove(list.SCAN.3)
gc()

mclapply(list.SCAN, function(x) arrange(x, Chr, Beg.Win))-> list.SCAN.3

#check

list.SCAN.3-> list.SCAN
remove(list.SCAN.3)
gc()

save(list.SCAN, file='list.SCAN.3.RData')
Store(list.SCAN)

#take simulation-based candidate windows.
mclapply(list.SCAN, function(x) x[which(x$P.val.NCDf0.5<(1/nsims)),])-> CANDf0.5; names(CANDf0.5)<-pops[1:7]; mclapply(CANDf0.5, function(x) arrange(x, Chr, Beg.Win))-> CANDf0.5
mclapply(list.SCAN, function(x) x[which(x$P.val.NCDf0.4<(1/nsims)),])-> CANDf0.4; names(CANDf0.4)<-pops[1:7]; mclapply(CANDf0.4, function(x) arrange(x, Chr, Beg.Win))-> CANDf0.4
mclapply(list.SCAN, function(x) x[which(x$P.val.NCDf0.3<(1/nsims)),])-> CANDf0.3; names(CANDf0.3)<-pops[1:7]; mclapply(CANDf0.3, function(x) arrange(x, Chr, Beg.Win))-> CANDf0.3
mclapply(list.SCAN, function(x) x[which(x$P.val.NCDf0.2<(1/nsims)),])-> CANDf0.2; names(CANDf0.2)<-pops[1:7]; mclapply(CANDf0.2, function(x) arrange(x, Chr, Beg.Win))-> CANDf0.2
mclapply(list.SCAN, function(x) x[which(x$P.val.NCDf0.1<(1/nsims)),])-> CANDf0.1; names(CANDf0.1)<-pops[1:7]; mclapply(CANDf0.1, function(x) arrange(x, Chr, Beg.Win))-> CANDf0.1

unlist(mclapply(CANDf0.5, function(x) nrow(x)))

#AWS  LWK  YRI  CEU  FIN  GBR  TSI 
#6088 6226 6854 6540 6074 6519 6403 

unlist(mclapply(CANDf0.4, function(x) nrow(x)))

#AWS  LWK  YRI  CEU  FIN  GBR  TSI 
#6596 6841 7420 6859 6587 6858 6661

unlist(mclapply(CANDf0.3, function(x) nrow(x)))

#AWS  LWK  YRI  CEU  FIN  GBR  TSI 
#6894 7364 8098 6145 6618 6173 5998 



Store(CANDf0.5); Store(CANDf0.4); Store(CANDf0.3); Store(CANDf0.2); Store(CANDf0.1);Store(list.SCAN) #now list.SCAN has everything I need.



#outlier windows


#nrow(list.SCAN[[1]]) #1657989 windows * 0.0005 == iapprox. 829

nrow(list.SCAN[[3]])*0.0005
#[1] 828.9945

top829f0.5<-mclapply(list.SCAN, function(x) arrange(x, Z.f0.5.P.val)[1:829,]) #top 829 windows ranked by tf=0.5
top829f0.4<-mclapply(list.SCAN, function(x) arrange(x, Z.f0.4.P.val)[1:829,]) #top 829 windows ranked by tf=0.4
top829f0.3<-mclapply(list.SCAN, function(x) arrange(x, Z.f0.3.P.val)[1:829,]) #top 829 windows ranked by tf=0.3
top829f0.2<-mclapply(list.SCAN, function(x) arrange(x, Z.f0.2.P.val)[1:829,]) #top 829 windows ranked by tf=0.2
top829f0.1<-mclapply(list.SCAN, function(x) arrange(x, Z.f0.1.P.val)[1:829,]) #top 829 windows ranked by tf=0.1

gc()

unlist(mclapply(top829f0.3, function(x) nrow(x)))
unlist(mclapply(top829f0.4, function(x) nrow(x)))
unlist(mclapply(top829f0.5, function(x) nrow(x)))


save(top829f0.5, file='top829f0.5.RData')

save(top829f0.4, file='top829f0.4.RData')

save(top829f0.3, file='top829f0.3.RData')

gc()



save(CANDf0.5, file='CANDf0.5.RData')

save(CANDf0.4, file='CANDf0.4.RData')

save(CANDf0.3, file='CANDf0.3.RData')


gc()

Store(top829f0.5)
Store(top829f0.4)
Store(top829f0.3)

#In  the next part we can start exploring these windows.
# The End


##*****************************************************************************

#bedtools in R
source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/read_scan_data/bedtools_inR.R')

#sort ensembl bed file

read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/final_encode.bed')-> bed2

gsub("chr", "", bed2$V1)-> bed2$V1

arrange(bed2, V1, V2, V3)-> bed2

paste0("chr", bed2$V1)-> bed2$V1

#top763 windows

####### MERGE ######## MERGE ######## MERGE ###########


mclapply(top829f0.5, function(x) arrange(x, Chr, Beg.Win))-> top829f0.5

mclapply(top829f0.4, function(x) arrange(x, Chr, Beg.Win))-> top829f0.4

mclapply(top829f0.3, function(x) arrange(x, Chr, Beg.Win))-> top829f0.3

#

mclapply(CANDf0.5, function(x) arrange(x, Chr, Beg.Win))-> CANDf0.5

mclapply(CANDf0.4, function(x) arrange(x, Chr, Beg.Win))-> CANDf0.4

mclapply(CANDf0.3, function(x) arrange(x, Chr, Beg.Win))-> CANDf0.3


for(i in 1:7){

with(top829f0.5[[i]], paste0("chr", Chr))-> top829f0.5[[i]]$Chr}

mclapply(1:7, function(x) bedTools.merge(bed1=select(top829f0.5[[x]], Chr, Beg.Win, End.Win, Win.ID)))-> merge.top829f0.5


for(i in 1:7){

with(CANDf0.5[[i]], paste0("chr", Chr))-> CANDf0.5[[i]]$Chr}

mclapply(1:7, function(x) bedTools.merge(bed1=select(CANDf0.5[[x]],Chr, Beg.Win, End.Win, Win.ID)))-> merge.CANDf0.5


for(i in 1:7){

with(top829f0.4[[i]], paste0("chr", Chr))-> top829f0.4[[i]]$Chr}

mclapply(1:7, function(x) bedTools.merge(bed1=select(top829f0.4[[x]],Chr, Beg.Win, End.Win, Win.ID)))-> merge.top829f0.4



for(i in 1:7){

with(CANDf0.4[[i]], paste0("chr", Chr))-> CANDf0.4[[i]]$Chr}

mclapply(1:7, function(x) bedTools.merge(bed1=select(CANDf0.4[[x]],Chr, Beg.Win, End.Win, Win.ID)))-> merge.CANDf0.4

#
for(i in 1:7){

with(top829f0.3[[i]], paste0("chr", Chr))-> top829f0.3[[i]]$Chr}

mclapply(1:7, function(x) bedTools.merge(bed1=select(top829f0.3[[x]], Chr, Beg.Win, End.Win, Win.ID)))-> merge.top829f0.3


for(i in 1:7){

with(CANDf0.3[[i]], paste0("chr", Chr))-> CANDf0.3[[i]]$Chr}

mclapply(1:7, function(x) bedTools.merge(bed1=select(CANDf0.3[[x]], Chr, Beg.Win, End.Win, Win.ID)))-> merge.CANDf0.3

####### MERGE ####### MERGE ######## MERGE #############


#### INTERSECT ##### INTERSECT ##### INTERSECT ########

mclapply(1:7, function(x) bedTools.2in(bed1=merge.top829f0.5[[x]], bed2=bed2))-> intersect.top829f0.5
mclapply(1:7, function(x) bedTools.2in(bed1=merge.top829f0.4[[x]], bed2=bed2))-> intersect.top829f0.4
mclapply(1:7, function(x) bedTools.2in(bed1=merge.top829f0.3[[x]], bed2=bed2))-> intersect.top829f0.3


mclapply(1:7, function(x) bedTools.2in(bed1=merge.CANDf0.5[[x]], bed2=bed2))-> intersect.CANDf0.5
mclapply(1:7, function(x) bedTools.2in(bed1=merge.CANDf0.4[[x]], bed2=bed2))-> intersect.CANDf0.4
mclapply(1:7, function(x) bedTools.2in(bed1=merge.CANDf0.3[[x]], bed2=bed2))-> intersect.CANDf0.3


Store(merge.CANDf0.5, merge.CANDf0.4, merge.CANDf0.3)
Store(intersect.CANDf0.5, intersect.CANDf0.4, intersect.CANDf0.3)
Store(merge.top829f0.5, merge.top829f0.4, merge.top829f0.3)
Store(intersect.top829f0.5, intersect.top829f0.4, intersect.top829f0.3)
############################################################
############################################################
#other scans
andres_AA<- read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/andres.2009.AA.bed')
andres_AA[1:15,]->andres_AA #remove triple entries
andres_EA<- read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/andres.2009.EA.bed')
andres_EA[1:31,]->andres_EA #remove triple entries
andres_AAandEA<- read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/andres.2009.AAandEA.bed')
andres_AAandEA[1:12,]->andres_AAandEA #remove triple entries (this is actually a problem with the inpur file)
DG_T2_YRI<-read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/DG.2014.T2.YRI.bed')
#this file has double entries for each line...
DG_T2_YRI[1:99,]->DG_T2_YRI
#DG_T2_YRI[-40,]->DG_T2_YRI
DG_T2_CEU<-read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/DG.2014.T2.CEU.bed')
DG_T1_CEU<-read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/DG.T1.CEU')
DG_T1_YRI<-read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/DG.T1.YRI')

names(andres_AA)<-c('chr','B', 'E', 'Name')
names(andres_EA)<-c('chr','B', 'E', 'Name')
names(andres_AAandEA)<-c('chr','B', 'E', 'Name')
names(DG_T2_YRI)<-c('chr','B', 'E', 'Name')
names(DG_T2_CEU)<-c('chr','B', 'E', 'Name')
all.andres<-unique(sort(as.character(rbind(andres_AA, andres_EA, andres_AAandEA)[,4])))
all.DG<-unique(c(as.character(rbind(DG_T2_CEU, DG_T2_YRI)[,4]), as.character(rbind(DG_T1_CEU, DG_T1_YRI)[,1])))
leffler<-c("FREM3", "HUS1", "MTRR","IGFBP7","PROKR2","ST3GAL1")
other.bal.sel<-unique(c(leffler, all.DG, all.andres))

Store(DG_T1_YRI, DG_T1_CEU, DG_T2_YRI, DG_T2_CEU, all.andres, leffler, all.DG, andres_AA, andres_EA, andres_AAandEA, other.bal.sel)

#################################################################


