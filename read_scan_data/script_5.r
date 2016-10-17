##############################################
#	handling extreme windows
#	Barbara Bitarello
#	Creation: 12.10.2016
#	Last modified: 12.10.2016
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

if(length(II)>0){
((list.SCAN[[j]][II,]$NCDf5-mean5.II)/sd5.II)->list.SCAN[[j]]$Dist.NCD.f0.5[II]
((list.SCAN[[j]][II,]$NCDf4-mean4.II)/sd4.II)->list.SCAN[[j]]$Dist.NCD.f0.4[II];((list.SCAN[[j]][II,]$NCDf3-mean3.II)/sd3.II)->list.SCAN[[j]]$Dist.NCD.f0.3[II]
((list.SCAN[[j]][II,]$NCDf2-mean2.II)/sd2.II)->list.SCAN[[j]]$Dist.NCD.f0.2[II];((list.SCAN[[j]][II,]$NCDf1-mean1.II)/sd1.II)->list.SCAN[[j]]$Dist.NCD.f0.1[II]} #there was an error in the dist f0.1 before...

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
((list.SCAN[[j]][I,]$NCDf3-mean3.bin)/sd3.bin)->list.SCAN[[j]]$Dist.NCD.f0.3[I]
((list.SCAN[[j]][I,]$NCDf2-mean2.bin)/sd2.bin)->list.SCAN[[j]]$Dist.NCD.f0.2[I]
((list.SCAN[[j]][I,]$NCDf1-mean1.bin)/sd1.bin)->list.SCAN[[j]]$Dist.NCD.f0.1[I]
}}}


#it works. now add this as as a " distance" collumn to the candidate windows data sets.
##############
#Now Europe#
##############

for (j in 4:7){ #Europe only
bin.list2[[j]][[length(bin.list2[[4]])]]->II
mean(l.bin.vec2.eu[[1]]$ncvFD_f0.5)->mean5.II;sd(l.bin.vec2.eu[[1]]$ncvFD_f0.5)-> sd5.II
mean(l.bin.vec2[[1]]$ncvFD_f0.4)->mean4.II;sd(l.bin.vec2[[1]]$ncvFD_f0.4)-> sd4.II
mean(l.bin.vec2[[1]]$ncvFD_f0.3)->mean3.II;sd(l.bin.vec2[[1]]$ncvFD_f0.3)-> sd3.II
mean(l.bin.vec2[[1]]$ncvFD_f0.2)->mean2.II;sd(l.bin.vec2[[1]]$ncvFD_f0.2)-> sd2.II
mean(l.bin.vec2[[1]]$ncvFD_f0.1)->mean1.II;sd(l.bin.vec2[[1]]$ncvFD_f0.1)-> sd1.II

if(length(II)>0){
((list.SCAN[[j]][II,]$NCDf5-mean5.II)/sd5.II)->list.SCAN[[j]]$Dist.NCD.f0.5[II]
((list.SCAN[[j]][II,]$NCDf4-mean4.II)/sd4.II)->list.SCAN[[j]]$Dist.NCD.f0.4[II]
((list.SCAN[[j]][II,]$NCDf3-mean3.II)/sd3.II)->list.SCAN[[j]]$Dist.NCD.f0.3[II]
((list.SCAN[[j]][II,]$NCDf2-mean2.II)/sd2.II)->list.SCAN[[j]]$Dist.NCD.f0.2[II]
((list.SCAN[[j]][II,]$NCDf1-mean1.II)/sd1.II)->list.SCAN[[j]]$Dist.NCD.f0.1[II]}

for (i in 1: (length(bin.list2[[4]])-1)){ #for all the bins except the last one (European pops)
I<-bin.list2[[j]][[i]]
mean(l.bin.vec1.eu[[i]]$ncvFD_f0.5)-> mean5.bin
mean(l.bin.vec1.eu[[i]]$ncvFD_f0.4)-> mean4.bin;mean(l.bin.vec1.eu[[i]]$ncvFD_f0.3)-> mean3.bin
mean(l.bin.vec1.eu[[i]]$ncvFD_f0.2)-> mean2.bin;mean(l.bin.vec1.eu[[i]]$ncvFD_f0.1)-> mean1.bin

sd(l.bin.vec1.eu[[i]]$ncvFD_f0.5)-> sd5.bin;sd(l.bin.vec1.eu[[i]]$ncvFD_f0.4)-> sd4.bin
sd(l.bin.vec1.eu[[i]]$ncvFD_f0.3)-> sd3.bin
sd(l.bin.vec1.eu[[i]]$ncvFD_f0.2)-> sd2.bin;sd(l.bin.vec1.eu[[i]]$ncvFD_f0.1)-> sd1.bin


if(length(I)>0){
((list.SCAN[[j]][I,]$NCDf5-mean5.bin)/sd5.bin)->list.SCAN[[j]]$Dist.NCD.f0.5[I]
((list.SCAN[[j]][I,]$NCDf4-mean4.bin)/sd4.bin)->list.SCAN[[j]]$Dist.NCD.f0.4[I]
((list.SCAN[[j]][I,]$NCDf3-mean3.bin)/sd3.bin)->list.SCAN[[j]]$Dist.NCD.f0.3[I]
((list.SCAN[[j]][I,]$NCDf2-mean2.bin)/sd2.bin)->list.SCAN[[j]]$Dist.NCD.f0.2[I]
((list.SCAN[[j]][I,]$NCDf1-mean1.bin)/sd1.bin)->list.SCAN[[j]]$Dist.NCD.f0.1[I]
}}}




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
mclapply(1:7, function(x) cbind(arrange(tmp2[[x]],Dist.NCD.f0.1), Z.f0.1.P.val=seq(1:nrow(tmp2[[x]]))/nrow(tmp2[[x]])))->list.SCAN
remove(tmp2);gc()


mclapply(list.SCAN, function(x) arrange(x, Chr, Beg.Win))-> list.SCAN
Store(list.SCAN)

#mclapply(1:7, function(x) with(list.SCAN[[x]], paste0(Chr, "|", Beg.Win, "|", End.Win)))-> Win.ID.scan
#take simulation-based candidate windows.
mclapply(list.SCAN, function(x) x[which(x$P.val.NCDf0.5<(1/nsims)),])-> CANDf0.5; names(CANDf0.5)<-pops[1:7]; mclapply(CANDf0.5, function(x) arrange(x, Chr, Beg.Win))-> CANDf0.5
mclapply(list.SCAN, function(x) x[which(x$P.val.NCDf0.4<(1/nsims)),])-> CANDf0.4; names(CANDf0.4)<-pops[1:7]; mclapply(CANDf0.4, function(x) arrange(x, Chr, Beg.Win))-> CANDf0.4
mclapply(list.SCAN, function(x) x[which(x$P.val.NCDf0.3<(1/nsims)),])-> CANDf0.3; names(CANDf0.3)<-pops[1:7]; mclapply(CANDf0.3, function(x) arrange(x, Chr, Beg.Win))-> CANDf0.3
mclapply(list.SCAN, function(x) x[which(x$P.val.NCDf0.2<(1/nsims)),])-> CANDf0.2; names(CANDf0.2)<-pops[1:7]; mclapply(CANDf0.2, function(x) arrange(x, Chr, Beg.Win))-> CANDf0.2
mclapply(list.SCAN, function(x) x[which(x$P.val.NCDf0.1<(1/nsims)),])-> CANDf0.1; names(CANDf0.1)<-pops[1:7]; mclapply(CANDf0.1, function(x) arrange(x, Chr, Beg.Win))-> CANDf0.1


Store(CANDf0.5); Store(CANDf0.4); Store(CANDf0.3); Store(CANDf0.2); Store(CANDf0.1);Store(list.SCAN) #now list.SCAN has everything I need.



#outlier windows


#nrow(list.SCAN[[1]]) #1525424 windows * 0.0005 == iapprox. 763

top763f0.5<-lapply(list.SCAN, function(x) arrange(x, Z.f0.5.P.val)[1:763,]) #top 763 windows ranked by feq=0.5
top763f0.4<-lapply(list.SCAN, function(x) arrange(x, Z.f0.4.P.val)[1:763,]) #top 763 windows ranked by feq=0.4
top763f0.3<-lapply(list.SCAN, function(x) arrange(x, Z.f0.3.P.val)[1:763,]) #top 763 windows ranked by feq=0.3
top763f0.2<-lapply(list.SCAN, function(x) arrange(x, Z.f0.2.P.val)[1:763,]) #top 763 windows ranked by feq=0.2
top763f0.1<-lapply(list.SCAN, function(x) arrange(x, Z.f0.1.P.val)[1:763,]) #top 763 windows ranked by feq=0.1

gc()

Store(top763f0.5); Store(top763f0.4); Store(top763.f0.3); Store(top763f0.2); Store(top763f0.1)


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

#first, sort

lapply(top763f0.5, function(x) arrange(x, Chr, Beg.Win, End.Win))-> top763f0.5

for(i in 1:7){

with(top763f0.5[[i]], paste0("chr", Chr))-> top763f0.5[[i]]$Chr}


for(i in 1:7){
arrange(top763f0.5[[i]], Beg.Win, End.Win)-> top763f0.5[[i]]}
mclapply(1:7, function(x) bedTools.merge(bed1=top763f0.5[[x]]))-> merge.top816f0.5
