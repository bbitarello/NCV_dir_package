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


for(i in 1:7){
m<-median(Union.CANDf0.5_0.4_0.3[[i]]$Nr.IS)
ggplot(Union.CANDf0.5_0.4_0.3[[i]], aes(x=Nr.IS)) + geom_density(colour='cornflowerblue') + geom_vline(xintercept=m)
ggsave(paste0('figures/Nr.IS.',pops[i],'.cand.pdf'))



m1<-median(list.SCAN[[i]]$Nr.IS)
ggplot(list.SCAN[[i]], aes(x=Nr.IS)) + geom_density(colour='darkgray') + geom_vline(xintercept=m1)
ggsave(paste0('figures/Nr.IS.genomic', pops[i], '.pdf'))



cat('Finished', pops[i], '\n')
}


for(i in 1:7){
m<-median(Union.top0.5_0.4_0.3[[i]]$Nr.IS)
ggplot(Union.top0.5_0.4_0.3[[i]], aes(x=Nr.IS)) + geom_density(colour='cornflowerblue') + geom_vline(xintercept=m)
ggsave(paste0('figures/Nr.IS.',pops[i],'.top829.pdf'))
cat('Finished', pops[i], '\n')
}


#
for(i in 1:7){
m<-median(Union.CANDf0.5_0.4_0.3[[i]]$PtoD)
ggplot(Union.CANDf0.5_0.4_0.3[[i]], aes(x=PtoD)) + geom_density(colour='cornflowerblue') + geom_vline(xintercept=m)
ggsave(paste0('figures/PtoD.',pops[i],'.cand.pdf'))


m1<-median(list.SCAN[[i]]$PtoD)
ggplot(list.SCAN[[i]], aes(x=PtoD)) + geom_density(colour='darkgray') + geom_vline(xintercept=m1)
ggsave(paste0('figures/PtoD.genomic', pops[i], '.pdf'))

cat('Finished', pops[i], '\n')
}


for(i in 1:7){
m<-median(Union.top0.5_0.4_0.3[[i]]$PtoD)
ggplot(Union.top0.5_0.4_0.3[[i]], aes(x=PtoD)) + geom_density(colour='cornflowerblue') + geom_vline(xintercept=m)
ggsave(paste0('figures/PtoD.',pops[i],'.top829.pdf'))

cat('Finished', pops[i], '\n')
}


for(i in 1:7){
m<-median(Union.CANDf0.5_0.4_0.3[[i]]$Nr.FDs)
ggplot(Union.CANDf0.5_0.4_0.3[[i]], aes(x=Nr.FDs)) + geom_density(colour='cornflowerblue') + geom_vline(xintercept=m)
ggsave(paste0('figures/Nr.FDs.',pops[i],'.cand.pdf'))


m1<-median(list.SCAN[[i]]$Nr.FDs)
ggplot(list.SCAN[[i]], aes(x=Nr.FDs)) + geom_density(colour='darkgray') + geom_vline(xintercept=m1)
ggsave(paste0('figures/Nr.FDs.genomic', pops[i], '.pdf'))

cat('Finished', pops[i], '\n')
}


for(i in 1:7){
m<-median(Union.top0.5_0.4_0.3[[i]]$Nr.FDs)
ggplot(Union.top0.5_0.4_0.3[[i]], aes(x=Nr.FDs)) + geom_density(colour='cornflowerblue') + geom_vline(xintercept=m)
ggsave(paste0('figures/Nr.FDs.',pops[i],'.top829.pdf'))

cat('Finished', pops[i], '\n')
}


####
#test
for(i in 1:7){
nrow(list.SCAN[[i]])-> n1;nrow(Union.CANDf0.5_0.4_0.3[[i]])-> n2;nrow(Union.top0.5_0.4_0.3[[i]])-> n3

cbind(rbind(select(list.SCAN[[i]], Nr.IS, Nr.FDs, Nr.SNPs), select(Union.CANDf0.5_0.4_0.3[[i]], Nr.IS, Nr.FDs, Nr.SNPs), select(Union.top0.5_0.4_0.3[[i]], Nr.IS, Nr.FDs, Nr.SNPs)), Type=c(rep('genomic', n1), rep('cand', n2), rep('Top829', n3)))-> temp

tmpname<-paste0('Nr.IS.',pops[[i]], '.pdf')
tmpname1<-paste0('Nr.FDs.',pops[[i]], '.pdf')
tmpname2<-paste0('Nr.SNPs.',pops[[i]], '.pdf')

ggplot(temp) + geom_density(aes(x = Nr.IS, colour = Type))
ggsave(paste0('figures/',tmpname))

ggplot(temp) + geom_density(aes(x = Nr.FDs, colour = Type))
ggsave(paste0('figures/',tmpname1))

ggplot(temp) + geom_density(aes(x = Nr.SNPs, colour = Type))
ggsave(paste0('figures/',tmpname2))


cat('Finished', pops[i], '\n')
}



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


