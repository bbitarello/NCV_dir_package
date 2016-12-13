########################################
#	Barbara Bitarello
#
#	Last modified: 07.12.2016
#	Assign tf to candidate windows
########################################


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
#

#find.gene(df, chr=6, name='HLA-B')
find.gene<-function(df=LWK.win,name1="HLA-B"){  #df can be changes for diff pops so that we can check p-value in each population!
as.numeric(as.character(dplyr::filter(hg19.coding.coords.bed, name %in% name1)$chr))->chr
which(unlist(mclapply(names.all.coding[[chr]], function(x) strsplit(x, ':', fixed=TRUE)[[1]][[1]]))==name1)->QUERY.POS
df[[chr]][[QUERY.POS]]-> QUERY.SUBSET
nrow(QUERY.SUBSET)->n1
return(list(query_subset=QUERY.SUBSET, query_pos=QUERY.POS, GENE=name1, number_windows=n1))
}  #currently this gives me correct results for all.coding, but not for the ALL.POPS.AF datasets. I am trying to fix this.
#
#
read.table('/mnt/sequencedb/PopGen/cesare/hg19/bedfiles/ensembl_genes_hg19.bed.gz')->hg19.coding.coords.bed
names(hg19.coding.coords.bed)<-c('chr', 'beg', 'end','name', 'type')



#
assign.tf<-function(df){
nrow(df)-> n
assigned.tf.per.window<-sapply(1:n, function(x) names(which.min(select(df[x,], Z.f0.5.P.val:Z.f0.3.P.val))))
min.p.value.per.window<-sapply(1:n, function(x) min(select(df[x,], Z.f0.5.P.val:Z.f0.3.P.val)))
cbind(assigned.ft=assigned.tf.per.window, min.p.value=min.p.value.per.window)-> res1
as.numeric(res1[,2])->res1[,2]
res1[which.min(res1[,2]),1]-> assigned.tf.gene
as.numeric(res1[which.min(res1[,2]),2])-> assigned.p.val.gene
return(list(assigned.per.window=res1,assigned.tf.per.gene= assigned.tf.gene, assigned.p.gene=assigned.p.val.gene))
}
####
####
Objects()

table(unlist(mclapply(1:nrow(Union.CANDf0.5_0.4_0.3[[2]]), function(x) as.numeric(gsub(".P.val", "",gsub("Z.f","",assign.tf(Union.CANDf0.5_0.4_0.3[[2]][x,])[[1]][[1]]))))))/nrow(Union.CANDf0.5_0.4_0.3[[2]])
#
#     0.3       0.4       0.5 
#0.5355097 0.1223061 0.3421842 

table(unlist(mclapply(1:nrow(Union.CANDf0.5_0.4_0.3[[3]]), function(x) as.numeric(gsub(".P.val", "",gsub("Z.f","",assign.tf(Union.CANDf0.5_0.4_0.3[[3]][x,])[[1]][[1]]))))))/nrow(Union.CANDf0.5_0.4_0.3[[3]])
#
#      0.3       0.4       0.5 
#0.4996324 0.1267724 0.3735952 

table(unlist(mclapply(1:nrow(Union.CANDf0.5_0.4_0.3[[6]]), function(x) as.numeric(gsub(".P.val", "",gsub("Z.f","",assign.tf(Union.CANDf0.5_0.4_0.3[[6]][x,])[[1]][[1]]))))))/nrow(Union.CANDf0.5_0.4_0.3[[6]])
#
#      0.3       0.4       0.5 
#0.4996324 0.1267724 0.3735952 

table(unlist(mclapply(1:nrow(Union.CANDf0.5_0.4_0.3[[7]]), function(x) as.numeric(gsub(".P.val", "",gsub("Z.f","",assign.tf(Union.CANDf0.5_0.4_0.3[[7]][x,])[[1]][[1]]))))))/nrow(Union.CANDf0.5_0.4_0.3[[7]])
#
#      0.3       0.4       0.5 
#0.4859944 0.1421030 0.3719026 



#but this doesn't make so much sense...i should do that only for windows in the intersection of tf.

table(unlist(mclapply(1:3715, function(x) as.numeric(gsub(".P.val","",gsub("Z.f","",assign.tf(subset(Union.CANDf0.5_0.4_0.3[[2]], Win.ID %in% intersect(intersect(CANDf0.5[[2]]$Win.ID, CANDf0.4[[2]]$Win.ID), CANDf0.3[[2]]$Win.ID))[x,])[[1]][[1]]))))))/3715


#
#      0.3       0.4       0.5 
#0.5138627 0.2080754 0.2780619 


#YRI
table(unlist(mclapply(1:4156, function(x) as.numeric(gsub(".P.val","",gsub("Z.f","",assign.tf(subset(Union.CANDf0.5_0.4_0.3[[3]], Win.ID %in% intersect(intersect(CANDf0.5[[3]]$Win.ID, CANDf0.4[[3]]$Win.ID), CANDf0.3[[3]]$Win.ID))[x,])[[1]][[1]]))))))/4156


#      0.3       0.4       0.5 
#0.5120308 0.2259384 0.2620308 


#GBR
table(unlist(mclapply(1:3426, function(x) as.numeric(gsub(".P.val","",gsub("Z.f","",assign.tf(subset(Union.CANDf0.5_0.4_0.3[[6]], Win.ID %in% intersect(intersect(CANDf0.5[[6]]$Win.ID, CANDf0.4[[6]]$Win.ID), CANDf0.3[[6]]$Win.ID))[x,])[[1]][[1]]))))))/3426

#     0.3       0.4       0.5 
#0.5683012 0.2180385 0.2136602 


table(unlist(mclapply(1:3405, function(x) as.numeric(gsub(".P.val","",gsub("Z.f","",assign.tf(subset(Union.CANDf0.5_0.4_0.3[[7]], Win.ID %in% intersect(intersect(CANDf0.5[[7]]$Win.ID, CANDf0.4[[7]]$Win.ID), CANDf0.3[[7]]$Win.ID))[x,])[[1]][[1]]))))))/3405

#      0.3       0.4       0.5 
#0.5406755 0.2414097 0.2179148 




#outliers

table(unlist(mclapply(1:421, function(x) as.numeric(gsub(".P.val","",gsub("Z.f","",assign.tf(subset(Union.top0.5_0.4_0.3[[2]], Win.ID %in% intersect(intersect(top829f0.5[[2]]$Win.ID, top829f0.4[[2]]$Win.ID), top829f0.3[[2]]$Win.ID))[x,])[[1]][[1]]))))))/421



#      0.3       0.4       0.5 
#0.5154394 0.2779097 0.2066508    #0.5 and 0.4 kind of inverted compared to signifcant windows.




#a different strategy:

test.col<-vector('list', 7)


for(i in 1:7){ #add the rest later

nrow(Union.top0.5_0.4_0.3[[i]])-> n
test.col[[i]]<-vector('list', n)

for (y in 1:n){
temp.df<-select(Union.top0.5_0.4_0.3[[i]][y,], Z.f0.5.P.val:Z.f0.3.P.val)

    colnames(temp.df)[which(temp.df == min(temp.df))]-> test.col[[i]][[y]]
    which(unlist(mclapply(test.col[[i]], function(x) length(x)==2)))->repl2
    
    if(length(repl2)>=1){
     for(j in 1: length(repl2)){
        paste0(test.col[[i]][[repl2[j]]][1], "|", test.col[[i]][[repl2[j]]][2])-> test.col[[i]][[repl2[j]]]}}
    which(unlist(lapply(test.col[[i]], function(x) length(x)==3)))-> repl3
    if(length(repl3)>=1){
        for(j in 1: length(repl3)){
            paste0(test.col[[i]][[repl3[j]]][1], "|", test.col[[i]][[repl3[j]]][2], "|", test.col[[i]][[repl3[j]]][3])-> test.col[[i]][[repl3[j]]]}}
}

cat(pops[i], ' done\n')

}


for ( i in 1:7){

setDT(cbind(Union.top0.5_0.4_0.3[[i]], AssignedTF=gsub(".P.val", "",gsub("Z.f","", unlist(test.col[[i]])))))-> Union.top0.5_0.4_0.3[[i]]
}

Store(Union.top0.5_0.4_0.3)

remove(test.col)
#significant


test.col2<-vector('list', 7)

for(i in 1:7){ #add the rest later

nrow(Union.CANDf0.5_0.4_0.3[[i]])-> n
test.col2[[i]]<-vector('list', n)

for (y in 1:n){
temp.df<-select(Union.CANDf0.5_0.4_0.3[[i]][y,], Z.f0.5.P.val:Z.f0.3.P.val)
    colnames(temp.df)[which(temp.df == min(temp.df))]-> test.col2[[i]][[y]]
    which(unlist(mclapply(test.col2[[i]], function(x) length(x)==2)))->repl2

    if(length(repl2)>=1){
     for(j in 1: length(repl2)){
        paste0(test.col2[[i]][[repl2[j]]][1], "|", test.col2[[i]][[repl2[j]]][2])-> test.col2[[i]][[repl2[j]]]}}
    which(unlist(lapply(test.col2[[i]], function(x) length(x)==3)))-> repl3
    if(length(repl3)>=1){
        for(j in 1: length(repl3)){
            paste0(test.col2[[i]][[repl3[j]]][1], "|", test.col2[[i]][[repl3[j]]][2], "|", test.col2[[i]][[repl3[j]]][3])-> test.col2[[i]][[repl3[j]]]}}
}

cat(pops[i], ' done\n')
}


for ( i in 1:7){

setDT(cbind(Union.CANDf0.5_0.4_0.3[[i]], AssignedTF=gsub(".P.val", "",gsub("Z.f","", unlist(test.col2[[i]])))))-> Union.CANDf0.5_0.4_0.3[[i]]
}

names(Union.top0.5_0.4_0.3)<- pops[1:7]
names(Union.CANDf0.5_0.4_0.3)<- pops[1:7]


Store(Union.CANDf0.5_0.4_0.3)
Store(Union.top0.5_0.4_0.3)


#select windows that are signiiicant for all tf per pop

intersect.all.tf<-vector('list', 7)
names(intersect.all.tf)<- pops[1:7]

for ( i in 1:7){
intersect(intersect(CANDf0.5[[i]]$Win.ID, CANDf0.4[[i]]$Win.ID), CANDf0.3[[i]]$Win.ID)-> intersect.all.tf[[i]][[1]]
intersect(intersect(top829f0.5[[i]]$Win.ID, top829f0.4[[i]]$Win.ID), top829f0.3[[i]]$Win.ID)-> intersect.all.tf[[i]][[2]]
}

sapply(1:7, function(x) table(select(filter(Union.top0.5_0.4_0.3[[x]], Win.ID %in% intersect.all.tf[[x]][[2]]), AssignedTF)))
sapply(1:7, function(x) table(select(filter(Union.CANDf0.5_0.4_0.3[[x]], Win.ID %in% intersect.all.tf[[x]][[1]]), AssignedTF)))




#now, to the assigned tf vlaue for the outlier shared and significant shared genes.


read.table('bedfiles/extreme.genes.bed', header=F)[,2]-> extreme.outliers
read.table('bedfiles/cand.extreme.genes.bed', header=F)-> extreme.cand


#first, get windoes that overlap them..

system.time(mclapply(extreme.outliers, function(x) find.gene(name1=x))-> df.outliers)
system.time(mclapply(extreme.cand, function(x) find.gene(name1=x))-> df.signif)
#next, assign. tf per pop

#summarise results

