##############################################################
#	Barbara Bitarello
#	Last modified:6.12.2016
#
#	Calculate some descriptive numbers for the paper.
##############################################################


#number of significant windows that overlap coding elements:



mclapply(1:7, function(i) setDT(rbind(intersect.CANDf0.5[[i]],intersect.CANDf0.4[[i]], intersect.CANDf0.3[[i]])))-> tmp


mclapply(1:7, function(i) setDT(rbind(intersect.top829f0.5[[i]],intersect.top829f0.4[[i]], intersect.top829f0.3[[i]])))-> tmp.1

mclapply(tmp, function(j) setDT(filter(j, V11=='protein_coding')))-> tmp2

mclapply(tmp.1, function(j) setDT(filter(j, V11=='protein_coding')))-> tmp.2

count2<-vector('list', 7)
count.2<-vector('list', 7)


for (x in 1:7){

length(unique(tmp2[[x]]$V4))-> count
length(unique(tmp.2[[x]]$V4))-> count.1

length(unique(unlist(sapply(1:count, function(z) strsplit(as.character(unique(tmp2[[x]]$V4)[z]), ",")))))-> count2[[x]]

length(unique(unlist(sapply(1:count.1, function(z) strsplit(as.character(unique(tmp.2[[x]]$V4)[z]), ",")))))-> count.2[[x]]

}


mean(as.vector(unlist(count2)/sapply(1:7, function(x) length(Union.CANDf0.5_0.4_0.3[[x]]$Win.ID)))[c(2,3,6,7)]) #40%

mean(as.vector(unlist(count.2)/sapply(1:7, function(x) length(Union.top0.5_0.4_0.3[[x]]$Win.ID)))[c(2,3,6,7)]) #40%

#[1] 0.3900329 0.3857228 0.3959262 0.4049743 0.4062533 0.4116164 0.4126266



#all scanned windows


unique(filter(intersect.scanned.windows, V11=='protein_coding')$V4)-> tempo1

system.time(length(unique(unlist(mclapply(1: length(tempo1), function(x) strsplit(as.character(tempo1[x]), ",")))))-> count3)

count3/nrow(list.SCAN[[2]]) #% of windows (scanned) that overlap protein coding genes. 77% (why is this so high??

