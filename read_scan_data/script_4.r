##############################################
#	handling extreme windows
#	Barbara Bitarello
#	Creation: 12.10.2016
#	Last modified: 12.10.2016
#
##############################################

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
