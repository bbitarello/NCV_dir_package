##########################################
#	Bárbara Bitarello
#
#	Last modified: 23.01/2017
###########################################
#Goals
#1. take outlier windows (merged version)
#2. take from the VCF files, regarding SNPs in those windows: snpID, chr, start(pos), end(pos), REF, ALT
#3. remove duplicates if any
#4. keep only chr, start, end, ref, alt
#5. run annovar for this file

####
#load packages
library(parallel)
library(dplyr)
library(SOAR)
library(SOAR)
Sys.setenv(R_LOCAL_CACHE="estsession")
library(data.table)


###
#create a function to get the SNPs

SNP_catch<-function(CHR, BEG, END, POP){#POP, a number from 1-7
header<-c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'Anc', 'AWS', 'LWK', 'YRI', 'CEU','FIN', 'GBR', 'TSI', 'CHB', 'CHS', 'JPT', 'MXL', 'CLM ','PUR')
n<-POP+6 #index the VCF columns, because in this parituclar file pop counts start in col 7.
PATH.FILE<-paste0('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/input_data/', CHR, '/', 'AC_13pops_',CHR,'.hg19.pantro2.Map50_100.TRF.SDs.tsv.gz')
a.file=tempfile(pattern = "", fileext = ".txt")
out=tempfile(pattern = "", fileext = ".txt")
options(scipen =99) # not to use scientific notation when writing out
gsub("chr", "", CHR)-> chr
command<-paste0('tabix ', PATH.FILE, ' ', chr, ':', BEG, '-', END, " >", out)

  cat(command,"\n")
  try(system(command))

  res=read.table(out,header=F)
  colnames(res)<-header
  which(res[,n]>0 & res[,n]<100) -> tmp

  res2<- res[tmp,c(1,2,3,4,5,6,n)]
return(res2)
}



##run the function

#e.g

mclapply(1:nrow(merge.top829f0.5[[2]]),function(x) SNP_catch(CHR=merge.top829f0.5[[2]][x,]$V1, BEG=merge.top829f0.5[[2]][x,]$V2, END=merge.top829f0.5[[2]][x,]$V3, POP=2))-> test

mclapply(1:nrow(merge.top829f0.4[[2]]),function(x) SNP_catch(CHR=merge.top829f0.4[[2]][x,]$V1, BEG=merge.top829f0.4[[2]][x,]$V2, END=merge.top829f0.4[[2]][x,]$V3, POP=2))-> test2

mclapply(1:nrow(merge.top829f0.3[[2]]),function(x) SNP_catch(CHR=merge.top829f0.3[[2]][x,]$V1, BEG=merge.top829f0.3[[2]][x,]$V2, END=merge.top829f0.3[[2]][x,]$V3, POP=2))-> test3

setDT(rbind(do.call(rbind, test), do.call(rbind, test2), do.call(rbind, test3)))-> test4

setkey(test4, ID)
unique(test4)-> test5

colnames(test5)[2]<-'BEG'

test6<-test5 %>% mutate(END=BEG) %>% select(CHROM, BEG, END, REF, ALT,ID)


write.table(test6, quote=F, sep="\t", col.names=F,row.names=F, file="bedfiles/annovar_outlier_LWK_input.txt")




mclapply(1:nrow(merge.top829f0.5[[6]]),function(x) SNP_catch(CHR=merge.top829f0.5[[6]][x,]$V1, BEG=merge.top829f0.5[[6]][x,]$V2, END=merge.top829f0.5[[6]][x,]$V3, POP=6))-> GBR.1

mclapply(1:nrow(merge.top829f0.4[[6]]),function(x) SNP_catch(CHR=merge.top829f0.4[[6]][x,]$V1, BEG=merge.top829f0.4[[6]][x,]$V2, END=merge.top829f0.4[[6]][x,]$V3, POP=6))-> GBR.2

mclapply(1:nrow(merge.top829f0.3[[6]]),function(x) SNP_catch(CHR=merge.top829f0.3[[6]][x,]$V1, BEG=merge.top829f0.3[[6]][x,]$V2, END=merge.top829f0.3[[6]][x,]$V3, POP=6))-> GBR.3

setDT(rbind(do.call(rbind, GBR.1), do.call(rbind, GBR.2), do.call(rbind, GBR.3)))-> GBR.4

setkey(GBR.4, ID)
unique(GBR.4)-> GBR.5

colnames(GBR.5)[2]<-'BEG'

GBR6<-GBR.5 %>% mutate(END=BEG) %>% select(CHROM, BEG, END, REF, ALT,ID)

write.table(GBR6, quote=F, sep="\t", col.names=F,row.names=F, file="bedfiles/annovar_outlier_GBR_input.txt")






mclapply(1:nrow(merge.top829f0.5[[3]]),function(x) SNP_catch(CHR=merge.top829f0.5[[3]][x,]$V1, BEG=merge.top829f0.5[[3]][x,]$V2, END=merge.top829f0.5[[3]][x,]$V3, POP=3))-> YRI.1

mclapply(1:nrow(merge.top829f0.4[[3]]),function(x) SNP_catch(CHR=merge.top829f0.4[[3]][x,]$V1, BEG=merge.top829f0.4[[3]][x,]$V2, END=merge.top829f0.4[[3]][x,]$V3, POP=3))-> YRI.2

mclapply(1:nrow(merge.top829f0.3[[3]]),function(x) SNP_catch(CHR=merge.top829f0.3[[3]][x,]$V1, BEG=merge.top829f0.3[[3]][x,]$V2, END=merge.top829f0.3[[3]][x,]$V3, POP=3))-> YRI.3

setDT(rbind(do.call(rbind, YRI.1), do.call(rbind, YRI.2), do.call(rbind,YRI.3)))-> YRI.4

setkey(YRI.4, ID)
unique(YRI.4)-> YRI.5

colnames(YRI.5)[2]<-'BEG'

YRI6<-YRI.5 %>% mutate(END=BEG) %>% select(CHROM, BEG, END, REF, ALT,ID)

write.table(YRI6, quote=F, sep="\t", col.names=F,row.names=F, file="bedfiles/annovar_outlier_YRI_input.txt")

#

mclapply(1:nrow(merge.top829f0.5[[7]]),function(x) SNP_catch(CHR=merge.top829f0.5[[7]][x,]$V1, BEG=merge.top829f0.5[[7]][x,]$V2, END=merge.top829f0.5[[7]][x,]$V3, POP=7))-> TSI.1

mclapply(1:nrow(merge.top829f0.4[[7]]),function(x) SNP_catch(CHR=merge.top829f0.4[[7]][x,]$V1, BEG=merge.top829f0.4[[7]][x,]$V2, END=merge.top829f0.4[[7]][x,]$V3, POP=7))-> TSI.2

mclapply(1:nrow(merge.top829f0.3[[7]]),function(x) SNP_catch(CHR=merge.top829f0.3[[7]][x,]$V1, BEG=merge.top829f0.3[[7]][x,]$V2, END=merge.top829f0.3[[7]][x,]$V3, POP=7))-> TSI.3

setDT(rbind(do.call(rbind, TSI.1), do.call(rbind, TSI.2), do.call(rbind,TSI.3)))-> TSI.4

setkey(TSI.4, ID)
unique(TSI.4)-> TSI.5

colnames(TSI.5)[2]<-'BEG'

TSI6<-TSI.5 %>% mutate(END=BEG) %>% select(CHROM, BEG, END, REF, ALT,ID)

write.table(TSI6, quote=F, sep="\t", col.names=F,row.names=F, file="bedfiles/annovar_outlier_TSI_input.txt")

#net, run annovar in the scratch folder
#e.g. command

# table_annovar.pl /mnt/sequencedb/PopGen/barbara/NCV_dir_package/read_scan_data/bedfiles/annovar_outlier_LWK_input.txt humandb/ -buildver hg19 -out myanno_LWK_outliers -remove -protocol wgEncodeGencodeBasicV19,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,1000g2015aug_eur,exac03,avsnp147,dbnsfp30a -operation g,r,r,f,f,f,f,f,f -nastring . -csvout


#read in annovar output

setDT(read.csv('/mnt/scratch/barbara/annovar/myanno_TSI_outliers.hg19_multianno.csv'))-> TSI_annovar
setDT(read.csv('/mnt/scratch/barbara/annovar/myanno_YRI_outliers.hg19_multianno.csv'))-> YRI_annovar
setDT(read.csv('/mnt/scratch/barbara/annovar/myanno_GBR_outliers.hg19_multianno.csv'))-> GBR_annovar
setDT(read.csv('/mnt/scratch/barbara/annovar/myanno_LWK_outliers.hg19_multianno.csv'))-> LWK_annovar



