##############################################
#	R function to plot SFS from VCF data
#
#	Barbara Bitarello
#	Last modified: 14.12. 2016
###############################################


SFS.function<-function(CHR,BEG, END, POP){#POP: a number form 1-7, where: 
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
	#if ANC==ALT, do 1-counts 
	toupper(res$Anc)-> res$Anc
	toupper(res$ALT)-> res$ALT
	subset(res, Anc %in% c("A", "C", "G", "T"))->res2
	temp<-which(res2$Anc==res2$ALT)
	100-res2[temp,n]-> res2[temp,n]
return(res2[,n])
}
#

SFS.function2<-function(CHR,BEG, END, POP){#POP: a number form 1-7, where: 
header<-c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'Anc', 'AWS', 'LWK', 'YRI', 'CEU','FIN', 'GBR', 'TSI', 'CHB', 'CHS', 'JPT', 'MXL', 'CLM ','PUR')
n<-POP+6 #index the VCF columns, because in this parituclar file pop counts start in col 7.
PATH.FILE<-paste0('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/input_data/chr', CHR, '/', 'AC_13pops_chr',CHR,'.hg19.pantro2.Map50_100.TRF.SDs.tsv.gz')
a.file=tempfile(pattern = "", fileext = ".txt")
out=tempfile(pattern = "", fileext = ".txt")
options(scipen =99) # not to use scientific notation when writing out
gsub("chr", "", CHR)-> chr
command<-paste0('tabix ', PATH.FILE, ' ', CHR, ':', BEG, '-', END, " >", out)

  cat(command,"\n")
  try(system(command))

  res=read.table(out,header=F)
  colnames(res)<-header
       toupper(res$Anc)-> res$Anc
        toupper(res$ALT)-> res$ALT
        subset(res, Anc %in% c("A", "C", "G", "T"))->res2
        temp<-which(res2$Anc==res2$ALT)
        100-res2[temp,n]-> res2[temp,n]
return(res2[,n])
}
#
