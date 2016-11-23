##############################################
#	R function to plot SFS from VCF data
#
#	Barbara Bitarello
#	Last modified: 
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
  #write bed formatted dataframes to tempfile
#  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
# command=paste(functionstring, "-i", a.file, "-nms", opt.string,">",out,sep=" ")
#	command=paste(functionstring, "-i", a.file, opt.string,">",out,sep=" ")

  cat(command,"\n")
  try(system(command))

  res=read.table(out,header=F)
  colnames(res)<-header
  #keep only Derived allele frequencies
 res[which(toupper(res$Anc) != toupper(res$ALT)),]-> res2

 # unlink(a.file)
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
  #write bed formatted dataframes to tempfile
#  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
# command=paste(functionstring, "-i", a.file, "-nms", opt.string,">",out,sep=" ")
#	command=paste(functionstring, "-i", a.file, opt.string,">",out,sep=" ")

  cat(command,"\n")
  try(system(command))

  res=read.table(out,header=F)
  colnames(res)<-header
  #keep only Derived allele frequencies
 res[which(toupper(res$Anc) != toupper(res$ALT)),]-> res2

 # unlink(a.file)
return(res2[,n])
}
#
