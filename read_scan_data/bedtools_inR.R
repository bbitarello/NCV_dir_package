#bedtools in R
#this is an edit from https://github.com/dariober/bioinformatics-cafe/blob/master/bedtools.R


bedTools.2in<-function(functionstring="bedtools intersect -wo",bed1,bed2,opt.string="")

{
  #create temp files
  a.file=tempfile(pattern = "", fileext = ".bed")
  b.file=tempfile(pattern = "", fileext = ".bed")
  out   =tempfile(pattern = "", fileext = ".bed")
  options(scipen =99) # not to use scientific notation when writing out
  #write bed formatted dataframes to tempfile
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
  # create the command string and call the command using system()
  command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  res=read.table(out,header=F)
  unlink(a.file);unlink(b.file);unlink(out);
  return(res)
}
 


bedTools.merge<-function(functionstring="mergeBed",bed1, opt.string="-c 4 -o collapse")
{

a.file=tempfile(pattern = "", fileext = ".bed")
 out   =tempfile(pattern = "", fileext = ".bed")
  options(scipen =99) # not to use scientific notation when writing out

  #write bed formatted dataframes to tempfile
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
# command=paste(functionstring, "-i", a.file, "-nms", opt.string,">",out,sep=" ")
	command=paste(functionstring, "-i", a.file, opt.string,">",out,sep=" ")

  cat(command,"\n")
  try(system(command))

  res=read.table(out,header=F)
  unlink(a.file)
return(res)
}

##test
