find.gene<-function(df=LWK.win,name1="HLA-B"){  #df can be changes for diff pops so that we can check p-value in each population!
unique(as.numeric(as.character(dplyr::filter(hg19.coding.coords.bed, name %in% name1)$chr)))->chr
which(unlist(mclapply(names.all.coding[[chr]], function(x) strsplit(x, ':', fixed=TRUE)[[1]][[1]]))==name1)->QUERY.POS
if(length(QUERY.POS)>1){
do.call('rbind', df[[chr]][QUERY.POS])-> QUERY.SUBSET}
else{
df[[chr]][[QUERY.POS]]-> QUERY.SUBSET}
nrow(QUERY.SUBSET)->n1
return(list(query_subset=QUERY.SUBSET, query_pos=QUERY.POS, GENE=name1, number_windows=n1))
}  #currently this gives me correct results for all.coding, but not for the ALL.POPS.AF datasets. I am trying to fix this.
#


read.table('/mnt/sequencedb/PopGen/cesare/hg19/bedfiles/ensembl_genes_hg19.bed.gz')->hg19.coding.coords.bed
names(hg19.coding.coords.bed)<-c('chr', 'beg', 'end','name', 'type')


