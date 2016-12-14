########################################################################################################
#	Barbara Bitarello
#	Last modified: 14.12.2016
#	Does: assigns tf values and gives the p-values for the most extreme window overlapping a gene
########################################################################################################


funcA<-function(x){
filter(hg19.coding.coords.bed, name==x)$chr-> chr
as.numeric(gsub(".P.val", "", gsub("Z.f","",assign.tf(find.gene(name1=x)[[1]])$assigned.tf.per.gene[[1]])))-> tf.LWK
as.numeric(gsub(".P.val", "", gsub("Z.f","",assign.tf(find.gene(name1=x)[[1]])$assigned.p.gene)))-> p.LWK
as.numeric(gsub(".P.val", "", gsub("Z.f","",assign.tf(find.gene(name1=x, df=YRI.win)[[1]])$assigned.tf.per.gene[[1]])))-> tf.YRI
as.numeric(gsub(".P.val", "", gsub("Z.f","",assign.tf(find.gene(name1=x, df=YRI.win)[[1]])$assigned.p.gene)))-> p.YRI
as.numeric(gsub(".P.val", "", gsub("Z.f","",assign.tf(find.gene(name1=x, df=GBR.win)[[1]])$assigned.tf.per.gene[[1]])))-> tf.GBR
as.numeric(gsub(".P.val", "", gsub("Z.f","",assign.tf(find.gene(name1=x, df=GBR.win)[[1]])$assigned.p.gene)))-> p.GBR
as.numeric(gsub(".P.val", "", gsub("Z.f","",assign.tf(find.gene(name1=x, df=TSI.win)[[1]])$assigned.tf.per.gene[[1]])))-> tf.TSI
as.numeric(gsub(".P.val", "", gsub("Z.f","",assign.tf(find.gene(name1=x, df=TSI.win)[[1]])$assigned.p.gene)))-> p.TSI
res<-data.frame(Chr=chr, Acronym=i, tf.LWK=as.numeric(tf.LWK), p.LWK=p.LWK, tf.YRI=tf.YRI, p.YRI=p.YRI, tf.GBR=tf.GBR, p.GBR=p.GBR, tf.TSI=tf.TSI, p.TSI=p.TSI,stringsAsFactors = FALSE)
return(res)
cat(x, ' done\n')
}

