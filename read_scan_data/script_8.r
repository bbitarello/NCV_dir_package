##################################
#	Make tables for paper
#
#	Barbara Bitarello
##################################

as.character(read.table('bedfiles/African.Genes.bed')$V1)-> Afr

res.afr<-vector('list', length(Afr))

for(i in Afr){

filter(hg19.coding.coords.bed, name==i)$chr-> chr

as.numeric(gsub(".P.val", "", gsub("Z.f","",assign.tf(find.gene(name1=i)[[1]])$assigned.tf.per.gene[[1]])))-> tf.LWK

as.numeric(gsub(".P.val", "", gsub("Z.f","",assign.tf(find.gene(name1=i)[[1]])$assigned.p.gene)))-> p.LWK

as.numeric(gsub(".P.val", "", gsub("Z.f","",assign.tf(find.gene(name1=i, df=YRI.win)[[1]])$assigned.tf.per.gene[[1]])))-> tf.YRI

as.numeric(gsub(".P.val", "", gsub("Z.f","",assign.tf(find.gene(name1=i, df=YRI.win)[[1]])$assigned.p.gene)))-> p.YRI

as.numeric(gsub(".P.val", "", gsub("Z.f","",assign.tf(find.gene(name1=i, df=GBR.win)[[1]])$assigned.tf.per.gene[[1]])))-> tf.GBR

as.numeric(gsub(".P.val", "", gsub("Z.f","",assign.tf(find.gene(name1=i, df=GBR.win)[[1]])$assigned.p.gene)))-> p.GBR


as.numeric(gsub(".P.val", "", gsub("Z.f","",assign.tf(find.gene(name1=i, df=TSI.win)[[1]])$assigned.tf.per.gene[[1]])))-> tf.TSI

as.numeric(gsub(".P.val", "", gsub("Z.f","",assign.tf(find.gene(name1=i, df=TSI.win)[[1]])$assigned.p.gene)))-> p.TSI

res.afr[[i]]<-data.frame(Chr=chr, Acronym=i, tf.LWK=as.numeric(tf.LWK), p.LWK=p.LWK, tf.YRI=tf.YRI, p.YRI=p.YRI, tf.GBR=tf.GBR, p.GBR=p.GBR, td.TSI=tf.TSI, p.TSI=p.TSI,stringsAsFactors = FALSE)
cat(i, 'done\n')
}
do.call('rbind',res.afr)-> res.afr2
write.table(res.afr2, file='bedfiles/Table_afr_manuscript.txt')
#


as.character(read.table('bedfiles/European.Genes.bed')$V1)-> Eur

res.eur<-vector('list', length(Eur))

for(i in Eur){

filter(hg19.coding.coords.bed, name==i)$chr-> chr

as.numeric(gsub(".P.val", "", gsub("Z.f","",assign.tf(find.gene(name1=i)[[1]])$assigned.tf.per.gene[[1]])))-> tf.LWK

as.numeric(gsub(".P.val", "", gsub("Z.f","",assign.tf(find.gene(name1=i)[[1]])$assigned.p.gene)))-> p.LWK

as.numeric(gsub(".P.val", "", gsub("Z.f","",assign.tf(find.gene(name1=i, df=YRI.win)[[1]])$assigned.tf.per.gene[[1]])))-> tf.YRI

as.numeric(gsub(".P.val", "", gsub("Z.f","",assign.tf(find.gene(name1=i, df=YRI.win)[[1]])$assigned.p.gene)))-> p.YRI

as.numeric(gsub(".P.val", "", gsub("Z.f","",assign.tf(find.gene(name1=i, df=GBR.win)[[1]])$assigned.tf.per.gene[[1]])))-> tf.GBR

as.numeric(gsub(".P.val", "", gsub("Z.f","",assign.tf(find.gene(name1=i, df=GBR.win)[[1]])$assigned.p.gene)))-> p.GBR

as.numeric(gsub(".P.val", "", gsub("Z.f","",assign.tf(find.gene(name1=i, df=TSI.win)[[1]])$assigned.tf.per.gene[[1]])))-> tf.TSI

as.numeric(gsub(".P.val", "", gsub("Z.f","",assign.tf(find.gene(name1=i, df=TSI.win)[[1]])$assigned.p.gene)))-> p.TSI

res.eur[[i]]<-data.frame(Chr=chr,Acronym=i, tf.LWK=as.numeric(tf.LWK), p.LWK=p.LWK, tf.YRI=tf.YRI, p.YRI=p.YRI, tf.GBR=tf.GBR, p.GBR=p.GBR, td.TSI=tf.TSI, p.TSI=p.TSI,stringsAsFactors = FALSE)
cat(i, 'done\n')
}

do.call('rbind',res.afr)-> res.eur2
write.table(res.eur2, file='bedfiles/Table_eur_manuscript.txt')

