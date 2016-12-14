##################################
#	Make tables for paper
#
#	Barbara Bitarello
#	Last modified: 13.12.2016	
##################################

library(SOAR)
Sys.setenv(R_LOCAL_CACHE="estsession")
source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/mclapply2.R')
source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/assign_tf_gene.R')

as.character(read.table('bedfiles/African.Genes.bed')$V1)-> Afr
mclapply2(Afr, function(y) funcA(y))-> res.afr
setDT(do.call('rbind',res.afr))-> res.afr2
write.table(res.afr2, file='bedfiles/Table_afr_manuscript.txt', row.names=F, quote=F)
Store(res.afr2)
#
as.character(read.table('bedfiles/European.Genes.bed')$V1)-> Eur
mclapply2(Eur, function(y) funcA(y))-> res.eur
setDT(do.call('rbind',res.eur))-> res.eur2
write.table(res.eur2, file='bedfiles/Table_eur_manuscript.txt')
Store(res.eur2)


#significant shared genes
as.character(read.table('bedfiles/cand.African.Genes.bed')$V1)-> Afr.cand
mclapply2(Afr.cand, function(y) funcA(y))-> res.afr.cand
setDT(do.call('rbind',res.afr.cand))-> res.afr.cand2
write.table(res.afr.cand2, file='bedfiles/Table_afr_cand_manuscript.txt')
Store(res.afr.cand2)
#
as.character(read.table('bedfiles/cand.European.Genes.bed')$V1)-> Eur.cand
mclapply2(Eur.cand, function(y) funcA(y))-> res.eur.cand
setDT(do.call('rbind',res.eur.cand))-> res.eur.cand2
write.table(res.eur.cand2, file='bedfiles/Table_eur_cand_manuscript.txt')
Store(res.eur.cand2)
