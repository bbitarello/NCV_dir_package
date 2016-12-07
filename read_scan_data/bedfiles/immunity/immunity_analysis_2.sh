###########################
# Immunity analysis
##########################
 
cd /mnt/sequencedb/PopGen/barbara/NCV_dir_package/read_scan_data/bedfiles/immunity
 
#make lists of genes per population

cat ../prot.cod.CANDf0.5.YRI.bed ../prot.cod.CANDf0.4.YRI.bed ../prot.cod.CANDf0.3.YRI.bed |sort|uniq > YRI.cand.all.bed
cat ../prot.cod.CANDf0.5.LWK.bed ../prot.cod.CANDf0.4.LWK.bed ../prot.cod.CANDf0.3.LWK.bed |sort|uniq > LWK.cand.all.bed
cat ../prot.cod.CANDf0.5.GBR.bed ../prot.cod.CANDf0.4.GBR.bed ../prot.cod.CANDf0.3.GBR.bed |sort|uniq > GBR.cand.all.bed
cat ../prot.cod.CANDf0.5.TSI.bed ../prot.cod.CANDf0.4.TSI.bed ../prot.cod.CANDf0.3.TSI.bed |sort|uniq > TSI.cand.all.bed
 

cat YRI.cand.all.bed TSI.cand.all.bed GBR.cand.all.bed LWK.cand.all.bed |sort|uniq > AllPops.cand.all.bed #2348 genes...


cat ../prot.cod.top829f0.5.YRI.bed ../prot.cod.top829f0.4.YRI.bed ../prot.cod.top829f0.3.YRI.bed |sort|uniq > YRI.outliers.all.bed
cat ../prot.cod.top829f0.5.LWK.bed ../prot.cod.top829f0.4.LWK.bed ../prot.cod.top829f0.3.LWK.bed |sort|uniq > LWK.outliers.all.bed
cat ../prot.cod.top829f0.5.GBR.bed ../prot.cod.top829f0.4.GBR.bed ../prot.cod.top829f0.3.GBR.bed |sort|uniq > GBR.outliers.all.bed
cat ../prot.cod.top829f0.5.TSI.bed ../prot.cod.top829f0.4.TSI.bed ../prot.cod.top829f0.3.TSI.bed |sort|uniq > TSI.outliers.all.bed
 

cat YRI.outliers.all.bed TSI.outliers.all.bed GBR.outliers.all.bed LWK.outliers.all.bed |sort|uniq > AllPops.outliers.all.bed #402 genes


 
#go to biomart ensembl hg19 and obtain GO terms for each gene
 
#download the file
 
GO_terms_all_outliers.txt
 

awk '$5!=""{print $1}' GO_terms_all_outliers.txt |sort|uniq -c|wc  #number of outlier genes with at least one GO term : 378

grep -f <(sed 's/ /_/g' 386_KeywordsSearchNameDescGO.txt) <(sed 's/ /_/g' GO_terms_all_outliers.txt)|awk '{print $1}'|sort|uniq > immune_outlier_genes.txt  #140 genes

#now the background genes
#scanned genes

bedtools intersect -a <(sed 's/^/chr/' <(sort -k 1,1 -k2,2n  ../background_windows.bed)) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > all.scanned.genes.bed

awk '$5!=""{print $1}' GO_term_all_scanned_genes.txt |sort|uniq -c|wc #17074 genes with at least one GO term

grep -f <(sed 's/ /_/g' 386_KeywordsSearchNameDescGO.txt) <(sed 's/ /_/g' GO_term_all_scanned_genes.txt)|awk '{print $1}'|sort|uniq > immune_all_genes.txt



awk '$5!=""{print $1}' GO_term_all_signif_genes.txt |sort|uniq -c|wc #2215 genes with at least one GO term

grep -f <(sed 's/ /_/g' 386_KeywordsSearchNameDescGO.txt) <(sed 's/ /_/g' GO_term_all_signif_genes.txt)|awk '{print $1}'|sort|uniq > immune_signif_genes.txt  #733
