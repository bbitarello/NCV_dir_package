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
 
 
grep -f <(sed 's/ /_/g' 386_KeywordsSearchNameDescGO.txt) <(sed 's/ /_/g' GO_terms_all_outliers.txt)|awk '{print $5}'  #wrong


#try this


grep -if <(grep -if <(sed 's/ /_/g' 386_KeywordsSearchNameDescGO.txt) <(sed  '/^$/d' <(awk '{print $6}' <(sed 's/ /_/g' GO_terms_all_outliers.txt)))) <(sed 's/ /_/g' GO_terms_all_outliers.txt)|awk '{print $1}'|sort|uniq|wc #126/378

grep -if <(grep -if <(sed 's/ /_/g' 386_KeywordsSearchNameDescGO.txt) <(sed  '/^$/d' <(awk '{print $6}' <(sed 's/ /_/g' GO_term_all_signif_genes.txt)))) <(sed 's/ /_/g' GO_term_all_signif_genes.txt)|awk '{print $1}'|sort|uniq|wc #694


grep -if <(grep -if <(sed 's/ /_/g' 386_KeywordsSearchNameDescGO.txt)  <(sed  '/^$/d' <(awk '{print $6}' <(sed 's/ /_/g' GO_term_all_scanned_genes.txt)))) <(sed 's/ /_/g' GO_term_all_scanned_genes.txt)|awk '{print $1}'|sort|uniq|wc #



#genes with at least one GO


awk '$5!=""{print $1}' GO_terms_all_outliers.txt |sort|uniq|wc #378 #126/378=33.3%
awk '$5!=""{print $1}' GO_term_all_signif_genes.txt |sort|uniq|wc #2215 #694/2215=31.3%
awk '$5!=""{print $1}' GO_term_all_scanned_genes.txt|sort|uniq|wc #17074

#now another analyses: take the significant GO categories from João's spreadsheet and check how many of them are immune-related.

grep -if <(sed 's/ /_/g' 386_KeywordsSearchNameDescGO.txt) <(sed 's/ /_/g' temp_GO_signif.txt|sort|uniq)

grep -if <(sed 's/ /_/g' 386_KeywordsSearchNameDescGO.txt) <(sed 's/ /_/g' temp_GO_outliers.txt|sort|uniq)
