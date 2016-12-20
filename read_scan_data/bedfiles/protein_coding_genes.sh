############################################################
#	Barbara Bitarello
#
#	Obtain list of protein-coding genes from bedfiles
###########################################################



bedtools intersect -a <(sed -e 's/^/chr/' Union.top829.0.5_0.4_0.3_LWK.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.Union.top829.0.5_0.4_0.3_LWK.bed

bedtools intersect -a <(sed -e 's/^/chr/' Union.top829.0.5_0.4_0.3_YRI.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.Union.top829.0.5_0.4_0.3_YRI.bed

bedtools intersect -a <(sed -e 's/^/chr/' Union.top829.0.5_0.4_0.3_GBR.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.Union.top829.0.5_0.4_0.3_GBR.bed

bedtools intersect -a <(sed -e 's/^/chr/' Union.top829.0.5_0.4_0.3_TSI.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.Union.top829.0.5_0.4_0.3_TSI.bed

bedtools intersect -a <(sed -e 's/^/chr/' Union.CANDf0.5_0.4_0.3_LWK.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.Union.CAND.0.5_0.4_0.3_LWK.bed

bedtools intersect -a <(sed -e 's/^/chr/' Union.CANDf0.5_0.4_0.3_YRI.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.Union.CAND.0.5_0.4_0.3_YRI.bed

bedtools intersect -a <(sed -e 's/^/chr/' Union.CANDf0.5_0.4_0.3_GBR.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.Union.CAND.0.5_0.4_0.3_GBR.bed

bedtools intersect -a <(sed -e 's/^/chr/' Union.CANDf0.5_0.4_0.3_TSI.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.Union.CAND.0.5_0.4_0.3_TSI.bed

#cand

bedtools intersect -a <(sed -e 's/^/chr/' CANDf0.5.LWK.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.CANDf0.5.LWK.bed

bedtools intersect -a <(sed -e 's/^/chr/' CANDf0.4.LWK.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.CANDf0.4.LWK.bed

bedtools intersect -a <(sed -e 's/^/chr/' CANDf0.3.LWK.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.CANDf0.3.LWK.bed

bedtools intersect -a <(sed -e 's/^/chr/' CANDf0.5.YRI.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.CANDf0.5.YRI.bed

bedtools intersect -a <(sed -e 's/^/chr/' CANDf0.4.YRI.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.CANDf0.4.YRI.bed

bedtools intersect -a <(sed -e 's/^/chr/' CANDf0.3.YRI.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.CANDf0.3.YRI.bed

bedtools intersect -a <(sed -e 's/^/chr/' CANDf0.5.GBR.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.CANDf0.5.GBR.bed

bedtools intersect -a <(sed -e 's/^/chr/' CANDf0.4.GBR.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.CANDf0.4.GBR.bed

bedtools intersect -a <(sed -e 's/^/chr/' CANDf0.3.GBR.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.CANDf0.3.GBR.bed

bedtools intersect -a <(sed -e 's/^/chr/' CANDf0.5.TSI.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.CANDf0.5.TSI.bed

bedtools intersect -a <(sed -e 's/^/chr/' CANDf0.4.TSI.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.CANDf0.4.TSI.bed

bedtools intersect -a <(sed -e 's/^/chr/' CANDf0.3.TSI.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.CANDf0.3.TSI.bed

#top

bedtools intersect -a <(sed -e 's/^/chr/' top829f0.5.LWK.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.top829f0.5.LWK.bed

bedtools intersect -a <(sed -e 's/^/chr/' top829f0.4.LWK.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.top829f0.4.LWK.bed

bedtools intersect -a <(sed -e 's/^/chr/' top829f0.3.LWK.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.top829f0.3.LWK.bed

bedtools intersect -a <(sed -e 's/^/chr/' top829f0.5.YRI.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.top829f0.5.YRI.bed

bedtools intersect -a <(sed -e 's/^/chr/' top829f0.4.YRI.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.top829f0.4.YRI.bed

bedtools intersect -a <(sed -e 's/^/chr/' top829f0.3.YRI.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.top829f0.3.YRI.bed

bedtools intersect -a <(sed -e 's/^/chr/' top829f0.5.GBR.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.top829f0.5.GBR.bed

bedtools intersect -a <(sed -e 's/^/chr/' top829f0.4.GBR.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.top829f0.4.GBR.bed

bedtools intersect -a <(sed -e 's/^/chr/' top829f0.3.GBR.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.top829f0.3.GBR.bed

bedtools intersect -a <(sed -e 's/^/chr/' top829f0.5.TSI.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.top829f0.5.TSI.bed

bedtools intersect -a <(sed -e 's/^/chr/' top829f0.4.TSI.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.top829f0.4.TSI.bed

bedtools intersect -a <(sed -e 's/^/chr/' top829f0.3.TSI.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.top829f0.3.TSI.bed


#restricted set of genes


sort <(cat prot.cod.top829f0.5.LWK.bed prot.cod.top829f0.4.LWK.bed prot.cod.top829f0.3.LWK.bed |sort|uniq) <(cat prot.cod.top829f0.5.YRI.bed prot.cod.top829f0.4.YRI.bed prot.cod.top829f0.3.YRI.bed |sort|uniq) |uniq -d > African.Genes.bed
#OR
grep -if prot.cod.Union.top829.0.5_0.4_0.3_LWK.bed prot.cod.Union.top829.0.5_0.4_0.3_YRI.bed > African.Genes.bed #either version yields the same 181 genes

sort <(cat prot.cod.top829f0.5.GBR.bed prot.cod.top829f0.4.GBR.bed prot.cod.top829f0.3.GBR.bed |sort|uniq) <(cat prot.cod.top829f0.5.TSI.bed prot.cod.top829f0.4.TSI.bed prot.cod.top829f0.3.TSI.bed |sort|uniq) |uniq -d > European.Genes.bed
#OR
grep -F -x -f prot.cod.Union.top829.0.5_0.4_0.3_GBR.bed prot.cod.Union.top829.0.5_0.4_0.3_TSI.bed > European.Genes.bed

grep -F -x -f African.Genes.bed European.Genes.bed  > afrANDeur.genes.bed #102

grep -F -x -v -f European.Genes.bed African.Genes.bed > just.Afr.genes.bed #79

grep -F -x -v -f African.Genes.bed European.Genes.bed  > just.Eur.genes.bed #84
 

#most exteme genes


sort afrANDeur.genes.bed just.Afr.genes.bed  just.Eur.genes.bed |uniq -c > extreme.genes.bed

### extreme set of significant genes

sort <(cat prot.cod.CANDf0.5.LWK.bed prot.cod.CANDf0.4.LWK.bed prot.cod.CANDf0.3.LWK.bed |sort|uniq) <(cat prot.cod.CANDf0.5.YRI.bed prot.cod.CANDf0.4.YRI.bed prot.cod.CANDf0.3.YRI.bed |sort|uniq) |uniq -d > cand.African.Genes.bed


sort <(cat prot.cod.CANDf0.5.GBR.bed prot.cod.CANDf0.4.GBR.bed prot.cod.CANDf0.3.GBR.bed |sort|uniq) <(cat prot.cod.CANDf0.5.TSI.bed prot.cod.CANDf0.4.TSI.bed prot.cod.CANDf0.3.TSI.bed |sort|uniq) |uniq -d > cand.European.Genes.bed


sort cand.African.Genes.bed cand.European.Genes.bed |uniq -d > cand.afrANDeur.genes.bed


grep -F -x -v -f cand.European.Genes.bed cand.African.Genes.bed > cand.just.Afr.genes.bed

grep -F -x -v -f cand.African.Genes.bed cand.European.Genes.bed  > cand.just.Eur.genes.bed
 

#most exteme genes


sort cand.afrANDeur.genes.bed cand.just.Afr.genes.bed  cand.just.Eur.genes.bed |uniq -c > cand.extreme.genes.bed


