############################################################
#	Barbara Bitarello
#
#	Obtain list of protein-coding genes from bedfiles
###########################################################



bedtools intersect -a <(sed -e 's/^/chr/' Union.top816.0.5_0.4_0.3_LWK.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.Union.top829.0.5_0.4_0.3_LWK.bed

bedtools intersect -a <(sed -e 's/^/chr/' Union.top816.0.5_0.4_0.3_YRI.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.Union.top829.0.5_0.4_0.3_YRI.bed

bedtools intersect -a <(sed -e 's/^/chr/' Union.top816.0.5_0.4_0.3_GBR.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.Union.top829.0.5_0.4_0.3_GBR.bed

bedtools intersect -a <(sed -e 's/^/chr/' Union.top816.0.5_0.4_0.3_TSI.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.Union.top829.0.5_0.4_0.3_TSI.bed

bedtools intersect -a <(sed -e 's/^/chr/' Union.CAND.0.5_0.4_0.3_LWK.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.Union.CAND.0.5_0.4_0.3_LWK.bed

bedtools intersect -a <(sed -e 's/^/chr/' Union.CAND.0.5_0.4_0.3_YRI.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.Union.CAND.0.5_0.4_0.3_YRI.bed

bedtools intersect -a <(sed -e 's/^/chr/' Union.CAND.0.5_0.4_0.3_GBR.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.Union.CAND.0.5_0.4_0.3_GBR.bed

bedtools intersect -a <(sed -e 's/^/chr/' Union.CAND.0.5_0.4_0.3_TSI.bed) -b /mnt/sequencedb/PopGen/barbara/scan_may_2014/10000_sims_per_bin/ensembl_hg19.bed -wo|awk '$9=="protein_coding"{print $8}'|sort|uniq > prot.cod.Union.CAND0.5_0.4_0.3_TSI.bed

