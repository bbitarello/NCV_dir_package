###################################################################
#	Barbara Bitarello
#
#	Last modified: 28.11.2016
#	Calcualte proportions of the genome in the scan, etc.
###################################################################


#positions that passed all filters (FD, TRF, SD, MAp50)

/mnt/sequencedb/PopGen/cesare/hg19/bedfiles/hg19.pantro2.Map50_100.TRF.SDs.bed.gz 

#this file explains how the above one was generated:

/mnt/sequencedb/PopGen/cesare/hg19/bedfiles/README 

# The new filtered data from 1000Genome phase1 are here:

/mnt/sequencedb/PopGen/cesare/1000G/phase1/vcfs/AC_13pops_chr9.hg19.pantro2.Map50_100.TRF.SDs.tsv.gz 





#number of positions in hg19 with orthology to chimp (chr1-22)

gunzip -c hg19.pantro2.bed.gz|awk '$a=$3-$2{print $a}'|awk '{s+=$1} END {print s}' #2571401270

#number of positions after all filters (chr1-22)

gunzip -c hg19.pantro2.Map50_100.TRF.SDs.bed.gz |awk '$a=$3-$2{print $a}'|awk '{s+=$1} END {print s}' #83% of all

#next: calculate how much was lost after the Nr.IS filter

#and check each filter separately: FD, TRF, Map50


#50-mer

gunzip -c /mnt/sequencedb/PopGen/cesare/hg19/mappability_tracks/wgEncodeCrgMapabilityAlign50mer_100.bed.gz|awk '$a=$3-$2{print $a}'|awk '{s+=$1} END {print s}' #2391379832 (92% of all)


#number of positions to remove with SF filter

cat /mnt/sequencedb/PopGen/joao/Exome_Project/Results/Filtering/SDs/Human/Human_MAPPED_COORDINATES.bed|awk '$1!="chrY"{print$0}'|awk '$1!="chrX"{print $0}'|awk '$a=$3-$2{print $a}'|awk '{s+=$1} END {print s}' #105361326 (96 is retained)



