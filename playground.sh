#####################################################
#	This is a playground, kind of obsolete.
#
#####################################################


#function test_func {
#i=$1
#TMPDIR=~/NCV_dir_package/scratch/temp_files/${CHROM}/bin${i}/
#POS=$(sed -n ${i}p ${POSITIONS})

# ~/NCV_dir_package/scripts/run_ncv_allpops_master_script.sh ${INPUT} ${POS} ${BP} ${SLIDE} ${TMPDIR} ${CHIMPfd} ${i}  &
#}

#export -f test_func

#for CHROM in {1..22}; do
#BP=3000
#SLIDE=1500
#LOGS=~/NCV_dir_package/scratch/logs/
#INPUT=~/NCV_dir_package/input_data/chr${CHROM}/AC_13pops.Map50_100.TRF.SDs.pantro2.tsv.gz
#CHIMPfd=~/NCV_dir_package/input_data/outgroup_files/fds.hg19_pantro2.${CHROM}.tsv.gz
#POSITIONS=~/NCV_dir_package/bins/scan_bin3Mb_chr${CHROM}.pos
#NBINS=$(wc -l ${POSITIONS} | cut -f 1 -d ' ') # the number of bins for the chromosome
#parallel --gnu --joblog ~/NCV_dir_package/scratch/logs/log test_func ::: `seq 1 $NBINS`
#done





