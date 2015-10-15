######################################################################################################
#   Bárbara D Bitarello
#    Created: 13.10.2015
#    Last modified: 15.10.2015
#    Description: this set of commands runs the NCV scan in parallel jobs using SGE (Sun Grid Engine)
#	          If it optimized for running at the MPI-EVA cluster (Leipzig)
#######################################################################################################


#calculate NCV in 948 parallel jobs (one for every 3Mbp)


###################################
# If cluster uses SGE
##################################
BP=3000
SLIDE=1500
LOGS=~/NCV_dir_package/scratch/logs/
    for CHROM in {1..22}; do
    INPUT=~/NCV_dir_package/input_data/chr${CHROM}/AC_13pops.Map50_100.TRF.SDs.pantro2.tsv.gz
        CHIMPfd=~/NCV_dir_package/input_data/outgroup_files/fds.hg19_pantro2.$CHROM.tsv.gz
        POSITIONS=~/NCV_dir_package/bins/scan_bin3Mb_chr${CHROM}.pos
            NBINS=$(wc -l ${POSITIONS} | cut -f 1 -d ' ') # the number of bins for the chromosome
            for i in `seq 1 ${NBINS}`; do
                  TMPDIR=/mnt/scratch/cee/bs_genomescan/ncv_run/chr${CHROM}/bin${i}/
                    TMPDIR=~/NCV_dir_package/scratch/temp_files/${CHROM}/bin${i}/
                               POS=$(sed -n ${i}p ${POSITIONS})

