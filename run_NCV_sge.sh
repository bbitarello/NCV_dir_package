######################################################################################################
#   Bárbara D Bitarello
#    Created: 13.10.2015
#    Last modified: 20.09.2016
#    Description: this set of commands runs the NCV scan in parallel jobs using SGE (Sun Grid Engine)
#	          If it optimized for running at the MPI-EVA cluster (Leipzig)
#######################################################################################################


#calculate NCV in 948 parallel jobs (one for every 3Mbp)


###################################
# If cluster uses SGE
##################################
BP=3000
SLIDE=1500
LOGS=/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scratch/logs
#LOGS=~/NCV_dir_package/scratch/logs/
    for CHROM in {1..22}; do
    INPUT=/mnt/sequencedb/PopGen/barbara/NCV_dir_package/input_data/chr${CHROM}/AC_13pops_chr${CHROM}.hg19.pantro2.Map50_100.TRF.SDs.tsv.gz
        CHIMPfd=/mnt/sequencedb/PopGen/barbara/NCV_dir_package/input_data/outgroup_files/fds.chr${CHROM}.hg19_pantro2.Map50_100.TRF.SDs.bed.gz
        POSITIONS=/mnt/sequencedb/PopGen/barbara/NCV_dir_package/bins/scan_bin3Mb_chr${CHROM}.pos
            NBINS=$(wc -l ${POSITIONS} | cut -f 1 -d ' ') # the number of bins for the chromosome
            for i in `seq 1 ${NBINS}`; do
                  TMPDIR=/mnt/sequencedb/PopGen/barbara/NCV_dir_package/tmp/chr${CHROM}/bin${i}/
                               POS=$(sed -n ${i}p ${POSITIONS})


qsub -e ${LOGS} -o ${LOGS} /mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/run_ncv_allpops_Rscript.sge ${INPUT} ${POS} ${BP} ${SLIDE} ${TMPDIR} ${CHIMPfd} ${i}
    done
done

