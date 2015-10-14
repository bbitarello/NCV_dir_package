########################################################################## 
#   Bárbara D Bitarello
#    Created: 13.10.2015
#    Last modified: 14.10.2015
#    Description: this set of commands runs the NCV scan in parallel jobs
##########################################################################


#calculate NCV in 948 parallel jobs (one for every 3Mbp)


###################################
# If cluster uses SGE
##################################
#BP=3000
#SLIDE=1500
#LOGS=~/NCV_dir_package/scratch/logs/
#    for CHROM in {1..22}; do
#    INPUT=~/NCV_dir_package/input_data/chr${CHROM}/AC_13pops.Map50_100.TRF.SDs.pantro2.tsv.gz
#        CHIMPfd=~/NCV_dir_package/input_data/outgroup_files/fds.hg19_pantro2.$CHROM.tsv.gz
#        POSITIONS=~/NCV_dir_package/bins/scan_bin3Mb_chr${CHROM}.pos
#            NBINS=$(wc -l ${POSITIONS} | cut -f 1 -d ' ') # the number of bins for the chromosome
#            for i in `seq 1 ${NBINS}`; do
            #       TMPDIR=/mnt/scratch/cee/bs_genomescan/ncv_run/chr${CHROM}/bin${i}/
#                    TMPDIR=~/NCV_dir_package/scratch/temp_files/${CHROM}/bin${i}/
#                               POS=$(sed -n ${i}p ${POSITIONS})
#                                  qsub -e ${LOGS} -o ${LOGS} run_ncv_allpops_Rscript.sge ${INPUT} ${POS} ${BP} ${SLIDE} ${TMPDIR} ${CHIMPfd} ${i}
#                done
#                done
#
# 
#############################
#If cluster does not use SGE
#############################

BP=3000
SLIDE=1500
LOGS=~/NCV_dir_package/scratch/logs/
    for CHROM in {1..22}; do    
    INPUT=~/NCV_dir_package/input_data/chr${CHROM}/AC_13pops.Map50_100.TRF.SDs.pantro2.tsv.gz
     CHIMPfd=~/NCV_dir_package/input_data/outgroup_files/fds.hg19_pantro2.$CHROM.tsv.gz
     POSITIONS=~/NCV_dir_package/bins/scan_bin3Mb_chr${CHROM}.pos
     NBINS=$(wc -l ${POSITIONS} | cut -f 1 -d ' ') # the number of bins for the chromosome
     for i in `seq 1 ${NBINS}`; do
          TMPDIR=~/NCV_dir_package/scratch/temp_files/${CHROM}/bin${i}/
            POS=$(sed -n ${i}p ${POSITIONS})

                ~/NCV_dir_package/scripts/run_ncv_allpops_master_script.sh ${INPUT} ${POS} ${BP} ${SLIDE} ${TMPDIR} ${CHIMPfd} ${i}  &   #hopefully the ampersand will do the trick for parallelizing.
done
done
#Without parallelizinh this is taking more than 50 to complete only a few jobs.
#The block above works, but the  commnads are not running in parallel

##########################################################################################
#using GNU parallel ##
#chech if parallel is isntaled with "which parallel"
#https://www.gnu.org/software/parallel/parallel_tutorial.html



function run_NCV {
CHROM=$1
BP=3000
SLIDE=1500
LOGS=~/NCV_dir_package/scratch/logs/
INPUT=~/NCV_dir_package/input_data/chr${CHROM}/AC_13pops.Map50_100.TRF.SDs.pantro2.tsv.gz
CHIMPfd=~/NCV_dir_package/input_data/outgroup_files/fds.hg19_pantro2.${CHROM}.tsv.gz
POSITIONS=~/NCV_dir_package/bins/scan_bin3Mb_chr${CHROM}.pos
NBINS=$(wc -l ${POSITIONS} | cut -f 1 -d ' ') # the number of bins for the chromosome

     for i in `seq 1 ${NBINS}`; do
	TMPDIR=~/NCV_dir_package/scratch/temp_files/${CHROM}/bin${i}/
            POS=$(sed -n ${i}p ${POSITIONS})

                ~/NCV_dir_package/scripts/run_ncv_allpops_master_script.sh ${INPUT} ${POS} ${BP} ${SLIDE} ${TMPDIR} ${CHIMPfd} ${i}  &   
done

}


export -f run_NCV


parallel --gnu --joblog ~/NCV_dir_package/scratch/logs/log run_NCV ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22

