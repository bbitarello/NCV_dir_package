#############################################################################################################################################
#	Author: Barbara Bitarello
#	Created: 15.10.2015
#	Last modified: 15.10.2015
#	Description: this bash function uses 'GNU parallel' to parallelize the jobs 
#		This current version does not allow specification of how much VMEM to use, etc, so it could still be improved, but it works.
#############################################################################################################################################
#Requirements (ask your system administrator to install these if  'which Rscript'and 'which parallel' do not return a path.
#Rscript
#GNU parallel

##########################################################################################
#using GNU parallel ##
#if you want to understand better and possibly make better use of Linux's parallelizing possibilities, check the link below
#https://www.gnu.org/software/parallel/parallel_tutorial.html

#this is a bash function 
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

#here we export the bash function
export -f run_NCV   

#here we run a separate job for each chromosome (but the bins within each chromosome will also parallelize)
parallel --gnu --joblog ~/NCV_dir_package/scratch/logs/log run_NCV ::: `seq 1 22`


