#Cesare de Filippo, MPI-EVA
# modified on 13-October-2015 by Barbara Bitarello 

RSCRIPT=~/NCV_dir_package/scripts/run_ncv_allpops_Rscript_v1.r
INPUT=${1} # the imput file as allele counts such as (AC_13pops.tsv.gz)
POSITIONS=${2} # the interval to use as CHROMOSOME, START and END positions
WINDOW=${3} # the window length in bp
SLIDE=${4}
TMPDIR=${5}
CHIMPfd=${6}
BIN=${7}
mkdir -p ${TMPDIR}
cd ${TMPDIR}

tabix -p vcf ${INPUT} ${POSITIONS} > tmp.ac 

${RSCRIPT} -i 'tmp.ac' -w ${WINDOW} -s ${SLIDE} -b ${BIN} -fd ${CHIMPfd}

cd ~/NCV_dir_package/
echo DONE

