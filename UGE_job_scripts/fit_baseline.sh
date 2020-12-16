#!/bin/bash
#$ -N fit_bl
#$ -S /bin/bash
#$ -cwd
#$ -o UGE$JOB_ID.o
#$ -j y
#$ -l m_mem_free=8G

# module load EBModules
# module load R/4.0.3-foss-2020a
export R_LIBS_USER=/grid/siepel/home_norepl/mo/R/x86_64-pc-linux-gnu-library/4.0

STATES=$1
LAMBDA=$2

echo "_START_$(date)"

./fit_baseline.R $STATES $LAMBDA

echo "_EXITSTAT_$?"
echo "_END_$(date)"
exit
