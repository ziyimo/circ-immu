#!/bin/bash
#$ -N fit_SIRS
#$ -S /bin/bash
#$ -cwd
#$ -o UGE$JOB_ID.o
#$ -j y
#$ -l m_mem_free=2G
#$ -pe threads 16
#$ -binding linear_per_task:1

# module load EBModules
# module load R/4.0.3-foss-2020a

export R_LIBS_USER=/grid/siepel/home_norepl/mo/R/x86_64-pc-linux-gnu-library/4.0

STATES=$1
R0MOD=$2
SHARE_INTRCPT=$3
LAMBDA=$4

echo "_START_$(date)"

Rscript fit_all.R $STATES $R0MOD $SHARE_INTRCPT $LAMBDA

echo "_EXITSTAT_$?"
echo "_END_$(date)"
exit
