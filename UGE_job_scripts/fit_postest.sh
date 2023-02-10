#!/bin/bash
#$ -N fit_SIRS
#$ -S /bin/bash
#$ -cwd
#$ -o UGE$JOB_ID.o
#$ -j y
#$ -l m_mem_free=2G
#$ -pe threads 32
#$ -binding linear_per_task:1

# module load EBModules
# module load R/4.0.3-foss-2020a

export R_LIBS_USER=/grid/siepel/home_norepl/mo/R/x86_64-pc-linux-gnu-library/4.0

REGION=$1
R0MOD=$2
#SHARE_INTRCPT=$3
LAMBDA=$3

echo "_START_$(date)"
echo "Model: ${REGION}/${R0MOD}; Lambda: ${LAMBDA}"

Rscript fit_postest.R $REGION $R0MOD $LAMBDA 30 # manually code the number of threads

echo "_EXITSTAT_$?"
echo "_END_$(date)"
exit
