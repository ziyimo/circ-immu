#!/bin/bash
#$ -N boot_SIRS
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

STATES=$1
BESTFIT=$2
TAG=$3
LAMBDA=$4

echo "_START_$(date)"
#echo "Model: ${R0MOD}; Lambda: ${LAMBDA}"

./fit_postest_boot.R $STATES $BESTFIT $TAG $LAMBDA 31 # manually code the number of threads

echo "_EXITSTAT_$?"
echo "_END_$(date)"
exit
