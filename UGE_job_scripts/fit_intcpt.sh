#!/bin/bash
#$ -N fit_intcpt
#$ -S /bin/bash
#$ -cwd
#$ -o UGE$JOB_ID.o
#$ -j y
#$ -l m_mem_free=2G
#$ -pe threads 44
#$ -binding linear_per_task:1

# module load EBModules
# module load R/4.0.3-foss-2020a

export R_LIBS_USER=/grid/siepel/home_norepl/mo/R/x86_64-pc-linux-gnu-library/4.0

STATES=$1
R0MOD=$2
INIT=$3
LAMBDA=$4

echo "_START_$(date)"

./fit_intcpt.R $STATES $R0MOD $INIT $LAMBDA 43 # manually code the number of threads

echo "_EXITSTAT_$?"
echo "_END_$(date)"
exit
