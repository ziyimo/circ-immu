#!/bin/bash
#$ -N fit_SEIH
#$ -S /bin/bash
#$ -cwd
#$ -o UGE$JOB_ID.o
#$ -j y
#$ -l m_mem_free=2G
#$ -pe threads 10
#$ -binding linear_per_task:1

# module load EBModules
# module load R/4.0.3-foss-2020a

export R_LIBS_USER=/grid/siepel/home_norepl/mo/R/x86_64-pc-linux-gnu-library/4.0

STATES=$1
R0MOD=$2

echo "_START_$(date)"
echo "Model: ${R0MOD}"

Rscript fit_hosp.R $STATES $R0MOD 8 # manually code the number of threads

echo "_EXITSTAT_$?"
echo "_END_$(date)"
exit
