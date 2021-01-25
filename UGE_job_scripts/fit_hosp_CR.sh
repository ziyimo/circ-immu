#!/bin/bash
#$ -N fit_SEIH
#$ -S /bin/bash
#$ -cwd
#$ -o UGE$JOB_ID.o
#$ -j y
#$ -l m_mem_free=1G
#$ -pe threads 5
#$ -binding linear_per_task:1

# module load EBModules
# module load R/4.0.3-foss-2020a

export R_LIBS_USER=/grid/siepel/home_norepl/mo/R/x86_64-pc-linux-gnu-library/4.0

R0MOD=$1
SEED=$2
#THD=$3

echo "_START_$(date)"
echo "Model: ${R0MOD}" #; Restart Threshold: ${THD}"

Rscript fit_CR.R $R0MOD $SEED

echo "_EXITSTAT_$?"
echo "_END_$(date)"
exit
