#!/bin/bash
#$ -N bootstrap_SEIH
#$ -S /bin/bash
#$ -cwd
#$ -o UGE$JOB_ID.o
#$ -j y
#$ -l m_mem_free=1G
#$ -pe threads 26
#$ -binding linear_per_task:1

# module load EBModules
# module load R/4.0.3-foss-2020a

export R_LIBS_USER=/grid/siepel/home_norepl/mo/R/x86_64-pc-linux-gnu-library/4.0

HOSPR=$1
IDTFR=$2

#SEED_SH=$3 # incb_prd, inf_prd, hpt_rate, hpt_prd, R0min, R0range, [R0_params]
#TUNE_MASK=$4 # prop_init, incb_prd, inf_prd, hpt_rate, hpt_prd, R0min, R0range

echo "_START_$(date)"
#echo "Model: ${R0MOD}"

./fit_hosp_boot.R states_49DC.tsv covid_hosp_fit/eta${HOSPR}_sd_iter.rds eta${HOSPR}_${IDTFR} 5:5:${HOSPR}:10:1:0.8 T:F:F:F:F:T:T 25

echo "_EXITSTAT_$?"
echo "_END_$(date)"
exit
