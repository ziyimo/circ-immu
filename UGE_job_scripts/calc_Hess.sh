#!/bin/bash
#$ -N calc_Hess
#$ -S /bin/bash
#$ -cwd
#$ -o UGE$JOB_ID.o
#$ -j y
#$ -l m_mem_free=2G
#$ -pe threads 26
#$ -binding linear_per_task:1

# module load EBModules
# module load R/4.0.3-foss-2020a

export R_LIBS_USER=/grid/siepel/home_norepl/mo/R/x86_64-pc-linux-gnu-library/4.0

MODE=$1
STATES=$2
R0MOD=$3
HYPRM=$4
NDEPS=$5
PARSCALE=$6

echo "_START_$(date)"
echo "Mode: ${MODE}; Model: ${R0MOD}; Lambda/hosp_rate: ${HYPRM}"

./hessian_${MODE}fit.R $STATES $R0MOD $HYPRM $NDEPS $PARSCALE 25 # manually code the number of threads

echo "_EXITSTAT_$?"
echo "_END_$(date)"
exit
