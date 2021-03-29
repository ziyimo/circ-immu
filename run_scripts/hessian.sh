#source ~/.bashrc
#conda activate r
for l in 0.005 0.001 0.01 0.006 0.007 0.008 0.009 0.004 0.003 0.002 0.02;do
    Rscript hessian_deathfit.R data/all_states.list sd fit_results/data/all_states.list_sd_02221502_iter11_iterative.fit3.rds $l 1 24
done
