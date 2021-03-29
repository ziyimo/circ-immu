#source ~/.bashrc
#conda activate r
Rscript fit_states_death.R data/all_states.list cos 0.001:1:0.8:180 24 &>cos.4.log
