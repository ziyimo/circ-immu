#source ~/.bashrc
#conda activate r
Rscript fit_states_death.R data/all_states.list hd 0.001:0.8:0.8:0:-80:0.32 32 &>hd.3.log
