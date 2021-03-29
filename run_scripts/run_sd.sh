#source ~/.bashrc
#conda activate r
Rscript fit_states_death.R data/all_states.list sd 0.001:0.8:1:-45:0.5:-13:0.32 24 &>sd.5.log
