#source ~/.bashrc
#conda activate r
Rscript fit_states_death.R data/all_states.list day 0.001:1:0.6:-80:0.4 24 &>day.4.log
