#source ~/.bashrc
#conda activate r
Rscript fit_cities_death_global_clean.R data/globalcity_clean.txt sd 0.001:1:1.2:-289:0.57:-104:0.379 32 &> sd_globalcityclean.2.log
