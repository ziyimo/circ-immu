#source ~/.bashrc
#conda activate r
#sleep 10m
Rscript fit_cities_global_clean_boot.R data/globalcity_clean.txt sd fit_results/data/globalcity_clean.txt_sd_03150903_iter10_iterative.city.glob.fit1.clean3.rds 18 &> cleancity.BS.rep3.log
