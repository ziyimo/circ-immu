# Prefix for bootstrap logs
BSPFX="cleancity.BS"
# Seed values from best sd fit
RDS="/grid/siepel/home_norepl/scheben/projects/virus/covid_branch/circ-immu/covid/fit/fit_results/data/globalcity_clean.txt_sd_03150903_iter10_iterative.city.glob.fit1.clean3.rds"
# Region list
LIST="/grid/siepel/home_norepl/scheben/projects/virus/covid_branch/circ-immu/covid/fit/data/globalcity_clean.txt"

# Get resampled param lists
#grep -A 1 -r Resamp_Iter ${BSPFX}*log| grep -v "\-\-"| paste - -| cut -d' ' -f5,8,11-| tr ' ' '\t'| awk '$1<$2'| cut -f3-| tr '\t' ':'| sed 's/:$//' > bootstrap_results_states.out

# Get initial infected per region
Rscript parse_seedvals.R ${RDS}| sed '/^$/d'| paste - -| grep "\\$"| sed 's/\$//'| sed 's/`/"/g'| sed 's/\[1\] //' > ${RDS%%.rds}.ssfit.txt.tmp
paste -d' ' <(cut -f2 ${RDS%%.rds}.ssfit.txt.tmp) <(cut -f1 ${RDS%%.rds}.ssfit.txt.tmp) > ${RDS%%.rds}.ssfit.txt
rm -f ${RDS%%.rds}.ssfit.txt.tmp


# Combine to generate cmd file
### SD ###
while read l; do
    #sed "s@^@Rscript plot_all_city_clean_permDST.R $LIST sd @" bootstrap_results_states.out| sed "s@$@ $l@" > makebs.cmd
    sed "s@^@Rscript plot_all_city_clean.R $LIST sd @" bootstrap_results_states.out| sed "s@\$@ $l@" >> makebs_sd.cmd
    sed "s@^@Rscript plot_all_city_clean_permST.R $LIST sd @" bootstrap_results_states.out| sed "s@\$@ $l@" >> makebs_ST.cmd
    sed "s@^@Rscript plot_all_city_clean_permDST.R $LIST sd @" bootstrap_results_states.out| sed "s@\$@ $l@" >> makebs_DST.cmd
done<${RDS%%.rds}.ssfit.txt

# Test breaker to manually check outputs
# Remove this after test
######
######
exit 1
# Generate the CSV outputs
bash makebs_sd.cmd
cd fit_results
#cat $LIST| while read m; do ls ${m}_*| while read l; do state=`echo "$l"| sed 's/_.*//'`; USCOUNTER=$(expr $USCOUNTER + 1); sed -i "s/^/${state},${USCOUNTER},/" $l;done;done
# Remove space from city names in filename
rename " " "_" *
# add region name to files
#for l in *csv; do name=`echo "$l" | sed 's/_R0_.*//'`; sed -i "s/^/${name},/" $l;sed -i 's/^/sd,/' $l;done
cat $LIST|sed 's/ /_/g' | while read m; do ls ${m}_*| while read l; do state=`echo "$l"| sed 's/_R0_.*//'`; USCOUNTER=$(expr $USCOUNTER + 1); sed -i "s/^/sd,${state},${USCOUNTER},/" $l;done;done
cat *csv| grep -v death > sd.out
rm -f *csv
cd ..

exit 1
### DST ###
# Generate the CSV outputs
bash makebs_DST.cmd
cd fit_results
# Remove space from city names in filename
rename " " "_" *
# add region name to files
cat $LIST| sed 's/ /_/g' |while read m; do ls ${m}_*| while read l; do state=`echo "$l"| sed 's/_R0_.*//'`; USCOUNTER=$(expr $USCOUNTER + 1); sed -i "s/^/sd_DST,${state},${USCOUNTER},/" $l;done;done
cat *csv| grep -v death > DST.out
rm -f *csv
cd ..
### ST ###
# Generate the CSV outputs
bash makebs_ST.cmd
cd fit_results
# Remove space from city names in filename
rename " " "_" *
# add region name to files
cat $LIST| sed 's/ /_/g' |while read m; do ls ${m}_*| while read l; do state=`echo "$l"| sed 's/_R0_.*//'`; USCOUNTER=$(expr $USCOUNTER + 1); sed -i "s/^/sd_ST,${state},${USCOUNTER},/" $l;done;done
cat *csv| grep -v death > ST.out
rm -f *csv
cd ..

# header for final
echo "model,state,bs,time,S,E,I,G,R,Rt,alpha,deathIncrease_cor,adj_G" > fit_results/sd_DST_ST.out
cat fit_results/ST.out fit_results/DST.out-fit_results/sd.out >> fit_results/sd_DST_ST.out






