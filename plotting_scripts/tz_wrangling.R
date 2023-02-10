#install.packages("lutz")

library(lutz)

state_pop_center <- read.table("state_lv_data/CenPop2010_Mean_ST.txt", header = TRUE, sep = ",")

no_juri <- nrow(state_pop_center)
off_center <- numeric(no_juri)

for (idx in seq(no_juri)){
  timezone <- tz_lookup_coords(state_pop_center$LATITUDE[idx], state_pop_center$LONGITUDE[idx], method = "accurate")
  tz_longitude <- 15*tz_offset("2021-01-01", timezone)$utc_offset_h
  off_center[idx] <- state_pop_center$LONGITUDE[idx]-tz_longitude
  
  cat(state_pop_center$STNAME[idx], timezone, tz_longitude, off_center[idx], "\n")
}

state_pop_center$OFFCENTER <- off_center

state_pop <- read.delim("state_lv_data/state_pop.tsv")
out_df <- merge(state_pop, state_pop_center, by.x = "state", by.y = "STNAME")
out_df <- subset(out_df, select = -c(STATEFP, POPULATION))

write.table(out_df, file = "state_lv_data/CenPop2010_Mean_ST_tzoffset.txt", quote = FALSE, sep = "\t", row.names = FALSE)
