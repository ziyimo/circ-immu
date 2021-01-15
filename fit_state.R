#!/usr/bin/env Rscript

library(DEoptim)

source("R0_mods.R")
source("SEIH_mod.R")

time_stamp <- format(Sys.time(), "%m%d%H")
args <- commandArgs(trailingOnly=TRUE)

state_code <- args[1]   # 2-letter state code
R0_mod <- args[2]       # functional form of R0

handle <- paste0(state_code, R0_mod, time_stamp)
sink(paste0("fit_results/", handle, ".log"))

## Load all state data
state_pop <- read.delim("state_lv_data/state_pop.tsv")
all_state_sun <- read.csv("state_lv_data/state_daily_sunrise_2019.csv")
all_state_day <- read.csv("state_lv_data/state_daytime_2019.csv")
all_state_hum <- read.csv("state_lv_data/state_humidity_2014_2018.csv")
all_state_hum <- all_state_hum[all_state_hum$X != "2016-02-29", ] # get rid of leap year Feb 29

cat(">>> Loading state data:", state_code, "<<<\n")
census_pop <- state_pop$pop[state_pop$code==state_code]

if (R0_mod %in% c("sd", "hsd")){
  sunob <- all_state_sun[[state_code]]/720 # 365 days, scaled to maximum 720 = 12 hours
}
if (R0_mod %in% c("day", "hd", "sd", "hsd")){
  dayob <- all_state_day[[state_code]]/1440
}
if (R0_mod %in% c("hum", "hd", "hsd")){
  climob <- all_state_hum[[state_code]] # multiple years
  climob <- matrix(climob, nrow = length(climob)/365, ncol = 365, byrow = TRUE)
  climob <- colMeans(climob) # 365 days
}

if (R0_mod == "cos"){
  varob <- list(seq(365))
} else if (R0_mod == "hum"){
  varob <- list(climob)
} else if (R0_mod == "day"){
  varob <- list(dayob)
} else if (R0_mod == "hd"){
  varob <- list(climob, dayob)
} else if (R0_mod == "sd"){
  varob <- list(sunob, dayob)
} else if (R0_mod == "hsd"){
  varob <- list(climob, sunob, dayob)
}

## Load flu data
covid_df <- readRDS("state_lv_data/state_hospitalization.rds")
covid_df$date <- as.numeric(as.Date(covid_df$date, format="%Y-%m-%d") - as.Date("2019-12-31", format="%Y-%m-%d"))

state_df <- subset(covid_df, state == state_code)
state_df <- state_df[order(state_df$date),]
state_df <- state_df[!is.na(state_df$hospitalizedCurrently), ]
state_df <- state_df[state_df$date <= 365, ]
state_hos <- subset(state_df, select = c("date", "hospitalizedCurrently"))

## Optimization
cat(">> Fitting:", state_code, "; R0 model:", R0_mod, "\n")

# seed_prop; hrate; R0_min; R0_range; [R0_prms]
prms_low = c(1e-5, 0.01, 0.3, 0.8, param_bounds[[R0_mod]]$low)
prms_high = c(0.1, 0.2, 3, 8, param_bounds[[R0_mod]]$high)
cat("seed_prop; hrate; R0_min; R0_range; [R0_prms]\n")
cat(paste(prms_low, sep = "\t"))
cat("\n")
cat(paste(prms_high, sep = "\t"))
cat("\n")

evo_optim <- DEoptim(pois_L, lower=prms_low, upper=prms_high,
                     control=DEoptim.control(trace = 5, reltol = 1e-6, itermax = 2000, steptol = 200),
                     R0_model = paste0("R0_", R0_mod), var_obs = varob,
                     hosp_df = state_hos, pop_size = census_pop)
saveRDS(evo_optim, file = paste0("fit_results/", handle, "_DE.rds"))
sink()
