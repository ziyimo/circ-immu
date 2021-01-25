#!/usr/bin/env Rscript

#.libPaths() # sanity check

source("R0_mods.R")
source("SEIH_mod.R")
var_cache <- ls()

time_stamp <- format(Sys.time(), "%m%d%H")
args <- commandArgs(trailingOnly=TRUE)

state_code <- args[1]
R0_mod <- args[2]        # R0 model, options: cos, hum, sun, hs
optimizer <- "pso"

library(optimizer, character.only = TRUE)
handle <- paste0(state_code, "_", R0_mod, "_", time_stamp, "_", optimizer)

## Load all state data
state_pop <- read.delim("state_lv_data/cr_pop.tsv")
all_state_sun <- read.csv("state_lv_data/census_reg_daily_sunrise_2019.csv")
#all_state_day <- read.csv("state_lv_data/state_daytime_2019.csv")
all_state_hum <- read.csv("state_lv_data/census_reg_humidity_2014_2018.csv")
all_state_hum <- all_state_hum[all_state_hum$date != "2016-02-29", ] # get rid of leap year Feb 29
## Load COVID hospitalization data
covid_df <- readRDS("state_lv_data/censusReg_hospitalization.rds")

cat(">>> Loading state data:", state_code, "<<<\n")
census_pop <- state_pop$pop[state_pop$Region == state_code]

if (R0_mod %in% c("sun", "hs", "sd", "hsd")){
  sunob <- all_state_sun[[state_code]]/720 # 365 days, scaled to maximum 720 = 12 hours
}
if (R0_mod %in% c("day", "hd", "sd", "hsd")){
  dayob <- all_state_day[[state_code]]/1440
}
if (R0_mod %in% c("hum", "hs", "hd", "hsd")){
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
} else if (R0_mod == "sun"){
  varob <- list(sunob)
} else if (R0_mod == "hs"){
  varob <- list(climob, sunob)
}

state_df <- subset(covid_df, Region == state_code)
state_df <- state_df[state_df$date <= 365, ]
state_hos <- subset(state_df, select = c("date", "hospitalizedCurrently"))

#state_data[[state_code]] <- list(idx=state_i, pop=census_pop, var=varob, epi=state_hos)


# seed_prop; hrate; R0_min; R0_range; [R0_prms]
prms_low = c(1e-5, 0.005, 0.6, 0.5, param_bounds[[R0_mod]]$low)
prms_high = c(0.05, 0.15, 1.6, 3, param_bounds[[R0_mod]]$high)
cat("seed_prop; hrate; R0_min; R0_range; [R0_prms]\n")
cat(paste(prms_low, sep = "\t"))
cat("\n")
cat(paste(prms_high, sep = "\t"))
cat("\n")

cat(">> Fitting with R0 model:", R0_mod, "\nLower:", prms_low, "\nUpper:", prms_high, "\n")
restart_thd <- c(1.1e-1, 1.1e-2, 1.1e-3, 1.1e-4)
seed <- rep(NA,length(prms_low))

negLL_wrapper <- function(norm_prms){ # everything else read as global variable
  orig_prms <- prms_low + norm_prms*(prms_high-prms_low)
  return(pois_L(orig_prms, paste0("R0_", R0_mod), varob, state_hos, census_pop))
}

for (thrhld in restart_thd){
  cat(">> Restart threshold:", thrhld, "\n")
  oo <- psoptim(seed, negLL_wrapper,
                lower=rep(0,length(prms_low)), upper=rep(1,length(prms_low)),
                control=list(trace=1, REPORT=5, maxit=10000, trace.stats=FALSE, maxit.stagnate=500,
                             max.restart=5, reltol=thrhld))
  seed <- oo$par
  cat("seed_prop; hrate; R0_min; R0_range; [R0_prms]\n")
  cat(prms_low + oo$par*(prms_high-prms_low))
  cat("\n")
}

saveRDS(oo, file = paste0("fit_results/", handle, ".rds"))

