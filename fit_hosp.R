#!/usr/bin/env Rscript

#.libPaths() # sanity check

source("R0_mods.R")
source("SEIH_mod.R")
var_cache <- ls()

time_stamp <- format(Sys.time(), "%m%d")
args <- commandArgs(trailingOnly=TRUE)

state_ls <- args[1]      # text file with a list of states to fit to
R0_mod <- args[2]        # R0 model, options: cos, hum, day, hd, sd, hsd
#share_intcpt <- args[3]  # "0" or "1"; deprecated
nthr <- as.numeric(args[3]) # limit no. of threads

if (nthr == 1){
  optMethod <- "pso"
  library("pso")
} else{
  optMethod <- "DEoptim"
  library("parallel")
  library("DEoptim")
}

handle <- paste0(state_ls, "_", R0_mod, "_", time_stamp, "_", optMethod)
#share_intcpt <- TRUE # as.logical(as.numeric(share_intcpt))
sink(paste0("fit_results/", handle, ".log"))

states <- read.delim(state_ls, header = FALSE)
no_states <- nrow(states)

## Load all state data
state_pop <- read.delim("state_lv_data/state_pop.tsv")
all_state_sun <- read.csv("state_lv_data/state_daily_sunrise_2019.csv")
all_state_day <- read.csv("state_lv_data/state_daytime_2019.csv")
all_state_hum <- read.csv("state_lv_data/state_humidity_2014_2018.csv")
all_state_hum <- all_state_hum[all_state_hum$X != "2016-02-29", ] # get rid of leap year Feb 29
## Load COVID hospitalization data
covid_df <- readRDS("state_lv_data/state_hospitalization.rds")
covid_df$date <- as.numeric(as.Date(covid_df$date, format="%Y-%m-%d") - as.Date("2019-12-31", format="%Y-%m-%d"))

state_data <- list()

for (state_i in seq(no_states)){
  state_code <- states$V1[state_i]
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
  
  state_df <- subset(covid_df, state == state_code)
  state_df <- state_df[order(state_df$date),]
  state_df <- state_df[!is.na(state_df$hospitalizedCurrently), ]
  state_df <- state_df[state_df$date <= 365, ]
  state_hos <- subset(state_df, select = c("date", "hospitalizedCurrently"))
  
  state_data[[state_code]] <- list(idx=state_i, pop=census_pop, var=varob, epi=state_hos)
}


prms_low = c(rep(1e-5, no_states), 0.01, 0.3, 0.8, param_bounds[[R0_mod]]$low)
prms_high = c(rep(0.1, no_states), 0.2, 3, 8, param_bounds[[R0_mod]]$high)

get_state_prms <- function(prms_vec, state_idx){
  return(c(prms_vec[state_idx], prms_vec[-1:-no_states]))
}

negLL_wrapper <- function(state_elem, p_vec){
  neg_LL <- pois_L(get_state_prms(p_vec, state_elem$idx), 
                     paste0("R0_", R0_mod),
                     state_elem$var,
                     state_elem$epi,
                     state_elem$pop)
  return(neg_LL)
}

all_state_negLL <- function(prms){ # everything else read as global variable
  states_negLL <- sapply(X=state_data, FUN=negLL_wrapper, p_vec=prms)
  return(sum(states_negLL))
}

cat(">> Fitting with R0 model:", R0_mod, "\nLower:", prms_low, "\nUpper:", prms_high, "\n")
if (optMethod == "DEoptim"){
  cat(detectCores(), "cores seen; limit to", nthr, "cores\n")
  clstr <- makeCluster(nthr, type = "FORK") # limits the number of cores to use
  oo <- DEoptim(all_state_negLL, lower=prms_low, upper=prms_high,
                       control=DEoptim.control(trace = 5, reltol = 1e-6, itermax = 2000, steptol = 200,
                                               cluster=clstr, packages=c("deSolve"),
                                               parVar=c(var_cache, "state_data", "negLL_wrapper", "get_state_prms", "R0_mod")))
  stopCluster(clstr)
} else if (optMethod == "pso") {
  oo <- psoptim(rep(NA,length(prms_low)), all_state_negLL, lower=prms_low, upper=prms_high,
                       control=list(trace=1, REPORT=5, maxit=5000, trace.stats=TRUE, maxit.stagnate=500))
}

cat(states$V1); cat("; hrate; R0_min; R0_range; [R0_prms]\n")
cat(oo$par)

saveRDS(oo, file = paste0("fit_results/", handle, ".rds"))
sink()
