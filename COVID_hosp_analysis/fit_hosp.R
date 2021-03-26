#!/usr/bin/env Rscript

#.libPaths() # sanity check
library("parallel")

source("R0_mods.R")
source("COVID_hosp_analysis/SEIH_mod.R")
var_cache <- ls()

time_stamp <- format(Sys.time(), "%m%d")
args <- commandArgs(trailingOnly=TRUE)

state_ls <- args[1]      # text file with a list of states to fit to
R0_mod <- args[2]        # R0 model, options: [see code]
macro_fixed <- as.numeric(strsplit(args[3], ":", fixed=TRUE)[[1]]) # only the fixed parameters!
macro_mask <- as.logical(strsplit(args[4], ":", fixed=TRUE)[[1]]) # a mask of TUNED parameters e.g. "T:F:F:T:F:T:T"
nthr <- as.numeric(args[5]) # limit no. of threads
optimizer <- "pso"

library(optimizer, character.only = TRUE)
handle <- paste0(state_ls, "_", R0_mod, "_", time_stamp, "_", optimizer)

states <- read.delim(state_ls, header = FALSE)
no_states <- nrow(states)

## Load all state data
state_pop <- read.delim("compiled_data/geo_clim_dem/state_pop.tsv")
all_state_sun <- read.csv("compiled_data/geo_clim_dem/state_daily_sunrise_2019.csv")
all_state_day <- read.csv("compiled_data/geo_clim_dem/state_daytime_2019.csv")
all_state_hum <- read.csv("compiled_data/geo_clim_dem/state_humidity_2014_2018.csv")
all_state_hum <- all_state_hum[all_state_hum$X != "2016-02-29", ] # get rid of leap year Feb 29
## Load COVID hospitalization data
covid_df <- readRDS("compiled_data/COVID_data/state_hospitalization.rds")
covid_df$date <- as.numeric(as.Date(covid_df$date, format="%Y-%m-%d") - as.Date("2019-12-31", format="%Y-%m-%d"))

state_data <- list()

for (state_i in seq(no_states)){
  state_code <- states$V1[state_i]
  #cat(">>> Loading state data:", state_code, "<<<\n")
  census_pop <- state_pop$pop[state_pop$code==state_code]
  
  if (R0_mod %in% c("sd", "hsd")){
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
  }
  
  state_df <- subset(covid_df, state == state_code)
  state_df <- state_df[order(state_df$date),]
  state_df <- state_df[!is.na(state_df$hospitalizedCurrently), ]
  state_df <- state_df[state_df$date <= 396, ] # take data till end of Jan 2021
  state_hos <- subset(state_df, select = c("date", "hospitalizedCurrently"))
  
  state_data[[state_code]] <- list(idx=state_i, pop=census_pop, var=varob, epi=state_hos)
}

cat("###Error func: Gaussian###\n")
negLL_wrapper <- function(state_elem, p_vec){
  neg_LL <- norm_L(get_state_prms(p_vec, state_elem$idx),
                   macro_fixed,
                   macro_mask,
                   paste0("R0_", R0_mod),
                   state_elem$var,
                   state_elem$epi,
                   state_elem$pop)
  return(neg_LL)
}

all_state_negLL <- function(norm_prms){ # everything else read as global variable
  orig_prms <- prms_low + norm_prms*(prms_high-prms_low)
  states_negLL <- parSapply(cl = clstr, X=state_data, FUN=negLL_wrapper, p_vec=orig_prms)
  return(sum(states_negLL))
}

## ROUND 1: state-specific I_init only
# shared, non-R0 model specific params - incb_prd, inf_prd, hpt_rate, hpt_prd, R0min, R0range
omnibus_low <- c(1, 1, 0.003, 1, 0.5, 0.2)
omnibus_high <- c(15, 15, 0.1, 15, 1.5, 2.5)

# I_init; [omnibus_prms]; [R0_prms]
prms_low <- c(rep(1e-5, no_states), omnibus_low[macro_mask[-1]], param_bounds[[R0_mod]]$low)
prms_high <- c(rep(0.1, no_states), omnibus_high[macro_mask[-1]], param_bounds[[R0_mod]]$high)

get_state_prms <- function(prms_vec, state_idx){
  return(c(prms_vec[state_idx], prms_vec[-1:-no_states]))
}

cat(detectCores(), "cores seen; limit to", nthr, "cores\n")
clstr <- makeCluster(nthr, type = "FORK") # limits the number of cores to use

cat(">> Fitting with R0 model:", R0_mod, "\n")
cat(">> Fixed params - ", parameter_names[!macro_mask], ":", macro_fixed, "\n")
cat(">> Fitted params - ", parameter_names[macro_mask], " + [R0 params]\n")
cat(">> Lower:", prms_low, "\n>> Upper:", prms_high, "\n")
restart_thd <- c(1.1e-2, 1.1e-3, 1.1e-4)

seed <- rep(NA,length(prms_low))
for (thrhld in restart_thd){
  cat(">> Restart threshold:", thrhld, "\n")
  oo <- psoptim(seed, all_state_negLL, lower=rep(0,length(prms_low)), upper=rep(1,length(prms_low)),
                control=list(trace=1, REPORT=5, maxit=10000, trace.stats=FALSE, maxit.stagnate=500,
                             max.restart=5, reltol=thrhld))
  seed <- oo$par
  cat(">>", R0_mod, "RND1:", oo$value, "\n")
  #cat(states$V1); cat("; hrate; R0_min; R0_range; [R0_prms]\n")
  cat(">> Fitted params - ", parameter_names[macro_mask], " + [R0 params]\n")
  cat(">> RND1:", prms_low + oo$par*(prms_high-prms_low), "\n")
}
stopCluster(clstr)
saveRDS(oo, file = paste0("fit_results/", handle, "_RND1.rds"))
