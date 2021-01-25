#!/usr/bin/env Rscript

#.libPaths() # sanity check
library("parallel")

source("R0_mods.R")
source("SEIH_mod.R")
var_cache <- ls()

time_stamp <- format(Sys.time(), "%m%d%H")
args <- commandArgs(trailingOnly=TRUE)

R0_mod <- args[1]        # R0 model, options: cos, hum, sun, hs
seed_prm <- as.numeric(strsplit(args[2], ":", fixed=TRUE)[[1]]) # seed, put "0:0" for default
nthr <- 4
optimizer <- "pso"

library(optimizer, character.only = TRUE)
handle <- paste0("CR_", R0_mod, "_", time_stamp, "_", optimizer)

## Load all state data
state_pop <- read.delim("state_lv_data/cr_pop.tsv")
all_state_sun <- read.csv("state_lv_data/census_reg_daily_sunrise_2019.csv")
all_state_day <- read.csv("state_lv_data/census_reg_daytime_2019.csv")
all_state_hum <- read.csv("state_lv_data/census_reg_humidity_2014_2018.csv")
all_state_hum <- all_state_hum[all_state_hum$date != "2016-02-29", ] # get rid of leap year Feb 29
## Load COVID hospitalization data
covid_df <- readRDS("state_lv_data/censusReg_hospitalization.rds")

state_data <- list()
no_states <- nrow(state_pop)

for (state_i in seq(no_states)){
  state_code <- state_pop$Region[state_i]
  cat(">>> Loading state data:", state_code, "<<<\n")
  census_pop <- state_pop$pop[state_i]
  
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
  
  state_data[[state_code]] <- list(idx=state_i, pop=census_pop, var=varob, epi=state_hos)
}

prms_low <- c(rep(1e-5, no_states), rep(0.003, no_states), 0.5, 0.2, param_bounds[[R0_mod]]$low)
prms_high <- c(rep(0.05, no_states), rep(0.13, no_states), 1.6, 2.5, param_bounds[[R0_mod]]$high)

get_state_prms <- function(prms_vec, state_idx){
  return(c(prms_vec[state_idx], prms_vec[no_states+state_idx], prms_vec[-1:-(2*no_states)]))
}

negLL_wrapper <- function(state_elem, p_vec){
  neg_LL <- norm_L(get_state_prms(p_vec, state_elem$idx), 
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

# fit <- readRDS("fit_results/states_trial8.tsv_day_0114_pso.rds")
# prms <- fit$par
# nthr <- 4

restart_thd <- c(1.1e-1, 1.1e-2, 1.1e-3, 1.1e-4)
cat(detectCores(), "cores seen; limit to", nthr, "cores\n")
clstr <- makeCluster(nthr, type = "FORK") # limits the number of cores to use
# system.time(states_negLLpar <- parSapply(cl = clstr, X=state_data, FUN=negLL_wrapper, p_vec=prms))
# system.time(states_negLLser <- sapply(X=state_data, FUN=negLL_wrapper, p_vec=prms))

if (sum(seed_prm) == 0){
  seed <- rep(NA,length(prms_low))
} else {
  seed <- (seed_prm-prms_low)/(prms_high-prms_low)
}

cat(">> Fitting with R0 model:", R0_mod, "\nLower:", prms_low, "\nUpper:", prms_high, "\nSeed:", prms_low + seed*(prms_high-prms_low), "\n")
for (thrhld in restart_thd){
  cat(">> Restart threshold:", thrhld, "\n")
  oo <- psoptim(seed, all_state_negLL, lower=rep(0,length(prms_low)), upper=rep(1,length(prms_low)),
                control=list(trace=1, REPORT=5, maxit=10000, trace.stats=FALSE, maxit.stagnate=500,
                             max.restart=5, reltol=thrhld))
  seed <- oo$par
  cat(state_pop$Region); cat("; hrate; R0_min; R0_range; [R0_prms]\n")
  cat(prms_low + oo$par*(prms_high-prms_low))
  cat("\n")
}
stopCluster(clstr)


if (R0_mod %in% c("sun", "sd", "hsd")){
  s0_idx <- ifelse(R0_mod %in% c("sun", "sd"), 2, 3)
  prms_low <- c(prms_low[1:(2*no_states)], rep(prms_low[2*no_states+2+s0_idx], no_states),
                prms_low[(2*no_states+1):(2*no_states+2)], prms_low[-1:(-2*no_states-2)][-s0_idx])
  prms_high <- c(prms_high[1:(2*no_states)], rep(prms_high[2*no_states+2+s0_idx], no_states),
                prms_high[(2*no_states+1):(2*no_states+2)], prms_high[-1:(-2*no_states-2)][-s0_idx])
  #if (sum(seed_prm) == 0){
  seed <- c(seed[1:(2*no_states)], rep(seed[2*no_states+2+s0_idx], no_states),
            seed[(2*no_states+1):(2*no_states+2)], seed[-1:(-2*no_states-2)][-s0_idx])
  #} else { seed <- (seed_prm-prms_low)/(prms_high-prms_low)}
  
  get_state_prms <- function(prms_vec, state_idx){
    return(c(prms_vec[state_idx], # state-specific seed_prop
             prms_vec[no_states+state_idx], # state-specific hosp rate
             prms_vec[(3*no_states+1):(3*no_states+1+s0_idx)], # shared: R0_min, R0_range, [R0 prms b4 s_0]
             prms_vec[2*no_states+state_idx], #)) # state_specific s_0
             prms_vec[-1:-(3*no_states+1+s0_idx)])) # shared: [R0 prms after s_0]
  }
  
  clstr <- makeCluster(nthr, type = "FORK") 
  cat(">> Refining R0 model:", R0_mod, "\nLower:", prms_low, "\nUpper:", prms_high, "\nSeed:", prms_low + seed*(prms_high-prms_low), "\n")
  
  for (thrhld in restart_thd){
    cat(">> Restart threshold:", thrhld, "\n")
    oo <- psoptim(seed, all_state_negLL, lower=rep(0,length(prms_low)), upper=rep(1,length(prms_low)),
                  control=list(trace=1, REPORT=5, maxit=10000, trace.stats=FALSE, maxit.stagnate=500,
                               max.restart=5, reltol=thrhld))
    seed <- oo$par
    cat(state_pop$Region); cat("; hrate; R0_min; R0_range; [R0_prms]\n")
    cat(prms_low + oo$par*(prms_high-prms_low))
    cat("\n")
  }
  stopCluster(clstr)
}

saveRDS(oo, file = paste0("fit_results/", handle, ".rds"))
#sink()
