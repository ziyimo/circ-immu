#!/usr/bin/env Rscript

#.libPaths() # sanity check
library("parallel")

source("R0_mods.R")
source("COVID_hosp_analysis/SEIH_mod.R")
var_cache <- ls()

time_stamp <- format(Sys.time(), "%m%d%H")
args <- commandArgs(trailingOnly=TRUE)

state_ls <- args[1]      # text file with a list of states to fit to
R0_mod <- args[2]        # R0 model, options: [see code]
prms_init <- as.numeric(strsplit(args[3], ":", fixed=TRUE)[[1]])
# xprop_initx, incb_prd, inf_prd, hpt_rate, hpt_prd, R0min, R0range, [R0_params]
macro_mask <- as.logical(strsplit(args[4], ":", fixed=TRUE)[[1]]) # a mask of TUNED parameters e.g. "T:F:F:T:F:T:T"
nthr <- as.numeric(args[5]) # limit no. of threads
optimizer <- "pso"

library(optimizer, character.only = TRUE)
handle <- paste0(c(time_stamp, R0_mod, parameter_names[macro_mask]), collapse = ".")

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
seed_ss <- list()
state2reg <- read.csv("compiled_data/geo_clim_dem/state2censusReg.csv")
cdv <- unique(state2reg$Division)
no_cdv <- length(cdv)
state2reg$CDV.Code <- match(state2reg$Division, cdv)

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
  state_df <- state_df[state_df$date <= 396, ] # take data till end of Jan 2021
  state_hos <- subset(state_df, select = c("date", "hospitalizedCurrently"))
  
  state_data[[state_code]] <- list(idx=state_i, pop=census_pop, var=varob, epi=state_hos, cdv=state2reg$CDV.Code[state2reg$State.Code == state_code])
}

# State specific params - prop_init
ss_low <- c(1e-5)
ss_high <- c(0.1)

# shared, non-R0 model specific params - incb_prd, inf_prd, hpt_rate, hpt_prd, R0min, R0range
omnibus_low <- c(1, 1, 0.003, 1, 0.8, 0.2)
omnibus_high <- c(15, 15, 0.1, 15, 1.5, 2.5)
no_omnibus <- sum(macro_mask[-1])

if (R0_mod %in% c("sun", "sd", "hsd")){
  s0_idx <- ifelse(R0_mod %in% c("sun", "sd"), 2, 3)
  # Shared params: [CD s_0], [omnibus params], [other R0 params]
  sh_low <- c(rep(0.4, no_cdv), omnibus_low[macro_mask[-1]], param_bounds[[R0_mod]]$low[-s0_idx])
  sh_high <- c(rep(0.7, no_cdv), omnibus_high[macro_mask[-1]], param_bounds[[R0_mod]]$high[-s0_idx])
  # constrain the range of s_0
  
  seed_sh <- c(rep(prms_init[6+s0_idx], no_cdv), prms_init[1:6][macro_mask[-1]], prms_init[-1:-6][-s0_idx])
} else {
  # Shared params: [omnibus params], [R0 params]
  sh_low <- c(omnibus_low[macro_mask[-1]], param_bounds[[R0_mod]]$low)
  sh_high <- c(omnibus_high[macro_mask[-1]], param_bounds[[R0_mod]]$high)
  
  seed_sh <- c(prms_init[1:6][macro_mask[-1]], prms_init[-1:-6])
}

cat(">> Seeded prms: [R0 prms] +",  parameter_names[-1][macro_mask[-1]], "\n>>", seed_sh, "\n")

macro_fixed <- prms_init[1:6][!macro_mask[-1]] # fixed parameters
cat(">> Fixed prms:",  parameter_names[!macro_mask], "\n>>", macro_fixed, "\n")

state_wrapper <- function(state_elem, ss_prms, sh_prms){
  #state_elem: element in the state_data list
  #ss_prms: vector of state specific parameters
  #sh_prms: vector of shared parameters 
  
  if (R0_mod %in% c("sun", "sd", "hsd")){
    state_prms <- c(ss_prms,
                    sh_prms[(no_cdv+1):(no_cdv+no_omnibus+s0_idx-1)],
                    sh_prms[state_elem$cdv],
                    sh_prms[-1:-(no_cdv+no_omnibus+s0_idx-1)])
  } else {
    state_prms <- c(ss_prms, sh_prms)
  }
  
  neg_LL <- norm_L(state_prms,
                   macro_fixed, # fixed parameters, defined globally
                   macro_mask,  # mask of tuned parameters, defined globally
                   paste0("R0_", R0_mod),
                   state_elem$var,
                   state_elem$epi,
                   state_elem$pop)
  return(neg_LL)
}

## Functions for fitting state specific parameters independently

fit_state <- function(state_dat, sh_vec){ # no seeding for 1D optimization
  oo <- optimize(ss_wrapper, c(0, 1), st_e=state_dat, sh_p=sh_vec)
  
  return(ss_low + oo$minimum*(ss_high-ss_low))
}

ss_wrapper <- function(ss_norm, st_e, sh_p){
  ss_p <- ss_low + ss_norm*(ss_high-ss_low)
  return(state_wrapper(st_e, ss_p, sh_p))
}


## Function for fitting shared parameters
all_state_negLL <- function(norm_sh_prms, ss_prms){
  orig_prms <- sh_low + norm_sh_prms*(sh_high-sh_low)
  states_negLL <- parSapply(cl = clstr, X=state_data, FUN=sh_wrapper, sh_orig=orig_prms, ss_ls=ss_prms)
  return(sum(states_negLL))
}

sh_wrapper <- function(st_e, sh_orig, ss_ls){
  ss_p <- ss_ls[[st_e$idx]]
  return(state_wrapper(st_e, ss_p, sh_orig))
}


cat(detectCores(), "cores seen; limit to", nthr, "cores\n")
clstr <- makeCluster(nthr, type = "FORK") # limits the number of cores to use
for (iter in seq(20)){
  cat(">> Iteration:", iter, "\n")
  ss_fit <- parLapply(cl = clstr, X = state_data, fun = fit_state, sh_vec = seed_sh)
  print(ss_fit)

  sh_fit <- psoptim((seed_sh-sh_low)/(sh_high-sh_low), all_state_negLL, ss_prms=ss_fit,
                    lower=rep(0,length(sh_low)), upper=rep(1,length(sh_low)),
                    control=list(trace=1, REPORT=5, maxit=10000, trace.stats=FALSE, maxit.stagnate=200))
  seed_sh <- sh_low + sh_fit$par*(sh_high-sh_low)

  cat(">>", R0_mod, "Iter", iter, "negLL:", sh_fit$value, "\n")
  cat(">> Fitted prms:",  parameter_names[macro_mask], "[R0_params] \n")
  cat(">> ", unlist(ss_fit), "\n>> ", seed_sh, "\n")
}

save(ss_fit, seed_sh, file = paste0("fit_results/", handle, "_iterative.rds"))
stopCluster(clstr)
