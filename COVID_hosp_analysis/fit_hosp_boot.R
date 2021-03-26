#!/usr/bin/env Rscript

#.libPaths() # sanity check
library("parallel")

source("R0_mods.R")
source("COVID_hosp_analysis/SEIH_mod.R")
var_cache <- ls()

args <- commandArgs(trailingOnly=TRUE)

state_ls <- args[1]      # text file with a list of states to fit to
best_fit <- args[2]      # path to the best fit parameters (to clamp I_init)
tag <- args[3]           # some tag to identify this bootstrap run
R0_mod <- "sd"
prms_init <- as.numeric(strsplit(args[4], ":", fixed=TRUE)[[1]])
# xprop_initx, incb_prd, inf_prd, hpt_rate, hpt_prd, R0min, R0range, [R0_params]
macro_mask <- as.logical(strsplit(args[5], ":", fixed=TRUE)[[1]]) # a mask of TUNED parameters e.g. "T:F:F:T:F:T:T"
nthr <- as.numeric(args[6]) # limit no. of threads
optimizer <- "pso"

library(optimizer, character.only = TRUE)
handle <- paste0(c(R0_mod, "boot", tag), collapse = ".")

states <- read.delim(state_ls, header = FALSE)
no_states <- nrow(states)

load(best_fit) # use 'ss_fit', ignore 'seed_sh'

## Load all state data
state_pop <- read.delim("compiled_data/geo_clim_dem/state_pop.tsv")
all_state_sun <- read.csv("compiled_data/geo_clim_dem/state_daily_sunrise_2019.csv")
all_state_day <- read.csv("compiled_data/geo_clim_dem/state_daytime_2019.csv")

## Load COVID hospitalization data
covid_df <- readRDS("compiled_data/COVID_data/state_hospitalization.rds")
covid_df$date <- as.numeric(as.Date(covid_df$date, format="%Y-%m-%d") - as.Date("2019-12-31", format="%Y-%m-%d"))

state_data <- list()
state2reg <- read.csv("compiled_data/geo_clim_dem/state2censusReg.csv")
cdv <- unique(state2reg$Division)
no_cdv <- length(cdv)
state2reg$CDV.Code <- match(state2reg$Division, cdv)

for (state_i in seq(no_states)){
  state_code <- states$V1[state_i]
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

# shared, non-R0 model specific params - incb_prd, inf_prd, hpt_rate, hpt_prd, R0min, R0range
omnibus_low <- c(1, 1, 0.003, 1, 0.8, 0.2)
omnibus_high <- c(15, 15, 0.1, 15, 1.5, 2.5)
no_omnibus <- sum(macro_mask[-1])

s0_idx <- 2
# Shared params: [CD s_0], [omnibus params], [other R0 params]
sh_low <- c(rep(0.4, no_cdv), omnibus_low[macro_mask[-1]], param_bounds[[R0_mod]]$low[-s0_idx])
sh_high <- c(rep(0.7, no_cdv), omnibus_high[macro_mask[-1]], param_bounds[[R0_mod]]$high[-s0_idx])
# constrain the range of s_0

seed_guide <- c(rep(0.5, no_cdv), seed_sh[-1:-no_cdv]) # seed all s_0 with 0.5

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

## Function for fitting shared parameters
all_state_negLL <- function(norm_sh_prms, ss_prms, boot_samps){
  orig_prms <- sh_low + norm_sh_prms*(sh_high-sh_low)
  states_negLL <- parSapply(cl = clstr, X=state_data, FUN=sh_wrapper, sh_orig=orig_prms, ss_ls=ss_prms)
  return(sum(states_negLL[boot_samps]))
}

sh_wrapper <- function(st_e, sh_orig, ss_ls){
  ss_p <- ss_ls[[st_e$idx]]
  return(state_wrapper(st_e, ss_p, sh_orig))
}

resamp_params <- matrix(nrow = 10, ncol = length(sh_low))

cat(detectCores(), "cores seen; limit to", nthr, "cores\n")
clstr <- makeCluster(nthr, type = "FORK") # limits the number of cores to use

for (iter in seq(10)){ # each run does 10 bootstrap resampling

  bssample <- sample(seq(no_states), size=no_states, replace = TRUE) # a bootstrap sample vector
  opt_negLL <- all_state_negLL((seed_sh-sh_low)/(sh_high-sh_low), ss_fit, bssample)
  cat(">>", R0_mod, "Resamp", iter, ":", bssample, "\n")
  cat(">> Opt param @ resamp:", opt_negLL, "\n")
  
  sh_fit <- psoptim((seed_guide-sh_low)/(sh_high-sh_low), all_state_negLL, ss_prms=ss_fit, boot_samps=bssample,
                    lower=rep(0,length(sh_low)), upper=rep(1,length(sh_low)),
                    control=list(trace=1, REPORT=5, maxit=10000, trace.stats=FALSE, maxit.stagnate=200))
  
  resamp_params[iter,] <- sh_low + sh_fit$par*(sh_high-sh_low)

  cat(">> Fitted prms:",  parameter_names[macro_mask], "[R0_params] \n")
  cat(">> Resamp_Iter", iter, "negLL:", sh_fit$value, "(ref: ", opt_negLL, ")\n")
  cat(">> Resamp_params:", resamp_params[iter,], "\n")
}

saveRDS(resamp_params, file = paste0("fit_results/", handle, ".rds"))
stopCluster(clstr)