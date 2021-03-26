#!/usr/bin/env Rscript

#.libPaths() # sanity check
library("parallel")

source("R0_mods.R")
source("SEIH_mod.R")

args <- commandArgs(trailingOnly=TRUE)

state_ls <- args[1]      # text file with a list of states to fit to
R0_mod <- args[2]        # R0 model, options: [see code]
hospr <- args[3]
discret <- args[4] # step for gradient approximation
grad2hess <- as.numeric(args[5]) # this is the parscale parameter in optimHess
nthr <- as.numeric(args[6]) # limit no. of threads

handle <- paste0("covid_hosp_fit/eta", hospr, "_", R0_mod, "_iter.rds")
hospr <- as.numeric(hospr)

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

load(handle) # ss_fit and seed_sh

macro_mask <- c(TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE)
no_omnibus <- sum(macro_mask[-1])
macro_fixed <- c(5, 5, hospr, 10)
cat(">> Fixed prms:",  parameter_names[!macro_mask], "\n>>", macro_fixed, "\n")

# shared, non-R0 model specific params - incb_prd, inf_prd, hpt_rate, hpt_prd, R0min, R0range
omnibus_low <- c(1, 1, 0.003, 1, 0.8, 0.2)
omnibus_high <- c(15, 15, 0.1, 15, 1.5, 1)

if (R0_mod %in% c("sun", "sd", "hsd")){
  s0_idx <- ifelse(R0_mod %in% c("sun", "sd"), 2, 3)
  # Shared params: [CD s_0], [omnibus params], [other R0 params]
  sh_low <- c(rep(0.4, no_cdv), omnibus_low[macro_mask[-1]], param_bounds[[R0_mod]]$low[-s0_idx])
  sh_high <- c(rep(0.7, no_cdv), omnibus_high[macro_mask[-1]], param_bounds[[R0_mod]]$high[-s0_idx])
  # constrain the range of s_0
} else {
  # Shared params: [omnibus params], [R0 params]
  sh_low <- c(omnibus_low[macro_mask[-1]], param_bounds[[R0_mod]]$low)
  sh_high <- c(omnibus_high[macro_mask[-1]], param_bounds[[R0_mod]]$high)
}

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
all_state_LL <- function(orig_sh_prms){
  states_negLL <- parSapply(cl = clstr, X=state_data, FUN=sh_wrapper,
                            sh_orig=orig_sh_prms, ss_ls=ss_fit)
  return(-sum(states_negLL))
}

sh_wrapper <- function(st_e, sh_orig, ss_ls){
  ss_p <- ss_ls[[st_e$idx]]
  return(state_wrapper(st_e, ss_p, sh_orig))
}

cat(detectCores(), "cores seen; limit to", nthr, "cores\n")
clstr <- makeCluster(nthr, type = "FORK") # limits the number of cores to use
# check function value at optimum
opt_val <- all_state_LL(seed_sh)
cat(handle, ":", opt_val, "\n")


prms_range <- sh_high-sh_low
hessian_step <- prms_range*as.numeric(discret)
cat(">>", hessian_step, "\n")
cat(">> grad2hess ratio:", grad2hess, "\n")
info_mtx <- -optimHess(seed_sh, all_state_LL,
                       control = list(parscale=rep(grad2hess, length(seed_sh)), ndeps=hessian_step))
# set gradient step relative the hessian step ^^, see https://docs.tibco.com/pub/enterprise-runtime-for-R/4.0.2/doc/html/Language_Reference/stats/optimHess.html

cat(info_mtx, "\n")
cov_mtx <- solve(info_mtx)
cat(">>", R0_mod, hospr, ";", grad2hess, discret, ":", diag(cov_mtx), "\n")

saveRDS(cov_mtx, file=paste0("covid_hosp_fit/", R0_mod, hospr, "_", grad2hess, "g", discret, "_covMtx.rds"))
stopCluster(clstr)