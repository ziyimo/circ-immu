#!/usr/bin/env Rscript

library(deSolve)
library(MASS)
source("R0_mods.R")
source("SEIH_mod.R")

args <- commandArgs(trailingOnly=TRUE)
hospr <- args[1]
#var_scale <- as.numeric(args[2]) # deprecated

states <- read.delim("states_49DC.tsv", header = FALSE)$V1
no_states <- length(states)
## Load all state data
state_pop <- read.delim("state_lv_data/state_pop.tsv")
all_state_sun <- read.csv("state_lv_data/state_daily_sunrise_2019.csv")
all_state_day <- read.csv("state_lv_data/state_daytime_2019.csv")
# all_state_hum <- read.csv("state_lv_data/state_humidity_2014_2018.csv")
# all_state_hum <- all_state_hum[all_state_hum$X != "2016-02-29", ] # get rid of leap year Feb 29

state2reg <- read.csv("state_lv_data/state2censusReg.csv")
state2reg <- subset(state2reg, select=-c(State))
CR <- c("Northeast", "South", "Midwest", "West")
cdv <- unique(state2reg$Division)
no_cdv <- length(cdv)
state2reg$CDV.Code <- match(state2reg$Division, cdv)

#param <- as.numeric(strsplit(param, " ", fixed=TRUE)[[1]])
load(paste0("covid_hosp_fit/eta", hospr, "_sd_iter.rds")) # ss_fit and seed_sh
resamp <- readRDS(paste0("covid_hosp_fit/sd.boot.eta", hospr, "_seeded.rds")) # load from bootstrapped parameter values
hospr <- as.numeric(hospr)

# cov_Mtx <- readRDS(paste0("covid_hosp_fit/sd", hospr, "_1g5e-3_covMtx.rds")) # read covariant matrix
# resamp <- mvrnorm(n=10, seed_sh, cov_Mtx*var_scale) # 10 replicates for now

prms_sets <- rbind(seed_sh, resamp)
no_reps <- nrow(prms_sets)

sim_SEIH <- function(prms, R0_func, var_ls, census_pop){
  prop_init <- prms[1]
  #hpt_rate <- hospr
  R0_min <- prms[2]
  R0_range <- prms[3]
  R0_prms <- prms[-1:-3]
  
  xstart <- c(S = census_pop*(1-prop_init), E = 0, I = census_pop*prop_init, H = 0, R = 0) # use actual state population
  times <- seq(73, 365+31) # from national emergency declaration into second half of 2021
  annual_R0 <- do.call(R0_func, c(var_ls, as.list(c(R0_prms, R0_min+R0_range, R0_min))))
  
  # some hard-coded parameters
  paras = list(sigma = 5,
               h = hospr, # defined as global variable!
               lambda = 5,
               gam = 5,
               k = 10,
               R0 = rep(annual_R0, times=2)) # pre-calculated R0 values
  
  predictions <- as.data.frame(ode(xstart, times, SEIH_R0, paras)) # run SIRS model
  
  annual_H1 <- sum(predictions$H[2:(396-72)] - 0.9*predictions$H[1:(396-73)]) # for k = 10
  annual_H2 <- sum(hospr/5*predictions$I[1:(396-72)]) # for lambda = 5
  ye_recvd <- predictions$R[396-72]/census_pop # proportion of recovered year end
  cat(annual_H1, annual_H2, ye_recvd, "\n") # sanity check
  
  return(list(hosp_traj=predictions$H, tot_hosp=annual_H1, year_end_Rem=ye_recvd))
}

####  CR level aggregation + counterfactual sim  ####

# DST mask
DST <- rep(0, times=365)
DST[69:306] <- 60/720 # for DST change in 2019

permDST <- rep(60/720, times=365)
permDST[69:306] <- 0

dst_trajs <- list()
nodst_trajs <- list()
permdst_trajs <- list()

for (item in states){ # record state lv. traj.
  dst_trajs[[item]] <- matrix(0, ncol = no_reps, nrow = length(seq(73, 365+31)))
  nodst_trajs[[item]] <- matrix(0, ncol = no_reps, nrow = length(seq(73, 365+31)))
  permdst_trajs[[item]] <- matrix(0, ncol = no_reps, nrow = length(seq(73, 365+31)))
}

nodst_rc <- data.frame(matrix(ncol = no_reps, nrow = no_states+1), row.names = c(states, "ALL"))
permdst_rc <- data.frame(matrix(ncol = no_reps, nrow = no_states+1), row.names = c(states, "ALL"))
year_end_R <- data.frame(matrix(ncol = no_reps, nrow = no_states), row.names = states) # proportion recovered

for (samp_i in seq(no_reps)){

  #annual_cnt <- data.frame(wdst=numeric(length(CR)), wodst=numeric(length(CR)), permdst=numeric(length(CR)), row.names = CR)
  annual_cnt <- data.frame(wdst=numeric(no_states), wodst=numeric(no_states), permdst=numeric(no_states), row.names = states)
  
  sh_params <- prms_sets[samp_i, ]
  for (state_i in seq(no_states)){
    state_eg <- states[state_i]
    census_pop <- state_pop$pop[state_pop$code==state_eg]
    sunob <- all_state_sun[[state_eg]]/720
    dayob <- all_state_day[[state_eg]]/1440
    # climob <- all_state_hum[[state_eg]] # multiple years
    # climob <- matrix(climob, nrow = length(climob)/365, ncol = 365, byrow = TRUE)
    # climob <- colMeans(climob) # 365 days
    
    s0_idx <- 2
    state_prms <- c(ss_fit[[state_eg]],
                    sh_params[(no_cdv+1):(no_cdv+2+s0_idx-1)],
                    sh_params[state2reg$CDV.Code[state2reg$State.Code == state_eg]],
                    sh_params[-1:-(no_cdv+2+s0_idx-1)])
    
    cat(state_eg, state_prms, "\n")
    
    state_CR <- state2reg$Region[state2reg$State.Code == state_eg]
    DST_sim <- sim_SEIH(state_prms, "R0_sd", list(sunob, dayob), census_pop)
    
    #dst_trajs[[state_CR]][, samp_i] <- dst_trajs[[state_CR]][, samp_i] + DST_sim$hosp_traj
    #annual_cnt[state_CR, "wdst"] <- annual_cnt[state_CR, "wdst"] + DST_sim$tot_hosp
    dst_trajs[[state_eg]][, samp_i] <- DST_sim$hosp_traj
    annual_cnt[state_eg, "wdst"] <- DST_sim$tot_hosp
    year_end_R[state_i, samp_i] <- DST_sim$year_end_Rem
    
    if (state_eg %in% c("AZ", "HI")){
      sansDST_sim <- DST_sim
      permDST_sim <- sim_SEIH(state_prms, "R0_sd", list(sunob+60/720, dayob), census_pop)
    } else {
      sansDST_sim <- sim_SEIH(state_prms, "R0_sd", list(sunob-DST, dayob), census_pop)
      permDST_sim <- sim_SEIH(state_prms, "R0_sd", list(sunob+permDST, dayob), census_pop)
    }
    #nodst_trajs[[state_CR]][, samp_i] <- nodst_trajs[[state_CR]][, samp_i] + sansDST_sim$hosp_traj
    #annual_cnt[state_CR, "wodst"] <- annual_cnt[state_CR, "wodst"] + sansDST_sim$tot_hosp
    nodst_trajs[[state_eg]][, samp_i] <- sansDST_sim$hosp_traj
    annual_cnt[state_eg, "wodst"] <-sansDST_sim$tot_hosp
    
    #permdst_trajs[[state_CR]][, samp_i] <- permdst_trajs[[state_CR]][, samp_i] + permDST_sim$hosp_traj
    #annual_cnt[state_CR, "permdst"] <- annual_cnt[state_CR, "permdst"] + permDST_sim$tot_hosp
    permdst_trajs[[state_eg]][, samp_i] <- permDST_sim$hosp_traj
    annual_cnt[state_eg, "permdst"] <- permDST_sim$tot_hosp
  }
  
  nodst_rc[1:no_states, samp_i] <- (annual_cnt$wodst - annual_cnt$wdst)/annual_cnt$wdst
  permdst_rc[1:no_states, samp_i] <- (annual_cnt$permdst - annual_cnt$wdst)/annual_cnt$wdst
  
  nodst_rc[(no_states+1), samp_i] <- (sum(annual_cnt$wodst) - sum(annual_cnt$wdst))/sum(annual_cnt$wdst)
  permdst_rc[(no_states+1), samp_i] <- (sum(annual_cnt$permdst) - sum(annual_cnt$wdst))/sum(annual_cnt$wdst)
  cat(">> Rep", samp_i, ":", nodst_rc[(no_states+1), samp_i], permdst_rc[(no_states+1), samp_i], "\n")
}

time_stamp <- format(Sys.time(), "%m%d%H")
#saveRDS(prms_sets, file = paste0("covid_hosp_fit/eta", hospr, "_var", var_scale,  "_", time_stamp, "_prmsamps.rds")) # deprecated for bootstrap method
save(dst_trajs, nodst_trajs, permdst_trajs, nodst_rc, permdst_rc, year_end_R,
     file = paste0("covid_hosp_fit/eta", hospr, "_bootsims_", time_stamp, ".RDa"))
