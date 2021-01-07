#!/usr/bin/env Rscript

library(deSolve)
source("SIRS_model.R")

args <- commandArgs(trailingOnly=TRUE)
state_ls <- "states_QC200.tsv"      # text file with a list of states to fit to
R0_mod <- "hs2d2"        # R0 model
R0_params <- as.numeric(strsplit(args[1], ":", fixed=TRUE)[[1]]) # do not include `c`
l_pnl <- args[2]         # lambda (for bookkeeping only)

handle <- paste0(state_ls, l_pnl, R0_mod, "_sims")

states <- read.delim(state_ls, header = FALSE)$V1

## Load all state data
state_pop <- read.delim("state_lv_data/state_pop.tsv")
all_state_sun <- read.csv("state_lv_data/state_daily_sunrise_2019.csv")
all_state_day <- read.csv("state_lv_data/state_daytime_2019.csv")
all_state_hum <- read.csv("state_lv_data/state_humidity_2014_2018.csv")
all_state_hum <- all_state_hum[all_state_hum$X != "2016-02-29", ] # get rid of leap year Feb 29

# DST mask
DST <- rep(0, times=365)
DST[69:306] <- 60/720 # for DST change in 2019

permDST <- rep(60/720, times=365)
permDST[69:306] <- 0

infection_dst <- data.frame(matrix(ncol = length(states), nrow = 364))
infection_nodst <- data.frame(matrix(ncol = length(states), nrow = 364))
infection_permdst <- data.frame(matrix(ncol = length(states), nrow = 364))
colnames(infection_dst) <- states
colnames(infection_nodst) <- states
colnames(infection_permdst) <- states

annual_cnt <- data.frame(wdst=numeric(length(states)), wodst=numeric(length(states)), permdst=numeric(length(states)), row.names = states)

sim_SIRS <- function(N_pop, annual_R0){
  
  xstart <- c(S = N_pop-1, I = 1, R = 0) # use actual state population
  # some hard-coded parameters
  paras = list(D = 5, 
               L = 40*7, 
               R0 = rep(head(annual_R0, 364), times = tot_years)[(offset+1):(364*tot_years)]) # pre-calculated R0 values
  predictions <- as.data.frame(ode(xstart, times, SIRS_R0, paras)) # run SIRS model
  
  ys <- (tot_years-1)*364-offset+1 # year start (last year)
  ye <- tot_years*364-offset       # year end
  
  annual_inf1 <- sum(predictions$I[ys:ye] - 0.8*predictions$I[(ys-1):(ye-1)]) # for D=5
  annual_inf2 <- sum(head(annual_R0, 364)/5*predictions$S[ys:ye]*predictions$I[ys:ye]/census_pop)
  cat(annual_inf1, annual_inf2, "\n") # sanity check
  
  return(list(inf_traj=predictions$I[ys:ye], tot_inf=annual_inf1))
}

for (state_i in seq(length(states))){
  state_code <- states[state_i]
  cat(">>> Simulating:", state_code, "<<<\n")
  census_pop <- state_pop$pop[state_pop$code==state_code]
  
  sunob <- all_state_sun[[state_code]]/720 # 365 days, scaled to maximum 720 = 12 hours
  dayob <- all_state_day[[state_code]]/1440

  climob <- all_state_hum[[state_code]] # multiple years
  climob <- matrix(climob, nrow = length(climob)/365, ncol = 365, byrow = TRUE)
  climob <- colMeans(climob) # 365 days
  
  DST_R0 <- do.call(paste0("R0_", R0_mod), c(list(climob, sunob, dayob), as.list(R0_params)))
  DST_sim <- sim_SIRS(census_pop, DST_R0)
  
  infection_dst[[state_code]] <- DST_sim$inf_traj
  annual_cnt[state_code, "wdst"] <- DST_sim$tot_inf
  
  if (state_code %in% c("AZ", "HI")){
    sansDST_R0 <- DST_R0
    permDST_R0 <- do.call(paste0("R0_", R0_mod), c(list(climob, sunob+60/720, dayob), as.list(R0_params)))
  } else {
    sansDST_R0 <- do.call(paste0("R0_", R0_mod), c(list(climob, sunob-DST, dayob), as.list(R0_params)))
    permDST_R0 <- do.call(paste0("R0_", R0_mod), c(list(climob, sunob+permDST, dayob), as.list(R0_params)))
  }
  sansDST_sim <- sim_SIRS(census_pop, sansDST_R0)
  infection_nodst[[state_code]] <- sansDST_sim$inf_traj
  annual_cnt[state_code, "wodst"] <- sansDST_sim$tot_inf
  
  permDST_sim <- sim_SIRS(census_pop, permDST_R0)
  infection_permdst[[state_code]] <- permDST_sim$inf_traj
  annual_cnt[state_code, "permdst"] <- permDST_sim$tot_inf
}

write.table(infection_dst, file = paste0("fit_results/", handle, "_dst.tsv"), quote=FALSE, sep="\t", row.names=FALSE)
write.table(infection_nodst, file = paste0("fit_results/", handle, "_sansdst.tsv"), quote=FALSE, sep="\t", row.names=FALSE)
write.table(infection_permdst, file = paste0("fit_results/", handle, "_permdst.tsv"), quote=FALSE, sep="\t", row.names=FALSE)
write.table(annual_cnt, file = paste0("fit_results/", handle, "_cnt.tsv"), quote=FALSE, sep="\t", row.names=TRUE)