#!/usr/bin/env Rscript

library(deSolve)
source("R0_mods.R")
source("SIRS_model.R")

args <- commandArgs(trailingOnly=TRUE)
state_ls <- "states_QC200.tsv"      # text file with a list of states to fit to
R0_mod <- "sd"        # R0 model
#R0_params <- as.numeric(strsplit(args[1], ":", fixed=TRUE)[[1]]) # do not include `c`
l_pnl <- args[1]         # lambda (for bookkeeping only)
mode <- args[2]

handle <- paste0(state_ls, l_pnl, R0_mod, "_", mode, "sims_trial")

fitted <- readRDS(Sys.glob(paste0("flu_fit/", paste(state_ls, l_pnl, R0_mod, "*.rds" ,sep="_"))))
prms_optim <- unname(fitted$optim$bestmem[-1]) # do not include `c`

if (mode == "hess"){
  library(MASS)
  ## Hessian ##
  cov_Mtx <- readRDS(paste0("flu_fit/", R0_mod, as.numeric(l_pnl), "_1g5e-3_covMtx.rds")) # read covariant matrix
  resamp <- mvrnorm(n=10, prms_optim, cov_Mtx) # 10 replicates for now
  saveRDS(rbind(prms_optim, resamp), file = paste0("flu_fit/", handle, "_prmsamps.rds"))
} else if (mode == "boot") {
  ## Bootstrap ##
  #resamp <- matrix(as.numeric(unlist(strsplit(prms_str, " +"))), byrow = TRUE, ncol = 4) # for premature checking
  resamp <- readRDS(paste0("flu_fit/sd_bootprms_", l_pnl, "_trial.rds")) # change handle
}

prms_sets <- rbind(prms_optim, resamp)
no_reps <- nrow(prms_sets)

states <- read.delim(state_ls, header = FALSE)$V1

## Load all state data
state_pop <- read.delim("state_lv_data/state_pop.tsv")
all_state_sun <- read.csv("state_lv_data/state_daily_sunrise_2019.csv")
all_state_day <- read.csv("state_lv_data/state_daytime_2019.csv")
# all_state_hum <- read.csv("state_lv_data/state_humidity_2014_2018.csv")
# all_state_hum <- all_state_hum[all_state_hum$X != "2016-02-29", ] # get rid of leap year Feb 29

# DST mask
DST <- rep(0, times=365)
DST[69:306] <- 60/720 # for DST change in 2019

permDST <- rep(60/720, times=365)
permDST[69:306] <- 0

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

dst_trajs <- matrix(ncol = no_reps, nrow = 364)
nodst_trajs <- matrix(ncol = no_reps, nrow = 364)
permdst_trajs <- matrix(ncol = no_reps, nrow = 364)

nodst_rc <- data.frame(matrix(ncol = no_reps, nrow = length(states)+1), row.names = c(states, "ALL"))
permdst_rc <- data.frame(matrix(ncol = no_reps, nrow = length(states)+1), row.names = c(states, "ALL"))

for (samp_i in seq(no_reps)){
  infection_dst <- data.frame(matrix(ncol = length(states), nrow = 364))
  infection_nodst <- data.frame(matrix(ncol = length(states), nrow = 364))
  infection_permdst <- data.frame(matrix(ncol = length(states), nrow = 364))
  colnames(infection_dst) <- states
  colnames(infection_nodst) <- states
  colnames(infection_permdst) <- states
  
  annual_cnt <- data.frame(wdst=numeric(length(states)), wodst=numeric(length(states)), permdst=numeric(length(states)), row.names = states)
  
  R0_params <- prms_sets[samp_i, ]
  for (state_i in seq(length(states))){
    state_code <- states[state_i]
    cat(">>> Simulating:", state_code, "<<<\n")
    census_pop <- state_pop$pop[state_pop$code==state_code]
    
    sunob <- all_state_sun[[state_code]]/720 # 365 days, scaled to maximum 720 = 12 hours
    dayob <- all_state_day[[state_code]]/1440
    
    # climob <- all_state_hum[[state_code]] # multiple years
    # climob <- matrix(climob, nrow = length(climob)/365, ncol = 365, byrow = TRUE)
    # climob <- colMeans(climob) # 365 days
    
    DST_R0 <- do.call(paste0("R0_", R0_mod), c(list(sunob, dayob), as.list(R0_params)))
    DST_sim <- sim_SIRS(census_pop, DST_R0)
    
    infection_dst[[state_code]] <- DST_sim$inf_traj
    annual_cnt[state_code, "wdst"] <- DST_sim$tot_inf
    
    if (state_code %in% c("AZ", "HI")){
      sansDST_R0 <- DST_R0
      permDST_R0 <- do.call(paste0("R0_", R0_mod), c(list(sunob+60/720, dayob), as.list(R0_params)))
    } else {
      sansDST_R0 <- do.call(paste0("R0_", R0_mod), c(list(sunob-DST, dayob), as.list(R0_params)))
      permDST_R0 <- do.call(paste0("R0_", R0_mod), c(list(sunob+permDST, dayob), as.list(R0_params)))
    }
    sansDST_sim <- sim_SIRS(census_pop, sansDST_R0)
    infection_nodst[[state_code]] <- sansDST_sim$inf_traj
    annual_cnt[state_code, "wodst"] <- sansDST_sim$tot_inf
    
    permDST_sim <- sim_SIRS(census_pop, permDST_R0)
    infection_permdst[[state_code]] <- permDST_sim$inf_traj
    annual_cnt[state_code, "permdst"] <- permDST_sim$tot_inf
  }
  dst_trajs[, samp_i] <- rowSums(infection_dst)
  nodst_trajs[, samp_i] <- rowSums(infection_nodst)
  permdst_trajs[, samp_i] <- rowSums(infection_permdst)
  
  nodst_rc[1:length(states), samp_i] <- (annual_cnt$wodst - annual_cnt$wdst)/annual_cnt$wdst
  permdst_rc[1:length(states), samp_i] <- (annual_cnt$permdst - annual_cnt$wdst)/annual_cnt$wdst
  
  nodst_rc[(length(states)+1), samp_i] <- (sum(annual_cnt$wodst) - sum(annual_cnt$wdst))/sum(annual_cnt$wdst)
  permdst_rc[(length(states)+1), samp_i] <- (sum(annual_cnt$permdst) - sum(annual_cnt$wdst))/sum(annual_cnt$wdst)
  
  cat(">> Rep", samp_i, ":", nodst_rc[(length(states)+1), samp_i], permdst_rc[(length(states)+1), samp_i], "\n")
}

write.table(dst_trajs, file = paste0("flu_fit/", handle, "_dstTraj.tsv"), quote=FALSE, sep="\t", row.names=FALSE)
write.table(nodst_trajs, file = paste0("flu_fit/", handle, "_sansdstTraj.tsv"), quote=FALSE, sep="\t", row.names=FALSE)
write.table(permdst_trajs, file = paste0("flu_fit/", handle, "_permdstTraj.tsv"), quote=FALSE, sep="\t", row.names=FALSE)
write.table(nodst_rc, file = paste0("flu_fit/", handle, "_sansdstRC.tsv"), quote=FALSE, sep="\t", row.names=TRUE)
write.table(permdst_rc, file = paste0("flu_fit/", handle, "_permdstRC.tsv"), quote=FALSE, sep="\t", row.names=TRUE)