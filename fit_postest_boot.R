#!/usr/bin/env Rscript

#.libPaths() # sanity check
library("parallel")
library("DEoptim")
library("binom")

source("R0_mods.R")
source("SIRS_model.R")
var_cache <- ls()

time_stamp <- format(Sys.time(), "%m%d")

args <- commandArgs(trailingOnly=TRUE)
  
state_ls <- args[1]      # text file with a list of states to fit to
best_fit <- args[2]   # best fit parameters (for sanity checking)
R0_mod <- "sd"        # R0 model, options: cos, hum, day, hd, sd, hsd
tag <- args[3]        # identifier of the bootstrap run
#share_intcpt <- args[3]  # "0" or "1"; deprecated
l_pnl <- args[4]         # lambda for penalized likelihood
nthr <- as.numeric(args[5]) # limit no. of threads

handle <- paste0(R0_mod, "_boot", l_pnl, "_", tag) #, "_", time_stamp)

l_pnl <- as.numeric(l_pnl)

#sink(paste0("fit_results/", handle, ".log"))

states <- read.delim(state_ls, header = FALSE)
no_states <- nrow(states)

## Load all state data
state_pop <- read.delim("state_lv_data/state_pop.tsv")
all_state_sun <- read.csv("state_lv_data/state_daily_sunrise_2019.csv")
all_state_day <- read.csv("state_lv_data/state_daytime_2019.csv")
all_state_hum <- read.csv("state_lv_data/state_humidity_2014_2018.csv")
all_state_hum <- all_state_hum[all_state_hum$X != "2016-02-29", ] # get rid of leap year Feb 29

bssample <- sort(sample(seq(no_states), size=no_states, replace = TRUE)) # a bootstrap sample vector
samp_states <- as.data.frame(table(bssample), stringsAsFactors = FALSE)
samp_states$bssample <- as.numeric(samp_states$bssample)

state_data <- list()

for (state_i in samp_states$bssample){
  state_code <- states$V1[state_i]
  cat(">>> Loading state data:", state_code, "<<<\n")
  census_pop <- state_pop$pop[state_pop$code==state_code]
  
  if (R0_mod %in% c("sun", "sd", "hsd")){
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
  } else if (R0_mod == "sun"){
    varob <- list(sunob)
  }
  
  epi_data <- load_state_epi(paste0("state_lv_data/Flu_data/flu_epi_", state_code, ".csv"))
  
  q_conf <- binom.confint(epi_data$k, epi_data$TT, conf.level = 1-1e-3, methods = c("exact"))
  q_99cap <- min(max(q_conf$upper), 1-1e-3)
  cat(">> q capped at", q_99cap, "\n")
  
  state_data[[state_code]] <- list(idx=state_i, pop=census_pop, var=varob, epi=epi_data, cap=q_99cap)
}

p_upper <- param_bounds[[R0_mod]]$high[-1]
p_lower <- param_bounds[[R0_mod]]$low[-1]

Lp_wrapper <- function(state_elem, p_vec){
  neg_LL <- binom_Lp(p_vec, 
                     paste0("R0_", R0_mod),
                     state_elem$var,
                     state_elem$epi,
                     state_elem$pop,
                     state_elem$cap,
                     l_pnl)
  return(neg_LL)
}

full_opt <- readRDS(best_fit)
param_opt <- unname(full_opt$optim$bestmem)
norm_const <- param_opt[1]

all_state_Lp <- function(boot_prms){ # everything else read as global variable
  prms <- c(norm_const, boot_prms) # does not include c for bootstrapping
  states_negLL <- sapply(X=state_data, FUN=Lp_wrapper, p_vec=prms)
  return(sum(states_negLL*samp_states$Freq)) # weigh by sample frequency
}

negLL_opt <- all_state_Lp(param_opt[-1])
cat(">> Resamp", tag, "with R0 model:", R0_mod, "; lambda =", l_pnl, "\nLower:", p_lower, "\nUpper:", p_upper, "\n")
cat(">> Samples:", samp_states$bssample, "\n")
cat(">> Frequency:", samp_states$Freq, "\n")
cat(">> Opt param @ resamp:", negLL_opt, "\n")

cat(detectCores(), "cores seen; limit to", nthr, "cores\n")
clstr <- makeCluster(nthr, type = "FORK") # limits the number of cores to use
evo_optim <- DEoptim(all_state_Lp, lower=p_lower, upper=p_upper,
                     control=DEoptim.control(trace = 5, reltol = 1e-6, itermax = 2000, steptol = 200,
                                             cluster=clstr, packages=c("deSolve"),
                                             parVar=c(var_cache, "state_data", "Lp_wrapper", "R0_mod", "l_pnl", "samp_states", "norm_const")))
stopCluster(clstr)

saveRDS(evo_optim, file = paste0("fit_results/", handle, ".rds"))
#sink()

if (FALSE){
  ## Combined bootstrapping results
  
  bts_files <- Sys.glob("fit_results/sd_boot1e1_noc*.rds")
  no_boots <- length(bts_files)
  boot_params <- matrix(0, nrow=0, ncol = 4) # adjust the no. of params
  
  for (iter in seq(no_boots)){
    optim_obj <- readRDS(bts_files[iter])
    boot_params <- rbind(boot_params, unname(optim_obj$optim$bestmem))
  }
  
  # some manual entries (close to finished runs)
  boot_params <- rbind(boot_params, c(-28.109160, 0.533124, -8.215940, 0.186342))
  boot_params <- rbind(boot_params, c(-18.537164, 0.572048, -3.905891, 0.006229))
  boot_params <- rbind(boot_params, c(-25.047151, 0.535022, -6.456461, 0.107673))
  
  saveRDS(boot_params, file = "flu_fit/sd_bootprms_1e1_Apr6.rds")
}