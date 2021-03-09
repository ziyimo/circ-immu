#!/usr/bin/env Rscript

#.libPaths() # sanity check
library("parallel")
library("binom")
#library("numDeriv")

source("R0_mods.R")
source("SIRS_model.R")

# var_cache <- ls()
# time_stamp <- format(Sys.time(), "%m%d")

if (interactive()){
  state_ls <- "states_QC200.tsv"
  R0_mod <- "sd" 
  l_pnl <- "1e1"
  nthr <- 4
} else {
  args <- commandArgs(trailingOnly=TRUE)
  
  state_ls <- args[1]      # text file with a list of states to fit to
  R0_mod <- args[2]        # R0 model, options: cos, hum, day, hd, sd, hsd
  #share_intcpt <- args[3]  # "0" or "1"; deprecated
  l_pnl <- args[3]         # lambda for penalized likelihood
  discret <- args[4] # step for gradient approximation, path to sampled parameters when below is -1
  grad2hess <- as.numeric(args[5]) # this is the parscale parameter in optimHess, put -1 for evaluating sampled parameters
  nthr <- as.numeric(args[6]) # limit no. of threads
}

handle <- Sys.glob(paste0("flu_fit/", paste(state_ls, l_pnl, R0_mod, "*.rds" ,sep="_")))
l_pnl <- as.numeric(l_pnl)

states <- read.delim(state_ls, header = FALSE)

## Load all state data
state_pop <- read.delim("state_lv_data/state_pop.tsv")
all_state_sun <- read.csv("state_lv_data/state_daily_sunrise_2019.csv")
all_state_day <- read.csv("state_lv_data/state_daytime_2019.csv")
all_state_hum <- read.csv("state_lv_data/state_humidity_2014_2018.csv")
all_state_hum <- all_state_hum[all_state_hum$X != "2016-02-29", ] # get rid of leap year Feb 29

state_data <- list()

for (state_i in seq(nrow(states))){
  state_code <- states$V1[state_i]
  #cat(">>> Loading state data:", state_code, "<<<\n")
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
  #cat(">> q capped at", q_99cap, "\n")
  
  state_data[[state_code]] <- list(idx=state_i, pop=census_pop, var=varob, epi=epi_data, cap=q_99cap)
}

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

all_state_Lp <- function(mod_prms){ # everything else read as global variable
  prms <- c(c_param, mod_prms)
  states_negLL <- parSapply(cl = clstr, X=state_data, FUN=Lp_wrapper, p_vec=prms)
  return(-sum(states_negLL)) # return positive log-likelihood
}

fitted <- readRDS(handle)
#fitted <- readRDS("flu_fit/states_QC200.tsv_1e1_sd_0111.rds")
fit_params <- unname(fitted$optim$bestmem)
c_param <- fit_params[1]
mod_param <- fit_params[-1]

cat(detectCores(), "cores seen; limit to", nthr, "cores\n")
clstr <- makeCluster(nthr, type = "FORK") # limits the number of cores to use
# check function value at optimum
opt_val <- all_state_Lp(mod_param)
cat(handle, ":", opt_val, fitted$optim$bestval, "\n")

if (FALSE){ # sanity check
  #load("covid_hosp_fit/eta.026_sd_iter.rds") # ss_fit and seed_sh
  cov_Mtx <- readRDS(paste0("flu_fit/", R0_mod, l_pnl, "_1g5e-3_covMtx.rds")) # read covariant matrix
  
  resamp <- mvrnorm(n=5, mod_param, cov_Mtx)
  
  for (i in seq(5)){
    cat(all_state_Lp(resamp[i,]), "\t")
  }
}

if (grad2hess == -1){
  par_samps <- readRDS(discret)
  
  for (samp_i in seq(nrow(par_samps))){
    cat(par_samps[samp_i,], ":", all_state_Lp(par_samps[samp_i,]), "\n")
  }
  
} else {
  prms_range <- param_bounds[[R0_mod]]$high - param_bounds[[R0_mod]]$low
  hessian_step <- prms_range[-1]*as.numeric(discret) # eliminate 'c' from range
  cat(">>", hessian_step, "\n")
  cat(">> grad2hess ratio:", grad2hess, "\n")
  #info_mtx <- hessian(all_state_Lp, fit_params)
  info_mtx <- -optimHess(mod_param, all_state_Lp,
                         control = list(parscale=rep(grad2hess, length(mod_param)), ndeps=hessian_step))
  cat(info_mtx, "\n")
  cov_mtx <- solve(info_mtx)
  cat(">>", R0_mod, l_pnl, ";", grad2hess, discret, ":", diag(cov_mtx), "\n")
  
  saveRDS(cov_mtx, file=paste0("flu_fit/", R0_mod, l_pnl, "_", grad2hess, "g", discret, "_covMtx.rds"))
}

stopCluster(clstr)
#sink()
