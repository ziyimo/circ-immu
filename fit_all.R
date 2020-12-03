#!/usr/bin/env Rscript

#.libPaths() # sanity check
library("parallel")
library("DEoptim")
library("binom")

source("SIRS_model.R")
var_cache <- ls()

time_stamp <- format(Sys.time(), "%m%d")

if (interactive()){
  state_ls <- "states_QC200.tsv"
  R0_mod <- "mixGE"
  share_intcpt <- "0" 
  l_pnl <- "1e2"
} else {
  args <- commandArgs(trailingOnly=TRUE)
  
  state_ls <- args[1]      # text file with a list of states to fit to
  R0_mod <- args[2]        # functional form of R0, options: exp, cdexp, bell, linEE, linGE, mixGE
  share_intcpt <- args[3]  # "0" or "1" 
  l_pnl <- args[4]         # lambda for penalized likelihood
}

handle <- paste0(state_ls, l_pnl, R0_mod, share_intcpt, "_", time_stamp)

share_intcpt <- as.logical(as.numeric(share_intcpt))
l_pnl <- as.numeric(l_pnl)

sink(paste0("fit_results/", handle, ".log"))

states <- read.delim(state_ls, header = FALSE)

## Load all state data
state_pop <- read.delim("state_lv_data/state_pop.tsv")
all_state_sun <- read.csv("state_lv_data/state_daily_sunrise_2019.csv")
all_state_hum <- read.csv("state_lv_data/state_humidity_2014_2018.csv")
all_state_hum <- all_state_hum[all_state_hum$X != "2016-02-29", ] # get rid of leap year Feb 29

state_data <- list()

for (state_i in seq(nrow(states))){
  state_code <- states$V1[state_i]
  cat(">>> Loading state data:", state_code, "<<<\n")
  census_pop <- state_pop$pop[state_pop$code==state_code]
  if (R0_mod %in% c("cdexp", "bell", "linEE", "linGE", "mixGE")){
    sunob <- all_state_sun[[state_code]]/720 # 365 days, scaled to maximum 720 = 12 hours
  }
  if (R0_mod %in% c("exp", "linEE", "linGE", "mixGE")){
    climob <- all_state_hum[[state_code]] # multiple years
    climob <- matrix(climob, nrow = length(climob)/365, ncol = 365, byrow = TRUE)
    climob <- colMeans(climob) # 365 days
  }
  
  if (R0_mod == "exp"){
    varob <- list(climob)
  } else if (R0_mod %in% c("cdexp", "bell")){
    varob <- list(sunob)
  } else if (R0_mod %in% c("linEE", "linGE", "mixGE")){
    varob <- list(sunob, climob)
  }
  
  epi_data <- load_state_epi(paste0("state_lv_data/Flu_data/flu_epi_", state_code, ".csv"))
  
  q_conf <- binom.confint(epi_data$k, epi_data$TT, conf.level = 1-1e-3, methods = c("exact"))
  q_99cap <- max(q_conf$upper)
  cat(">> q capped at", q_99cap, "\n")
  
  state_data[[state_code]] <- list(idx=state_i, pop=census_pop, var=varob, epi=epi_data, cap=q_99cap)
}

if (share_intcpt){ # all states share same parameter set
  p_upper <- param_bounds[[R0_mod]]$high
  p_lower <- param_bounds[[R0_mod]]$low
  
  get_state_prms <- function(prms_vec, state_idx){
    return(prms_vec)
  }
} else {
  p_upper <- c(param_bounds[[R0_mod]]$high[-3], rep(1, nrow(states)))
  p_lower <- c(param_bounds[[R0_mod]]$low[-3], rep(0, nrow(states)))
  
  if (R0_mod %in% c("cdexp", "bell")){
    get_state_prms <- function(prms_vec, state_idx){
      return(c(prms_vec[1:2], prms_vec[2+state_idx]))
    }
  } else if (R0_mod %in% c("linEE", "linGE")){
    get_state_prms <- function(prms_vec, state_idx){
      return(c(prms_vec[1:2], prms_vec[3+state_idx], prms_vec[3]))
    }
  } else if (R0_mod == "mixGE"){
    get_state_prms <- function(prms_vec, state_idx){
      return(c(prms_vec[1:2], prms_vec[4+state_idx], prms_vec[3:4]))
    }
  }
}

Lp_wrapper <- function(state_elem, p_vec){
  neg_LL <- binom_Lp(get_state_prms(p_vec, state_elem$idx), 
                     paste0("R0_", R0_mod),
                     state_elem$var,
                     state_elem$epi,
                     state_elem$pop,
                     state_elem$cap,
                     l_pnl)
  return(neg_LL)
}

all_state_Lp <- function(prms){ # everything else read as global variable
  states_negLL <- sapply(X=state_data, FUN=Lp_wrapper, p_vec=prms)
  return(sum(states_negLL))
}

cat(">> Fitting with R0 model:", R0_mod, "; lambda =", l_pnl, "\nLower:", p_lower, "\nUpper:", p_upper, "\n")
cat(detectCores(), "cores seen\n")
clstr <- makeCluster(15, type = "FORK") # hard-code the number of cores to use
evo_optim <- DEoptim(all_state_Lp, lower=p_lower, upper=p_upper,
                     control=DEoptim.control(trace = 5, reltol = 1e-6, itermax = 1000, steptol = 200,
                                             parallelType=1, cluster=clstr, packages=c("deSolve"),
                                             parVar=c(var_cache, "state_data", "Lp_wrapper", "get_state_prms", "R0_mod", "l_pnl")))

saveRDS(evo_optim, file = paste0("fit_results/", handle, ".rds"))
sink()
