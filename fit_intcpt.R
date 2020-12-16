#!/usr/bin/env Rscript

#.libPaths() # sanity check
library("binom")
source("SIRS_model.R")

grad_step <- 1e-5 # step size for finite-difference approx. to the gradient, default: 1e-3
var_cache <- ls()
time_stamp <- format(Sys.time(), "%m%d")

args <- commandArgs(trailingOnly=TRUE)

state_ls <- args[1]      # text file with a list of states to fit to
R0_mod <- args[2]        # functional form of R0, options: cdexp, bell, linEE, linGE, mixGE
init <- as.numeric(strsplit(args[3], ":", fixed=TRUE)[[1]])          # initial values, separated by ":" 
l_pnl <- args[4]         # lambda for penalized likelihood
nthr <- as.numeric(args[5]) # limit no. of threads, put `0` for decoupled fitting

handle <- paste0(state_ls, l_pnl, R0_mod, "_intcpt", time_stamp)
l_pnl <- as.numeric(l_pnl)

sink(paste0("fit_results/", handle, ".log"))

states <- read.delim(state_ls, header = FALSE, stringsAsFactors = FALSE)

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
  
  if (R0_mod %in% c("cdexp", "bell")){
    varob <- list(sunob)
  } else if (R0_mod %in% c("linEE", "linGE", "mixGE")){
    varob <- list(sunob, climob)
  }
  
  epi_data <- load_state_epi(paste0("state_lv_data/Flu_data/flu_epi_", state_code, ".csv"))
  
  q_conf <- binom.confint(epi_data$k, epi_data$TT, conf.level = 1-1e-3, methods = c("exact"))
  q_99cap <- min(max(q_conf$upper), 1-1e-3)
  cat(">> q capped at", q_99cap, "\n")
  
  state_data[[state_code]] <- list(idx=state_i, pop=census_pop, var=varob, epi=epi_data, cap=q_99cap)
}

if (nthr == 0){
  Lp_wrapper <- function(s_0, state_elem){
    neg_LL <- binom_Lp(append(init[-3], s_0, after = 2),
                       paste0("R0_", R0_mod),
                       state_elem$var,
                       state_elem$epi,
                       state_elem$pop,
                       state_elem$cap,
                       l_pnl)
    return(neg_LL)
  }
  find_s0 <- function(XX) optimize(Lp_wrapper, c(0, 1), state_elem = XX)
  
  cat(">> Fitting s_0 decoupled; lambda =", l_pnl, "\n")
  s0_fit <- lapply(X=state_data, FUN=find_s0)
  s0_df <- do.call(rbind.data.frame, s0_fit)
  
  write.table(s0_df, file = paste0("fit_results/", handle, ".tsv"), quote=FALSE, sep="\t", row.names=FALSE)
  cat("Sum of likelihood:", sum(s0_df$objective))
  
} else {
  library("optimParallel")
  
  p_upper <- c(param_bounds[[R0_mod]]$high[-3], rep(1, nrow(states)))
  p_lower <- c(param_bounds[[R0_mod]]$low[-3], rep(0, nrow(states)))
  p_init <- c(init[-3], rep(init[3], nrow(states)))
  
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

  Lp_wrapper <- function(state_elem, p_vec){
    neg_LL <- binom_Lp(get_state_prms(p_vec, state_elem$idx), 
                       paste0("R0_", R0_mod),
                       state_elem$var,
                       state_elem$epi,
                       state_elem$pop,
                       state_elem$cap,
                       l_pnl)
    # if (is.infinite(neg_LL)){
    #   cat("Infinite negLL: state-", state_elem$idx, "; param set:", p_vec, "\n")
    # }
    return(neg_LL)
  }
  
  all_state_Lp <- function(prms){ # everything else read as global variable
    states_negLL <- sapply(X=state_data, FUN=Lp_wrapper, p_vec=prms)
    return(sum(states_negLL))
  }
  
  cat(">> Fitting with R0 model:", R0_mod, "; lambda =", l_pnl,
      "\nLower:", p_lower, "\nUpper:", p_upper, "\nInit:", p_init, "\n")
  cat(detectCores(), "cores seen; limit to", nthr, "cores\n")
  clstr <- makeCluster(nthr, type = "FORK") # hard-code the number of cores to use
  clusterExport(clstr, c(var_cache, "state_data", "Lp_wrapper", "get_state_prms", "R0_mod", "l_pnl"))
  
  soln <- optimParallel(p_init, all_state_Lp, lower = p_lower, upper = p_upper,
                        control = list(trace=3, parscale=(p_upper-p_lower),
                                       ndeps=rep(grad_step, length(p_init)),
                                       REPORT=5, factr=1e-6/.Machine$double.eps),
                        parallel = list(cl=clstr, forward=TRUE))
  #packages=c("deSolve")
  stopCluster(clstr)
  saveRDS(soln, file = paste0("fit_results/", handle, ".rds"))
}

sink()
