#!/usr/bin/env Rscript

#.libPaths() # sanity check
library("parallel")
library("DEoptim")
library("binom")

source("R0_mods.R")
source("flu_analysis/SIRS_model.R")
var_cache <- ls()

time_stamp <- format(Sys.time(), "%m%d")
args <- commandArgs(trailingOnly=TRUE)

state_ls <- args[1]      # text file with a list of states to fit to
R0_mod <- args[2]        # R0 model, options: cos, hum, day, hd, sd, hsd
l_pnl <- args[3]         # lambda for penalized likelihood
outdir <- args[4]
nthr <- as.numeric(args[5]) # limit no. of threads

handle <- paste0(state_ls, "_", l_pnl, "_", R0_mod, "_", time_stamp)

share_intcpt <- TRUE
l_pnl <- as.numeric(l_pnl)

sink(paste0(outdir, "/", handle, ".log"))

states <- read.delim(state_ls, header = FALSE)

## Load all state data
state_pop <- read.delim("compiled_data/geo_clim_dem/state_pop.tsv")
all_state_sun <- read.csv("compiled_data/geo_clim_dem/state_daily_sunrise_2019.csv")
all_state_day <- read.csv("compiled_data/geo_clim_dem/state_daytime_2019.csv")
all_state_hum <- read.csv("compiled_data/geo_clim_dem/state_humidity_2014_2018.csv")
all_state_hum <- all_state_hum[all_state_hum$X != "2016-02-29", ] # get rid of leap year Feb 29

state_data <- list()

for (state_i in seq(nrow(states))){
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
  
  epi_data <- load_state_epi(paste0("compiled_data/Flu_data/flu_epi_", state_code, ".csv"))
  
  q_conf <- binom.confint(epi_data$k, epi_data$TT, conf.level = 1-1e-3, methods = c("exact"))
  q_99cap <- min(max(q_conf$upper), 1-1e-3)
  cat(">> q capped at", q_99cap, "\n")
  
  state_data[[state_code]] <- list(idx=state_i, pop=census_pop, var=varob, epi=epi_data, cap=q_99cap)
}

p_upper <- param_bounds[[R0_mod]]$high
p_lower <- param_bounds[[R0_mod]]$low

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

all_state_Lp <- function(prms){ # everything else read as global variable
  states_negLL <- sapply(X=state_data, FUN=Lp_wrapper, p_vec=prms)
  return(sum(states_negLL))
}

cat(">> Fitting with R0 model:", R0_mod, "; lambda =", l_pnl, "\nLower:", p_lower, "\nUpper:", p_upper, "\n")
cat(detectCores(), "cores seen; limit to", nthr, "cores\n")
clstr <- makeCluster(nthr, type = "FORK") # limits the number of cores to use
evo_optim <- DEoptim(all_state_Lp, lower=p_lower, upper=p_upper,
                     control=DEoptim.control(trace = 5, reltol = 1e-6, itermax = 2000, steptol = 200,
                                             cluster=clstr, packages=c("deSolve"),
                                             parVar=c(var_cache, "state_data", "Lp_wrapper", "R0_mod", "l_pnl")))
stopCluster(clstr)

saveRDS(evo_optim, file = paste0(outdir, "/", handle, ".rds"))
sink()
