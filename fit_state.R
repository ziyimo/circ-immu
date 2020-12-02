#!/usr/bin/env Rscript

library(DEoptim)
library(binom)

source("SIRS_model.R")

# grad_step <- 1e-5 # step size for finite-difference approx. to the gradient
# T_start <- 10 # starting temp for simulated annealing, default is 10

time_stamp <- format(Sys.time(), "%m%d%H")
args <- commandArgs(trailingOnly=TRUE)

state_code <- args[1]   # 2-letter state code
# fit_var <- args[2]      # `sun` or `cli` or `both`, implied by the R0 model
R0_mod <- args[2]       # functional form of R0, options: exp, cdexp, bell, linEE, linGE, mixGE
optim_arg <- "new"       # argument regarding optimization: 
l_pnl <- as.numeric(args[3]) # lambda for penalized likelihood

#suffix <- ifelse(optim_arg=="GS", "GS", "DE")
handle <- paste0(state_code, l_pnl, R0_mod, time_stamp)
sink(paste0("fit_results/", handle, ".log"))

## Load population
state_pop <- read.delim("state_lv_data/state_pop.tsv")
census_pop <- state_pop$pop[state_pop$code==state_code]

if (R0_mod %in% c("cdexp", "bell", "linEE", "linGE", "mixGE")){
  ## Load sunrise data
  all_state_sun <- read.csv("state_lv_data/state_daily_sunrise_2019.csv")
  sunob <- all_state_sun[[state_code]]/720 # 365 days, scaled to maximum 720 = 12 hours
}

if (R0_mod %in% c("exp", "linEE", "linGE", "mixGE")){
  ## Load humidity data
  all_state_hum <- read.csv("state_lv_data/state_humidity_2014_2018.csv")
  all_state_hum <- all_state_hum[all_state_hum$X != "2016-02-29", ] # get rid of leap year Feb 29
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

## Load flu data
epi_data <- load_state_epi(paste0("state_lv_data/Flu_data/flu_epi_", state_code, ".csv"))

# Calculate q_cap
q_conf <- binom.confint(epi_data$k, epi_data$TT, conf.level = 1-1e-3, methods = c("exact"))
q_99cap <- max(q_conf$upper)
cat(">> q capped at", q_99cap, "\n")

## Optimization
cat(">> Fitting:", state_code, "; R0 model:", R0_mod, "; lambda=", l_pnl, "\n")

if (optim_arg == "new"){
  cat(">> Run DEoptim from scratch\n")
  evo_optim <- DEoptim(binom_Lp, lower=param_bounds[[R0_mod]]$low, upper=param_bounds[[R0_mod]]$high,
                         control=DEoptim.control(trace = 5, reltol = 1e-5, itermax = 500, steptol = 100),
                         R0_model = paste0("R0_", R0_mod), var_obs = varob,
                         epi_df = epi_data, pop_size = census_pop, q_cap = q_99cap, lambda = l_pnl)
  saveRDS(evo_optim, file = paste0("fit_results/", handle, "_DE.rds"))
}

#saveRDS(fit_result, file = paste0("fit_results/", state_code, l_pnl, fit_var, time_stamp, suffix, ".rds"))
sink()
