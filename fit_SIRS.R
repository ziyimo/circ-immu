#!/usr/bin/env Rscript

library(DEoptim)
library(binom)
library(deSolve)
source("SIRS_model.R")

time_stamp <- format(Sys.time(), "%m%d-%H%M")
args <- commandArgs(trailingOnly=TRUE)

state_code <- args[1]   # 2-letter state code
fit_var <- args[2]      # `sun` or `cli` or `both`
#prms_init <- args[3]    # supply a path to DEoptim .rds file, or `new` to run DEoptim from scratch
optim_method <- args[3] # `DE` (genetic), `LBFGSB`, `GS` (grid search, for sun and cli only)
#optim_method <- "hyb" # DE followed by LBFGSB

## Load population
state_pop <- read.delim("state_lv_data/state_pop.tsv")
census_pop <- state_pop$pop[state_pop$code==state_code]

if (fit_var == "sun" | fit_var == "both"){
  ## Load sunrise data
  all_state_sun <- read.csv("state_lv_data/state_daily_sunrise_2019.csv")
  sunob <- all_state_sun[[state_code]]/720 # 365 days, scaled to maximum 720 = 12 hours
}

if (fit_var == "cli" | fit_var == "both"){
  ## Load humidity data
  all_state_hum <- read.csv("state_lv_data/state_humidity_2014_2018.csv")
  all_state_hum <- all_state_hum[all_state_hum$X != "2016-02-29", ] # get rid of leap year Feb 29
  climob <- all_state_hum[[state_code]] # multiple years
  climob <- matrix(climob, nrow = length(climob)/365, ncol = 365, byrow = TRUE)
  climob <- colMeans(climob) # 365 days
}

## Load flu data
epiob <- read.csv(paste0("state_lv_data/Flu_data/flu_epi_", state_code, ".csv")) # blank field automatically NA
epiob <- epiob[epiob$WEEK <= 52, ] # cap year to 52 weeks

# calculate relative date to the beginning of the 1st year in the data
y0 <- epiob$YEAR[1]
rel_date <- (epiob$YEAR-y0)*364 + epiob$WEEK*7 - 3 # center on Thursday

k <- epiob$TOTAL.A + epiob$TOTAL.B
TT <- epiob$TOTAL.SPECIMENS
pi <- epiob$X.UNWEIGHTED.ILI/100

epi_data <- data.frame("rel_date"=rel_date, "k"=k, "TT"=TT, "pi"=pi)

missing <- is.na(TT) | (TT == 0)
# Calculate q_cap
q_99 <- binom.confint(k[!missing], TT[!missing], conf.level = 0.99, methods = c("exact"))
q_99cap <- max(q_99$upper)

sink(paste0("fit_results/", state_code, fit_var, optim_method, time_stamp, ".log"))
cat("q capped at", q_99cap, "\n")

## Optimization
cat("Fitting:", state_code, fit_var, "with", optim_method, "\n")

if (fit_var == "both"){
  if (optim_method == "DE"){
    fit_result <- DEoptim(binom2_L, lower=c(0, 0, -3000, 0), upper=c(100, 1, 0, 0.022),
                           control=DEoptim.control(trace = 5, reltol = 1e-5),
                           var1_obs = sunob, var2_obs = climob,
                           epi_df = epi_data, pop_size = census_pop, q_cap = q_99cap)
    #saveRDS(evo_optim, file = paste0("fit_results/", state_code, fit_var, time_stamp, "_DE.rds"))
  } else if (optim_method == "LBFGSB") {
    # evo_optim <- readRDS(prms_init)
    # cat("Initial param vals:", evo_optim[["optim"]][["bestmem"]], "; negLL:", evo_optim[["optim"]][["bestval"]], "\n")
    fit_result <- optim(c(10, 0.5, -100, 0.01), binom2_L,
                        var1_obs = sunob, var2_obs = climob,
                        epi_df = epi_data, pop_size = census_pop, q_cap = q_99cap,
                        lower=c(0, 0, -3000, 0), upper=c(100, 1, 0, 0.022),
                        control = list(trace=3, parscale=c(100, 1, 3000, 0.022), REPORT=5), method = "L-BFGS-B")
  }
} else {
  if (fit_var == "cli"){
    param_low <- c(-3000, 0)
    param_high <- c(0, 0.022)
    param_init <- c(-100, 0.01)
    param_scale <- c(3000, 0.022)
    feed_var <- climob
  } else if (fit_var == "sun"){
    param_low <- c(0, 0)
    param_high <- c(100, 1)
    param_init <- c(10, 0.5)
    param_scale <- c(100, 1)
    feed_var <- sunob
  }
  if (optim_method == "DE"){
    fit_result <- DEoptim(binom1_L, lower=param_low, upper=param_high,
                           control=DEoptim.control(trace = 5, reltol = 1e-5),
                           var_obs = feed_var, epi_df = epi_data, pop_size = census_pop, q_cap = q_99cap)
    #saveRDS(evo_optim, file = paste0("fit_results/", state_code, fit_var, time_stamp, "_DE.rds"))
  } else if (optim_method == "LBFGSB") {
    # evo_optim <- readRDS(prms_init)
    # cat("Initial param vals:", evo_optim[["optim"]][["bestmem"]], "; negLL:", evo_optim[["optim"]][["bestval"]], "\n")
    fit_result <- optim(param_init, binom1_L,
                        var_obs = feed_var, epi_df = epi_data, pop_size = census_pop, q_cap = q_99cap,
                        lower=param_low, upper=param_high,
                        control = list(trace=3, parscale=param_scale, REPORT=5), method = "L-BFGS-B")
  } else if (optim_method == "GS") {
    dim_size <- 50
    k_range <- seq(param_low[1], param_high[1], length.out = dim_size)
    q0_range <- seq(param_low[2], param_high[2], length.out = dim_size)
    LL_surf <- matrix(nrow = dim_size, ncol = dim_size)
    
    for (i in 1:dim_size){
      for (j in 1:dim_size){
        LL_surf[i, j] <-  binom1_L(c(k_range[i], q0_range[j]), feed_var, epi_data, census_pop, q_99cap)
        cat(k_range[i], q0_range[j], ":", LL_surf[i, j], "\n")
      }
    }
    fit_result <- list(k_range, q0_range, LL_surf)
  }
}

saveRDS(fit_result, file = paste0("fit_results/", state_code, fit_var, optim_method, time_stamp, ".rds"))
sink()
