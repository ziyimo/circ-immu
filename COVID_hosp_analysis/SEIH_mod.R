#!/usr/bin/env Rscript

library(deSolve)

param_bounds <- list()
param_bounds[["cos"]] <- list(low = c(1), high = c(364))
param_bounds[["hum"]] <- list(low = c(-300), high = c(0))
param_bounds[["day"]] <- list(low = c(-500, 0), high = c(0, 1))
param_bounds[["hd"]] <- list(low = c(-300, -500, 0), high = c(0, 0, 1))
param_bounds[["sd"]] <- list(low = c(-500, 0, -500, 0), high = c(0, 1, 0, 1))
param_bounds[["hsd"]] <- list(low = c(-300, -500, 0, -500, 0), high = c(0, 0, 1, 0, 1))

SEIH_R0 <- function(time, state, theta){
  ## Parameters:
  R0_vec <- theta[["R0"]] # vector pre-calculated R0 values
  sigma <- theta[["sigma"]] # Incubation period
  h <- theta[["h"]] # hospitalization rate
  lambda <- theta[["lambda"]] # time from symptom onset to hospitalization
  gam <- theta[["gam"]] # infectious period
  k <- theta[["k"]] # length of hospitalization
  
  R0_t <- R0_vec[time] # get R0 at time t
  
  ## States:
  S <- state["S"] # suceptible
  E <- state["E"] # exposed
  I <- state["I"] # infected
  H <- state["H"] # hospitalized
  R <- state["R"] # removed
  N <- S + E + I + H + R # total population
  
  # Ordinary differential equations
  beta = R0_t/gam
  
  dS <- -beta * S * I/N
  dE <- beta * S * I/N - E/sigma
  dI <- E/sigma - h/lambda*I - (1-h)/gam*I
  dH <- h/lambda*I - H/k
  dR <- (1-h)/gam*I + H/k
  
  return(list(c(dS, dE, dI, dH, dR)))
}

parameter_names <- c("prop_init", "incb_prd", "inf_prd", "hpt_rate", "hpt_prd", "R0min", "R0range")
run_SEIH <- function(R0_func, var_ls, census_pop, prop_init, incb_prd, inf_prd, hpt_rate, hpt_prd, R0min, R0range, R0_params){
    
  xstart <- c(S = census_pop*(1-prop_init), E = 0, I = census_pop*prop_init, H = 0, R = 0) # use actual state population
  times <- seq(73, 396) # from national emergency declaration to 01/31/2021
  annual_R0 <- do.call(R0_func, c(var_ls, as.list(c(R0_params, R0min+R0range, R0min))))
  
  paras = list(sigma = incb_prd,
               h = hpt_rate,
               lambda = inf_prd,
               gam = inf_prd,
               k = hpt_prd,
               R0 = rep(annual_R0, times=2)) # pre-calculated R0 values
  
  predictions <- as.data.frame(ode(xstart, times, SEIH_R0, paras)) # run SIRS model
  return(predictions)
}

norm_L <- function(prms, fixed_prms, tune_mask, R0_model, var_obs, hosp_df, pop_size){
  # 7 parameters: prop_init, incb_prd, inf_prd, hpt_rate, hpt_prd, R0min, R0range
  no_tune <- sum(tune_mask)
  tune_prms <- prms[1:no_tune]
  R0_prms <- prms[-1:-no_tune]
  
  all_prms <- numeric(7)
  all_prms[tune_mask] <- tune_prms
  all_prms[!tune_mask] <- fixed_prms
  
  mod_pred <- do.call(run_SEIH, c(list(R0_model, var_obs, pop_size), as.list(all_prms), list(R0_prms)))
  result <- merge(mod_pred, hosp_df, by.x = "time", by.y = "date")
  neg_log_L <- -sum(dnorm(result$hospitalizedCurrently, mean = result$H, log = TRUE))
  return(neg_log_L)
}