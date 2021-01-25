#!/usr/bin/env Rscript

library(deSolve)

param_bounds[["cos"]] <- list(low = c(1), high = c(364))
param_bounds[["hum"]] <- list(low = c(-300), high = c(0))
param_bounds[["sun"]] <- list(low = c(-500, 0), high = c(0, 1))
param_bounds[["hs"]] <- list(low = c(-300, -500, 0), high = c(0, 0, 1))

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


run_SEIH <- function(census_pop, prop_init, hpt_rate, R0_func, R0_params, var_ls){
    
  xstart <- c(S = census_pop*(1-prop_init), E = 0, I = census_pop*prop_init, H = 0, R = 0) # use actual state population
  times <- seq(73, 365) # from national emergency declaration to end of the year
  annual_R0 <- do.call(R0_func, c(var_ls, as.list(R0_params)))
  
  # some hard-coded parameters
  paras = list(sigma = 5,
               h = hpt_rate,
               lambda = 5,
               gam = 5,
               k = 10,
               R0 = annual_R0) # pre-calculated R0 values
  
  predictions <- as.data.frame(ode(xstart, times, SEIH_R0, paras)) # run SIRS model
  return(predictions)
}

pois_L <- function(prms, R0_model, var_obs, hosp_df, pop_size){
  seed_prop <- prms[1]
  hrate <- prms[2]
  R0_min <- prms[3]
  R0_range <- prms[4]
  R0_prms <- prms[-1:-4]
  
  mod_pred <- run_SEIH(pop_size, seed_prop, hrate, R0_model, c(R0_prms, R0_min+R0_range, R0_min), var_obs)
  result <- merge(mod_pred, hosp_df, by.x = "time", by.y = "date")
  neg_log_L <- -sum(dpois(result$hospitalizedCurrently, result$H, log = TRUE))
  return(neg_log_L)
}

norm_L <- function(prms, R0_model, var_obs, hosp_df, pop_size){
  seed_prop <- prms[1]
  hrate <- prms[2]
  R0_min <- prms[3]
  R0_range <- prms[4]
  R0_prms <- prms[-1:-4]
  
  mod_pred <- run_SEIH(pop_size, seed_prop, hrate, R0_model, c(R0_prms, R0_min+R0_range, R0_min), var_obs)
  result <- merge(mod_pred, hosp_df, by.x = "time", by.y = "date")
  neg_log_L <- -sum(dnorm(result$hospitalizedCurrently, mean = result$H, log = TRUE))
  return(neg_log_L)
}