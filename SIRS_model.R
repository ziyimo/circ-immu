#!/usr/bin/env Rscript

library("deSolve")

param_bounds[["cos"]] <- list(low = c(0, 1), high = c(1, 364))
param_bounds[["hum"]] <- list(low = c(0, -300), high = c(1, 0))
param_bounds[["day"]] <- list(low = c(0, -100, 0), high = c(1, 0, 1))
param_bounds[["hd"]] <- list(low = c(0, -300, -100, 0), high = c(1, 0, 0, 1))
param_bounds[["sd"]] <- list(low = c(0, -100, 0, -100, 0), high = c(1, 0, 1, 0, 1))
param_bounds[["hsd"]] <- list(low = c(0, -300, -100, 0, -100, 0), high = c(1, 0, 0, 1, 0, 1))

param_bounds[["sun"]] <- list(low = c(0, -100, 0), high = c(1, 0, 1))

# SIRS model
SIRS_R0 <- function(time, state, theta){
  ## Parameters:
  R0_vec <- theta[["R0"]] # vector pre-calculated R0 values
  D <- theta[["D"]] # scalar: Infectious perioud
  L <- theta[["L"]] # scalar: Duration of immunity
  
  R0_t <- R0_vec[time] # get R0 at time t
  
  ## States:
  S <- state["S"] # suceptible
  I <- state["I"] # infected
  R <- state["R"] # recovered
  N <- S + I + R # total population
  
  # Ordinary differential equations
  beta = R0_t/D
  
  dS <- (R/L) -beta * S * I/N
  dI <- beta * S * I/N - (I/D)
  dR <- (I/D) - (R/L)
  
  return(list(c(dS, dI, dR)))
}

# Global Var (for efficiency)
offset <- 0 # specify when to seed the single infection
tot_years <- 50 # no. of years to run SIRS model till steady state
times <- seq(1, 364 * tot_years, by = 1)[1:(364*tot_years-offset)]

# SIRS model prediction given variables
SIRSvar_pred <- function(R0_func, R0_params, var_ls, census_pop){
  # R0_func: string of function name
  # R0_params: vector of parameters
  # var_ls: list of vector(s) of variables (sunrise, humidity)
  
  xstart <- c(S = census_pop-1, I = 1, R = 0) # use actual state population
  annual_R0 <- do.call(R0_func, c(var_ls, as.list(R0_params)))
  
  # some hard-coded parameters
  paras = list(D = 5, 
               L = 40*7, 
               R0 = rep(head(annual_R0, 364), times = tot_years)[(offset+1):(364*tot_years)]) # pre-calculated R0 values
  
  predictions <- as.data.frame(ode(xstart, times, SIRS_R0, paras)) # run SIRS model
  raw_p <- predictions$I/census_pop
  return(raw_p)
  }

# Generate q (binom probability) from raw values of p
p2q <- function(mod_pred, epi_df, c_scaling){ # epi_df should be filtered ahead of time!
  
  data_years <- ceiling(max(epi_df$rel_date)/364)
  mod_pred <- mod_pred[((tot_years-data_years)*364-offset+1):(tot_years*364-offset)] # take only the number of years corresponding to data
  
  p_weekly <- mod_pred[epi_df$rel_date] # model prediction on particular days
  q <- c_scaling*p_weekly/epi_df$pi

  return(q)
}

## Wrapper for (penalized) likelihood function, separate `c`, for plotting likelihood surface purpose
binom_L <- function(prms, R0_model, var_obs, epi_df, pop_size, q_cap, c = 1, lambda = 0){
  state_p <- SIRSvar_pred(R0_model, prms, var_obs, pop_size)
  state_q <- p2q(state_p, epi_df, c)
  q_adj <- pmin(state_q, q_cap) # cap the value of q, this is q'
  neg_log_L <- -sum(dbinom(epi_df$k, size=epi_df$TT, prob=q_adj, log=TRUE)) + lambda*sum(state_q-q_adj)
  return(neg_log_L)
}

## Wrappers for penalized likelihood function, includes `c` for fitting
binom_Lp <- function(prms, R0_model, var_obs, epi_df, pop_size, q_cap, lambda = 0){
  state_p <- SIRSvar_pred(R0_model, prms[-1], var_obs, pop_size)
  state_q <- p2q(state_p, epi_df, prms[1])
  q_adj <- pmin(state_q, q_cap) # cap the value of q, this is q'
  neg_log_L <- -sum(dbinom(epi_df$k, size=epi_df$TT, prob=q_adj, log=TRUE)) + lambda*sum(state_q-q_adj)
  return(neg_log_L)
}

## Function for loading epidemiological data from csv file
load_state_epi <- function(file_path){
  epiob <- read.csv(file_path) # blank field automatically NA
  epiob <- epiob[epiob$WEEK <= 52, ] # cap year to 52 weeks
  
  # calculate relative date to the beginning of the 1st year in the data
  y0 <- epiob$YEAR[1]
  rel_date <- (epiob$YEAR-y0)*364 + epiob$WEEK*7 - 3 # center on Thursday
  
  k <- epiob$TOTAL.A + epiob$TOTAL.B
  TT <- epiob$TOTAL.SPECIMENS
  pi <- epiob$X.UNWEIGHTED.ILI/100
  mask <- is.na(TT) | (TT == 0) | (pi == 0) # missing data: no. of test is 0 or NA, or no sympotomatic patient
  cat(sum(mask), "weeks masked", "\n")
  
  epi_df <- data.frame("rel_date"=rel_date, "k"=k, "TT"=TT, "pi"=pi)
  epi_df <- epi_df[!mask, ]
  
  return(epi_df)
}
