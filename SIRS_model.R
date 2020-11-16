#!/usr/bin/env Rscript

library(deSolve)

# R0 models
R0_sig <- function(q_t, k, q_0, R0base, R0min){ # single variable sigmoid
  R0 = (R0base-R0min)/(1+exp(-k*(q_t-q_0)))+R0min
  return(R0)
}

R0_sig2 <- function(q_t, r_t, k1, q_0, k2, r_0, R0base, R0min){ # two variable sigmoid
  R0 = (R0base-R0min)/(1+exp(-k1*(q_t-q_0)-k2*(r_t-r_0)))+R0min
  return(R0)
}

######## Plot model #######
if(interactive()){
  library(ggplot2)
  library(gridExtra)
  
  clean <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                 text = element_text(size=20))
  
  sun_range = seq(0, 1, by=0.01)
  cli_range = seq(0, 0.03, by=0.0003)
  
  p1 <- ggplot() + 
    geom_line(data = data.frame("hum"= cli_range, "R0"= R0_sig(cli_range, -100, 0.025, 2, 1.2)), 
              aes(x = hum, y=R0), color = "#E69F00") +
    geom_line(data = data.frame("hum"= cli_range, "R0"= R0_sig(cli_range, -1500, 0.015, 2, 1.2)), 
              aes(x = hum, y=R0), color = "#0072B2") +
    geom_line(data = data.frame("hum"= cli_range, "R0"= R0_sig(cli_range, -300, 0.01, 2, 1.2)), 
              aes(x = hum, y=R0), color = "#009E73") + clean
  
  p2 <- ggplot() + 
    geom_line(data = data.frame("sun"= sun_range, "R0"= R0_sig(sun_range, 100, 0.7, 2, 1.2)), 
              aes(x = sun, y=R0), color = "#E69F00") +
    geom_line(data = data.frame("sun"= sun_range, "R0"= R0_sig(sun_range, 10, 0.5, 2, 1.2)), 
              aes(x = sun, y=R0), color = "#0072B2") +
    geom_line(data = data.frame("sun"= sun_range, "R0"= R0_sig(sun_range, 30, 0.5, 2, 1.2)), 
              aes(x = sun, y=R0), color = "#009E73") + clean
  grid.arrange(p1, p2, ncol = 1) 
  
  sun_cli_df <- expand.grid(sun_range, cli_range)
  colnames(sun_cli_df) <- c("sun", "hum")
  
  sun_cli_df$R0 <- R0_sig2(sun_cli_df$sun, sun_cli_df$hum, 3, 0.5, -100, 0.01, 2, 1.2)
  p1 <- ggplot(sun_cli_df, aes(sun, hum, z = R0)) + geom_contour_filled() + clean
  
  sun_cli_df$R0 <- R0_sig2(sun_cli_df$sun, sun_cli_df$hum, 10, 0.5, -100, 0.01, 2, 1.2)
  p2 <- ggplot(sun_cli_df, aes(sun, hum, z = R0)) + geom_contour_filled() + clean
  
  sun_cli_df$R0 <- R0_sig2(sun_cli_df$sun, sun_cli_df$hum, 3, 0.5, -300, 0.01, 2, 1.2)
  p3 <- ggplot(sun_cli_df, aes(sun, hum, z = R0)) + geom_contour_filled() + clean
  
  sun_cli_df$R0 <- R0_sig2(sun_cli_df$sun, sun_cli_df$hum, 10, 0.5, -300, 0.01, 2, 1.2)
  p4 <- ggplot(sun_cli_df, aes(sun, hum, z = R0)) + geom_contour_filled() + clean
  
  grid.arrange(p1, p2, p3, p4, ncol = 2)
  grid.arrange(p2, p3, ncol = 1)
}
###########################

# SIRS model
SIRS_R0sig <- function(time, state, theta){
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

# single variable model prediction
SIRS1var_pred <- function(sig_params, var_observed, census_pop){
  k_step = sig_params[1]
  q_center = sig_params[2]
  
  # some hard-coded parameters
  offset = 0 # specify when to seed the single infection
  tot_years = 50
  xstart = c(S = census_pop-1, I = 1, R = 0) # use actual state population
  times = seq(1, 364 * tot_years, by = 1)[1:(364*tot_years-offset)]
  
  q_ob = rep(head(var_observed, 364), times = tot_years)[(offset+1):(364*tot_years)] # replicate variable over multiple years
  paras = list(D = 5, 
               L = 40*7, 
               R0 = R0_sig(q_ob, k_step, q_center, R0base = 2, R0min = 1.2)) # pre-calculate R0 values
  
  predictions <- as.data.frame(ode(xstart, times, SIRS_R0sig, paras)) # run SIRS model
  
  raw_p <- predictions$I/census_pop
  return(raw_p)
  }

# 2 variable model prediction
SIRS2var_pred <- function(sig2_params, var1_obs, var2_obs, census_pop){
  k1_step = sig2_params[1]
  q_center = sig2_params[2]
  k2_step = sig2_params[3]
  r_center = sig2_params[4]
  
  # some hard-coded parameters
  offset = 0 # specify when to seed the single infection
  tot_years = 50
  xstart = c(S = census_pop-1, I = 1, R = 0) # use actual state population
  times = seq(1, 364 * tot_years, by = 1)[1:(364*tot_years-offset)]
  
  q_ob = rep(head(var1_obs, 364), times = tot_years)[(offset+1):(364*tot_years)]
  r_ob = rep(head(var2_obs, 364), times = tot_years)[(offset+1):(364*tot_years)]
    
  paras = list(D = 5, 
               L = 40*7, 
               R0 = R0_sig2(q_ob, r_ob,
                            k1_step, q_center, k2_step, r_center,
                            R0base = 2, R0min = 1.2)) # pre-calculate R0 values
  
  predictions <- as.data.frame(ode(xstart, times, SIRS_R0sig, paras)) # run SIRS model
  
  raw_p <- predictions$I/census_pop
  return(raw_p)
}

# Generate q (binom probability) from raw values of p
p2q <- function(mod_pred, epi_df, q_cap, c_scaling){
  
  offset = 0
  tot_years = 50
  data_years <- ceiling(max(epi_df$rel_date)/364)
  model_range <- seq((tot_years-data_years)*364-offset+1, tot_years*364-offset)
  mod_pred <- mod_pred[model_range] # take only the number of years corresponding to data
  mod_pred <- mod_pred*c_scaling
  
  p_weekly <- mod_pred[epi_df$rel_date] # model prediction on particular days
  mask <- is.na(epi_df$TT) | (epi_df$TT == 0) # drop na or no test
  p_weekly <- p_weekly[!mask]
  epi_weekly <- epi_df[!mask, ]
  
  q <- p_weekly/epi_weekly$pi
  q_adj <- pmin(q, q_cap) # cap the value of q
  
  return(list(q_adj, epi_weekly))
}

## Wrappers for likelihood function
binom2_L <- function(prms, var1_obs, var2_obs, epi_df, pop_size, q_cap, c = 1){ # 2 variable
  state_p <- SIRS2var_pred(prms, var1_obs, var2_obs, pop_size)
  mod_dat <- p2q(state_p, epi_df, q_cap, c)
  q_adj <- mod_dat[[1]]
  epi_weekly <- mod_dat[[2]]
  neg_log_L <- -sum(dbinom(epi_weekly$k, size=epi_weekly$TT, prob=q_adj, log=TRUE))
  return(neg_log_L)
}

binom1_L <- function(prms, var_obs, epi_df, pop_size, q_cap, c = 1){ # single variable
  state_p <- SIRS1var_pred(prms, var_obs, pop_size)
  mod_dat <- p2q(state_p, epi_df, q_cap, c)
  q_adj <- mod_dat[[1]]
  epi_weekly <- mod_dat[[2]]
  neg_log_L <- -sum(dbinom(epi_weekly$k, size=epi_weekly$TT, prob=q_adj, log=TRUE))
  return(neg_log_L)
}
