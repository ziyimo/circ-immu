#!/usr/bin/env Rscript

library(deSolve)

# R0 models
R0_sig <- function(q_t, k, b, R0base, R0min){ # single variable sigmoid
  R0 = (R0base-R0min)/(1+exp(-k*q_t+b))+R0min
  return(R0)
}

R0_sig2 <- function(q_t, r_t, k1, k2, b, R0base, R0min){ # two variable sigmoid
  R0 = (R0base-R0min)/(1+exp(-k1*q_t-k2*r_t+b))+R0min
  return(R0)
}

######## Plot R0 models #######
if(FALSE){
  library(ggplot2)
  library(gridExtra)
  
  clean <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                 text = element_text(size=20))
  
  sun_range = seq(0, 1, by=0.01)
  cli_range = seq(0, 0.03, by=0.0003)
  
  p1 <- ggplot() + 
    geom_line(data = data.frame("hum"= cli_range, "R0"= R0_sig(cli_range, -100, -2.5, 2, 1.2)), 
              aes(x = hum, y=R0), color = "#E69F00") +
    geom_line(data = data.frame("hum"= cli_range, "R0"= R0_sig(cli_range, -1500, -15, 2, 1.2)), 
              aes(x = hum, y=R0), color = "#0072B2") +
    geom_line(data = data.frame("hum"= cli_range, "R0"= R0_sig(cli_range, -300, -3, 2, 1.2)), 
              aes(x = hum, y=R0), color = "#009E73") + clean
  
  p2 <- ggplot() + 
    geom_line(data = data.frame("sun"= sun_range, "R0"= R0_sig(sun_range, 100, 70, 2, 1.2)), 
              aes(x = sun, y=R0), color = "#E69F00") +
    geom_line(data = data.frame("sun"= sun_range, "R0"= R0_sig(sun_range, 10, 5, 2, 1.2)), 
              aes(x = sun, y=R0), color = "#0072B2") +
    geom_line(data = data.frame("sun"= sun_range, "R0"= R0_sig(sun_range, 30, 15, 2, 1.2)), 
              aes(x = sun, y=R0), color = "#009E73") + clean
  grid.arrange(p1, p2, ncol = 1) 
  
  sun_cli_df <- expand.grid(sun_range, cli_range)
  colnames(sun_cli_df) <- c("sun", "hum")
  
  sun_cli_df$R0 <- R0_sig2(sun_cli_df$sun, sun_cli_df$hum, 3, -100, 0.5, 2, 1.2)
  p1 <- ggplot(sun_cli_df, aes(sun, hum, z = R0)) + geom_contour_filled() + clean
  
  sun_cli_df$R0 <- R0_sig2(sun_cli_df$sun, sun_cli_df$hum, 10, -100, 4, 2, 1.2)
  p2 <- ggplot(sun_cli_df, aes(sun, hum, z = R0)) + geom_contour_filled() + clean
  
  sun_cli_df$R0 <- R0_sig2(sun_cli_df$sun, sun_cli_df$hum, 3, -300, -1.5, 2, 1.2)
  p3 <- ggplot(sun_cli_df, aes(sun, hum, z = R0)) + geom_contour_filled() + clean
  
  sun_cli_df$R0 <- R0_sig2(sun_cli_df$sun, sun_cli_df$hum, 10, -300, 2, 2, 1.2)
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

# Global Var (for efficiency)
offset <- 0 # specify when to seed the single infection
tot_years <- 50 # no. of years to run SIRS model till steady state
times <- seq(1, 364 * tot_years, by = 1)[1:(364*tot_years-offset)]

# single variable model prediction
SIRS1var_pred <- function(sig_params, var1, census_pop){
  k_step <- sig_params[1]
  b_intcpt <- sig_params[2]
  
  xstart <- c(S = census_pop-1, I = 1, R = 0) # use actual state population
  annual_R0 <- R0_sig(var1, k_step, b_intcpt, R0base = 2, R0min = 1.2)
  # some hard-coded parameters
  paras = list(D = 5, 
               L = 40*7, 
               R0 = rep(head(annual_R0, 364), times = tot_years)[(offset+1):(364*tot_years)]) # pre-calculated R0 values
  
  predictions <- as.data.frame(ode(xstart, times, SIRS_R0sig, paras)) # run SIRS model
  raw_p <- predictions$I/census_pop
  return(raw_p)
  }

# 2 variable model prediction
SIRS2var_pred <- function(sig2_params, var1, var2, census_pop){
  k1_step <- sig2_params[1]
  k2_step <- sig2_params[2]
  b_intcpt <- sig2_params[3]
  
  xstart <- c(S = census_pop-1, I = 1, R = 0) # use actual state population
  annual_R0 <- R0_sig2(var1, var2, k1_step, k2_step, b_intcpt, R0base = 2, R0min = 1.2)
  # some hard-coded parameters
  paras = list(D = 5, 
               L = 40*7, 
               R0 = rep(head(annual_R0, 364), times = tot_years)[(offset+1):(364*tot_years)]) # pre-calculated R0 values
  
  predictions <- as.data.frame(ode(xstart, times, SIRS_R0sig, paras)) # run SIRS model
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
binom1_L <- function(prms, var_obs, epi_df, pop_size, q_cap, c = 1, lambda = 0){ # single variable
  state_p <- SIRS1var_pred(prms, var_obs, pop_size)
  state_q <- p2q(state_p, epi_df, c)
  q_adj <- pmin(state_q, q_cap) # cap the value of q, this is q'
  neg_log_L <- -sum(dbinom(epi_df$k, size=epi_df$TT, prob=q_adj, log=TRUE)) + lambda*sum(state_q-q_adj)
  return(neg_log_L)
}

## Wrappers for penalized likelihood function, includes `c` for fitting
binom2_Lp <- function(prms, var1_obs, var2_obs, epi_df, pop_size, q_cap, lambda = 0){ # 2 variable
  state_p <- SIRS2var_pred(prms[1:3], var1_obs, var2_obs, pop_size)
  state_q <- p2q(state_p, epi_df, prms[4])
  q_adj <- pmin(state_q, q_cap) # cap the value of q, this is q'
  neg_log_L <- -sum(dbinom(epi_df$k, size=epi_df$TT, prob=q_adj, log=TRUE)) + lambda*sum(state_q-q_adj)
  return(neg_log_L)
}

binom1_Lp <- function(prms, var_obs, epi_df, pop_size, q_cap, lambda = 0){ # single variable
  state_p <- SIRS1var_pred(prms[1:2], var_obs, pop_size)
  state_q <- p2q(state_p, epi_df, prms[3])
  q_adj <- pmin(state_q, q_cap) # cap the value of q, this is q'
  neg_log_L <- -sum(dbinom(epi_df$k, size=epi_df$TT, prob=q_adj, log=TRUE)) + lambda*sum(state_q-q_adj)
  return(neg_log_L)
}