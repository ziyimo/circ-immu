#!/usr/bin/env Rscript

library("deSolve")

param_bounds <- list() # boundaries of parameters, first element is "c"

# sinusoidal R0 model
R0_cos <- function(t, phi, R0base = 2, R0min = 1.2){ # a sinusoidal baseline model
  R0 <- (R0base - R0min)/2*cos(2*pi/364*(t - phi)) + (R0base + R0min)/2
  return(R0)
}
param_bounds[["cos"]] <- list(low = c(0, 1), high = c(1, 364))

# single var R0 models

R0_hum <- function(h_t, alpha, R0base = 2, R0min = 1.2){ # Baker et al. exponential decay
  R0 <- exp(alpha*h_t + log(R0base - R0min)) + R0min
  return(R0)
}
param_bounds[["hum"]] <- list(low = c(0, -300), high = c(1, 0))

R0_day <- R0_hum # alias for exponential decay modeling of daytime
param_bounds[["day"]] <- list(low = c(0, -100), high = c(1, 0))
  
# composite R0 models

R0_hd <- function(h_t, d_t, alpha_1, alpha_2, R0base = 2, R0min = 1.2){ # humidity + daytime
  R0 <- exp(alpha_1*h_t + alpha_2*d_t + log(R0base - R0min)) + R0min
  return(R0)
}
param_bounds[["hd"]] <- list(low = c(0, -300, -100), high = c(1, 0, 0))

R0_sd <- function(s_t, d_t, alpha_1, alpha_2, R0base = 2, R0min = 1.2){ # sunrise + daytime
  R0 <- exp(alpha_1*(1-s_t) + alpha_2*d_t + log(R0base - R0min)) + R0min
  return(R0)
}
param_bounds[["sd"]] <- list(low = c(0, -100, -100), high = c(1, 0, 0))

R0_hsd <- function(h_t, s_t, d_t, alpha_1, alpha_2, alpha_3, R0base = 2, R0min = 1.2){ # humidity + sunrise + daytime
  R0 <- exp(alpha_1*h_t + alpha_2*(1-s_t) + alpha_3*d_t + log(R0base - R0min)) + R0min
  return(R0)
}
param_bounds[["hsd"]] <- list(low = c(0, -300, -100, -100), high = c(1, 0, 0, 0))


######## Plot R0 models #######
if(FALSE){
  library(ggplot2)
  library(gridExtra)
  
  plot(seq(364), R0_cos(seq(364), 120))
  
  clean <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                 text = element_text(size=20))
  
  sun_range = seq(0, 1, by=0.01)
  day_range = seq(0, 1, by=0.01)
  cli_range = seq(0, 0.03, by=0.0003)
  
  ggplot() + 
    geom_line(data = data.frame("hum"= cli_range, "R0"= R0_hum(cli_range, -10)), 
              aes(x = hum, y=R0), color = "#E69F00") +
    geom_line(data = data.frame("hum"= cli_range, "R0"= R0_hum(cli_range, -100)), 
              aes(x = hum, y=R0), color = "#0072B2") +
    geom_line(data = data.frame("hum"= cli_range, "R0"= R0_hum(cli_range, -300)), 
              aes(x = hum, y=R0), color = "#009E73") + clean
  
  ggplot() + 
    geom_line(data = data.frame("day"= day_range, "R0"= R0_day(sun_range, -1)), 
              aes(x = day, y=R0), color = "#E69F00") +
    geom_line(data = data.frame("day"= day_range, "R0"= R0_day(sun_range, -10)), 
              aes(x = day, y=R0), color = "#0072B2") +
    geom_line(data = data.frame("day"= day_range, "R0"= R0_day(sun_range, -100)), 
              aes(x = day, y=R0), color = "#009E73") + clean
  
  #grid.arrange(p1, p2, ncol = 1)
  
  hum_day_df <- expand.grid(cli_range, day_range)
  colnames(hum_day_df) <- c("hum", "day")
  
  sun_day_df <- expand.grid(sun_range, day_range)
  colnames(sun_day_df) <- c("sun", "day")
  
  hum_day_df$R0 <- R0_hd(hum_day_df$hum, hum_day_df$day, -100, -1)
  p1 <- ggplot(hum_day_df, aes(hum, day, z = R0)) + geom_contour_filled() + clean
  
  hum_day_df$R0 <- R0_hd(hum_day_df$hum, hum_day_df$day, -30, -5)
  p2 <- ggplot(hum_day_df, aes(hum, day, z = R0)) + geom_contour_filled() + clean
  
  sun_day_df$R0 <- R0_sd(sun_day_df$sun, sun_day_df$day, -100, -3)
  p3 <- ggplot(sun_day_df, aes(sun, day, z = R0)) + geom_contour_filled() + clean
  
  sun_day_df$R0 <- R0_sd(sun_day_df$sun, sun_day_df$day, -100, -100)
  p4 <- ggplot(sun_day_df, aes(sun, day, z = R0)) + geom_contour_filled() + clean
  
  grid.arrange(p1, p2, p3, p4, ncol = 2)
  grid.arrange(p2, p3, ncol = 1)
}
###########################

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
