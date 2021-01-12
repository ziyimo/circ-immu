#!/usr/bin/env Rscript

suppressWarnings(suppressMessages(library("deSolve")))

# Boundaries for inital infected/recovered
# Note inf are multiplied by 100, rec are multiplied by 1000
inf_low <- 1
inf_high <- 400000
rec_low <- 1
rec_high <- 4000000

param_bounds <- list() # boundaries of parameters, first element is "c"

# sinusoidal R0 model
R0_cos <- function(t, phi, R0base = 2.5, R0min = 0.8){ # a sinusoidal baseline model
  R0 <- (R0base - R0min)/2*cos(2*pi/364*(t - phi)) + (R0base + R0min)/2
  return(R0)
}
param_bounds[["cos"]] <- list(low = c(inf_low, rec_low, 1), high = c(inf_high, rec_high, 364))

# single var R0 models
R0_bell <- function(q_t, alpha, q_0, R0base = 2.5, R0min = 0.8){
  #R0 <- exp(-((q_t-mu)^2)/(2*sigma^2) + log(R0base - R0min)) + R0min
  R0 <- exp(alpha*(q_t-q_0)^2 + log(R0base - R0min)) + R0min # reparameterized
  return(R0)
}
param_bounds[["bell"]] <- list(low = c(inf_low, rec_low, -1000, 0), high = c(inf_high, rec_high, 1, 0, 1))

R0_hum <- R0_bell
param_bounds[["hum"]] <- list(low = c(inf_low, rec_low, -1e5, 0), high = c(inf_high, rec_high, 0, 0.021))
# > max(all_state_hum)
# [1] 0.02134766

R0_day <- R0_bell
param_bounds[["day"]] <- list(low = c(inf_low, rec_low, -100, 0), high = c(inf_high, rec_high, 0, 1))

# composite R0 models

R0_linGG <- function(q_t, r_t, alpha_1, q_0, alpha_2, r_0, R0base = 2.5, R0min = 0.8){
  R0 <- exp(alpha_1*(q_t-q_0)^2 + alpha_2*(r_t-r_0)^2 + log(R0base - R0min)) + R0min
  return(R0)
}

R0_hd <- R0_linGG
param_bounds[["hd"]] <- list(low = c(inf_low, rec_low, -1e5, 0, -100, 0), high = c(inf_high, rec_high, 0, 0.021, 0, 1))

R0_sd <- R0_linGG
param_bounds[["sd"]] <- list(low = c(inf_low, rec_low, -100, 0, -100, 0), high = c(inf_high, rec_high, 0, 1, 0, 1))

R0_hsd <- function(h_t, s_t, d_t, alpha_1, h_0, alpha_2, s_0, alpha_3, d_0, R0base = 2.5, R0min = 0.8){
  R0 <- exp(alpha_1*(h_t-h_0)^2 + alpha_2*(s_t-s_0)^2 + alpha_3*(d_t-d_0)^2 + log(R0base - R0min)) + R0min
  return(R0)
}
param_bounds[["hsd"]] <- list(low = c(inf_low, rec_low, -1e5, 0, -100, 0, -100, 0), high = c(inf_high, rec_high, 0, 0.021, 0, 1, 0, 1))


###########################

SEIRD_R0 <- function(time, state, theta){
  ## Parameters:
  R0_vec <- theta[["R0"]] # vector pre-calculated R0 values
  delta <- theta[["delta"]]
  alpha <- theta[["alpha"]]
  gamma <- theta[["gamma"]] # infectious period
  rho <- theta[["rho"]]

  R0_t <- R0_vec[time] # get R0 at time t
  
  ## States:
  S <- state["S"] # suceptible
  E <- state["E"] # suceptible
  I <- state["I"] # infected
  R <- state["R"] # recovered
  N <- S + E + I + R
  
  # Include deaths
  #D <-  state["D"]
  #N <- S + E + I + R + D
  
  # Ordinary differential equations
  beta = R0_t*gamma
  dS = -beta * S * I / N
  dE = beta * S * I / N - delta * E  
  dI = delta * E - gamma * I
  dR =  gamma * I
  
  return(list(c(dS, dE, dI, dR)))
  
  # Include deaths
  #dI = delta * E - (1 - alpha) * gamma * I - alpha * rho * I
  #dR = (1 - alpha) * gamma * I 
  #dD = alpha * rho * I
  #return(list(c(dS, dE, dI, dR, dD)))
}

###########################

SEIRDvar_pred <- function(R0_func, R0_params, var_ls, census_pop){
  # R0_func: string of function name
  # R0_params: vector of parameters
  # var_ls: list of vector(s) of variables (sunrise, humidity)
  # Set initial infected and recovered
  inf_init <- R0_params[1]
  inf_recov <- R0_params[2]
  # fixed initial infected and recovered (LA, 1. June)
  #inf_init <- 107900
  #inf_recov <- 625614
  
  # extract the params used to calculate R0 by removing first two
  R0_params <- R0_params[-c(1,2)]
  c_total_days <- length(var_ls[[1]]) # specify when to seed the single infection
  c_times <- seq(1, c_total_days, by = 1)[1:(c_total_days)]
  xstart <- c(S = census_pop-(inf_init+inf_recov),E = 0, I = inf_init, R = inf_recov) # use actual state population
  
  annual_R0 <- do.call(R0_func, c(var_ls, as.list(R0_params)))
    # some hard-coded parameters
  paras = list(alpha = 0.01,  # 1% death rate
               #L = 280, # duration of immunity - not used for COVID
               rho = 1 / 9,  # 9 days from infection until death
               delta = 1 / 5,  # incubation period of five days
               gamma = 1 / 5, # infection lasts 5 days
               R0 = annual_R0) # pre-calculated R0 values
  
  predictions <- as.data.frame(ode(xstart, c_times, SEIRD_R0, paras)) # run SIRS model
  #raw_p <- predictions$I/census_pop
  return(predictions)
}

## Wrappers for SEIRD likelihood function
binom_seird <- function(prms, R0_model, var_obs, obs_deaths, pop_size){
  
  # set start date of model and filter var_obs accordingly
  cityname <- head(obs_deaths$Name,1)
  cityname <- sub(" ", "",cityname)
  obs_deaths <- obs_deaths$Cases_New
  lambda <- 1
  alpha <- 0.01
  alpha_vec <- rep(alpha,length(obs_deaths))
  rho <- 0.1
  # predicted infections
  predictions <- SEIRDvar_pred(R0_model, prms, var_obs, pop_size)
  
  # Log intermediate resilts and prediction table
  #outlog <- paste(paste(prms,collapse = '_'),"_log_pred.csv", sep = '')
  #write.table(predictions, file = outlog, sep = ",", quote = FALSE)
  
  # Link infections to the deaths 9 days later for fitting
  #obs_deaths<- head(obs_deaths,-9)
  #predictions <- tail(predictions, -9)
  
  # Average time from infection to death is 10 days
  rho <- 10
  alpha <- 0.01
  alpha_vec <- rep(alpha,length(obs_deaths))
  predictions$obs_I <- obs_deaths
  # cap minimum predicted infection at value of obs_deaths*100*rho
  predictions$adj_I <- ifelse(predictions$obs_I > predictions$I,predictions$obs_I*100*rho, predictions$I)
  predictions$obs_I <- predictions$obs_I*100*rho
  lambda = 1
  # Penalize capped infection counts
  neg_log_L <- -sum(dbinom(as.integer(obs_deaths), size = as.integer(0.1*predictions$adj_I), prob = alpha_vec)) + lambda*sum((predictions$adj_I-predictions$I)/predictions$adj_I)
  return(neg_log_L)
}

binom_seird_p <- function(prms, R0_model, var_obs, obs_deaths, pop_size){
  # function is a copy of binom_seird but only outputs plot
  suppressWarnings(suppressMessages(require(cowplot)))
  suppressWarnings(suppressMessages(require(reshape2)))
  suppressWarnings(suppressMessages(require(ggplot2)))
  # set start date of model and filter var_obs accordingly
  cityname <- head(obs_deaths$Name,1)
  cityname <- sub(" ", "",cityname)
  obs_deaths <- obs_deaths$Cases_New
  lambda <- 1
  alpha <- 0.01
  alpha_vec <- rep(alpha,length(obs_deaths))
  rho <- 0.1
  # predicted infections
  predictions <- SEIRDvar_pred(R0_model, prms, var_obs, pop_size)
  
  # Log intermediate resilts and prediction table
  #outlog <- paste(paste(prms,collapse = '_'),"_log_pred.csv", sep = '')
  #write.table(predictions, file = outlog, sep = ",", quote = FALSE)
  
  # Link infections to the deaths 9 days later for fitting
  #obs_deaths<- head(obs_deaths,-9)
  #predictions <- tail(predictions, -9)
  
  # Average time from infection to death is 10 days
  rho <- 10
  alpha <- 0.01
  alpha_vec <- rep(alpha,length(obs_deaths))
  predictions$obs_I <- obs_deaths
  # cap minimum predicted infection at value of obs_deaths*100*rho
  predictions$adj_I <- ifelse(predictions$obs_I > predictions$I,predictions$obs_I*100*rho, predictions$I)
  predictions$obs_I <- predictions$obs_I*100*rho
  lambda = 1
  # Penalize capped infection counts
  neg_log_L <- -sum(dbinom(as.integer(obs_deaths), size = as.integer(0.1*predictions$adj_I), prob = alpha_vec)) + lambda*sum((predictions$adj_I-predictions$I)/predictions$adj_I)
  # Plot Fitting results and R0 estimates
  
  xm <- melt(predictions , id.vars=c("time"))
  xm <- subset(xm, variable == "adj_I"| variable == "I" | variable == "obs_I")
  p1 <- ggplot(xm, aes(x=time,color=variable,y=value)) + geom_line(position=position_jitter(w=1, h=1),size=2,alpha = 0.4 ) +
    ylab("People") +
    xlab("Day")+
    theme(legend.position="top",
          plot.title = element_text(size = 12),
          legend.title=element_text(size=8),
          legend.text=element_text(size=8),
          axis.text=element_text(size=12),
          axis.title=element_text(size=12)) +
    ggtitle(paste("-loglik=",neg_log_L, sep = " "))
  
  R0_params <- prms[-c(1,2)]
  annual_R0 <- do.call(R0_model, c(var_obs, as.list(R0_params)))
  annual_R0 <- as.data.frame(annual_R0)
  p2 <- ggplot(annual_R0, aes(x=seq_along(annual_R0), y=annual_R0)) + geom_point()+
          theme(plot.title = element_text(size = 12),
          axis.text=element_text(size=12),
          axis.title=element_text(size=12))
  outplot <- paste("fit_results/plot_",cityname,"_",R0_model,"_",paste(prms,collapse = '_'),"_",neg_log_L,".pdf", sep = '')
  plot_grid(p1, p2,align = 'v',ncol = 1)
  ggsave(outplot)
}
