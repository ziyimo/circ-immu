#!/usr/bin/env Rscript

suppressWarnings(suppressMessages(library("deSolve")))

param_bounds <- list() # boundaries of parameters, first element is "c"

# sinusoidal R0 model
R0_cos <- function(t, phi, R0base = 2, R0min = 1.2){ # a sinusoidal baseline model
  R0 <- (R0base - R0min)/2*cos(2*pi/364*(t - phi)) + (R0base + R0min)/2
  return(R0)
}

# single var R0 models
R0_hum <- function(h_t, alpha, R0base = 2, R0min = 1.2){ # Baker et al. exponential decay
  R0 <- exp(alpha*h_t + log(R0base - R0min)) + R0min
  return(R0)
}

R0_day <- function(d_t, alpha, d_0, R0base = 2, R0min = 1.2){
  R0 <- exp(alpha*(d_t-d_0)^2 + log(R0base - R0min)) + R0min
  return(R0)
}

R0_sun <- R0_day

# composite R0 models
R0_hd <- function(h_t, d_t, alpha_1, alpha_2, d_0, R0base = 2, R0min = 1.2){ # humidity + daytime
  R0 <- exp(alpha_1*h_t + alpha_2*(d_t-d_0)^2 + log(R0base - R0min)) + R0min
  return(R0)
}

R0_hs <- R0_hd

R0_sd <- function(s_t, d_t, alpha_1, s_0, alpha_2, d_0, R0base = 2, R0min = 1.2){
  R0 <- exp(alpha_1*(s_t-s_0)^2 + alpha_2*(d_t-d_0)^2 + log(R0base - R0min)) + R0min
  return(R0)
}

R0_hsd <- function(h_t, s_t, d_t, alpha_1, alpha_2, s_0, alpha_3, d_0, R0base = 2, R0min = 1.2){
  R0 <- exp(alpha_1*h_t + alpha_2*(s_t-s_0)^2 + alpha_3*(d_t-d_0)^2 + log(R0base - R0min)) + R0min
  return(R0)
}


#param_bounds[["cos"]] <- list(low = c(inf_low, rec_low,alpha_low, 1,R0range_low,R0min_low), high = c(inf_high, rec_high,alpha_high, 364,R0range_high,R0min_high))
#param_bounds[["bell"]] <- list(low = c(inf_low, rec_low, alpha_low, -1000, 0, R0range_low,R0min_low), high = c(inf_high, rec_high, alpha_high, 1, 0, 1,R0range_high,R0min_high))
#param_bounds[["day"]] <- list(low = c(inf_low, rec_low,alpha_low, -100, 0,R0range_low,R0min_low), high = c(inf_high, rec_high, alpha_high, 0, 1,R0range_high,R0min_high))
#param_bounds[["hd"]] <- list(low = c(inf_low, rec_low,alpha_low, -1e5, 0, -100, 0,R0range_low,R0min_low), high = c(inf_high, rec_high,alpha_high, 0, 0.021, 0, 1,R0range_high,R0min_high))
#param_bounds[["sd"]] <- list(low = c(inf_low, rec_low,alpha_low, -100, 0, -100, 0,R0range_low,R0min_low), high = c(inf_high, rec_high,alpha_high, 0, 1, 0, 1,R0range_high,R0min_high))
#param_bounds[["hsd"]] <- list(low = c(inf_low, rec_low,alpha_low, -1e5, 0, -100, 0, -100, 0,R0range_low,R0min_low), high = c(inf_high, rec_high,alpha_high, 0, 0.021, 0, 1, 0, 1,R0range_high,R0min_high))

param_bounds[["cos"]] <- list(low = c(1), high = c(364))
param_bounds[["hum"]] <- list(low = c(-300), high = c(0))
param_bounds[["sun"]] <- list(low = c(-500, 0), high = c(0, 1))
param_bounds[["hs"]] <- list(low = c(-300, -500, 0), high = c(0, 0, 1))

param_bounds[["day"]] <- list(low = c(-500, 0), high = c(0, 1))
param_bounds[["hd"]] <- list(low = c(-300, -500, 0), high = c(0, 0, 1))
param_bounds[["sd"]] <- list(low = c(-500, 0, -500, 0), high = c(0, 1, 0, 1))
param_bounds[["hsd"]] <- list(low = c(-300, -500, 0, -500, 0), high = c(0, 0, 1, 0, 1))


###########################

SEIRD_R0 <- function(time, state, theta){
  ## Parameters:
  R0_vec <- theta[["R0"]] # vector pre-calculated R0 values
  delta <- theta[["delta"]]
  alpha <- theta[["alpha"]]
  gamma <- theta[["gamma"]] # infectious period
  rho <- theta[["rho"]]
  psi <- theta[["psi"]]
  
  R0_t <- R0_vec[time] # get R0 at time t
  
  ## States:
  S <- state["S"] # suceptible
  E <- state["E"] # suceptible
  I <- state["I"] # infected
  G <- state["G"] # resolving
  R <- state["R"] # recovered
  N <- S + E + I + G + R
  
  # Ordinary differential equations
  beta = R0_t*gamma
  dS = -beta * S * I / N
  dE = beta * S * I / N - delta * E  
  dI = delta * E - gamma * I 
  dG =  gamma * I - G*psi
  dR =  G*psi
  return(list(c(dS, dE, dI, dG, dR)))
  
}

###########################

SEIRDvar_pred <- function(R0_func, R0_params, var_ls, census_pop,obs_deaths){
  # R0_func: string of function name
  # R0_params: vector of parameters
  # var_ls: list of vector(s) of variables (sunrise, humidity)
  # Set initial infected and recovered
  inf_init <- R0_params[1] * census_pop
  inf_recov <- R0_params[2] * census_pop
  param_alpha = R0_params[3]
  R0_min <- R0_params[4]
  R0_range <- R0_params[5]
  # extract the params used to calculate R0 by removing first five
  R0_params <- R0_params[-c(1,2,3,4,5)]
  
  c_total_days <- length(obs_deaths) # specify when to seed the single infection
  c_times <- seq(1, c_total_days, by = 1)[1:(c_total_days)]
  xstart <- c(S = census_pop-(inf_init+inf_recov),E = 0, I = inf_init, G = 0, R = inf_recov) # use actual state population
  annual_R0 <- do.call(R0_func, c(var_ls, as.list(c(R0_params, R0_min+R0_range, R0_min))))
  if(R0_func == "R0_cos"){
	  annual_R0 <- annual_R0[1:c_total_days]
  }

  # some hard-coded parameters
  paras = list(alpha = param_alpha,  # 1% death rate
               #L = 280, # duration of immunity - not used for COVID
               rho = 1 / 9,  # 9 days from infection until death
               delta = 1 / 5,  # incubation period of five days
               gamma = 1 / 5, # infection lasts 5 days
               psi = 1 / 10, # resolving takes 10 days
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
  alpha <- prms[3]
  alpha_vec <- rep(alpha,length(obs_deaths))
  rho <- 0.1
  # predicted infections
  predictions <- SEIRDvar_pred(R0_model, prms, var_obs, pop_size,obs_deaths)
  alpha_vec <- rep(alpha,length(obs_deaths))
  predictions$obs_D <- obs_deaths
  # cap minimum predicted infection at deaths+1
  predictions$adj_G <- ifelse(predictions$obs_D > predictions$G,predictions$obs_D+1, predictions$G)
  lambda = 1
  # Penalize capped infection counts
  neg_log_L <- -sum(dbinom(as.integer(obs_deaths), size = as.integer(predictions$adj_G), prob = alpha_vec, log = TRUE)) + lambda*sum(predictions$adj_G-predictions$G)
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
  #alpha <- 0.01
  alpha <- prms[3]
  alpha_vec <- rep(alpha,length(obs_deaths))
  rho <- 0.1
  # predicted infections
  predictions <- SEIRDvar_pred(R0_model, prms, var_obs, pop_size,obs_deaths)
  
  # Average time from infection to death is 10 days
  alpha_vec <- rep(alpha,length(obs_deaths))
  predictions$obs_D <- obs_deaths
  # cap minimum predicted infection at value of obs_deaths*100*rho
  predictions$adj_G <- ifelse(predictions$obs_D > predictions$G,predictions$obs_D+1, predictions$G)
  predictions$obs_G <- predictions$obs_D/alpha
  lambda = 1
  # Penalize capped infection counts
  neg_log_L <- -sum(dbinom(as.integer(obs_deaths), size = as.integer(predictions$adj_G), prob = alpha_vec, log = TRUE)) + lambda*sum(predictions$adj_G-predictions$G)
  # Plot Fitting results and R0 estimates
  xm <- melt(predictions , id.vars=c("time"))
  xm <- subset(xm, variable == "adj_G"| variable == "G" | variable == "obs_G")
  p1 <- ggplot(xm, aes(x=time,color=variable,y=value)) + geom_line(position=position_jitter(w=1, h=1),size=2,alpha = 0.4 ) +
    ylab("People") +
    xlab("Day")+
    theme(legend.position="top",
          plot.title = element_text(size = 12),
          legend.title=element_text(size=8),
          legend.text=element_text(size=8),
          axis.text=element_text(size=12),
          axis.title=element_text(size=12)) +
    ggtitle(paste("-loglik=-",neg_log_L, sep = " "))
  
  R0_min <- prms[4]
  R0_range <- prms[5]
  R0_params <- prms[-c(1,2,3,4,5)]
  annual_R0 <- do.call(R0_model, c(var_obs, as.list(c(R0_params, R0_min+R0_range, R0_min))))
  if(R0_model == "R0_cos"){
	  annual_R0 <- annual_R0[1:length(obs_deaths)]
  }
  # Add R0 to predictions for log
  annual_R0 <- as.data.frame(annual_R0)
  predictions$R0 <- annual_R0$annual_R0
  p2 <- ggplot(annual_R0, aes(x=seq_along(annual_R0), y=annual_R0)) + geom_point()+
          theme(plot.title = element_text(size = 12),
          axis.text=element_text(size=12),
          axis.title=element_text(size=12))
  prms <- round(prms,digits = 3)
  outplot <- paste("fit_results/",cityname,"_",R0_model,"_",paste(prms,collapse = '_'),"_",neg_log_L,".pdf", sep = '')
  plot_grid(p1, p2,align = 'v',ncol = 1)
  ggsave(outplot)
  outlog <- paste("fit_results/",cityname,"_",R0_model,"_",paste(prms,collapse = '_'),"_log_pred.csv", sep = '')
  write.table(predictions, file = outlog, sep = ",", quote = FALSE)
}

