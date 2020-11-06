#!/usr/bin/env Rscript

library(deSolve)
library(DEoptim)

R0_sig <- function(q_t, k, q_0, R0base, R0min){
  R0 = (R0base-R0min)/(1+exp(-k*(q_t-q_0)))+R0min
  return(R0)
}

# SIRS model
SIRS_R0sig <- function(time, state, theta) {
  ## Parameters:
  varob <- theta[["varob"]] # observed variable, either sun or climate
  D <- theta[["D"]]
  L <- theta[["L"]]
  R0base <- theta[["R0max"]]
  R0min <- theta[["R0min"]]
  k <- theta[["k_step"]] # steepness
  q_0 <- theta[["q_center"]] # center of the sigmoid
  
  # get humidity or sunrise at time t
  varobt <- varob[time]
  
  ## States
  S <- state["S"]
  I <- state["I"]
  R <- state["R"]
  # Total population
  N <- S + I + R
  
  # Ordinary differential equations
  beta = R0_sig(varobt, k, q_0, R0base, R0min)/D
  
  dS <- (R/L) -beta * S * I/N
  dI <- beta * S * I/N - (I/D)
  dR <- (I/D) - (R/L)
  
  return(list(c(dS, dI, dR)))
}

binom_likelihood <- function(sig_params, var_observed, epi_df, census_pop, p_hat_range){
  ## sig_params = c(k_step, q_center) ##
  ## var_observed: climate/sunrise data (annual) ##
  ## epi_df: dataframe of epidemiological data ##
  ## p_hat_range = c(p_hat_max, p_hat_min) ##
  
  offset = 0 # specify when to seed the single infection
  tot_years = 50
  xstart = c(S = census_pop-1, I = 1, R = 0) # use actual state population
  times = seq(1, 364 * tot_years, by = 1)[1:(364*tot_years-offset)]
  
  # hard code some parameters
  SIRS_params = list(D = 5,
                     L = 40*7,
                     R0max = 2,
                     R0min = 1.2,
                     k_step = sig_params[1],
                     q_center = sig_params[2],
                     varob = rep(head(var_observed, 364), times = tot_years)[(offset+1):(364*tot_years)]
                     )
  predictions <- as.data.frame(ode(xstart, times, SIRS_R0sig, SIRS_params))
  state_p <- predictions$I/census_pop

  data_years <- ceiling(max(epi_df$rel_date)/364)
  model_range = seq((tot_years-data_years)*364-offset+1, tot_years*364-offset)
  model_p <- state_p[model_range] # take only the number of years corresponding to data
  
  # apply min-max scaling for model predicted probability
  scaler <- sum(p_hat_range)/(max(model_p)+min(model_p))
  model_p <- model_p*scaler
  
  # model prediction on particular days
  p_weekly <- model_p[rel_date]
  # Drop data point violating model assumption and na fields
  bad <- p_weekly > epi_df$pi
  mask <- is.na(epi_df$TT) | (epi_df$TT == 0)
  drop <- bad | mask
  epi_weekly <- epi_df[!drop, ]
  p_weekly <- p_weekly[!drop]
  q <- p_weekly/epi_weekly$pi
  
  # Calculate negative log likelihood
  neg_log_L <- -sum(dbinom(epi_weekly$k, size=epi_weekly$TT, prob=q, log=TRUE))
  return(neg_log_L)  
}

args <- commandArgs(trailingOnly=TRUE)
# Data preprocessing

state_code <- args[1]
#state_code <- "TX"

## Load population
state_pop <- read.delim("state_lv_data/state_pop.tsv")
census_pop <- state_pop$pop[state_pop$code==state_code]

## Load sunrise data
all_state_sun <- read.csv("state_lv_data/state_daily_sunrise_2019.csv")
sunob <- all_state_sun[[state_code]]/720 # 365 days, scaled to maximum 720 = 12 hours

## Load flu data
epiob <- read.csv(paste0("state_lv_data/Flu_data/flu_epi_", state_code, ".csv")) # blank field automatically NA
epiob <- epiob[epiob$WEEK <= 52, ] # cap year to 52 weeks

# calculate relative date to the beginning of the 1st year in the data
y0 <- epiob$YEAR[1]
rel_date <- (epiob$YEAR-y0)*364 + epiob$WEEK*7 - 3

k <- epiob$TOTAL.A + epiob$TOTAL.B
TT <- epiob$TOTAL.SPECIMENS
pi <- epiob$X.UNWEIGHTED.ILI/100

epi_data <- data.frame("rel_date"=rel_date, "k"=k, "TT"=TT, "pi"=pi)

# calculate scaling factors
p_df <- data.frame("p_hat"=pi*k/TT, "week"=epiob$WEEK)
meanP_week <- aggregate(p_df, by=list(epiob$WEEK), FUN=mean, na.rm=TRUE)
#qplot(meanP_week$week, meanP_week$p_hat)
p_hat_max <- max(meanP_week$p_hat, na.rm = TRUE)
p_hat_min <- min(meanP_week$p_hat, na.rm = TRUE)

## Optimization

sink(paste0(state_code, "_sunOptim.log"))
evo_optim <- DEoptim(binom_likelihood, lower=c(0, 0), upper=c(100, 1),
        control=DEoptim.control(trace = 5, reltol = 1e-5),
        var_observed = sunob, epi_df = epi_data, census_pop = census_pop, p_hat_range = c(p_hat_max, p_hat_min))

saveRDS(evo_optim, file = paste0(state_code, "_DEoptim.rds"))

NM_optim <- optim(evo_optim[["optim"]][["bestmem"]], binom_likelihood,
      var_observed = sunob, epi_df = epi_data, census_pop = census_pop, p_hat_range = c(p_hat_max, p_hat_min),
      control = list(trace = 1), method = "Nelder-Mead")

saveRDS(NM_optim, file = paste0(state_code, "_NMoptim.rds"))
sink()