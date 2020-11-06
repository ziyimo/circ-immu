#!/usr/bin/env Rscript

library(deSolve)

# R0 dependency
R0_exp <- function(q_t, alpha, R0base, R0min){
  R0 = exp(alpha*q_t + log(R0base - R0min)) + R0min
  return(R0)
}

# SIRS model
SIRS_R0exp <- function(time, state, theta) {
  ## Parameters:
  varob <- theta[["varob"]] # observed variable, either sun or climate
  D <- theta[["D"]]
  L <- theta[["L"]]
  R0base <- theta[["R0max"]]
  R0min <- theta[["R0min"]]
  alpha <- theta[["alpha"]] # dependency, "a"
  
  # get humidity or sunrise at time t
  varobt <- varob[time]
  
  ## States
  S <- state["S"]
  I <- state["I"]
  R <- state["R"]
  # Total population
  N <- S + I + R
  
  # Ordinary differential equations
  beta = R0_exp(varobt, alpha, R0base, R0min)/D
  
  dS <- (R/L) - beta * S * I/N
  dI <- beta * S * I/N - (I/D)
  dR <- (I/D) - (R/L)
  
  return(list(c(dS, dI, dR)))
}

binom_L <- function(exp_para, var_observed, epi_df, census_pop, p_hat_range){
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
                     alpha = exp_para,
                     varob = rep(head(var_observed, 364), times = tot_years)[(offset+1):(364*tot_years)]
  )
  predictions <- as.data.frame(ode(xstart, times, SIRS_R0exp, SIRS_params))
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

state_code <- args[1]
#state_code <- "AZ"
prm_range <- seq(0, -360)

## Load population
state_pop <- read.delim("state_lv_data/state_pop.tsv")
census_pop <- state_pop$pop[state_pop$code==state_code]

## Load humidity data
all_state_hum <- read.csv("state_lv_data/state_humidity_2014_2018.csv")
all_state_hum <- all_state_hum[all_state_hum$X != "2016-02-29", ] # get rid of leap year Feb 29
climob <- all_state_hum[[state_code]] # multiple years
climob <- matrix(climob, nrow = length(climob)/365, ncol = 365, byrow = TRUE)
climob <- colMeans(climob) # 365 days

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
#plot(meanP_week$week, meanP_week$p_hat)
p_hat_max <- max(meanP_week$p_hat, na.rm = TRUE)
p_hat_min <- min(meanP_week$p_hat, na.rm = TRUE)

neg_ll <- rep(NA, length(prm_range))

sink(paste0(state_code, "_cliOptim.log"))
for (prm_idx in 1:length(prm_range)) {
  neg_ll[prm_idx] <- binom_L(prm_range[prm_idx], climob, epi_data, census_pop, c(p_hat_max, p_hat_min))
  cat(state_code, prm_range[prm_idx], neg_ll[prm_idx], "\n")
}
sink()

write.table(cbind(prm_range, neg_ll), file = paste0(state_code, "_cli_negLL.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE)