#!/usr/bin/env Rscript

library(ggplot2)
library(gridExtra)
library(deSolve)

#args <- commandArgs(trailingOnly=TRUE)

state_code <- "AZ"

## Load population
state_pop <- read.delim("state_lv_data/state_pop.tsv")
census_pop <- state_pop$pop[state_pop$code==state_code]

## Load sunrise data
all_state_sun <- read.csv("state_lv_data/state_daily_sunrise_2019.csv")
sunob <- all_state_sun[[state_code]]/720 # 365 days, scaled to maximum 720 = 12 hours

## Load humidity data
all_state_hum <- read.csv("state_lv_data/state_humidity_2014_2018.csv")
all_state_hum <- all_state_hum[all_state_hum$X != "2016-02-29", ] # get rid of leap year Feb 29
climob <- all_state_hum[[state_code]] # multiple years
climob <- matrix(climob, nrow = length(climob)/365, ncol = 365, byrow = TRUE)
climob <- colMeans(climob) # 365 days

p1 <- qplot(seq(365), -sunob)
p2 <- qplot(seq(365), climob)
grid.arrange(p1, p2, ncol = 1)

## Load flu data
epiob <- read.csv(paste0("state_lv_data/Flu_data/flu_epi_", state_code, ".csv")) # blank field automatically NA
epiob <- epiob[epiob$WEEK <= 52, ] # cap year to 52 weeks

# calculate relative date to the beginning of the 1st year in the data
y0 <- epiob$YEAR[1]
rel_date <- (epiob$YEAR-y0)*364 + epiob$WEEK*7 - 3
no_years <- ceiling(max(rel_date)/365)

k <- epiob$TOTAL.A + epiob$TOTAL.B
TT <- epiob$TOTAL.SPECIMENS
pi <- epiob$X.UNWEIGHTED.ILI/100
mask <- is.na(TT) | (TT == 0) # missing data: no. of test is 0 or NA

# calculate scaling factors
p_df <- data.frame("p_hat"=pi*k/TT, "week"=epiob$WEEK)
meanP_week <- aggregate(p_df, by=list(epiob$WEEK), FUN=mean, na.rm=TRUE)
qplot(meanP_week$week, meanP_week$p_hat)
p_hat_max <- max(meanP_week$p_hat, na.rm = TRUE)
p_hat_min <- min(meanP_week$p_hat, na.rm = TRUE)

# R0 dependency

R0_exp <- function(q_t, alpha, R0base, R0min){
  R0 = exp(alpha*q_t + log(R0base - R0min)) + R0min
  return(R0)
}

R0_sig <- function(q_t, k, q_0, R0base, R0min){
  R0 = (R0base-R0min)/(1+exp(-k*(q_t-q_0)))+R0min
  return(R0)
}

x_range = seq(0, 1, by=0.01)
p1 <- ggplot() + 
  geom_line(data = data.frame("X"= x_range, "R0"= R0_exp(x_range, -1, 2.5, 1.5)), 
            aes(x = X, y=R0), color = "blue") +
  geom_line(data = data.frame("X"= x_range, "R0"= R0_exp(x_range, -10, 2.5, 1.5)), 
            aes(x = X, y=R0), color = "red") +
  geom_line(data = data.frame("X"= x_range, "R0"= R0_exp(x_range, -50, 2.5, 1.5)), 
            aes(x = X, y=R0), color = "green")
p2 <- ggplot() + 
  geom_line(data = data.frame("X"= x_range, "R0"= R0_sig(x_range, 1, 0.5, 2.5, 1.5)), 
            aes(x = X, y=R0), color = "blue") +
  geom_line(data = data.frame("X"= x_range, "R0"= R0_sig(x_range, 10, 0.5, 2.5, 1.5)), 
            aes(x = X, y=R0), color = "red") +
  geom_line(data = data.frame("X"= x_range, "R0"= R0_sig(x_range, 100, 0.5, 2.5, 1.5)), 
            aes(x = X, y=R0), color = "green")
grid.arrange(p1, p2, ncol = 1)

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

## Test that steady state reached regardless of when the 1st infection is seeded
offset = 0
xstart = c(S = census_pop-1, I = 1, R = 0)
times = seq(1, 364 * 50, by = 1)[1:(364*50-offset)]

paras = list(D = 5, 
             L = 40*7, 
             R0max = 2.2,
             R0min = 1.2,
             k_step = 10,
             q_center = 0.5,
             varob = rep(head(sunob, 364), times = 50)[(offset+1):(364*50)]) # sunob or climob

predictions <- as.data.frame(ode(xstart, times, SIRS_R0sig, paras))
state_p <- predictions$I/census_pop

plot_range = seq((50-no_years)*364-offset+1, 50*364-offset)
#qplot(1:length(plot_range), state_p[plot_range])
x_range = seq(length(plot_range))
model_p = state_p[plot_range]

# Min-max scaling
#model_p <- (model_p-min(model_p))/(max(model_p) - min(model_p))*(p_hat_max - p_hat_min)+p_hat_min
scaler = (p_hat_max-p_hat_min)/(max(model_p)-min(model_p))
model_p <- model_p*scaler

ggplot() + 
  geom_line(data = data.frame("X"= x_range, "p"= model_p), 
            aes(x = X, y=p), color = "blue") +
  geom_point(data = data.frame("X"= rel_date, "p_hat"= p_df$p_hat), 
            aes(x = X, y=p_hat), color = "red")
## PASSED


ggplot() + 
  geom_line(data = data.frame("X"= rel_date, "pi"= pi), 
            aes(x = X, y=pi), color = "blue") +
  geom_point(data = data.frame("X"= rel_date, "p_hat"= p_df$p_hat), 
             aes(x = X, y=p_hat), color = "red") +
  geom_point(data = data.frame("X"= rel_date, "p"= p_weekly), 
             aes(x = X, y=p), color = "orange")
# model prediction on particular days
p_weekly <- model_p[rel_date]
# Drop data point violating model assumption and na fields
bad <- p_weekly > pi
drop <- bad | mask


p_weekly <- p_weekly[!drop]
pi_weekly <- pi[!drop]
k_weekly <- k[!drop]
TT_weekly <- TT[!drop]
q <- p_weekly/pi_weekly

ggplot() + 
  geom_point(data = data.frame("X"= rel_date[!drop], "p_hat"= k_weekly/TT_weekly), 
            aes(x = X, y=p_hat), color = "blue") +
  geom_line(data = data.frame("X"= rel_date[!drop], "q"= q), 
             aes(x = X, y=q), color = "orange")

# Calc. binomial probability
k_T_q <- cbind(k_weekly, TT_weekly, q)
#neg_log_L <- -sum(apply(k_T_q, 1, function(xx) dbinom(xx[1], size=xx[2], prob=xx[3], log=TRUE)))
neg_log_L <- -sum(dbinom(k_weekly, size=TT_weekly, prob=q, log=TRUE))
