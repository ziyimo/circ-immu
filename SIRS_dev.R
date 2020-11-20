#!/usr/bin/env Rscript

library(ggplot2)
library(gridExtra)

#library(deSolve)
library(binom)

source("SIRS_model.R")

#args <- commandArgs(trailingOnly=TRUE)

state_code <- "NY"

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
#no_years <- ceiling(max(rel_date)/365)

k <- epiob$TOTAL.A + epiob$TOTAL.B
TT <- epiob$TOTAL.SPECIMENS
pi <- epiob$X.UNWEIGHTED.ILI/100
mask <- is.na(TT) | (TT == 0) | (pi == 0) # missing data: no. of test is 0 or NA, or no sympotomatic patient
cat(sum(mask), "weeks masked", "\n")

epi_data <- data.frame("rel_date"=rel_date, "k"=k, "TT"=TT, "pi"=pi)
epi_data <- epi_data[!mask, ]

# # calculate scaling factors
# p_df <- data.frame("p_hat"=pi*k/TT, "week"=epiob$WEEK)
# meanP_week <- aggregate(p_df, by=list(epiob$WEEK), FUN=mean, na.rm=TRUE)
# qplot(meanP_week$week, meanP_week$p_hat)
# p_hat_max <- max(meanP_week$p_hat, na.rm = TRUE)
# p_hat_min <- min(meanP_week$p_hat, na.rm = TRUE)
# cat("p_hat_max:", p_hat_max, "\n")
# 
# R0max_paras <- list(D = 5, L = 40*7, R0 = rep(2, 364*50)) # use 50 years as burn in
# R0max_pred <- as.data.frame(ode(c(S = census_pop-1, I = 1, R = 0),
#                                 seq(364*50),
#                                 SIRS_R0sig,
#                                 R0max_paras))
# 
# R0max_p <- R0max_pred$I/census_pop
# 
# plot_range = seq(45*364+1, 50*364) # take last 5 years
# #plot(1:length(plot_range), R0max_p[plot_range])
# p_max <- max(R0max_p[plot_range])
# cat("p_R0max:", p_max, "\n")
# 
# scaler = min(1, p_hat_max/p_max)
# cat("c:", scaler, "\n")

# Calculate q_cap
q_conf <- binom.confint(k, TT, conf.level = 1-1e-3, methods = c("exact"))
q_conf$date <- rel_date
ggplot() +
  geom_line(data = q_conf, aes(x = date, y=mean), color = "blue") +
  geom_ribbon(data = q_conf, aes(x=date, ymin=lower, ymax=upper), fill="blue", alpha=0.5) +
  geom_line(data = data.frame("X"= epi_data$rel_date, "pi"= epi_data$pi),
            aes(x = X, y=pi), color = "red") +
  geom_line(data = data.frame("X"= epi_data$rel_date, "p_hat"= k/TT*pi),
            aes(x = X, y=p_hat), color = "orange")

q_cap <- max(q_conf$upper)

#### Test that steady state reached regardless of when the 1st infection is seeded #####
## PASSED

state_p <- SIRS1var_pred(c(16.32242, 16.32242*0.5535876), sunob, census_pop)
state_p <- SIRS1var_pred(c(-272.7839, -272.7839*0.01), climob, census_pop)
state_p <- SIRS2var_pred(c(17.01981, -15.6381, 17.01981*0.5493769-15.6381*0.01599056), sunob, climob, census_pop)

plot_range = seq(0*364+1, 10*364) # take last 5 years
plot(1:length(plot_range), state_p[plot_range])


# Legacy
#low_q <- epi_df$k/epi_df$TT < 0.01
#mask <- is.na(epi_df$TT) | (epi_df$TT == 0)
#drop <- low_q | mask
# asc_corr <- epi_weekly$pi
# asc_corr[low_q] = 1

state_q <- p2q(state_p, epi_data, c_scaling = 0.6)
q_adj <- pmin(state_q, q_cap)

# q_adj <- mod_dat[[1]]
# epi_weekly <- mod_dat[[2]]
# 
# p1 <- ggplot() +
#   geom_point(data = data.frame("X"= epi_data$rel_date, "p_hat"= epi_data$k/epi_data$TT*epi_data$pi),
#              aes(x = X, y=p_hat), color = "red") +
#   geom_line(data = data.frame("X"= epi_data$rel_date, "p"= p_weekly),
#             aes(x = X, y=p), color = "orange") +
#   geom_line(data = data.frame("X"= epi_data$rel_date, "pi"= epi_weekly$pi),
#             aes(x = X, y=pi), color = "red")
# 
# p2 <- ggplot() +
#   geom_point(data = data.frame("X"= epi_weekly$rel_date, "q_hat"= epi_weekly$k/epi_weekly$TT),
#             aes(x = X, y=q_hat), color = "blue") +
#   geom_line(data = data.frame("X"= epi_weekly$rel_date, "pi"= epi_weekly$pi),
#             aes(x = X, y=pi), color = "red") +
#   geom_line(data = data.frame("X"= epi_weekly$rel_date, "p"= p_weekly),
#             aes(x = X, y=p), color = "orange")
# grid.arrange(p1, p2, ncol = 1)

p_low <- qbinom(0.025, size=epi_weekly$TT, prob=q_adj)/epi_weekly$TT
p_high <- qbinom(0.975, size=epi_weekly$TT, prob=q_adj)/epi_weekly$TT

ggplot() +
  geom_point(data = data.frame("X"= epi_data$rel_date, "q_hat"= epi_data$k/epi_data$TT), 
             aes(x = X, y=q_hat), color = "red") +
  geom_line(data = data.frame("X"= epi_data$rel_date, "q"= state_q),
          aes(x = X, y=q), color = "blue") +
  geom_line(data = data.frame("X"= epi_data$rel_date, "q_prime"= q_adj),
            aes(x = X, y=q_prime), color = "green")
  # geom_ribbon(data = data.frame("X"= epi_weekly$rel_date, "low"= p_low, "high"= p_high),
  #             aes(x=X, ymin=low, ymax=high), fill="blue", alpha=0.5)

ggplot() +
  geom_line(data = data.frame("X"= epi_data$rel_date, "pi"= epi_data$pi),
          aes(x = X, y=pi), color = "red")

# Calc. binomial probability
#k_T_q <- cbind(k_weekly, TT_weekly, q)
#neg_log_L <- -sum(apply(k_T_q, 1, function(xx) dbinom(xx[1], size=xx[2], prob=xx[3], log=TRUE)))
-sum(dbinom(epi_data$k, size=epi_data$TT, prob=q_adj, log=TRUE)) + 1e2*sum(state_q-q_adj)
