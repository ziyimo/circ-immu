#!/usr/bin/env Rscript

library(ggplot2)
library(gridExtra)

#library(deSolve)
library("binom")

source("R0_mods.R")
source("SIRS_model.R")

####### check missing data #######
state_pop <- read.delim("state_lv_data/state_pop.tsv")
state_pop$missing <- rep_len(NaN, nrow(state_pop))

for (state_idx in seq(nrow(state_pop))){
  ## Load flu data
  epiob <- read.csv(paste0("state_lv_data/Flu_data/flu_epi_", state_pop$code[state_idx], ".csv")) # blank field automatically NA
  epiob <- epiob[epiob$WEEK <= 52, ] # cap year to 52 weeks
  
  k <- epiob$TOTAL.A + epiob$TOTAL.B
  TT <- epiob$TOTAL.SPECIMENS
  pi <- epiob$X.UNWEIGHTED.ILI/100
  mask <- is.na(TT) | (TT == 0) | (pi == 0) # missing data: no. of test is 0 or NA, or no sympotomatic patient
  cat(state_pop$code[state_idx], sum(mask), "weeks masked", "\n")
  state_pop$missing[state_idx] <- sum(mask)
}
state_pop <- state_pop[order(state_pop$missing, decreasing = TRUE), ]
write.table(state_pop, "/Users/mo/Downloads/missing_data.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
hist(state_pop$missing, breaks = 20)
QC <- state_pop$missing < 200
View(state_pop[QC,])
write.table(state_pop$code[QC], "states_QC200.tsv", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
###################################

## Load all state data
state_pop <- read.delim("state_lv_data/state_pop.tsv")
all_state_sun <- read.csv("state_lv_data/state_daily_sunrise_2019.csv")
all_state_day <- read.csv("state_lv_data/state_daytime_2019.csv")
all_state_hum <- read.csv("state_lv_data/state_humidity_2014_2018.csv")
all_state_hum <- all_state_hum[all_state_hum$X != "2016-02-29", ] # get rid of leap year Feb 29

######## Humidity model sanity check #######

#shaman_states <- c("AZ", "IL", "NY", "WA")
shaman_states <- read.delim("states_QC200.tsv", header = FALSE)$V1

hum_fit <- readRDS("flu_fit/states_QC200.tsv_1e1_hum_1221.rds")
hum_alpha <- unname(hum_fit$optim$bestmem[2])

dat_I <- c()
mod_I <- c()

for (state_code in shaman_states){
  ## Load population
  census_pop <- state_pop$pop[state_pop$code==state_code]
  ## Load sunrise data
  #sunob <- all_state_sun[[state_code]]/720 # 365 days, scaled to maximum 720 = 12 hours
  
  ## Load humidity data
  climob <- all_state_hum[[state_code]] # multiple years
  climob <- matrix(climob, nrow = length(climob)/365, ncol = 365, byrow = TRUE)
  climob <- colMeans(climob) # 365 days
  
  epiob <- read.csv(paste0("state_lv_data/Flu_data/flu_epi_", state_code, ".csv")) # blank field automatically NA
  epiob <- epiob[epiob$WEEK <= 52, ] # cap year to 52 weeks
  mask <- is.na(epiob$TOTAL.SPECIMENS) | (epiob$TOTAL.SPECIMENS == 0) | (epiob$X.UNWEIGHTED.ILI == 0) # missing data: no. of test is 0 or NA, or no sympotomatic patient
  cat(">>", state_code, sum(mask), "weeks masked", "\n")
  epiob <- epiob[!mask, ]
  
  k <- epiob$TOTAL.A + epiob$TOTAL.B
  TT <- epiob$TOTAL.SPECIMENS
  pi <- epiob$X.UNWEIGHTED.ILI/100
  epiob$est_I <- k/TT*pi
  avg_inf <- aggregate(est_I ~ WEEK, data = epiob, mean)
  
  hum_pred <- SIRSvar_pred("R0_hum", hum_alpha, list(climob), census_pop)
  yearly_pred <- hum_pred[((tot_years-1)*364+1):(tot_years*364)][(avg_inf$WEEK*7-3)] # center on Thursday
  
  cat(">>", state_code, "r =", cor(avg_inf$est_I, yearly_pred), "\n")
  # p1 <- qplot(seq(52) , avg_inf$est_I)
  # p2 <- qplot(seq(52), yearly_pred)
  # grid.arrange(p1, p2, ncol = 1)
  
  dat_I <- c(dat_I, avg_inf$est_I)
  mod_I <- c(mod_I, yearly_pred)
}

cor(dat_I, mod_I, method = "spearman")

#############################################

##### calculate estimator of infection rate #####
infections_df <- data.frame(YEAR=rep(2010:2020, each=52), WEEK=rep(seq(52), times=11))[c(-seq(39), -(560:572)),]
pop_df <- data.frame(YEAR=rep(2010:2020, each=52), WEEK=rep(seq(52), times=11))[c(-seq(39), -(560:572)),]

states <- read.delim("states_QC200.tsv", header = FALSE)
state_pop <- read.delim("state_lv_data/state_pop.tsv")
#c <- 0.239310

for (state_code in states$V1){
  census_pop <- state_pop$pop[state_pop$code==state_code]
  epiob <- read.csv(paste0("state_lv_data/Flu_data/flu_epi_", state_code, ".csv")) # blank field automatically NA
  epiob <- epiob[epiob$WEEK <= 52, ] # cap year to 52 weeks
  mask <- is.na(epiob$TOTAL.SPECIMENS) | (epiob$TOTAL.SPECIMENS == 0) | (epiob$X.UNWEIGHTED.ILI == 0) # missing data: no. of test is 0 or NA, or no sympotomatic patient
  cat(sum(mask), "weeks masked", "\n")
  epiob <- epiob[!mask, ]
  # calculate relative date to the beginning of the 1st year in the data
  #rel_date <- (epiob$YEAR-y0)*364 + epiob$WEEK*7 - 3 # center on Thursday
  
  k <- epiob$TOTAL.A + epiob$TOTAL.B
  TT <- epiob$TOTAL.SPECIMENS
  pi <- epiob$X.UNWEIGHTED.ILI/100
  
  state_inf <- data.frame(YEAR=epiob$YEAR, WEEK=epiob$WEEK)
  state_inf[[state_code]] <- k/TT*pi*census_pop
  state_tot <- data.frame(YEAR=epiob$YEAR, WEEK=epiob$WEEK)
  state_tot[[state_code]] <- census_pop
  
  infections_df <- merge(infections_df, state_inf, all.x = TRUE)
  pop_df <- merge(pop_df, state_tot, all.x = TRUE)
}

write.table(infections_df, "flu_fit/p_estimator.tsv", quote = FALSE, sep="\t", row.names = FALSE)
write.table(pop_df, "flu_fit/denom_estimator.tsv", quote = FALSE, sep="\t", row.names = FALSE)
################################################

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
q_conf <- binom.confint(epi_data$k, epi_data$TT, conf.level = 1-1e-3, methods = c("exact"))
q_conf$date <- epi_data$rel_date
ggplot() +
  geom_line(data = q_conf, aes(x = date, y=mean), color = "blue") +
  geom_ribbon(data = q_conf, aes(x=date, ymin=lower, ymax=upper), fill="blue", alpha=0.5) #+
  # geom_line(data = data.frame("X"= epi_data$rel_date, "pi"= epi_data$pi),
  #           aes(x = X, y=pi), color = "red") +
  # geom_line(data = data.frame("X"= epi_data$rel_date, "p_hat"= k/TT*pi),
  #           aes(x = X, y=p_hat), color = "orange")

q_cap <- min(max(q_conf$upper), 1-1e-3)

### Troubleshoot const R0 ###
const_R0 <- 1.2

xstart <- c(S = census_pop-1, I = 1, R = 0) # use actual state population
# some hard-coded parameters
paras = list(D = 5, 
             L = 40*7, 
             R0 = rep(const_R0, times = 364*tot_years)[(offset+1):(364*tot_years)]) # pre-calculated R0 values

predictions <- as.data.frame(ode(xstart, times, SIRS_R0, paras)) # run SIRS model
raw_p <- predictions$I/census_pop
state_q <- p2q(raw_p, epi_data, 1)
q_adj <- pmin(state_q, q_cap) # cap the value of q, this is q'
neg_log_L <- -sum(dbinom(epi_data$k, size=epi_data$TT, prob=q_adj, log=TRUE)) + 10*sum(state_q-q_adj)


#############################

#### Test that steady state reached regardless of when the 1st infection is seeded #####

clean <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"),
               text = element_text(size=20))

# Global Var (for efficiency)
offset <- 180 # specify when to seed the single infection
tot_years <- 10 # no. of years to run SIRS model till steady state
times <- seq(1, 364 * tot_years, by = 1)[1:(364*tot_years-offset)]

# SIRS model prediction given variables
#SIRSvar_pred <- function(R0_func, R0_params, var_ls, census_pop){
# R0_func: string of function name
# R0_params: vector of parameters
# var_ls: list of vector(s) of variables (sunrise, humidity)

xstart <- c(S = census_pop-1, I = 1, R = 0) # use actual state population
annual_R0 <- do.call("R0_hum", c(list(climob), as.list(c(-60))))

# some hard-coded parameters
paras = list(D = 5,
             L = 40*7,
             R0 = rep(head(annual_R0, 364), times = tot_years)[(offset+1):(364*tot_years)]) # pre-calculated R0 values

predictions2 <- as.data.frame(ode(xstart, times, SIRS_R0, paras)) # run SIRS model
predictions2$time <- predictions2$time+offset

p1 <- ggplot() +
  geom_line(data = predictions1, aes(x = time, y=S), color = "black") +
  geom_line(data = predictions1, aes(x = time, y=I), color = "red") +
  geom_line(data = predictions1, aes(x = time, y=R), color = "blue") +
  xlim(0, tot_years*364)

#p2 <- 
ggplot() +
  geom_line(data = predictions2, aes(x = time, y=S), color = "black") +
  geom_line(data = predictions2, aes(x = time, y=I), color = "red") +
  geom_line(data = predictions2, aes(x = time, y=R), color = "blue") + 
  xlim((tot_years-10)*364, tot_years*364) + clean

grid.arrange(p1, p2, ncol = 1)


##### PASSED #####

state_p <- SIRSvar_pred("R0_exp", c(-100), list(climob), census_pop)

state_p <- SIRSvar_pred("R0_cdexp", c(-10, 0.6), list(sunob), census_pop)
state_p <- SIRSvar_pred("R0_bell", c(-50, 0.6), list(sunob), census_pop)

state_p <- SIRSvar_pred("R0_linEE" ,c(-10, 0.6, -50), list(sunob, climob), census_pop)
state_p <- SIRSvar_pred("R0_linGE" ,c(-50, 0.6, -50), list(sunob, climob), census_pop)
state_p <- SIRSvar_pred("R0_mixGE" ,c(-10, 0.6, -50, 0.75), list(sunob, climob), census_pop)

plot_range = seq(45*364+1, 50*364) # take last 5 years
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

### Experiment R0 models ###

fitted <- readRDS("fit_results/states_QC200.tsv_1e1_sun_0130.rds")
s_params <- unname(fitted$optim$bestmem[-1]) # do not include `c`
s_c <- fitted$optim$bestmem[1]

fitted <- readRDS("flu_fit/states_QC200.tsv_1e1_day_0111.rds")
d_params <- unname(fitted$optim$bestmem[-1]) # do not include `c`
d_c <- fitted$optim$bestmem[1]

fitted <- readRDS("flu_fit/states_QC200.tsv_1e1_sd_0111.rds")
sd_params <- unname(fitted$optim$bestmem[-1]) # do not include `c`
sd_c <- fitted$optim$bestmem[1]

all_state_sun <- read.csv("state_lv_data/state_daily_sunrise_2019.csv")
all_state_day <- read.csv("state_lv_data/state_daytime_2019.csv")

stateOI <- "AL"
census_pop <- state_pop$pop[state_pop$code==stateOI]
sunob <- all_state_sun[[stateOI]]/720 # 365 days, scaled to maximum 720 = 12 hours
dayob <- all_state_day[[stateOI]]/1440

## Plot R0 values
#source("R0_mods.R")
p1 <- ggplot() + 
  geom_line(data = data.frame("date"= seq(365), "R0"= do.call("R0_sun", c(list(sunob), as.list(s_params)))), 
            aes(x = date, y=R0), color = "#F0E442") + clean

p2 <- 
  ggplot() + 
    # geom_line(data = data.frame("date"= seq(365), "R0"= do.call("R0_sun", c(list(sunob), as.list(s_params)))), 
    #           aes(x = date, y=R0), color = "#F0E442") + 
  geom_line(data = data.frame("date"= seq(365), "R0"= do.call("R0_day", c(list(dayob), as.list(d_params)))), 
            aes(x = date, y=R0), color = "#E69F00") +
  geom_line(data = data.frame("date"= seq(365), "R0"= do.call("R0_sd", c(list(sunob, dayob), as.list(sd_params)))), 
            aes(x = date, y=R0), color = "#D55E00") + clean

grid.arrange(p1, p2, ncol = 2)

epiob <- read.csv(paste0("state_lv_data/Flu_data/flu_epi_", stateOI, ".csv")) # blank field automatically NA
epiob <- epiob[epiob$WEEK <= 52, ] # cap year to 52 weeks
mask <- is.na(epiob$TOTAL.SPECIMENS) | (epiob$TOTAL.SPECIMENS == 0) | (epiob$X.UNWEIGHTED.ILI == 0) # missing data: no. of test is 0 or NA, or no sympotomatic patient
cat(">>", stateOI, sum(mask), "weeks masked", "\n")
epiob <- epiob[!mask, ]

k <- epiob$TOTAL.A + epiob$TOTAL.B
TT <- epiob$TOTAL.SPECIMENS
pi <- epiob$X.UNWEIGHTED.ILI/100
epiob$est_I <- k/TT*pi
avg_inf <- aggregate(est_I ~ WEEK, data = epiob, median)

s_pred <- SIRSvar_pred("R0_sun", s_params, list(sunob), census_pop)
s_pred <- s_pred[((tot_years-1)*364+1):(tot_years*364)][(avg_inf$WEEK*7-3)] # center on Thursday

d_pred <- SIRSvar_pred("R0_day", d_params, list(dayob), census_pop)
d_pred <- d_pred[((tot_years-1)*364+1):(tot_years*364)][(avg_inf$WEEK*7-3)] # center on Thursday

sd_pred <- SIRSvar_pred("R0_sd", sd_params, list(sunob, dayob), census_pop)
sd_pred <- sd_pred[((tot_years-1)*364+1):(tot_years*364)][(avg_inf$WEEK*7-3)] # center on Thursday

ggplot(data = data.frame(week=seq(52), obs=avg_inf$est_I, s=s_pred*s_c, d=d_pred*d_c, sd=sd_pred*sd_c)) +
  geom_point(aes(x = week, y=obs), color = "#000000") +
  geom_line(aes(x = week, y=s), color = "#F0E442") +
  geom_line(aes(x = week, y=d), color = "#E69F00") +
  geom_line(aes(x = week, y=sd), color = "#D55E00") + clean
  


library(ggplot2)
library(gridExtra)

clean <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"),
               text = element_text(size=20))

sun_range = seq(0, 1, by=0.01)
day_range = seq(0, 1, by=0.01)
cli_range = seq(0, 0.03, by=0.0003)

R0_expd <- function(x_t, alpha, R0base = 2, R0min = 1.2){
  R0 <- exp(alpha*x_t + log(R0base - R0min)) + R0min
  return(R0)
}

R0_expg <- function(x_t, alpha, R0base = 2, R0min = 1.2){
  R0 <- exp(alpha*(1-x_t) + log(R0base - R0min)) + R0min
  return(R0)
}

R0_bellg <- function(x_t, alpha, R0base = 2, R0min = 1.2){
  R0 <- exp(alpha*(1-x_t)^2 + log(R0base - R0min)) + R0min
  return(R0)
}

ggplot() + 
  geom_line(data = data.frame("hum"= cli_range, "R0"= R0_expd(cli_range, -10)), 
            aes(x = hum, y=R0), color = "#E69F00") +
  geom_line(data = data.frame("hum"= cli_range, "R0"= R0_expd(cli_range, -100)), 
            aes(x = hum, y=R0), color = "#0072B2") +
  geom_line(data = data.frame("hum"= cli_range, "R0"= R0_expd(cli_range, -300)), 
            aes(x = hum, y=R0), color = "#009E73") + clean

ggplot() + 
  geom_line(data = data.frame("day"= day_range, "R0"= R0_expd(day_range, -1)), 
            aes(x = day, y=R0), color = "#E69F00") +
  geom_line(data = data.frame("day"= day_range, "R0"= R0_expd(day_range, -2.5)), 
            aes(x = day, y=R0), color = "#0072B2") +
  geom_line(data = data.frame("day"= day_range, "R0"= R0_expd(day_range, -5)), 
            aes(x = day, y=R0), color = "#009E73") + clean

p1 <- ggplot() + 
  geom_line(data = data.frame("day"= sun_range, "R0"= R0_expg(sun_range, -1)), 
            aes(x = day, y=R0), color = "#E69F00") +
  geom_line(data = data.frame("day"= sun_range, "R0"= R0_expg(sun_range, -3)), 
            aes(x = day, y=R0), color = "#0072B2") +
  geom_line(data = data.frame("day"= sun_range, "R0"= R0_expg(sun_range, -6)), 
            aes(x = day, y=R0), color = "#009E73") + clean

p2 <- ggplot() + 
  geom_line(data = data.frame("day"= sun_range, "R0"= R0_bellg(sun_range, -1)), 
            aes(x = day, y=R0), color = "#E69F00") +
  geom_line(data = data.frame("day"= sun_range, "R0"= R0_bellg(sun_range, -3)), 
            aes(x = day, y=R0), color = "#0072B2") +
  geom_line(data = data.frame("day"= sun_range, "R0"= R0_bellg(sun_range, -8)), 
            aes(x = day, y=R0), color = "#009E73") + clean

grid.arrange(p1, p2, ncol = 1)


all_state_sun <- read.csv("state_lv_data/state_daily_sunrise_2019.csv")
all_state_day <- read.csv("state_lv_data/state_daytime_2019.csv")

hist(as.vector(t(all_state_day[c(-1, -9)]/1440)))
hist(as.vector(t(all_state_sun[c(-1, -9)]/720)))
