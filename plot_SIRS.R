#!/usr/bin/env Rscript

library(deSolve)
library(DEoptim)
library(ggplot2)
library(gridExtra)

######## Plot model #######

sun_range = seq(0, 1, by=0.01)
cli_range = seq(0, 0.03, by=0.0003)

p1 <- ggplot() + 
  geom_line(data = data.frame("hum"= cli_range, "R0"= R0_exp(cli_range, -10, 2, 1.2)), 
            aes(x = hum, y=R0), color = "blue") +
  geom_line(data = data.frame("hum"= cli_range, "R0"= R0_exp(cli_range, -100, 2, 1.2)), 
            aes(x = hum, y=R0), color = "red") +
  geom_line(data = data.frame("hum"= cli_range, "R0"= R0_exp(cli_range, -300, 2, 1.2)), 
            aes(x = hum, y=R0), color = "green")

p2 <- ggplot() + 
  geom_line(data = data.frame("sun"= sun_range, "R0"= R0_sig(sun_range, 10, 0.7, 2, 1.2)), 
            aes(x = sun, y=R0), color = "blue") +
  geom_line(data = data.frame("sun"= sun_range, "R0"= R0_sig(sun_range, 10, 0.5, 2, 1.2)), 
            aes(x = sun, y=R0), color = "red") +
  geom_line(data = data.frame("sun"= sun_range, "R0"= R0_sig(sun_range, 100, 0.5, 2, 1.2)), 
            aes(x = sun, y=R0), color = "green")
grid.arrange(p1, p2, ncol = 1)


###########################

state_code <- "TX"

## Load fitting results

DE_fit <- readRDS(paste0(state_code, "_DEoptim.rds"))
cat(DE_fit$optim$bestmem, ";", DE_fit$optim$bestval)

NM_fit <- readRDS(paste0(state_code, "_NMoptim.rds"))
cat(NM_fit$par, ";", NM_fit$value)

cli_negLL <- read.table(paste0(state_code, "_cli_negLL.tsv"), sep = "\t", col.names = c("alpha", "negLL"))
qplot(cli_negLL$alpha, cli_negLL$negLL)
cli_fit <-cli_negLL[which.min(cli_negLL$negLL), ]


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

################# Test binomial likelihood calc. (sunrise) #################
sig_params = DE_fit$optim$bestmem
var_observed = sunob
epi_df = epi_data
census_pop = census_pop
p_hat_range = c(p_hat_max, p_hat_min)

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

################# END OF SUNRISE #################

################# binomial likelihood (climate) #################
exp_para = cli_fit$alpha
var_observed = climob
epi_df = epi_data
census_pop = census_pop
p_hat_range = c(p_hat_max, p_hat_min)

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

################# END OF CLIMATE #################

state_p <- predictions$I/census_pop

data_years <- ceiling(max(epi_df$rel_date)/364)
model_range = seq((tot_years-data_years)*364-offset+1, tot_years*364-offset)
model_p <- state_p[model_range] # take only the number of years corresponding to data

# apply min-max scaling for model predicted probability
#scaler <- p_hat_range/(max(model_p)-min(model_p)) # if R0 is constant, min-max model_p are equal!!
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

q_low <- qbinom(0.025, size=epi_weekly$TT, prob=q)/epi_weekly$TT
q_high <- qbinom(0.975, size=epi_weekly$TT, prob=q)/epi_weekly$TT

ggplot() + 
  geom_line(data = data.frame("X"= epi_weekly$rel_date, "q"= q), 
            aes(x = X, y=q), color = "blue") +
  geom_ribbon(data = data.frame("X"= epi_weekly$rel_date, "low"= q_low, "high"= q_high),
              aes(x=X, ymin=low, ymax=high), fill="blue", alpha=0.2) +
  geom_point(data = data.frame("X"= epi_weekly$rel_date, "q_hat"= epi_weekly$k/epi_weekly$TT), 
             aes(x = X, y=q_hat), color = "red", size=0.5)
  # geom_line(data = data.frame("X"= epi_weekly$rel_date, "pi"= epi_weekly$pi), 
  #           aes(x = X, y=pi), color = "orange") +
  # geom_line(data = data.frame("X"= epi_weekly$rel_date, "p"= p_weekly), 
  #           aes(x = X, y=p), color = "green") +

