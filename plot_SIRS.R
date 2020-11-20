#!/usr/bin/env Rscript

library(deSolve)
library(DEoptim)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(binom)

source("SIRS_model.R")

######## Plot likelihood surface #########

LL_surf <- readRDS("fit_results/TX0sun111917GS.rds")
LL_mtx <- LL_surf[[3]]
colnames(LL_mtx) <- LL_surf[[2]]
rownames(LL_mtx) <- LL_surf[[1]]

long_LL <- melt(LL_mtx)
colnames(long_LL) <- c("k", "b", "R0")
long_LL[which.min(long_LL$R0), ]

ggplot(long_LL, aes(k, b)) +
  geom_raster(aes(fill = R0), interpolate = FALSE) + 
  scale_fill_distiller(palette = "Spectral", direction = -1, trans = "log")

##########################################

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

missing <- is.na(TT) | (TT == 0)
# Calculate q_cap
q_99 <- binom.confint(k[!missing], TT[!missing], conf.level = 0.99, methods = c("exact"))
q_99cap <- max(q_99$upper)

## Load fitting results

DE_fit <- readRDS("fit_results/NY0.6both111700_DE.rds")
cat(DE_fit$optim$bestmem, ";", DE_fit$optim$bestval)

SANN_fit <- readRDS(paste0("fit_results/TX1sun111615SANN.rds"))
cat(SANN_fit$par, ";", SANN_fit$value)

# choose 1
best <- DE_fit$optim$bestmem
best <- SANN_fit$par
# choose 1
state_p <- SIRS1var_pred(best, sunob, census_pop)           # sunrise model
state_p <- SIRS1var_pred(best, climob, census_pop)          # climate model
state_p <- SIRS2var_pred(best, sunob, climob, census_pop)   # combinaed model

mod_dat <- p2q(state_p, epi_data, q_99cap, 0.6) # <= remember to change scaling here
q_adj <- mod_dat[[1]]
epi_weekly <- mod_dat[[2]]

neg_log_L <- -sum(dbinom(epi_weekly$k, size=epi_weekly$TT, prob=q_adj, log=TRUE))

## Plot in p space
data_years <- ceiling(max(epi_weekly$rel_date)/364)
model_range <- seq((50-data_years)*364+1, 50*364)
p_adj <- state_p[model_range] # take only the number of years corresponding to data
p_adj <- p_adj*1

ggplot() + 
  geom_line(data = data.frame("X"= seq(length(p_adj)), "p"= p_adj), 
            aes(x = X, y=p), color = "blue") +
  geom_point(data = data.frame("X"= epi_weekly$rel_date, "p_hat"= epi_weekly$k/epi_weekly$TT*epi_weekly$pi), 
             aes(x = X, y=p_hat), color = "red", size=0.5)

## Plot in q space
q_low <- qbinom(0.025, size=epi_weekly$TT, prob=q_adj)/epi_weekly$TT
q_high <- qbinom(0.975, size=epi_weekly$TT, prob=q_adj)/epi_weekly$TT

ggplot() + 
  geom_line(data = data.frame("X"= epi_weekly$rel_date, "q"= q_adj), 
            aes(x = X, y=q), color = "blue") +
  geom_ribbon(data = data.frame("X"= epi_weekly$rel_date, "low"= q_low, "high"= q_high),
              aes(x=X, ymin=low, ymax=high), fill="blue", alpha=0.2) +
  geom_point(data = data.frame("X"= epi_weekly$rel_date, "q_hat"= epi_weekly$k/epi_weekly$TT), 
             aes(x = X, y=q_hat), color = "red", size=0.5)

