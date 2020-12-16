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

### Load data
state_pop <- read.delim("state_lv_data/state_pop.tsv")
all_state_sun <- read.csv("state_lv_data/state_daily_sunrise_2019.csv")
all_state_hum <- read.csv("state_lv_data/state_humidity_2014_2018.csv")
all_state_hum <- all_state_hum[all_state_hum$X != "2016-02-29", ] # get rid of leap year Feb 29

state_code <- "AZ"

## Load population
census_pop <- state_pop$pop[state_pop$code==state_code]
## Load sunrise data
sunob <- all_state_sun[[state_code]]/720 # 365 days, scaled to maximum 720 = 12 hours
## Load humidity data
climob <- all_state_hum[[state_code]] # multiple years
climob <- matrix(climob, nrow = length(climob)/365, ncol = 365, byrow = TRUE)
climob <- colMeans(climob) # 365 days

## Load flu data
epi_data <- load_state_epi(paste0("state_lv_data/Flu_data/flu_epi_", state_code, ".csv"))

## Calculate q_cap
q_conf <- binom.confint(epi_data$k, epi_data$TT, conf.level = 1-1e-3, methods = c("exact"))
q_99cap <- max(q_conf$upper)
cat(">> q capped at", q_99cap, "\n")

## Run model

R0_mod <- "linGE"
lambda <- 1e3

if (R0_mod == "exp"){
  varob <- list(climob)
} else if (R0_mod %in% c("cdexp", "bell")){
  varob <- list(sunob)
} else if (R0_mod %in% c("linEE", "linGE", "mixGE")){
  varob <- list(sunob, climob)
}

# DE_fit <- readRDS("fit_results/CA1both112008_DE.rds")
# cat(DE_fit$optim$bestmem, ";", DE_fit$optim$bestval)
# 
# SANN_fit <- readRDS(paste0("fit_results/TX1sun111615SANN.rds"))
# cat(SANN_fit$par, ";", SANN_fit$value)

best <- c(0.223, -7.251, 0.913, -19.430) #DE_fit$optim$bestmem, SANN_fit$par

state_p <- SIRSvar_pred(paste0("R0_", R0_mod), best[-1], varob, census_pop)
state_q <- p2q(state_p, epi_data, best[1])
q_adj <- pmin(state_q, q_99cap)

neg_log_L <- -sum(dbinom(epi_data$k, size=epi_data$TT, prob=q_adj, log=TRUE)) + lambda*sum(state_q-q_adj)
cat(neg_log_L)

## Plot in p space
data_years <- ceiling(max(epi_data$rel_date)/364)
model_range <- seq((50-data_years)*364+1, 50*364)
p_adj <- state_p[model_range] # take only the number of years corresponding to data
p_adj <- p_adj*best[1]

ggplot() + 
  geom_line(data = data.frame("X"= seq(length(p_adj)), "p"= p_adj), 
            aes(x = X, y=p), color = "blue") +
  geom_point(data = data.frame("X"= epi_data$rel_date, "p_hat"= epi_data$k/epi_data$TT*epi_data$pi), 
             aes(x = X, y=p_hat), color = "red", size=0.5)

## Plot in q space
q_low <- qbinom(0.025, size=epi_data$TT, prob=q_adj)/epi_data$TT
q_high <- qbinom(0.975, size=epi_data$TT, prob=q_adj)/epi_data$TT

ggplot() + 
  geom_line(data = data.frame("X"= epi_data$rel_date, "q"= q_adj), 
            aes(x = X, y=q), color = "blue") +
  geom_ribbon(data = data.frame("X"= epi_data$rel_date, "low"= q_low, "high"= q_high),
              aes(x=X, ymin=low, ymax=high), fill="blue", alpha=0.2) +
  geom_point(data = data.frame("X"= epi_data$rel_date, "q_hat"= epi_data$k/epi_data$TT), 
             aes(x = X, y=q_hat), color = "red", size=0.5)

