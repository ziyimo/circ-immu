#!/usr/bin/env Rscript

library(deSolve)
source("R0_mods.R")
source("SIRS_model.R")

state_ls <- "US_states.tsv"      # text file with a list of states to fit to

fitted <- readRDS("fit_results/US_1e1_sd_0203.rds")
sd_optim <- unname(fitted$optim$bestmem)

fitted <- readRDS("fit_results/US_1e1_day_0202.rds")
d_optim <- unname(fitted$optim$bestmem)

fitted <- readRDS("fit_results/US_1e1_tmp_0202.rds")
t_optim <- unname(fitted$optim$bestmem)

states <- read.delim(state_ls, header = FALSE)$V1

## Load all state data
state_pop <- read.delim("env_covar/US_pop.tsv")
all_state_sun <- read.csv("env_covar/US_sunrise_2019.csv")
all_state_day <- read.csv("env_covar/US_daytime_2019.csv")
all_state_tmp <- read.csv(paste0("env_covar/US_temperature_2010_18.csv"))



sim_SIRS <- function(N_pop, annual_R0){
  
  xstart <- c(S = N_pop-1, I = 1, R = 0) # use actual state population
  # some hard-coded parameters
  paras = list(D = 5, 
               L = 40*7, 
               R0 = rep(head(annual_R0, 364), times = tot_years)[(offset+1):(364*tot_years)]) # pre-calculated R0 values
  predictions <- as.data.frame(ode(xstart, times, SIRS_R0, paras)) # run SIRS model
  
  ys <- (tot_years-1)*364-offset+1 # year start (last year)
  ye <- tot_years*364-offset       # year end
  
  annual_inf1 <- sum(predictions$I[ys:ye] - 0.8*predictions$I[(ys-1):(ye-1)]) # for D=5
  annual_inf2 <- sum(head(annual_R0, 364)/5*predictions$S[ys:ye]*predictions$I[ys:ye]/census_pop)
  cat(annual_inf1, annual_inf2, "\n") # sanity check
  
  return(list(inf_traj=predictions$I[ys:ye], tot_inf=annual_inf1))
}


infection_sd <- data.frame(matrix(ncol = length(states), nrow = 364))
infection_d <- data.frame(matrix(ncol = length(states), nrow = 364))
infection_t <- data.frame(matrix(ncol = length(states), nrow = 364))
colnames(infection_sd) <- states
colnames(infection_d) <- states
colnames(infection_t) <- states

for (state_i in seq(length(states))){
  state_code <- states[state_i]
  cat(">>> Simulating:", state_code, "<<<\n")
  census_pop <- state_pop$pop[state_pop$code==state_code]
  
  sunob <- all_state_sun[[state_code]]/720 # 365 days, scaled to maximum 720 = 12 hours
  dayob <- all_state_day[[state_code]]/1440
  
  tmpob <- all_state_tmp[[state_code]]
  tmpob <- matrix(tmpob, nrow = length(tmpob)/365, ncol = 365, byrow = TRUE)
  tmpob <- colMeans(tmpob) # 365 days
  
  sd_R0 <- do.call("R0_sd", c(list(sunob, dayob), as.list(sd_optim[-1])))
  sd_sim <- sim_SIRS(census_pop, sd_R0)
  infection_sd[[state_code]] <- sd_sim$inf_traj
  
  d_R0 <- do.call("R0_day", c(list(dayob), as.list(d_optim[-1])))
  d_sim <- sim_SIRS(census_pop, d_R0)
  infection_d[[state_code]] <- d_sim$inf_traj
  
  t_R0 <- do.call("R0_tmp", c(list(tmpob), as.list(t_optim[-1])))
  t_sim <- sim_SIRS(census_pop, t_R0)
  infection_t[[state_code]] <- t_sim$inf_traj

}

write.table(infection_sd, file = "flu_fit/annual_inf_sd.tsv", quote=FALSE, sep="\t", row.names=FALSE)
write.table(infection_d, file = "flu_fit/annual_inf_d.tsv", quote=FALSE, sep="\t", row.names=FALSE)
write.table(infection_t, file = "flu_fit/annual_inf_t.tsv", quote=FALSE, sep="\t", row.names=FALSE)
