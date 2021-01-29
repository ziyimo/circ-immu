#!/usr/bin/env Rscript

library("binom")

source("R0_mods.R")
source("SIRS_model.R")

args <- commandArgs(trailingOnly=TRUE)
state_ls <- "states_QC200.tsv"      # text file with a list of states to fit to
#prms <- as.numeric(strsplit(args[3], ":", fixed=TRUE)[[1]])          # fitted params, separated by ":" 
l_pnl <- args[1]         # lambda for penalized likelihood

fit_prms <- list()
R0_mods <- c("cos", "hum", "day", "sd", "hd", "hsd")        # functional form of R0
no_prms <- c(2, 2, 3, 5, 4, 6)

calc_BIC <- function(negLL, k, n=523*40){
  return(k*log(n)+2*negLL)
}

for (idx in seq(length(R0_mods))){
  R0m <- R0_mods[idx]
  fitted <- readRDS(Sys.glob(paste0("flu_fit/", paste(state_ls, l_pnl, R0m, "*.rds" ,sep="_"))))
  neg_LL <- fitted$optim$bestval
  best_prms <- fitted$optim$bestmem
  
  fit_prms[[R0m]] <- unname(best_prms)
  cat(R0m, neg_LL, calc_BIC(neg_LL, no_prms[idx]), best_prms, "\n", sep = "\t")
}

#l_pnl <- as.numeric(l_pnl)

states <- read.delim(state_ls, header = FALSE, stringsAsFactors = FALSE)$V1
neg_LL_df <- data.frame(matrix(ncol = length(states), nrow = length(R0_mods)))
colnames(neg_LL_df) <- states

## Load all state data
state_pop <- read.delim("state_lv_data/state_pop.tsv")
all_state_sun <- read.csv("state_lv_data/state_daily_sunrise_2019.csv")
all_state_day <- read.csv("state_lv_data/state_daytime_2019.csv")
all_state_hum <- read.csv("state_lv_data/state_humidity_2014_2018.csv")
all_state_hum <- all_state_hum[all_state_hum$X != "2016-02-29", ] # get rid of leap year Feb 29
#state_data <- list()

for (state_i in seq(length(states))){
  state_code <- states[state_i]
  #cat(">>> Loading state data:", state_code, "<<<\n")
  census_pop <- state_pop$pop[state_pop$code==state_code]
  
  sunob <- all_state_sun[[state_code]]/720 # 365 days, scaled to maximum 720 = 12 hours
  dayob <- all_state_day[[state_code]]/1440
  climob <- all_state_hum[[state_code]] # multiple years
  climob <- matrix(climob, nrow = length(climob)/365, ncol = 365, byrow = TRUE)
  climob <- colMeans(climob) # 365 days
  
  epi_data <- load_state_epi(paste0("state_lv_data/Flu_data/flu_epi_", state_code, ".csv"))
  
  q_conf <- binom.confint(epi_data$k, epi_data$TT, conf.level = 1-1e-3, methods = c("exact"))
  q_99cap <- min(max(q_conf$upper), 1-1e-3)
  #cat(">> q capped at", q_99cap, "\n")
  
  #state_data[[state_code]] <- list(idx=state_i, pop=census_pop, var=varob, epi=epi_data, cap=q_99cap)
  neg_LL_df[[state_code]][1] <- binom_Lp(fit_prms$cos, "R0_cos", list(seq(365)), epi_data, census_pop, q_99cap, as.numeric(l_pnl))
  neg_LL_df[[state_code]][2] <- binom_Lp(fit_prms$hum, "R0_hum", list(climob), epi_data, census_pop, q_99cap, as.numeric(l_pnl))
  neg_LL_df[[state_code]][3] <- binom_Lp(fit_prms$day, "R0_day", list(dayob), epi_data, census_pop, q_99cap, as.numeric(l_pnl))
  neg_LL_df[[state_code]][4] <- binom_Lp(fit_prms$sd, "R0_sd", list(sunob, dayob), epi_data, census_pop, q_99cap, as.numeric(l_pnl))
  neg_LL_df[[state_code]][5] <- binom_Lp(fit_prms$hd, "R0_hd", list(climob, dayob), epi_data, census_pop, q_99cap, as.numeric(l_pnl))
  neg_LL_df[[state_code]][6] <- binom_Lp(fit_prms$hsd, "R0_hsd", list(climob, sunob, dayob), epi_data, census_pop, q_99cap, as.numeric(l_pnl))
}


write.table(neg_LL_df, file = paste0("flu_fit/", "state_fits", "_", l_pnl, ".tsv"), quote=FALSE, sep="\t", row.names=R0_mods)
cat("Sum up to:", rowSums(neg_LL_df))
