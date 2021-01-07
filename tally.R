#!/usr/bin/env Rscript

library("binom")
source("SIRS_model.R")

args <- commandArgs(trailingOnly=TRUE)
state_ls <- args[1]      # text file with a list of states to fit to
R0_mod <- args[2]        # functional form of R0
prms <- as.numeric(strsplit(args[3], ":", fixed=TRUE)[[1]])          # fitted params, separated by ":" 
l_pnl <- args[4]         # lambda for penalized likelihood

handle <- paste0(state_ls, l_pnl, R0_mod, "_stateFit")
l_pnl <- as.numeric(l_pnl)

sink(paste0("fit_results/", handle, ".log"))

states <- read.delim(state_ls, header = FALSE, stringsAsFactors = FALSE)

## Load all state data
state_pop <- read.delim("state_lv_data/state_pop.tsv")
all_state_sun <- read.csv("state_lv_data/state_daily_sunrise_2019.csv")
all_state_day <- read.csv("state_lv_data/state_daytime_2019.csv")
all_state_hum <- read.csv("state_lv_data/state_humidity_2014_2018.csv")
all_state_hum <- all_state_hum[all_state_hum$X != "2016-02-29", ] # get rid of leap year Feb 29

state_data <- list()

for (state_i in seq(nrow(states))){
  state_code <- states$V1[state_i]
  cat(">>> Loading state data:", state_code, "<<<\n")
  census_pop <- state_pop$pop[state_pop$code==state_code]
  
  if (R0_mod %in% c("s2d2", "hs2d2")){
    sunob <- all_state_sun[[state_code]]/720 # 365 days, scaled to maximum 720 = 12 hours
  }
  if (R0_mod %in% c("day2", "hd2", "s2d2", "hs2d2")){
    dayob <- all_state_day[[state_code]]/1440
  }
  if (R0_mod %in% c("hum", "hd2", "hs2d2")){
    climob <- all_state_hum[[state_code]] # multiple years
    climob <- matrix(climob, nrow = length(climob)/365, ncol = 365, byrow = TRUE)
    climob <- colMeans(climob) # 365 days
  }
  
  if (R0_mod == "cos"){
    varob <- list(seq(365))
  } else if (R0_mod == "hum"){
    varob <- list(climob)
  } else if (R0_mod == "day2"){
    varob <- list(dayob)
  } else if (R0_mod == "hd2"){
    varob <- list(climob, dayob)
  } else if (R0_mod == "s2d2"){
    varob <- list(sunob, dayob)
  } else if (R0_mod == "hs2d2"){
    varob <- list(climob, sunob, dayob)
  }
  
  epi_data <- load_state_epi(paste0("state_lv_data/Flu_data/flu_epi_", state_code, ".csv"))
  
  q_conf <- binom.confint(epi_data$k, epi_data$TT, conf.level = 1-1e-3, methods = c("exact"))
  q_99cap <- min(max(q_conf$upper), 1-1e-3)
  cat(">> q capped at", q_99cap, "\n")
  
  state_data[[state_code]] <- list(idx=state_i, pop=census_pop, var=varob, epi=epi_data, cap=q_99cap)
}

get_negLL <- function(state_elem){
  neg_LL <- binom_Lp(prms,
                     paste0("R0_", R0_mod),
                     state_elem$var,
                     state_elem$epi,
                     state_elem$pop,
                     state_elem$cap,
                     l_pnl)
  return(neg_LL)
}

states$negLL <- sapply(X=state_data, FUN=get_negLL)

write.table(states, file = paste0("fit_results/", handle, ".tsv"), quote=FALSE, sep="\t", row.names=FALSE)
cat("Sum up to:", sum(states$negLL))
sink()