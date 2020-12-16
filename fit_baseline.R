#!/usr/bin/env Rscript

#.libPaths() # sanity check
library("binom")
source("SIRS_model.R")

time_stamp <- format(Sys.time(), "%m%d%H")

args <- commandArgs(trailingOnly=TRUE)
state_ls <- args[1]      # text file with a list of states to fit to
R0_lower <- args[2]
R0_upper <- 2
l_pnl <- args[3]         # lambda for penalized likelihood

handle <- paste0(state_ls, l_pnl, "const", R0_lower, "_", time_stamp)
R0_lower <- as.numeric(R0_lower)
l_pnl <- as.numeric(l_pnl)

sink(paste0("fit_results/", handle, ".log"))

states <- read.delim(state_ls, header = FALSE, stringsAsFactors = FALSE)$V1

## Load all state data
state_pop <- read.delim("state_lv_data/state_pop.tsv")

state_data <- list()

for (state_i in seq(length(states))){
  state_code <- states[state_i]
  cat(">>> Loading state data:", state_code, "<<<\n")
  census_pop <- state_pop$pop[state_pop$code==state_code]
  
  epi_data <- load_state_epi(paste0("state_lv_data/Flu_data/flu_epi_", state_code, ".csv"))
  
  q_conf <- binom.confint(epi_data$k, epi_data$TT, conf.level = 1-1e-3, methods = c("exact"))
  q_99cap <- min(max(q_conf$upper), 1-1e-3)
  cat(">> q capped at", q_99cap, "\n")
  
  state_data[[state_code]] <- list(idx=state_i, pop=census_pop, epi=epi_data, cap=q_99cap)
}

Lp_const <- function(const_R0, state_elem){
  
  xstart <- c(S = state_elem$pop-1, I = 1, R = 0) # use actual state population
  # some hard-coded parameters
  paras = list(D = 5, 
               L = 40*7, 
               R0 = rep(const_R0, times = 364*tot_years)[(offset+1):(364*tot_years)]) # pre-calculated R0 values
  
  predictions <- as.data.frame(ode(xstart, times, SIRS_R0, paras)) # run SIRS model
  raw_p <- predictions$I/state_elem$pop
  state_q <- p2q(raw_p, state_elem$epi, 1)
  q_adj <- pmin(state_q, state_elem$cap) # cap the value of q, this is q'
  neg_log_L <- -sum(dbinom(state_elem$epi$k, size=state_elem$epi$TT, prob=q_adj, log=TRUE)) + l_pnl*sum(state_q-q_adj)
  return(neg_log_L)
}

findR0 <- function(XX) optimize(Lp_const, c(R0_lower, R0_upper), state_elem = XX)

cat(">> Fitting the baseline model; lambda =", l_pnl, "\n")

R0_fit <- lapply(X=state_data, FUN=findR0)

baseline_df <- data.frame(state=character(length(states)),
                          R0_const=numeric(length(states)),
                          negLL=numeric(length(states)))

for (state_i in seq(length(states))){
  state_code <- states[state_i]
  baseline_df$state[state_i] <- state_code
  baseline_df$R0_const[state_i] <- R0_fit[[state_code]]$minimum
  baseline_df$negLL[state_i] <- R0_fit[[state_code]]$objective
}

write.table(baseline_df, file = paste0("fit_results/", handle, ".tsv"), quote=FALSE, sep="\t", row.names=FALSE)
cat("Baseline likelihood:", sum(baseline_df$negLL))
sink()
