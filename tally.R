#!/usr/bin/env Rscript

library(DEoptim)

states <- c("AZ", "CA", "IL", "NY", "TX", "WA")

fit_df <- data.frame(sun_nLL=numeric(length(states)), sun_k=numeric(length(states)), sun_x0=numeric(length(states)),
                     cli_nLL=numeric(length(states)), cli_alpha=numeric(length(states)), row.names = states)

for (state_idx in seq(length(states))) {
  
  state_code <- states[state_idx]
  DE_fit <- readRDS(paste0("fit_results/", state_code, "_DEoptim.rds"))
  cat(state_code, "/sun, DE:", DE_fit$optim$bestmem, ";", DE_fit$optim$bestval, "\n")
  fit_df$sun_nLL[state_idx] <- DE_fit$optim$bestval
  fit_df$sun_k[state_idx] <- DE_fit$optim$bestmem[1]
  fit_df$sun_x0[state_idx] <- DE_fit$optim$bestmem[2]
  
  NM_fit <- readRDS(paste0("fit_results/", state_code, "_NMoptim.rds"))
  cat(state_code, "/sun, NM:", NM_fit$par, ";", NM_fit$value, "\n")
  
  cli_negLL <- read.table(paste0("fit_results/", state_code, "_cli_negLL.tsv"), sep = "\t", col.names = c("alpha", "negLL"))
  #qplot(cli_negLL$alpha, cli_negLL$negLL)
  cli_fit <-cli_negLL[which.min(cli_negLL$negLL), ]
  cat(state_code, "/cli:", cli_fit$alpha, ";", cli_fit$negLL, "\n")
  fit_df$cli_nLL[state_idx] <- cli_fit$negLL
  fit_df$cli_alpha[state_idx] <- cli_fit$alpha
    
}

write.csv(fit_df, file="fit_results/prelim.csv")
