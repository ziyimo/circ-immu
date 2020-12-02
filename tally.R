#!/usr/bin/env Rscript

library(stringr)
library(DEoptim)

fit_files <- list.files("fit_results", pattern = "*\\.rds")

fit_df <- NULL

for (file_name in fit_files) {
  
  code_method <- str_extract_all(file_name, "[A-Z]+")[[1]]
  lambda_val <- as.numeric(str_extract(file_name, "[0-9]+"))
  variable <- str_extract(file_name, "[a-z]+")
  
  if (code_method[2] == "DE"){
    DE_fit <- readRDS(paste0("fit_results/", file_name))
    negLL <- DE_fit$optim$bestval
    bestParam <- DE_fit$optim$bestmem
  } else if (code_method[2] != "GS"){
    optim_fit <- readRDS(paste0("fit_results/", file_name))
    bestParam <- optim_fit$par
    negLL <- optim_fit$value
  } else{
    next
  }
  
  if (variable == "both"){
    fit_df <- rbind(fit_df, c(code_method, lambda_val, variable, negLL, bestParam))
  } else if (variable == "cli"){
    fit_df <- rbind(fit_df, c(code_method, lambda_val, variable, negLL, "NA", bestParam))
  } else if (variable == "sun"){
    fit_df <- rbind(fit_df, c(code_method, lambda_val, variable, negLL, bestParam[1], "NA", bestParam[2:3]))
  }
}

fit_df <- as.data.frame(fit_df)
colnames(fit_df) <- c("state", "method", "lambda", "var", "negLL", "k_sun", "k_cli", "b", "c")

write.csv(fit_df, file="fit_results/prelim.csv")
