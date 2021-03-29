#!/usr/bin/env Rscript

suppressWarnings(suppressMessages(library("pso", character.only = TRUE)))
suppressWarnings(suppressMessages(library("parallel")))
suppressWarnings(suppressMessages(library("binom")))
source("Covid_state_joint.R")

### Parse user args
args <- commandArgs(trailingOnly=TRUE)
state_ls <- args[1]  # 2-letter state code
R0_mod <- args[2] # functional form of R0, options: day,hum,hd,sd,hsd
seed_sh <- as.numeric(strsplit(args[3], ":", fixed=TRUE)[[1]])  # h_rate, R0min, R0range, [R0_params]
nthr <- args[4]
# text file with a list of states to fit to

### Logging
time_stamp <- format(Sys.time(), "%m%d%H%m")

handle <- paste0(state_ls, "_", R0_mod,"_", time_stamp)
handle <- sub(" ", "", handle)
cat(paste0("fit_results/", handle, ".log"))

# Optimizer
optimizer <- "pso"

### Load population size data
states <- read.delim(state_ls, header = FALSE, stringsAsFactors=FALSE)
no_cities <- nrow(states)

## Load all state data
state_pop <- read.delim("state_lv_data/state_pop.tsv")
all_state_sun <- read.csv("state_lv_data/state_daily_sunrise_2019.csv")
all_state_day <- read.csv("state_lv_data/state_daytime_2019.csv")
all_state_hum <- read.csv("state_lv_data/state_humidity_2014_2018.csv")
all_state_hum <- all_state_hum[all_state_hum$X != "2016-02-29", ] # get rid of leap year Feb 29
## Load COVID hospitalization data
#template <- readRDS("data/state_hospitalization.rds")

covid_df <- read.csv("data/states_death.csv", header = TRUE)
covid_df$date <- as.numeric(as.Date(covid_df$date, format="%Y-%m-%d") - as.Date("2019-12-31", format="%Y-%m-%d"))
state2reg <- read.csv("state_lv_data/state2censusReg.csv")
cdv <- unique(state2reg$Division)
no_cdv <- length(cdv)
state2reg$CDV.Code <- match(state2reg$Division, cdv)

city_data <- list()
seed_ss <- list()

for (state_i in seq(no_cities)){
  state_code <- states$V1[state_i]
  cat(">>> Loading state data:", state_code, "<<<\n")
  census_pop <- state_pop$pop[state_pop$code==state_code]
  
  if (R0_mod %in% c("sun", "hs", "sd", "hsd")){
    sunob <- all_state_sun[[state_code]]/720 # 365 days, scaled to maximum 720 = 12 hours
  }
  if (R0_mod %in% c("day", "hd", "sd", "hsd")){
    dayob <- all_state_day[[state_code]]/1440
  }
  if (R0_mod %in% c("hum", "hs", "hd", "hsd")){
    climob <- all_state_hum[[state_code]] # multiple years
    climob <- matrix(climob, nrow = length(climob)/365, ncol = 365, byrow = TRUE)
    climob <- colMeans(climob) # 365 days
  }
  
  if (R0_mod == "cos"){
    varob <- list(seq(365))
  } else if (R0_mod == "hum"){
    varob <- list(climob)
  } else if (R0_mod == "day"){
    varob <- list(dayob)
  } else if (R0_mod == "hd"){
    varob <- list(climob, dayob)
  } else if (R0_mod == "sd"){
    varob <- list(sunob, dayob)
  } else if (R0_mod == "hsd"){
    varob <- list(climob, sunob, dayob)
  } else if (R0_mod == "sun"){
    varob <- list(sunob)
  } else if (R0_mod == "hs"){
    varob <- list(climob, sunob)
  }
  
  state_df <- subset(covid_df, state == state_code)
  state_df <- state_df[order(state_df$date),]
  state_df <- state_df[!is.na(state_df$deathIncrease_cor), ]
  state_df <- state_df[state_df$date <= 396, ]
  state_deaths <- subset(state_df, select = c("date", "deathIncrease_cor"))
  
  city_data[[state_code]] <- list(idx=state_i, pop=census_pop, var=varob, epi=state_deaths, cdv=state2reg$CDV.Code[state2reg$State.Code == state_code])
}


# State specific params - I_init
ss_low <- c(1e-5)
ss_high <- c(0.1)

if (R0_mod %in% c("sun", "sd", "hsd")){
  s0_idx <- ifelse(R0_mod %in% c("sun", "sd"), 2, 3)
  # Shared params: [CD s_0], hosp_rate, R0min, R0range, [other R0 params]
  sh_low <- c(rep(0.4, no_cdv), 0.001, 0.8, 0.6, param_bounds[[R0_mod]]$low[-s0_idx])
  sh_high <- c(rep(0.7, no_cdv), 0.001, 1.2, 1.3, param_bounds[[R0_mod]]$high[-s0_idx])
  
  seed_sh <- c(rep(seed_sh[3+s0_idx], no_cdv), seed_sh[-(3+s0_idx)])
} else {
  # Shared params: h_rate, R0min, R0range, [R0 params]
  sh_low <- c(0.001, 0.8, 1.2, param_bounds[[R0_mod]]$low)
  sh_high <- c(0.001, 1.2, 1.6, param_bounds[[R0_mod]]$high)
}

state_wrapper <- function(state_elem, ss_prms, sh_prms){
  if (R0_mod %in% c("sun", "sd", "hsd")){
    # Init, alpha,R0min,R0range,s_0,other
    state_prms <- c(ss_prms, sh_prms[(no_cdv+1):(no_cdv+3+s0_idx-1)], sh_prms[state_elem$cdv], sh_prms[-1:-(no_cdv+3+s0_idx-1)])
    cat(state_prms,"\n")
  } else {
    state_prms <- c(ss_prms, sh_prms)
  }
  # inf_init,alpha,R0min,R0range,s,s0,d,d0
  neg_LL <- binom_seird(state_prms, 
                   paste0("R0_", R0_mod),
                   state_elem$var,
                   state_elem$epi,
                   state_elem$pop)
  
  return(neg_LL)  
}

fit_state <- function(state_dat, sh_vec){ # no seeding for 1D optimization
  oo <- optimize(ss_wrapper, c(0, 1), st_e=state_dat, sh_p=sh_vec)
  
  return(ss_low + oo$minimum*(ss_high-ss_low))
}


ss_wrapper <- function(ss_norm, st_e, sh_p){
  ss_p <- ss_low + ss_norm*(ss_high-ss_low)
  return(state_wrapper(st_e, ss_p, sh_p))
}

all_state_negLL_test <- function(norm_sh_prms, ss_prms){
  orig_prms <- sh_low + norm_sh_prms*(sh_high-sh_low)
  states_negLL <- sh_wrapper_test(st_e=city_data$CA, sh_orig=orig_prms, ss_ls=ss_prms)
  return(sum(states_negLL))
}

## Function for fitting shared parameters
all_state_negLL <- function(norm_sh_prms, ss_prms){
  orig_prms <- sh_low + norm_sh_prms*(sh_high-sh_low)
  states_negLL <- parSapply(cl = clstr, X=city_data, FUN=sh_wrapper, sh_orig=orig_prms, ss_ls=ss_prms)
  return(sum(states_negLL))
}

sh_wrapper_test <- function(st_e, sh_orig, ss_ls){
  ss_p <- ss_ls
  return(state_wrapper(st_e, ss_p, sh_orig))
}

sh_wrapper <- function(st_e, sh_orig, ss_ls){
  ss_p <- ss_ls[[st_e$idx]]
  return(state_wrapper(st_e, ss_p, sh_orig))
}


cat(detectCores(), "cores seen; limit to", nthr, "cores\n")
clstr <- makeCluster(nthr, type = "FORK") # limits the number of cores to use
for (iter in seq(20)){
  cat(">> Iteration:", iter, "\n")
  ss_fit <- parLapply(cl = clstr, X = city_data, fun = fit_state, sh_vec = seed_sh)
  print(ss_fit)
  cat(seed_sh,"\n")
  sh_fit <- psoptim((seed_sh-sh_low)/(sh_high-sh_low), all_state_negLL, ss_prms=ss_fit,
                    lower=rep(0,length(sh_low)), upper=rep(1,length(sh_low)),
                    control=list(trace=1, REPORT=5, maxit=5000, trace.stats=FALSE, maxit.stagnate=200))
  seed_sh <- sh_low + sh_fit$par*(sh_high-sh_low)

  cat(">>", R0_mod, "Iter", iter, "negLL:", sh_fit$value, "\n")
  cat(">> ", unlist(ss_fit), "\n>> ", seed_sh, "\n")
  save(ss_fit, seed_sh, file = paste0("fit_results/", handle,"_iter",iter, "_iterative.fit3.rds"))
}

save(ss_fit, seed_sh, file = paste0("fit_results/", handle, "_iterative.fit3.rds"))
stopCluster(clstr)
