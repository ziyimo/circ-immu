#!/usr/bin/env Rscript
#setwd("C:\\Users\\Armin\\ws\\projects\\virus\\covid_branch\\circ-immu\\iter_test")

suppressWarnings(suppressMessages(library("pso", character.only = TRUE)))
suppressWarnings(suppressMessages(library("parallel")))
suppressWarnings(suppressMessages(library("binom")))
source("Covid_state_joint.R")

### Parse user args
args <- commandArgs(trailingOnly=TRUE)
#argsLen <- length(args)
state_ls <- args[1]  # 2-letter state code
R0_mod <- args[2] # functional form of R0, options: day,hum,hd,sd,hsd
#h_rate, R0min, R0range,alpha_1, s_0, alpha_2, d_0
#seedstr <- "0.01:1.2:1:-200:0.5:-300:0.8"
nthr <- 1
# text file with a list of states to fit to
seed_sh <- as.numeric(strsplit(args[3], ":", fixed=TRUE)[[1]])  # h_rate, R0min, R0range, [R0_params]
init_inf <- as.numeric(args[4])
target_state <- args[5]


### Logging
time_stamp <- format(Sys.time(), "%m%d%H%m")

handle <- paste0(state_ls, "_", R0_mod,"_", time_stamp)
handle <- sub(" ", "", handle)
#sink(paste0("fit_results2/", handle, ".log"))
#/cat(paste0("fit_results/", handle, ".log"))

# Optimizer
optimizer <- "pso"
#restart_thd <- as.numeric(args[5]) # restart threshold for pso optimizer
#library(optimizer, character.only = TRUE)

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
  #cat(">>> Loading state data:", state_code, "<<<\n")
  census_pop <- state_pop$pop[state_pop$code==state_code]
  
  if (R0_mod %in% c("sun", "hs", "sd", "hsd")){
    # in US DST starts on 10 March (day69) and ends on 3. Nov (day 307)
    # Add 60 min to day 1-68
    s1 <- head(all_state_sun[[state_code]],68)
    # Get unchanged days 69-306
    s2 <-  head(tail(all_state_sun[[state_code]],297),238) - 60
    # Add 60 min to 307-365
    s3 <- tail(all_state_sun[[state_code]],59)
    # permanent DST
    s <- c(s1,s2,s3)
    sunob <- s/720 # 365 days, scaled to maximum 720 = 12 hours
    # To abolish DST
    # Just subtract 60 from s2 and do leave s1 and s3 unchanged
    #sunob <- all_state_sun[[state_code]]/720 # 365 days, scaled to maximum 720 = 12 hours
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
if (R0_mod %in% c("sun", "sd", "hsd")){
  s0_idx <- ifelse(R0_mod %in% c("sun", "sd"), 2, 3)
  # Shared params: [CD s_0], hosp_rate, R0min, R0range, [other R0 params]
  # constrain the range of s_0
  
} else {
  # Shared params: h_rate, R0min, R0range, [R0 params]
  sh_low <- c(0.002, 0.8, 0.2, param_bounds[[R0_mod]]$low)
  sh_high <- c(0.01, 1.5, 2.5, param_bounds[[R0_mod]]$high)
}

### TEST

#state_dat <- city_data$WY
#sh_vec <- seed_sh
#too <- optimize(ss_wrapper, c(0, 1), st_e=state_dat, sh_p=sh_vec)
#ss_prms <- ss_low + too$minimum*(ss_high-ss_low)
#ss_ls <- ss_prms
#ss_p <- ss_ls
#st_e <- state_dat
#orig_prms <- sh_vec
#sh_orig <- orig_prms
#state_wrapper(st_e, ss_p, sh_orig)
#ss_fit <- ss_ls
#psoptim((seed_sh-sh_low)/(sh_high-sh_low), all_state_negLL_test, ss_prms=ss_fit,
#        lower=rep(0,length(sh_low)), upper=rep(1,length(sh_low)),
#        control=list(trace=1, REPORT=5, maxit=10000, trace.stats=FALSE, maxit.stagnate=200))



#state_wrapper(city_data$CA, c(1.661004e-05), c(rep(0.5,9),0.1,1.2,2.5,-20,-30,0.5770212))

#state_wrapper(city_data$CA, c(1.661004e-05), c(0.4671362,0.5365699,0.5053176,0.6486928,0.5940593,0.4785501,0.5320025, 0.5327572, 0.6067473, 0.009778883, 1.036737, 2.067482, -500, -341.6809, 0.6026052))
### TEST

state_wrapper <- function(state_elem, ss_prms, sh_prms, target_state){
  #state_elem: element in the state_data list
  #ss_prms: vector of state specific parameters
  #sh_prms: vector of shared parameters
  
  #sink("tmp_negLL.log",append=TRUE)
  #cat("Full parms","\n",sh_prms,"\n")
  #cat("alpha: ", sh_prms[(no_cdv+1):(no_cdv+3+s0_idx-1)],"\n")
  #cat("s0: ", sh_prms[state_elem$cdv],"\n")
  #cat("other: ",sh_prms[-1:-(no_cdv+3+s0_idx-1)],"\n")
  #print(state_elem)
  if (R0_mod %in% c("sun", "sd", "hsd")){
    # Init, alpha,R0min,R0range,s_0,other
    state_prms <- c(ss_prms, sh_prms[(no_cdv+1):(no_cdv+3+s0_idx-1)], sh_prms[state_elem$cdv], sh_prms[-1:-(no_cdv+3+s0_idx-1)])
    cat(state_prms,"\n")
  } else {
    state_prms <- c(ss_prms, sh_prms)
  }
  # inf_init,alpha,R0min,R0range,s,s0,d,d0
  # 0.002915468 0.005872427 1.257415 2.604619 -147.4783 0.4605626 -97.57329 0.03931712 
  neg_LL <- binom_seird_p(state_prms, 
                   paste0("R0_", R0_mod),
                   state_elem$var,
                   state_elem$epi,
                   state_elem$pop,
                   target_state)
  
  #print(neg_LL)
  #sink()
  return(neg_LL)  
}

state_wrapper(get(target_state, city_data),init_inf,seed_sh, target_state)

