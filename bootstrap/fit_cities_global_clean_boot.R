#!/usr/bin/env Rscript
#setwd("C:\\Users\\Armin\\ws\\projects\\virus\\covid_branch\\circ-immu\\iter_test")
suppressWarnings(suppressMessages(library("pso", character.only = TRUE)))
suppressWarnings(suppressMessages(library("parallel")))
suppressWarnings(suppressMessages(library("binom")))
source("Covid_state_joint.R")

args <- commandArgs(trailingOnly=TRUE)
city_ls <- args[1]  # 2-letter state code
R0_mod <- args[2] # functional form of R0, options: day,hum,hd,sd,hsd
#*# `seed_sh` needs an update, should be loaded from the best fitting results
seed_file <- args[3]  # h_rate, R0min, R0range, [R0_params]
nthr <- args[4]
death_rate <- 0.001

load(seed_file)
#
### Logging
time_stamp <- format(Sys.time(), "%m%d%H%m")

handle <- paste0(city_ls, "_", R0_mod,"_", time_stamp)
handle <- sub(" ", "", handle)
#sink(paste0("fit_results2/", handle, ".log"))
cat(paste0("fit_results/", handle, ".log"))

# Optimizer
optimizer <- "pso"
#restart_thd <- as.numeric(args[5]) # restart threshold for pso optimizer
#library(optimizer, character.only = TRUE)

### Load population size data
cities <- read.delim(city_ls, header = FALSE,stringsAsFactors=F)
no_cities <- nrow(cities)
pops<- read.csv("data/covid_cities_populations.csv", sep = ",", stringsAsFactors = FALSE,header = TRUE, check.names = FALSE)

#NEW
state2reg <- read.csv("data/city2censusReg_clean.csv")
cdv <- unique(state2reg$cregion)
city <- unique(state2reg$city)
no_cdv <- length(cdv)
state2reg$CDV.Code <- match(state2reg$city, city)
# END NEW

boundaries <- read.csv("data/cities_max_min_s0.csv", sep = ",", stringsAsFactors = FALSE,header = TRUE, check.names = FALSE)
state2reg <- merge(state2reg,boundaries, by=c("city"))
state2reg <- state2reg[order(state2reg$CDV.Code),]
state2reg$min_s0 <- ceiling(state2reg$min_s0*100) / 100
state2reg$max_s0 <- floor(state2reg$max_s0*100) / 100

### Load covid data
dt <- read.csv("data/covid_cities_corrected.csv", sep = ",", stringsAsFactors = FALSE,header = TRUE, check.names = FALSE)
dt$Date <- as.Date(dt$Date)
cities_df <- subset(dt, Type == "Deaths")
# Convert to Day of Year format
cities_df$Date <- as.numeric(as.Date(cities_df$Date, format="%Y-%m-%d") - as.Date("2019-12-31", format="%Y-%m-%d"))


### Load environmental data

sun <- read.csv("data/cities_sunrise.csv",sep = ",", stringsAsFactors = FALSE,header = FALSE, na.strings="",check.names = FALSE)
hum <- read.csv("data/cities_humidity.csv",sep = ",", stringsAsFactors = FALSE,header = FALSE, na.strings="",check.names = FALSE)
day <- read.csv("data/cities_daylight.csv",sep = ",", stringsAsFactors = FALSE,header = FALSE, na.strings="",check.names = FALSE)

# Process sun data
names(sun) <- c("Date","City","Sunrise")
sun <- sun[with(sun, order(City, Date)),]
sun$Date <- as.Date(sun$Date, "%Y-%m-%d")
sun$Date <- as.numeric(strftime(sun$Date , format = "%j"))
#city_sun <- subset(sun, City == cityname)
sun$Sunrise <- sun$Sunrise/720
# Process hum data
names(hum) <- c("Date","City","Humidity")
hum <- hum[with(hum, order(City, Date)),]
hum$Date <- as.Date(hum$Date, "%Y-%m-%d")
hum$Date <- as.numeric(strftime(hum$Date , format = "%j"))
# Process day data
names(day) <- c("Date","City","Day")
day <- day[with(day, order(City, Date)),]
day$Date <- as.Date(day$Date, "%Y-%m-%d")
day$Date <- as.numeric(strftime(day$Date , format = "%j"))
day$Day <- day$Day/1440

city_data <- list()
seed_ss <- list()

for (city_i in seq(no_cities)){
  city_code <- cities$V1[city_i]
  cat(">>> Loading state data:", city_code, "<<<\n")
  census_pop <- pops$Population[pops$Name==city_code]
  
  # Get first case and last day
  city_df <- subset(cities_df, Name == city_code)
  # Get start day of pandemic
  #first_case <- head(subset(city_df, Cases_New > 0),1)
  # First exposure is 15 days before first death
  #first_case_t <- first_case$Date
  # Last day of data
  #cat(max(city_df$Date))
  #last_day <- ifelse(max(city_df$Date) > 396,396,max(city_df$Date))
  # This assumes first infection is always after 15. Jan to work
  #city_df <- subset(city_df , Date >= first_case_t & Date <= last_day)
  # remove unneeded
  keeps <- c("Date", "Cases_New")
  city_df <- city_df[keeps]
  names(city_df) <- c("date","deathIncrease_cor")
  
  if (R0_mod %in% c("sun", "hs", "sd", "hsd")){
    city_sun <- subset(sun, City == city_code)
    #sunob <- subset(city_sun, Date >= first_case_t & Date <= last_day)
    sunob  <- city_sun$Sunrise
  }
  if (R0_mod %in% c("day", "hd", "sd", "hsd")){
    city_day <- subset(day, City == city_code)
    #dayob <- subset(city_day, Date >= first_case_t & Date <= last_day)
    dayob  <- city_day$Day
  }
  if (R0_mod %in% c("hum", "hs", "hd", "hsd")){
    city_hum <- subset(hum, City == city_code)
    #humob <- subset(city_hum, Date >= first_case_t & Date <= last_day)
    humob  <- city_hum$Humidity
  }
  
  if (R0_mod == "cos"){
    # baseline R0 outputs 364 days
    # other models output t days based on death data
    # need to add merge step to reconcile in SEIRDvar_pred()
    varob <- list(seq(365))
    #varob <- list(seq(length(city_df$Date)))

  } else if (R0_mod == "hum"){
    varob <- list(humob)
  } else if (R0_mod == "day"){
    varob <- list(dayob)
  } else if (R0_mod == "hd"){
    varob <- list(humob, dayob)
  } else if (R0_mod == "sd"){
    varob <- list(sunob, dayob)
  } else if (R0_mod == "hsd"){
    varob <- list(humob, sunob, dayob)
  } else if (R0_mod == "sun"){
    varob <- list(sunob)
  } else if (R0_mod == "hs"){
    varob <- list(humob, sunob)
  }
  
  # TO DO
  # 1) Add ordering by date city to city_df
  # 2) When combining with predictions use merge on date
  #    This ensures that if a date is missing, it won't cause a row mismatch

  #state_df <- subset(covid_df, state == state_code)
  #state_df <- state_df[order(state_df$date),]
  #state_df <- state_df[!is.na(state_df$hospitalizedCurrently), ]
  #state_df <- state_df[state_df$date <= 365, ]
  #state_hos <- subset(state_df, select = c("date", "hospitalizedCurrently"))
  
  # city df contains cols: Date,Type,Name,Cases_New
  #city_data[[city_code]] <- list(idx=city_i, pop=census_pop, var=varob, epi=city_df)
  city_data[[city_code]] <- list(idx=city_i, pop=census_pop, var=varob, epi=city_df, cdv=state2reg$CDV.Code[state2reg$city == city_code])
}

if (R0_mod %in% c("sun", "sd", "hsd")){
  s0_idx <- ifelse(R0_mod %in% c("sun", "sd"), 2, 3)
  # Shared params: [CD s_0], hosp_rate, R0min, R0range, [other R0 params]
  
  #seed_sh <- c(rep(seed_sh[2+s0_idx], no_cdv), seed_sh[-(2+s0_idx)]) #*# this becomes deprecated
  sh_low <- c(as.numeric(state2reg$min_s0), 0.001, 0.8, 0.6, param_bounds[[R0_mod]]$low[-s0_idx])
  sh_high <- c(as.numeric(state2reg$max_s0), 0.001, 1.2, 1.3, param_bounds[[R0_mod]]$high[-s0_idx])
  # constrain the range of s_0
  #seed_sh <- c(rep(seed_sh[3+s0_idx], no_cities), seed_sh[-(3+s0_idx)])

  #*# `seed_sh` (should be loaded from file already) is the best fit shared parameters, but we will not use that to seed the bootstrap runs (this parameter name is somewhat of a misnomer, but I am horrible and didn't bother changing it)
  #*# Instead we introduce a different seed `seed_guide`, this is basically the best fit parameter, but with all the s_0 values set at 0.5
  #*# I found that in reality this works pretty well, better than not seeding the runs at all
  bad_alpha_pos <- no_cities + 1
  #seed_sh <- seed_sh[-bad_alpha_pos]
  seed_guide <- c(rep(0.59, no_cities), seed_sh[-1:-no_cities]) #*# seed all s_0 with 0.5
  #seed_guide <- seed_guide[-bad_alpha_pos]

} else {
  #*# this is probably also deprecated, we only care about bootstrapping for sd models
  # Shared params: h_rate, R0min, R0range, [R0 params]
  sh_low <- c(0.8, 1.2, param_bounds[[R0_mod]]$low)
  sh_high <- c(1.2, 1.6, param_bounds[[R0_mod]]$high)
}
print(sh_low)
print(sh_high)
print(seed_sh)
print(seed_guide)

state_wrapper <- function(state_elem, ss_prms, sh_prms){
  #state_elem: element in the state_data list
  #ss_prms: vector of state specific parameters
  #sh_prms: vector of shared parameters
  
  #sink("tmp_negLL.log",append=TRUE)
  #cat("Full parms","\n",sh_prms,"\n")
  #cat("alpha: ", sh_prms[(no_cdv+1):(no_cdv+3+s0_idx-1)],"\n")
  #cat("s0: ", sh_prms[state_elem$cdv],"\n")
  #cat("other: ",sh_prms[-1:-(no_cdv+3+s0_idx-1)],"\n")
  if (R0_mod %in% c("sun", "sd", "hsd")){
    # Init, alpha,R0min,R0range,s_0,other
    #state_prms <- c(ss_prms,death_rate, sh_prms[(no_cities+1):(no_cities+s0_idx+1)], sh_prms[state_elem$cdv], sh_prms[-1:-(no_cities+s0_idx+1)])
    state_prms <- c(ss_prms,0.001, sh_prms[(no_cities+2):(no_cities+3+s0_idx-1)], sh_prms[state_elem$cdv], sh_prms[-1:-(no_cities+3+s0_idx-1)])
    #cat(state_prms,"\n")
  } else {
    state_prms <- c(ss_prms,sh_prms)
  }

  #print(s0_idx)
  #print(sh_prms)
  #print(state_elem$pop)
  # inf_init,alpha,R0min,R0range,s,s0,d,d0
  # 0.002915468 0.005872427 1.257415 2.604619 -147.4783 0.4605626 -97.57329 0.03931712 
  neg_LL <- binom_seird(state_prms, 
                        paste0("R0_", R0_mod),
                        state_elem$var,
                        state_elem$epi,
                        state_elem$pop)
  
  #print(neg_LL)
  #sink()
  return(neg_LL)  
}
#*# probably don't need the functions for fitting initial infections anymore
# fit_state <- function(state_dat, sh_vec){ # no seeding for 1D optimization
#   oo <- optimize(ss_wrapper, c(0, 1), st_e=state_dat, sh_p=sh_vec)
  
#   return(ss_low + oo$minimum*(ss_high-ss_low))
# }

# ss_wrapper <- function(ss_norm, st_e, sh_p){
#   ss_p <- ss_low + ss_norm*(ss_high-ss_low)
#   return(state_wrapper(st_e, ss_p, sh_p))
# }

## Function for fitting shared parameters
#*# this is the key change to the script, modifying likelihood calculation
all_state_negLL <- function(norm_sh_prms, ss_prms, boot_samps){ #*# extra argument 'boot_samps' to the likelihood function
  #*# We calculate the per-state likelihood as usual, but return the likelihood of bootstrap samples
  #*# 'boot_samps' is a vector of indices 
  orig_prms <- sh_low + norm_sh_prms*(sh_high-sh_low)
# set deaths
  states_negLL <- parSapply(cl = clstr, X=city_data, FUN=sh_wrapper, sh_orig=orig_prms, ss_ls=ss_prms)
  return(sum(states_negLL[boot_samps])) #*# return only the likelihood calculations for the bootstrap samples
}


sh_wrapper <- function(st_e, sh_orig, ss_ls){
  ss_p <- ss_ls[[st_e$idx]]
  return(state_wrapper(st_e, ss_p, sh_orig))
}

no_bsit <- 10 #*# manually set the number of bootstrap sampling
resamp_params <- matrix(nrow = no_bsit, ncol = length(sh_low)) #*# pre-allocate a matrix of bootstrapped parameters

cat(detectCores(), "cores seen; limit to", nthr, "cores\n")
clstr <- makeCluster(nthr, type = "FORK") # limits the number of cores to use

for (iter in seq(no_bsit)){ #*# sequentially run bootstrap iterations, we can submit multiple jobs in parallel, so far I have been doing 10 jobs, each with 10 bootstrap runs

  #*# we don't fit the initial infections anymore, they are fixed
  #ss_fit <- parLapply(cl = clstr, X = city_data, fun = fit_state, sh_vec = seed_sh)
  print(seed_sh)
  bssample <- sample(seq(no_cities), size=no_cities, replace = TRUE) #*# sampling with replacement the state indices
  opt_negLL <- all_state_negLL((seed_sh-sh_low)/(sh_high-sh_low), ss_fit, bssample) #*#
  #*# sanity check: the likelihood of the "best-fit" parameters given the bootstrap sample
  cat(">>", R0_mod, "Resamp", iter, ":", bssample, "\n") #*#
  cat(">> Opt param @ resamp:", opt_negLL, "\n") #*#

  #*# note we are seeding with 'seed_guide'
  sh_fit <- psoptim((seed_guide-sh_low)/(sh_high-sh_low), all_state_negLL, ss_prms=ss_fit, boot_samps=bssample, #*# the additional bootstrap sample parameter
                    lower=rep(0,length(sh_low)), upper=rep(1,length(sh_low)),
                    control=list(trace=1, REPORT=5, maxit=5000, trace.stats=FALSE, maxit.stagnate=200))
  resamp_params[iter,] <- sh_low + sh_fit$par*(sh_high-sh_low) #*# store in the matrix directly
  #*# bunch of updated bookkeeping/sanity checks
  cat(">> Resamp_Iter", iter, "negLL:", sh_fit$value, "(ref: ", opt_negLL, ")\n") #*# this is for checking that the parameters fitted *specifically* to the bootstrapped sample indeed result in a better fit (lower negative log likelihood) that the original best-fit parameters
  cat(">> Resamp_params:", resamp_params[iter,], "\n")
}

save(resamp_params, file = paste0("fit_results/", handle,"_cleancit_BS.rds")) #*# whatever name
stopCluster(clstr)


