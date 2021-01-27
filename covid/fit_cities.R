#!/usr/bin/env Rscript

suppressWarnings(suppressMessages(library("pso", character.only = TRUE)))
suppressWarnings(suppressMessages(library("parallel")))
suppressWarnings(suppressMessages(library("binom")))
source("Covid_model_joint.R")

### Parse user args
args <- commandArgs(trailingOnly=TRUE)
#argsLen <- length(args)
city_ls <- args[1]    # 2-letter state code
R0_mod <- args[2] # functional form of R0, options: day,hum,hd,sd,hsd
nthr <- args[3]
# text file with a list of states to fit to

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

### Load covid data
dt <- read.csv("data/covid_cities_corrected.csv", sep = ",", stringsAsFactors = FALSE,header = TRUE, check.names = FALSE)
dt$Date <- as.Date(dt$Date)
cities_df <- subset(dt, Type == "Deaths" & format.Date(Date, "%Y")=="2020")
# Remove leap year Feb 29
cities_df <- cities_df[cities_df[["Date"]] != "2020-02-29", ]
# Convert to Day of Year format
cities_df$Date <- as.numeric(strftime(cities_df$Date, format = "%j"))
# Adjust Day of Year to remove leap day
cities_df$Date <- cities_df$Date - ifelse(cities_df$Date >= 60, 1, 0)

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

for (city_i in seq(no_cities)){
  city_code <- cities$V1[city_i]
  cat(">>> Loading state data:", city_code, "<<<\n")
  census_pop <- pops$Population[pops$Name==city_code]
  
  # Get first case and last day
  city_df <- subset(cities_df, Name == city_code)
  # Get start day of pandemic
  first_case <- head(subset(city_df, Cases_New > 0),1)
  # First exposure is 15 days before first death
  first_case_t <- first_case$Date
  # Last day of data
  last_day <- ifelse(max(city_df$Date) > 358,358,max(city_df$Date))
  # This assumes first infection is always after 15. Jan to work
  city_df <- subset(city_df , Date >= first_case_t & Date <= last_day)
  
  
  if (R0_mod %in% c("sun", "hs", "sd", "hsd")){
    city_sun <- subset(sun, City == city_code)
    sunob <- subset(city_sun, Date >= first_case_t & Date <= last_day)
    sunob  <- sunob$Sunrise
  }
  if (R0_mod %in% c("day", "hd", "sd", "hsd")){
    city_day <- subset(day, City == city_code)
    dayob <- subset(city_day, Date >= first_case_t & Date <= last_day)
    dayob  <- dayob$Day
  }
  if (R0_mod %in% c("hum", "hs", "hd", "hsd")){
    city_hum <- subset(hum, City == city_code)
    humob <- subset(city_hum, Date >= first_case_t & Date <= last_day)
    humob  <- humob$Humidity
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
  city_data[[city_code]] <- list(idx=city_i, pop=census_pop, var=varob, epi=city_df)
}

# Boundaries for inital infected/recovered
inf_low <- 0.00001
inf_high <- 0.05
rec_low <- 0
rec_high <- 0
alpha_low <- 0.01
alpha_high <- 0.01
R0range_low <- 0.2
R0range_high <- 5
R0min_low <- 0.6
R0min_high <- 1.6

# Positional params
# Lower bounds
prms_low <- c(rep(inf_low, no_cities), #1: seed_prop - starting infected prop 
              rep(rec_low, no_cities), #2: seed_prop - starting recovered prop 
              rep(alpha_low, no_cities), #3: alpha - fatality rate
              R0min_low, #4: R0_min
              R0range_low, #5: R0_range
              param_bounds[[R0_mod]]$low) #6: Model prms
# Upper bounds
prms_high <- c(rep(inf_high, no_cities), 
               rep(rec_high, no_cities),
               rep(alpha_high, no_cities),
               R0min_high, 
               R0range_high, 
               param_bounds[[R0_mod]]$high)

## PSO optimization all cities
# Extract each cities param set from super param set
# Use param positions to extract, based on order in prms_low/prms_high
get_city_prms <- function(prms_vec, city_idx){
  return(c(prms_vec[city_idx], # city infection rate
           prms_vec[no_cities+city_idx], # city recovered rate
           prms_vec[(no_cities*2)+city_idx], # city death rate
           prms_vec[-1:-(3*no_cities)])) # All other R0 and model params
}


negLL_wrapper <- function(city_elem, p_vec){
  neg_LL <- binom_seird(get_city_prms(p_vec, city_elem$idx), 
                   paste0("R0_", R0_mod),
                   city_elem$var,
                   city_elem$epi,
                   city_elem$pop)
  # Plot small fraction of results to help monitor fitting progress
  if(runif(1) > 0.9995){ 
      binom_seird_p(get_city_prms(p_vec, city_elem$idx), 
                    paste0("R0_", R0_mod),
                    city_elem$var,
                    city_elem$epi,
                    city_elem$pop)
  }
  return(neg_LL)
}
cat(detectCores(), "cores seen; limit to", nthr, "cores\n")
clstr <- makeCluster(nthr, type = "FORK") # limits the number of cores to use

all_state_negLL <- function(norm_prms){ # everything else read as global variable
  orig_prms <- prms_low + norm_prms*(prms_high-prms_low)
  states_negLL <- parSapply(cl = clstr, X=city_data, FUN=negLL_wrapper, p_vec=orig_prms)
  return(sum(states_negLL))
}

cat(">> Fitting with R0 model:", R0_mod, "\nLower:", prms_low, "\nUpper:", prms_high, "\n")
restart_thd <- c(1.1e-2, 1.1e-3, 1.1e-4)
seed <- rep(NA,length(prms_low))

for (thrhld in restart_thd){
  cat(">> Restart threshold:", thrhld, "\n")
  # optimizer only takes in dimensions equal to num prms
  # just uses 0-1 for each prm
  oo <- psoptim(seed, all_state_negLL, lower=rep(0,length(prms_low)), upper=rep(1,length(prms_low)),
                control=list(trace=1, REPORT=5, maxit=10000, trace.stats=FALSE, maxit.stagnate=500,
                             max.restart=5, reltol=thrhld))
  seed <- oo$par
  cat(cities$V1); cat("; hrate; R0_min; R0_range; [R0_prms]\n")
  cat(prms_low + oo$par*(prms_high-prms_low))
  cat("\n")
}
stopCluster(clstr)
saveRDS(oo, file = paste0("fit_results/", handle, ".run1.rds"))

### Fit with different s_0
if (R0_mod %in% c("sun", "sd", "hsd")){
  s0_idx <- ifelse(R0_mod %in% c("sun", "sd"), 2, 3)
  prms_low <- c(prms_low[1:(3*no_cities)], # initial city specific prms
                rep(prms_low[3*no_cities+2+s0_idx],no_cities), # new city specific s0
                prms_low[(3*no_cities+1):(3*no_cities+2)],
                prms_low[-1:(-3*no_cities-2)][-s0_idx])
  prms_high <- c(prms_high[1:(3*no_cities)], 
                 rep(prms_high[3*no_cities+2+s0_idx], no_cities),
                 prms_high[(3*no_cities+1):(3*no_cities+2)], 
                 prms_high[-1:(-3*no_cities-2)][-s0_idx])
  seed <- c(seed[1:(3*no_cities)], 
            rep(seed[3*no_cities+2+s0_idx], no_cities),
            seed[(3*no_cities+1):(3*no_cities+2)], 
            seed[-1:(-3*no_cities-2)][-s0_idx])
  
  get_city_prms <- function(prms_vec, city_idx){
    return(c(prms_vec[city_idx], # state-specific seed_prop
             prms_vec[no_cities+city_idx], # state-specific recov seed
             prms_vec[no_cities*2+city_idx], # state-specific death rate
             prms_vec[(4*no_cities+1):(4*no_cities+1+s0_idx)], # shared: R0_min, R0_range, [R0 prms b4 s_0]
             prms_vec[3*no_cities+city_idx], #)) # state_specific s_0
             prms_vec[-1:-(4*no_cities+1+s0_idx)])) # shared: [R0 prms after s_0]
  }
  
  clstr <- makeCluster(nthr, type = "FORK") 
  cat(">> Refining R0 model:", R0_mod, "\nLower:", prms_low, "\nUpper:", prms_high, "\n")
  
  for (thrhld in restart_thd){
    cat(">> Restart threshold:", thrhld, "\n")
    oo <- psoptim(seed, all_state_negLL, lower=rep(0,length(prms_low)), upper=rep(1,length(prms_low)),
                  control=list(trace=1, REPORT=5, maxit=10000, trace.stats=FALSE, maxit.stagnate=500,
                               max.restart=5, reltol=thrhld))
    seed <- oo$par
    cat(cities$V1); cat("; hrate; R0_min; R0_range; [R0_prms]\n")
    cat(prms_low + oo$par*(prms_high-prms_low))
    cat("\n")
  }
  stopCluster(clstr)
}

saveRDS(oo, file = paste0("fit_results/", handle, ".rds"))
