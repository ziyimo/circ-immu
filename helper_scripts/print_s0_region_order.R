#!/usr/bin/env Rscript
#setwd("C:\\Users\\Armin\\ws\\projects\\virus\\covid_branch\\circ-immu\\iter_test")

suppressWarnings(suppressMessages(library("pso", character.only = TRUE)))
suppressWarnings(suppressMessages(library("parallel")))
suppressWarnings(suppressMessages(library("binom")))
source("Covid_state_joint.R")

### Parse user args
args <- commandArgs(trailingOnly=TRUE)
#argsLen <- length(args)
city_ls <- args[1]  # 2-letter state code
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

handle <- paste0(city_ls, "_", R0_mod,"_", time_stamp)
handle <- sub(" ", "", handle)
#sink(paste0("fit_results2/", handle, ".log"))

# Optimizer
optimizer <- "pso"
#restart_thd <- as.numeric(args[5]) # restart threshold for pso optimizer
#library(optimizer, character.only = TRUE)

### Load population size data
cities <- read.delim(city_ls, header = FALSE,stringsAsFactors=F)
no_cities <- nrow(cities)
pops<- read.csv("data/covid_cities_populations.csv", sep = ",", stringsAsFactors = FALSE,header = TRUE, check.names = FALSE)

#NEW
state2reg <- read.csv("data/city2censusReg.csv")
cdv <- unique(state2reg$region)
no_cdv <- length(cdv)
state2reg$CDV.Code <- match(state2reg$region, cdv)
# END NEW

### Load covid data
dt <- read.csv("data/covid_cities_corrected.csv", sep = ",", stringsAsFactors = FALSE,header = TRUE, check.names = FALSE)
dt$Date <- as.Date(dt$Date)
#cities_df <- subset(dt, Type == "Deaths" & format.Date(Date, "%Y")=="2020")
cities_df <- subset(dt, Type == "Deaths")
# Remove leap year Feb 29
#cities_df <- cities_df[cities_df[["Date"]] != "2020-02-29", ]
# Convert to Day of Year format
cities_df$Date <- as.numeric(as.Date(cities_df$Date, format="%Y-%m-%d") - as.Date("2019-12-31", format="%Y-%m-%d"))
#cities_df$Date <- as.numeric(strftime(cities_df$Date, format = "%j"))
# Adjust Day of Year to remove leap day
#cities_df$Date <- cities_df$Date - ifelse(cities_df$Date >= 60, 1, 0)

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
  census_pop <- pops$Population[pops$Name==city_code]
  
  # Get first case and last day
  city_df <- subset(cities_df, Name == city_code)
  # Get start day of pandemic
  first_case <- head(subset(city_df, Cases_New > 0),1)
  # First exposure is 15 days before first death
  first_case_t <- first_case$Date
  # Last day of data
  last_day <- ifelse(max(city_df$Date) > 396,396,max(city_df$Date))
  # This assumes first infection is always after 15. Jan to work
  city_df <- subset(city_df , Date >= first_case_t & Date <= last_day)
  # remove unneeded
  keeps <- c("Date", "Cases_New")
  city_df <- city_df[keeps]
  names(city_df) <- c("date","deathIncrease_cor")
  
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
  #city_data[[city_code]] <- list(idx=city_i, pop=census_pop, var=varob, epi=city_df)
  city_data[[city_code]] <- list(idx=city_i, pop=census_pop, var=varob, epi=city_df, cdv=state2reg$CDV.Code[state2reg$city == city_code])
}


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
    #cat(state_prms,"\n")
  } else {
    state_prms <- c(ss_prms, sh_prms)
  }
  #print(state_prms)
  #print(state_elem$var)
  #print(state_elem$epi)
  #print(state_elem$pop)
  #print(target_state)
  cat(target_state,"\t",ss_prms,"\t",sh_prms[state_elem$cdv],"\n") 
}

#cat(target_state,"\n")
#cat(init_inf,"\n")
#print(seed_sh)
#cat(R0_mod)
state_wrapper(get(target_state, city_data),init_inf,seed_sh, target_state)

# s0 value is state_elem$cdv
