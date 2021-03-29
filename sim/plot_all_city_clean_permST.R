#!/usr/bin/env Rscript
#setwd("C:\\Users\\Armin\\ws\\projects\\virus\\covid_branch\\circ-immu\\iter_test")
suppressWarnings(suppressMessages(library("pso", character.only = TRUE)))
suppressWarnings(suppressMessages(library("parallel")))
suppressWarnings(suppressMessages(library("binom")))
source("Covid_state_joint.R")


### Parse user args
args <- commandArgs(trailingOnly=TRUE)
#argsLen <- length(args)
city_ls <- args[1]    # 2-letter state code
R0_mod <- args[2] # functional form of R0, options: day,hum,hd,sd,hsd
nthr <- 1
# text file with a list of states to fit to
seedstr <- args[3]
seed_sh <- as.numeric(strsplit(seedstr, ":", fixed=TRUE)[[1]])  # h_rate, R0min, R0range, [R0_params]
init_inf <- as.numeric(args[4])
target_state <- args[5]

### Logging
time_stamp <- format(Sys.time(), "%m%d%H%m")

handle <- paste0(city_ls, "_", R0_mod,"_", time_stamp)
handle <- sub(" ", "", handle)
#sink(paste0("fit_results2/", handle, ".log"))
#cat(paste0("fit_results/", handle, ".log"))

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
'%ni%' <- Negate('%in%')
for (city_i in seq(no_cities)){
  city_code <- cities$V1[city_i]
  #cat(">>> Loading state data:", city_code, "<<<\n")
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
    permcit <- c(
    "Moscow",
    "Osaka",
    "Hiroshima",
    "Phoenix",
    "Lima",
    "Delhi",
    "Tokyo",
    "Sao Paulo",
    "Rio de Janeiro",
    "Tucson",
    "Honolulu"
    )
    # Times are for 2018
    #DST 84-301 +1
    eurcit <- c(
    "Birmingham",
    "London",
    "Cologne",
    "Hamburg",
    "Stuttgart",
    "Berlin",
    "Rome",
    "Madrid",
    "Stockholm")
    #DST 133-224 -1
    chile <- c("Santiago")
    #DST 91-301 +1
    mexico <- c("Mexico City")
    city_sun <- subset(sun, City == city_code)
    if (city_code %ni% permcit){
        if (city_code %in% eurcit){
            start = 83
            end = 300
            s1 <- head(city_sun$Sunrise,start)
            s2 <-  head(tail(city_sun$Sunrise,(365-start)),end - start) - (60 / 720)
            s3 <- tail(city_sun$Sunrise,(365-end))
        } else if (city_code %in% chile){
            start = 132
            end = 223
            s1 <- head(city_sun$Sunrise,start)
            s2 <-  head(tail(city_sun$Sunrise,(365-start)),end - start) + (60 / 720)
            s3 <- tail(city_sun$Sunrise,(365-end))

        } else if (city_code %in% mexico){
            start = 90
            end = 300
            s1 <- head(city_sun$Sunrise,start)
            s2 <-  head(tail(city_sun$Sunrise,(365-start)),end - start) - (60 / 720)
            s3 <- tail(city_sun$Sunrise,(365-end))

        } else if (city_code %in% c("Melbourne")){
            start = 90
            end = 279
            s1 <- head(city_sun$Sunrise,start) 
            s2 <-  head(tail(city_sun$Sunrise,(365-start)),end - start) + (60 / 720)
            s3 <- tail(city_sun$Sunrise,(365-end))
        
        } else {
            start = 69
            end = 307
            s1 <- head(city_sun$Sunrise,start)
            s2 <-  head(tail(city_sun$Sunrise,(365-start)),end - start) - (60 / 720)
            s3 <- tail(city_sun$Sunrise,(365-end))
        }

    # permanent DST
    s <- c(s1,s2,s3)
    #sunob <- s/720 # 365 days, scaled to maximum 720 = 12 hours
    city_sun$Sunrise <- s
    #write.table(city_sun, "myDF.csv", sep = ",", append = T)
    #sunob <- subset(city_sun, Date >= first_case_t & Date <= last_day)
    sunob  <- city_sun$Sunrise


    # Do not modify
    # DST 89-299 +1
    # Euro DST
    # US DST
    # Other DST
    } else {
        sunob  <- city_sun$Sunrise
        #write.table(city_sun, "myDF.csv", sep = ",", append = T)
    }


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

#print(city_data)
# State specific params - I_init
ss_low <- c(1e-5)
ss_high <- c(0.1)

if (R0_mod %in% c("sun", "sd", "hsd")){
  s0_idx <- ifelse(R0_mod %in% c("sun", "sd"), 2, 3)
  # Shared params: [CD s_0], hosp_rate, R0min, R0range, [other R0 params]
  sh_low <- c(as.numeric(state2reg$min_s0), 0.001, 0.8, 0.6, param_bounds[[R0_mod]]$low[-s0_idx])
  sh_high <- c(as.numeric(state2reg$max_s0), 0.001, 1.2, 1.3, param_bounds[[R0_mod]]$high[-s0_idx])
  # constrain the range of s_0
  
} else {
  # Shared params: h_rate, R0min, R0range, [R0 params]
  sh_low <- c(0.001, 0.8, 1.2, param_bounds[[R0_mod]]$low)
  sh_high <- c(0.001, 1.2, 1.6, param_bounds[[R0_mod]]$high)
}

state_wrapper <- function(state_elem, ss_prms, sh_prms, target_state){
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
    state_prms <- c(ss_prms, sh_prms[(no_cities+1):(no_cities+3+s0_idx-1)], sh_prms[state_elem$cdv], sh_prms[-1:-(no_cities+3+s0_idx-1)])
    #cat(state_prms,"\n")
  } else {
    state_prms <- c(ss_prms, sh_prms)
  }

  #sink("tmp.log")
  #print(no_cdv)
  #print(s0_idx)
  #print(sh_prms)
  #print(state_elem$pop)
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

