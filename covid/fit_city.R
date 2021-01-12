#!/usr/bin/env Rscript

suppressWarnings(suppressMessages(library(DEoptim)))
suppressWarnings(suppressMessages(library(binom)))
source("Covid_model.R")

####################################################
# To do to improve covid fitting 18.12.2020        #
####################################################

# For deaths condition death rate on vector of variable empirical death rate per day
# Calculate this as smoothed cases / deaths per day
# Add dampening factor to R0 to account for social distancing - this may not be necessary if truncating data to leave out first peak

####################################################
# To do to improve covid fitting 01.01.2021        #
####################################################

# Fit city-specific scaling factor for death rate


##############################################################################


#setwd("C:\\Users\\Armin\\ws\\projects\\virus\\push_dir\\circ-immu\\covid")

### Parse user args
args <- commandArgs(trailingOnly=TRUE)
argsLen <- length(args)
cityname <- args[1]   # 2-letter state code
R0_mod <- args[2]       # functional form of R0, options: day,hum,hd,sd,hsd

if (argsLen == 2){
  optim_arg <- "optim"
  cat(">> Optimization run started \n")
} else if (argsLen == 3){
  optim_arg <- "plot"
  prm_ls <- args[3] # comma sep list of params e.g. "15.4,499,-73.14,0.789"
  prms <- as.numeric(unlist(strsplit(prm_ls, split=",")))
  cat(">> Plotting run started \n")
} else {
  stop('Please supply two args for optim or three args for plot')
}

### Logging
time_stamp <- format(Sys.time(), "%m%d%H")

handle <- paste0(cityname, R0_mod, time_stamp)
sink(paste0("fit_results/", handle, ".log"))

### Load population size data

pops<- read.csv("data/covid_cities_populations.csv", sep = ",", stringsAsFactors = FALSE,header = TRUE, check.names = FALSE)
census_pop = subset(pops, Name==cityname)$Population

### Load covid data

dt <- read.csv("data/covid_cities_corrected.csv", sep = ",", stringsAsFactors = FALSE,header = TRUE, check.names = FALSE)
dt$Date <- as.Date(dt$Date)
# Keep only 2020 for selected city
city_df <- subset(dt, Type == "Deaths" & Name == cityname & format.Date(Date, "%Y")=="2020")
# Remove leap year Feb 29
city_df <- city_df[city_df[["Date"]] != "2020-02-29", ]
# Convert to Day of Year format
city_df$Date <- as.numeric(strftime(city_df$Date, format = "%j"))
# Adjust Day of Year to remove leap day
city_df$Date <- city_df$Date - ifelse(city_df$Date >= 60, 1, 0)
# Get start day of pandemic
#first_case <- head(subset(city_df, Cases_New > 0),1)
# First exposure is 15 days before first death
#first_case_t <- first_case$Date - 15
# Use 1st July as start of modelling
first_case_t <- 183
# Last day of data
last_day <- tail(city_df$Date,1)
# This assumes first infection is always after 15. Jan to work
city_df <- subset(city_df , Date >= first_case_t)

### Load environmental data

sun <- read.csv("data/cities_sunrise.csv",sep = ",", stringsAsFactors = FALSE,header = FALSE, na.strings="",check.names = FALSE)
hum <- read.csv("data/cities_daylight.csv",sep = ",", stringsAsFactors = FALSE,header = FALSE, na.strings="",check.names = FALSE)
day <- read.csv("data/cities_daylight.csv",sep = ",", stringsAsFactors = FALSE,header = FALSE, na.strings="",check.names = FALSE)

# Process sun data
names(sun) <- c("Date","City","Sunrise")
sun$Date <- as.Date(sun$Date, "%Y-%m-%d")
sun$Date <- as.numeric(strftime(sun$Date , format = "%j"))
city_sun <- subset(sun, City == cityname)
city_sun$Sunrise <- city_sun$Sunrise/720
#sunob <- city_sun[, -c(2)]
sunob <- subset(city_sun, Date >= first_case_t & Date <= last_day)
sunob  <- sunob$Sunrise
# Process hum data
names(hum) <- c("Date","City","Humidity")
hum$Date <- as.Date(hum$Date, "%Y-%m-%d")
hum$Date <- as.numeric(strftime(hum$Date , format = "%j"))
city_hum <- subset(hum, City == cityname)
climob <- subset(city_hum, Date >= first_case_t & Date <= last_day)
climob  <- climob$Humidity
# Process day data
names(day) <- c("Date","City","Day")
day$Date <- as.Date(day$Date, "%Y-%m-%d")
day$Date <- as.numeric(strftime(day$Date , format = "%j"))
city_day <- subset(day, City == cityname)
city_day$Day <- city_day$Day/1440
dayob <- subset(city_day, Date >= first_case_t & Date <= last_day)
dayob  <- dayob$Day

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
}

## Optimization
cat(">> Fitting:", cityname, "; R0 model:", R0_mod, "\n")

if (optim_arg == "optim"){
  cat(">> Run DEoptim from scratch\n")
  evo_optim <- DEoptim(binom_seird, lower=param_bounds[[R0_mod]]$low, upper=param_bounds[[R0_mod]]$high,
                       control=DEoptim.control(trace = 5, reltol = 1e-5, itermax = 500, steptol = 100),
                       R0_model = paste0("R0_", R0_mod), var_obs = varob,
                       obs_deaths = city_df, pop_size = census_pop)
  saveRDS(evo_optim, file = paste0("fit_results/", handle, "_DE.rds"))
}
sink()

if (optim_arg == "plot"){
  binom_seird_p(prms = prms, R0_model = paste0("R0_", R0_mod), var_obs = varob,
                obs_deaths = city_df, pop_size = census_pop)
}

