library(ggplot2)
library(gridExtra)

clean <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"),
               text = element_text(size=20))

library(deSolve)

source("R0_mods.R")
source("SEIH_mod.R")

## Load all state data
state_pop <- read.delim("state_lv_data/state_pop.tsv")
all_state_sun <- read.csv("state_lv_data/state_daily_sunrise_2019.csv")
all_state_day <- read.csv("state_lv_data/state_daytime_2019.csv")
all_state_hum <- read.csv("state_lv_data/state_humidity_2014_2018.csv")
all_state_hum <- all_state_hum[all_state_hum$X != "2016-02-29", ] # get rid of leap year Feb 29

########
# covid_df <- read.csv("state_lv_data/covid_tracking_proj.csv")
# covid_df <- subset(covid_df, select = c(date, state, hospitalized, hospitalizedCumulative, hospitalizedCurrently, hospitalizedIncrease))
# saveRDS(covid_df, "state_lv_data/state_hospitalization.rds")

covid_df <- readRDS("state_lv_data/state_hospitalization.rds")

covid_df$date <- as.numeric(as.Date(covid_df$date, format="%Y-%m-%d") - as.Date("2019-12-31", format="%Y-%m-%d"))
#all_days <- covid_df$date[!is.na(covid_df$hospitalizedCurrently)]
#table(all_days) # Mar 13, day 73, national emergency

state_eg <- "TX"
census_pop <- state_pop$pop[state_pop$code==state_eg]
sunob <- all_state_sun[[state_eg]]/720
dayob <- all_state_day[[state_eg]]/1440
climob <- all_state_hum[[state_eg]] # multiple years
climob <- matrix(climob, nrow = length(climob)/365, ncol = 365, byrow = TRUE)
climob <- colMeans(climob) # 365 days

state_df <- subset(covid_df, state == state_eg)
#rownames(state_df) <- state_df$date
state_df <- state_df[order(state_df$date),]
state_df <- state_df[!is.na(state_df$hospitalizedCurrently), ]
state_df <- state_df[state_df$date <= 365, ]
state_hos <- subset(state_df, select = c("date", "hospitalizedCurrently"))
########

run_instance <- function(prms, R0_model, var_obs){
  seed_prop <- prms[1]
  hrate <- prms[2]
  R0_min <- prms[3]
  R0_range <- prms[4]
  R0_prms <- prms[-1:-4]
  
  R0_vals <- do.call(R0_model, c(var_obs, as.list(c(R0_prms, R0_min+R0_range, R0_min))))
  mod_pred <- run_SEIH(census_pop, seed_prop, hrate, R0_model, c(R0_prms, R0_min+R0_range, R0_min), var_obs)
  result <- merge(mod_pred, state_hos, by.x = "time", by.y = "date")
  cat(-sum(dpois(result$hospitalizedCurrently, result$H, log = TRUE)))
  
  p1 <- ggplot(data = data.frame("time"= 73:365, "R0"= R0_vals[73:365])) +
    geom_line(aes(x = time, y=R0), color = "black") + clean
  
  p4 <- ggplot(data = mod_pred) + 
    geom_line(aes(x = time, y=E), color = "#CC79A7") +
    geom_line(aes(x = time, y=I), color = "#D55E00") +
    geom_line(aes(x = time, y=H), color = "#0072B2") + clean
  
  p2 <- ggplot(data = mod_pred) + 
    geom_line(aes(x = time, y=S), color = "#000000") +
    geom_line(aes(x = time, y=E), color = "#CC79A7") +
    geom_line(aes(x = time, y=I), color = "#D55E00") +
    geom_line(aes(x = time, y=H), color = "#0072B2") + #clean
    geom_line(aes(x = time, y=R), color = "#009E73") + clean
  
  p3 <- ggplot(data = result) + 
    geom_line(aes(x = time, y=hospitalizedCurrently), color = "#D55E00") +
    geom_line(aes(x = time, y=H), color = "#0072B2") + clean
  
  grid.arrange(p1, p2, p3, p4, ncol = 2)
  #return(neg_log_L)
}

result <- readRDS("fit_results/states_trial8.tsv_sd_0114_pso.rds")
opt_prms <- result$par
fit_prms <- c(opt_prms[5], opt_prms[-1:-8])

fit_prms <- c(0.006788, 0.026612, 0.970366, 3.029638, -268.792228, 0.556643, -14.620948, 0.120207)
run_instance(fit_prms, "R0_sd", list(sunob, dayob))
fit_prms <- c(0.004630, 0.012029, 1.075839, 1.23119, -62.059135, 0.379025)
run_instance(fit_prms, "R0_day", list(dayob))
# fit_prms <- c(0.012906, 0.010619, 1.109247, 0.943467, 319.942628)
# run_instance(fit_prms, "R0_cos", list(73:365))
