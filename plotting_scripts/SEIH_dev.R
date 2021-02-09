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

#### Data curation ####
covid_df <- read.csv("state_lv_data/covid_tracking_proj.csv")
covid_df <- subset(covid_df, select = c(date, state, hospitalized, hospitalizedCumulative, hospitalizedCurrently, hospitalizedIncrease))
saveRDS(covid_df, "state_lv_data/state_hospitalization.rds")
#######################

covid_df <- readRDS("state_lv_data/state_hospitalization.rds")

covid_df$date <- as.numeric(as.Date(covid_df$date, format="%Y-%m-%d") - as.Date("2019-12-31", format="%Y-%m-%d"))
#all_days <- covid_df$date[!is.na(covid_df$hospitalizedCurrently)]
#table(all_days) # Mar 13, day 73, national emergency

#### Curate census region data ####
state2reg <- read.csv("state_lv_data/state2censusReg.csv")
state2reg <- subset(state2reg, select=-c(State))
covid_df <- merge(covid_df, state2reg, by.x = "state", by.y = "State.Code")

hosp_census_reg <- aggregate(hospitalizedCurrently ~ date+Region, covid_df, sum)
reg_mask <- hosp_census_reg$Region == "South"
plot(hosp_census_reg$date[reg_mask], hosp_census_reg$hospitalizedCurrently[reg_mask])
saveRDS(hosp_census_reg, "state_lv_data/censusReg_hospitalization.rds")

state_pop <- merge(state_pop, state2reg, by.x = "code", by.y = "State.Code")
cr_pop <- aggregate(pop ~ Region, state_pop, sum)
write.table(cr_pop, file = "state_lv_data/cr_pop.tsv", quote=FALSE, sep="\t", row.names = FALSE)

p1 <- qplot(seq(365), state_data$Northeast$var[[2]]) + clean
p2 <- qplot(seq(365), state_data$Midwest$var[[2]]) + clean
p3 <- qplot(seq(365), state_data$South$var[[2]]) + clean
p4 <- qplot(seq(365), state_data$West$var[[2]]) + clean
grid.arrange(p1, p2, p3, p4, ncol = 2)

p1 <- qplot(state_data$Northeast$epi$date, state_data$Northeast$epi$hospitalizedCurrently) + clean
p2 <- qplot(state_data$Midwest$epi$date, state_data$Midwest$epi$hospitalizedCurrently) + clean
p3 <- qplot(state_data$South$epi$date, state_data$South$epi$hospitalizedCurrently) + clean
p4 <- qplot(state_data$West$epi$date, state_data$West$epi$hospitalizedCurrently) + clean
grid.arrange(p1, p2, p3, p4, ncol = 2)

##################################

states <- read.delim("states_49DC.tsv", header = FALSE)$V1
no_states <- length(states)
state2reg <- read.csv("state_lv_data/state2censusReg.csv")
state2reg <- subset(state2reg, select=-c(State))
CR <- c("Northeast", "South", "Midwest", "West")
cdv <- unique(state2reg$Division)
no_cdv <- length(cdv)
state2reg$CDV.Code <- match(state2reg$Division, cdv)

CR_hosp <- readRDS("state_lv_data/censusReg_hospitalization.rds")

mod_hosp <- data.frame(matrix(0, ncol = 4, nrow = length(seq(73, 365+31))))
colnames(mod_hosp) <- CR
load("fit_results/states_49DC.tsv_sun_020420_pso_iterative.rds") # ss_fit and seed_sh
load("covid_hosp_fit/states_49DC.tsv_day_020418_pso_iterative.rds")

sim_SEIH <- function(prms, R0_func, var_ls, census_pop){
  prop_init <- prms[1]
  hpt_rate <- prms[2]
  R0_min <- prms[3]
  R0_range <- prms[4]
  R0_prms <- prms[-1:-4]
  
  xstart <- c(S = census_pop*(1-prop_init), E = 0, I = census_pop*prop_init, H = 0, R = 0) # use actual state population
  times <- seq(73, 365+31) # from national emergency declaration into second half of 2021
  annual_R0 <- do.call(R0_func, c(var_ls, as.list(c(R0_prms, R0_min+R0_range, R0_min))))
  
  # some hard-coded parameters
  paras = list(sigma = 5,
               h = hpt_rate,
               lambda = 5,
               gam = 5,
               k = 10,
               R0 = rep(annual_R0, times=2)) # pre-calculated R0 values
  
  predictions <- as.data.frame(ode(xstart, times, SEIH_R0, paras)) # run SIRS model
  
  
  annual_H1 <- sum(predictions$H[74:365] - 0.9*predictions$H[73:364]) # for k = 10
  annual_H2 <- sum(hpt_rate/5*predictions$I[74:365]) # for lambda = 5
  cat(annual_H1, annual_H2, "\n") # sanity check
  
  return(list(hosp_traj=predictions$H, tot_hosp=annual_H1))
}

for (state_i in seq(no_states)){
  state_eg <- states[state_i]
  census_pop <- state_pop$pop[state_pop$code==state_eg]
  #sunob <- all_state_sun[[state_eg]]/720
  dayob <- all_state_day[[state_eg]]/1440
  # climob <- all_state_hum[[state_eg]] # multiple years
  # climob <- matrix(climob, nrow = length(climob)/365, ncol = 365, byrow = TRUE)
  # climob <- colMeans(climob) # 365 days
  
  #state_prms <- c(ss_fit[[state_eg]][1], seed_sh, ss_fit[[state_eg]][2])
  state_prms <- c(ss_fit[[state_eg]], seed_sh)
  
  cat(state_eg, state_prms, "\n")
  
  state_CR <- state2reg$Region[state2reg$State.Code == state_eg]
  #DST_sim <- sim_SEIH(state_prms, "R0_sun", list(sunob), census_pop)
  DST_sim <- sim_SEIH(state_prms, "R0_day", list(dayob), census_pop)
  
  mod_hosp[[state_CR]] <- mod_hosp[[state_CR]] + DST_sim$hosp_traj
}

mod_hosp$date <- seq(73, 365+31)

p1 <- 
  ggplot() + 
  geom_line(data = mod_hosp, aes(x = date, y=Northeast), color = "#009E73", size=1.2) +
  geom_line(data = subset(CR_hosp, Region == "Northeast"), aes(x = date, y=hospitalizedCurrently), color = "#000000", size=1.2) + clean

p2 <- ggplot() + 
  geom_line(data = mod_hosp, aes(x = date, y=Midwest), color = "#009E73", size=1.2) +
  geom_line(data = subset(CR_hosp, Region == "Midwest"), aes(x = date, y=hospitalizedCurrently), color = "#000000", size=1.2) + clean

p3 <- ggplot() + 
  geom_line(data = mod_hosp, aes(x = date, y=South), color = "#009E73", size=1.2) +
  geom_line(data = subset(CR_hosp, Region == "South"), aes(x = date, y=hospitalizedCurrently), color = "#000000", size=1.2) + clean

p4 <- ggplot() + 
  geom_line(data = mod_hosp, aes(x = date, y=West), color = "#009E73", size=1.2) +
  geom_line(data = subset(CR_hosp, Region == "West"), aes(x = date, y=hospitalizedCurrently), color = "#000000", size=1.2) + clean

grid.arrange(p1, p2, p3, p4, ncol = 2)

########################

state_eg <- "MA"
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
state_df <- state_df[state_df$date <= 396, ]
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
  result <- merge(mod_pred, state_hos, by.x = "time", by.y = "date", all.x=TRUE)
  #cat(-sum(dpois(result$hospitalizedCurrently, result$H, log = TRUE)))
  
  p1 <- ggplot(data = data.frame("time"= 73:(365+72), "R0"= c(R0_vals[73:365], R0_vals[1:72]))) +
    geom_line(aes(x = time, y=R0), color = "black") + clean
  
  # p4 <- ggplot(data = mod_pred) + 
  #   geom_line(aes(x = time, y=E), color = "#CC79A7") +
  #   geom_line(aes(x = time, y=I), color = "#D55E00") +
  #   geom_line(aes(x = time, y=H), color = "#0072B2") + clean
  # 
  # p2 <- ggplot(data = mod_pred) + 
  #   geom_line(aes(x = time, y=S), color = "#000000") +
  #   geom_line(aes(x = time, y=E), color = "#CC79A7") +
  #   geom_line(aes(x = time, y=I), color = "#D55E00") +
  #   geom_line(aes(x = time, y=H), color = "#0072B2") + #clean
  #   geom_line(aes(x = time, y=R), color = "#009E73") + clean
  
  p3 <- ggplot(data = result) + 
    geom_line(aes(x = time, y=hospitalizedCurrently), color = "#D55E00") +
    geom_line(aes(x = time, y=H), color = "#0072B2") + 
    xlim(73, 365+72) + clean
  
  grid.arrange(p1, p3, ncol = 1)
  #return(neg_log_L)
}

# state_code <- read.delim("states_49DC.tsv", header = FALSE)
# no_states <- nrow(state_code)
# prm_names <- c(state_code$V1, "hrate", "R0_min", "R0_range")
# 
# examine <- data.frame(cos=numeric(length(prm_names)),
#                       day=numeric(length(prm_names)),
#                       sd=numeric(length(prm_names)),
#                       sd_s0=numeric(length(prm_names)),
#                       hd=numeric(length(prm_names)),
#                       hsd=numeric(length(prm_names)),
#                       hsd_s0=numeric(length(prm_names)),
#                       row.names = prm_names)
# 
# for (R0_mod in c("cos", "day", "hd")){
#   prms_low = c(rep(1e-5, no_states), 0.01, 0.3, 0.8, param_bounds[[R0_mod]]$low)
#   prms_high = c(rep(0.1, no_states), 0.2, 3, 8, param_bounds[[R0_mod]]$high)
#   result <- readRDS(paste0("fit_results/states_49DC.tsv_", R0_mod, "_011817_pso.rds"))
#   opt_prms <- result$par
#   orig_prms <- prms_low + opt_prms*(prms_high-prms_low)
#   
#   examine[[R0_mod]] <- orig_prms[1:(no_states+3)]
#   cat(paste0(R0_mod, ":")) ; cat(orig_prms[-1:-(no_states+3)]); cat("\n")
# }
# 
# for (R0_mod in c("sd", "hsd")){
#   s0_idx <- ifelse(R0_mod == "sd", 2, 3)
#   prms_low <- c(rep(1e-5, no_states), rep(param_bounds[[R0_mod]]$low[s0_idx], no_states),
#                 0.01, 0.3, 0.8, param_bounds[[R0_mod]]$low[-s0_idx])
#   prms_high = c(rep(0.1, no_states), rep(param_bounds[[R0_mod]]$high[s0_idx], no_states),
#                 0.2, 3, 8, param_bounds[[R0_mod]]$high[-s0_idx])
#   
#   result <- readRDS(paste0("fit_results/states_49DC.tsv_", R0_mod, "_011817_pso.rds"))
#   opt_prms <- result$par
#   orig_prms <- prms_low + opt_prms*(prms_high-prms_low)
#   
#   examine[[R0_mod]] <- c(orig_prms[1:no_states], orig_prms[(2*no_states+1):(2*no_states+3)])
#   examine[[paste0(R0_mod, "_s0")]][1:no_states] <- c(orig_prms[(no_states+1):(2*no_states)])
#   cat(paste0(R0_mod, ":")) ; cat(orig_prms[-1:-(2*no_states+3)]); cat("\n")
# }

load("fit_results/states_49DC.tsv_sun_020420_pso_iterative.rds") # ss_fit and seed_sh
fit_prms <- c(ss_fit[[state_eg]][1], seed_sh, ss_fit[[state_eg]][2])
#fit_prms <- c(0.0001930657, 0.0386439073, 1.185171, 2.106715, -340.1579, 0.2060147761)
run_instance(fit_prms, "R0_sun", list(sunob))

load("covid_hosp_fit/states_49DC.tsv_day_020418_pso_iterative.rds")
fit_prms <- c(ss_fit[[state_eg]], seed_sh)
run_instance(fit_prms, "R0_day", list(dayob))

load("covid_hosp_fit/states_49DC.tsv_sd_020521_pso_iterative.rds")
fit_prms <- c(ss_fit[[state_eg]],
                seed_sh[(no_cdv+1):(no_cdv+4)],
                seed_sh[state2reg$CDV.Code[state2reg$State.Code == state_eg]],
                seed_sh[-1:-(no_cdv+4)])
run_instance(fit_prms, "R0_sd", list(sunob, dayob))
