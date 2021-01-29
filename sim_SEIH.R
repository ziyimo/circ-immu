library(ggplot2)
library(gridExtra)

clean <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"),
               text = element_text(size=20))

library(deSolve)

source("R0_mods.R")
source("SEIH_mod.R")

states <- read.delim("states_49DC.tsv", header = FALSE)$V1
no_states <- length(states)
## Load all state data
state_pop <- read.delim("state_lv_data/state_pop.tsv")
all_state_sun <- read.csv("state_lv_data/state_daily_sunrise_2019.csv")
all_state_day <- read.csv("state_lv_data/state_daytime_2019.csv")
all_state_hum <- read.csv("state_lv_data/state_humidity_2014_2018.csv")
all_state_hum <- all_state_hum[all_state_hum$X != "2016-02-29", ] # get rid of leap year Feb 29

covid_df <- readRDS("state_lv_data/state_hospitalization.rds")

covid_df$date <- as.numeric(as.Date(covid_df$date, format="%Y-%m-%d") - as.Date("2019-12-31", format="%Y-%m-%d"))

state2reg <- read.csv("state_lv_data/state2censusReg.csv")
state2reg <- subset(state2reg, select=-c(State))
CR <- c("Northeast", "South", "Midwest", "West")

param <- "0.002253395 0.004099292 0.001187666 0.003824723 0.002925207 0.01087801 0.009133795 0.03703352 0.003161158 0.002423917 0.01176747 0.001649684 0.005277838 0.002426359 0.001268024 0.002893894 0.002026028 0.004181413 0.0002852294 0.007146568 0.0139065 0.003157778 0.002429602 0.002498922 0.002114651 0.0006085891 0.0009459835 0.002762578 3.914838e-05 0.04092271 0.002191752 0.04971995 0.001879071 0.0002826201 0.002100939 0.001301621 0.000311214 0.005456829 0.004271333 0.003132766 0.001622551 0.002092081 0.002301657 0.001167583 0.001981641 0.00411364 0.0005887581 0.0002794262 0.000741401 0.01209923 0.4930169 0.5976622 0.4473113 0.4974304 0.4866465 0.5768411 0.5598271 0.6745206 0.5258702 0.5212076 0.1102569 0.4531817 0.5406531 0.539495 0.4670216 0.509153 0.5088651 0.5492728 0.2937096 0.6648065 0.5661439 0.5549252 0.4995165 0.5195733 0.5132375 0.5170742 0.4918857 0.5512954 0.9235613 0.5619274 0.5096342 0.5948038 0.4532768 0.8939135 0.5221441 0.5135969 0.1430238 0.6521453 0.6553496 0.5031417 0.5364245 0.5074912 0.5393127 0.4279327 0.2610154 0.5003523 0.3169015 0.4298926 0.4233851 0.8696278 0.008829422 1.188416 1.635429 -0.00612538 -63.09184 -132.4911 0.4073375"
param <- as.numeric(strsplit(param, " ", fixed=TRUE)[[1]])

# DST mask
DST <- rep(0, times=365)
DST[69:306] <- 60/720 # for DST change in 2019

permDST <- rep(60/720, times=365)
permDST[69:306] <- 0

hosp_dst <- data.frame(matrix(0, ncol = 4, nrow = length(seq(73, 365+180))))
hosp_nodst <- data.frame(matrix(0, ncol = 4, nrow = length(seq(73, 365+180))))
hosp_permdst <- data.frame(matrix(0, ncol = 4, nrow = length(seq(73, 365+180))))
colnames(hosp_dst) <- CR
colnames(hosp_nodst) <- CR
colnames(hosp_permdst) <- CR

annual_cnt <- data.frame(wdst=numeric(length(CR)), wodst=numeric(length(CR)), permdst=numeric(length(CR)), row.names = CR)

sim_SEIH <- function(prms, R0_func, var_ls, census_pop){
  prop_init <- prms[1]
  hpt_rate <- prms[2]
  R0_min <- prms[3]
  R0_range <- prms[4]
  R0_prms <- prms[-1:-4]
  
  xstart <- c(S = census_pop*(1-prop_init), E = 0, I = census_pop*prop_init, H = 0, R = 0) # use actual state population
  times <- seq(73, 365+180) # from national emergency declaration into second half of 2021
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
  sunob <- all_state_sun[[state_eg]]/720
  dayob <- all_state_day[[state_eg]]/1440
  climob <- all_state_hum[[state_eg]] # multiple years
  climob <- matrix(climob, nrow = length(climob)/365, ncol = 365, byrow = TRUE)
  climob <- colMeans(climob) # 365 days
  
  s0_idx <- 3
  state_prms <- c(param[state_i], # state-specific seed_prop
                  param[(2*no_states+1):(2*no_states+2+s0_idx)], # shared: R0_min, R0_range, [R0 prms b4 s_0]
                  param[no_states+state_i], #)) # state_specific s_0
                  param[-1:-(2*no_states+2+s0_idx)])
  cat(state_eg, state_prms, "\n")
  
  state_CR <- state2reg$Region[state2reg$State.Code == state_eg]
  DST_sim <- sim_SEIH(state_prms, "R0_hsd", list(climob, sunob, dayob), census_pop)
  
  hosp_dst[[state_CR]] <- hosp_dst[[state_CR]] + DST_sim$hosp_traj
  annual_cnt[state_CR, "wdst"] <- annual_cnt[state_CR, "wdst"] + DST_sim$tot_hosp
  
  if (state_eg %in% c("AZ", "HI")){
    sansDST_sim <- DST_sim
    permDST_sim <- sim_SEIH(state_prms, "R0_hsd", list(climob, sunob+60/720, dayob), census_pop)
  } else {
    sansDST_sim <- sim_SEIH(state_prms, "R0_hsd", list(climob, sunob-DST, dayob), census_pop)
    permDST_sim <- sim_SEIH(state_prms, "R0_hsd", list(climob, sunob+permDST, dayob), census_pop)
  }
  hosp_nodst[[state_CR]] <- hosp_nodst[[state_CR]] + sansDST_sim$hosp_traj
  annual_cnt[state_CR, "wodst"] <- annual_cnt[state_CR, "wodst"] + sansDST_sim$tot_hosp
  
  hosp_permdst[[state_CR]] <- hosp_permdst[[state_CR]] + permDST_sim$hosp_traj
  annual_cnt[state_CR, "permdst"] <- annual_cnt[state_CR, "permdst"] + permDST_sim$tot_hosp
}

CR_hosp <- readRDS("state_lv_data/censusReg_hospitalization.rds")
#CR_hosp <- CR_hosp[CR_hosp$date <= 365, ]
hosp_dst$date <- seq(73, 365+180)
hosp_nodst$date <- seq(73, 365+180)
hosp_permdst$date <- seq(73, 365+180)

p1 <- 
  ggplot() + 
  geom_line(data = hosp_dst, aes(x = date, y=Northeast), color = "#009E73", size=1.2) +
  geom_line(data = hosp_nodst, aes(x = date, y=Northeast), color = "#0072B2", linetype="dashed", size=1) +
  geom_line(data = hosp_permdst, aes(x = date, y=Northeast), color = "#D55E00", linetype="dashed", size=1) +
  geom_line(data = subset(CR_hosp, Region == "Northeast"), aes(x = date, y=hospitalizedCurrently), color = "#000000", size=1.2) + clean

p2 <- ggplot() + 
  geom_line(data = hosp_dst, aes(x = date, y=Midwest), color = "#009E73", size=1.2) +
  geom_line(data = hosp_nodst, aes(x = date, y=Midwest), color = "#0072B2", linetype="dashed", size=1) +
  geom_line(data = hosp_permdst, aes(x = date, y=Midwest), color = "#D55E00", linetype="dashed", size=1) +
  geom_line(data = subset(CR_hosp, Region == "Midwest"), aes(x = date, y=hospitalizedCurrently), color = "#000000", size=1.2) + clean

p3 <- ggplot() + 
  geom_line(data = hosp_dst, aes(x = date, y=South), color = "#009E73", size=1.2) +
  geom_line(data = hosp_nodst, aes(x = date, y=South), color = "#0072B2", linetype="dashed", size=1) +
  geom_line(data = hosp_permdst, aes(x = date, y=South), color = "#D55E00", linetype="dashed", size=1) +
  geom_line(data = subset(CR_hosp, Region == "South"), aes(x = date, y=hospitalizedCurrently), color = "#000000", size=1.2) + clean

p4 <- ggplot() + 
  geom_line(data = hosp_dst, aes(x = date, y=West), color = "#009E73", size=1.2) +
  geom_line(data = hosp_nodst, aes(x = date, y=West), color = "#0072B2", linetype="dashed", size=1) +
  geom_line(data = hosp_permdst, aes(x = date, y=West), color = "#D55E00", linetype="dashed", size=1) +
  geom_line(data = subset(CR_hosp, Region == "West"), aes(x = date, y=hospitalizedCurrently), color = "#000000", size=1.2) + clean

grid.arrange(p1, p2, p3, p4, ncol = 2)

