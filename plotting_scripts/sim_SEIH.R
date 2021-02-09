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

state2reg <- read.csv("state_lv_data/state2censusReg.csv")
state2reg <- subset(state2reg, select=-c(State))
CR <- c("Northeast", "South", "Midwest", "West")
cdv <- unique(state2reg$Division)
no_cdv <- length(cdv)
state2reg$CDV.Code <- match(state2reg$Division, cdv)

#param <- as.numeric(strsplit(param, " ", fixed=TRUE)[[1]])
load("covid_hosp_fit/states_49DC.tsv_sd_020521_pso_iterative.rds") # ss_fit and seed_sh

# DST mask
DST <- rep(0, times=365)
DST[69:306] <- 60/720 # for DST change in 2019

permDST <- rep(60/720, times=365)
permDST[69:306] <- 0


hosp_dst <- data.frame(matrix(0, ncol = 4, nrow = length(seq(73, 365+31))))
hosp_nodst <- data.frame(matrix(0, ncol = 4, nrow = length(seq(73, 365+31))))
hosp_permdst <- data.frame(matrix(0, ncol = 4, nrow = length(seq(73, 365+31))))
colnames(hosp_dst) <- CR
colnames(hosp_nodst) <- CR
colnames(hosp_permdst) <- CR

annual_cnt <- data.frame(wdst=numeric(length(CR)), wodst=numeric(length(CR)), permdst=numeric(length(CR)), row.names = CR)
#annual_cnt <- data.frame(wdst=numeric(no_states), wodst=numeric(no_states), permdst=numeric(no_states), row.names = states)

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
  sunob <- all_state_sun[[state_eg]]/720
  dayob <- all_state_day[[state_eg]]/1440
  # climob <- all_state_hum[[state_eg]] # multiple years
  # climob <- matrix(climob, nrow = length(climob)/365, ncol = 365, byrow = TRUE)
  # climob <- colMeans(climob) # 365 days
  
  s0_idx <- 2
  state_prms <- c(ss_fit[[state_eg]],
                  seed_sh[(no_cdv+1):(no_cdv+3+s0_idx-1)],
                  seed_sh[state2reg$CDV.Code[state2reg$State.Code == state_eg]],
                  seed_sh[-1:-(no_cdv+3+s0_idx-1)])
  
  cat(state_eg, state_prms, "\n")
  
  state_CR <- state2reg$Region[state2reg$State.Code == state_eg]
  DST_sim <- sim_SEIH(state_prms, "R0_sd", list(sunob, dayob), census_pop)
  
  hosp_dst[[state_CR]] <- hosp_dst[[state_CR]] + DST_sim$hosp_traj
  annual_cnt[state_CR, "wdst"] <- annual_cnt[state_CR, "wdst"] + DST_sim$tot_hosp
  #annual_cnt[state_eg, "wdst"] <- DST_sim$tot_hosp
  
  if (state_eg %in% c("AZ", "HI")){
    sansDST_sim <- DST_sim
    permDST_sim <- sim_SEIH(state_prms, "R0_sd", list(sunob+60/720, dayob), census_pop)
  } else {
    sansDST_sim <- sim_SEIH(state_prms, "R0_sd", list(sunob-DST, dayob), census_pop)
    permDST_sim <- sim_SEIH(state_prms, "R0_sd", list(sunob+permDST, dayob), census_pop)
  }
  hosp_nodst[[state_CR]] <- hosp_nodst[[state_CR]] + sansDST_sim$hosp_traj
  #annual_cnt[state_eg, "wodst"] <-sansDST_sim$tot_hosp
  annual_cnt[state_CR, "wodst"] <- annual_cnt[state_CR, "wodst"] + sansDST_sim$tot_hosp
  
  hosp_permdst[[state_CR]] <- hosp_permdst[[state_CR]] + permDST_sim$hosp_traj
  #annual_cnt[state_eg, "permdst"] <- permDST_sim$tot_hosp
  annual_cnt[state_CR, "permdst"] <- annual_cnt[state_CR, "permdst"] + permDST_sim$tot_hosp
}

CR_hosp <- readRDS("state_lv_data/censusReg_hospitalization.rds")
#CR_hosp <- CR_hosp[CR_hosp$date <= 365, ]
hosp_dst$date <- seq(73, 365+31)
hosp_nodst$date <- seq(73, 365+31)
hosp_permdst$date <- seq(73, 365+31)

hosp_dst <- hosp_dst[1:(396-72),]
hosp_nodst <- hosp_nodst[1:(396-72),]
hosp_permdst <- hosp_permdst[1:(396-72),]

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

if (FALSE){
  covid_df <- readRDS("state_lv_data/state_hospitalization.rds")
  covid_df$date <- as.numeric(as.Date(covid_df$date, format="%Y-%m-%d") - as.Date("2019-12-31", format="%Y-%m-%d"))
  
  
  state_eg <- "MA"
  state_i <- match(state_eg, states)
  census_pop <- state_pop$pop[state_pop$code==state_eg]
  sunob <- all_state_sun[[state_eg]]/720
  dayob <- all_state_day[[state_eg]]/1440
  climob <- all_state_hum[[state_eg]] # multiple years
  climob <- matrix(climob, nrow = length(climob)/365, ncol = 365, byrow = TRUE)
  climob <- colMeans(climob) # 365 days
  
  state_df <- subset(covid_df, state == state_eg)
  state_df <- state_df[order(state_df$date),]
  state_df <- state_df[!is.na(state_df$hospitalizedCurrently), ]
  state_df <- state_df[state_df$date <= 365, ]
  state_hos <- subset(state_df, select = c("date", "hospitalizedCurrently"))
  
  # s0_idx <- 3
  # state_prms <- c(param[state_i], # state-specific seed_prop
  #                 param[(2*no_states+1):(2*no_states+2+s0_idx)], # shared: R0_min, R0_range, [R0 prms b4 s_0]
  #                 param[no_states+state_i], #)) # state_specific s_0
  #                 param[-1:-(2*no_states+2+s0_idx)])
  state_prms <- c(param[state_i], param[-1:-no_states])
  cat(state_eg, state_prms, "\n")
  
  # state_hd_prms <- c(hd_param[state_i], hd_param[-1:-no_states])
  # cat(state_eg, state_hd_prms, "\n")
  
  DST_sim <- sim_SEIH(state_prms, "R0_hsd", list(climob, sunob, dayob), census_pop)
  
  #hd_sim <- sim_SEIH(state_hd_prms, "R0_hd", list(climob, dayob), census_pop)
  
  if (state_eg %in% c("AZ", "HI")){
    sansDST_sim <- DST_sim
    permDST_sim <- sim_SEIH(state_prms, "R0_hsd", list(climob, sunob+60/720, dayob), census_pop)
  } else {
    sansDST_sim <- sim_SEIH(state_prms, "R0_hsd", list(climob, sunob-DST, dayob), census_pop)
    permDST_sim <- sim_SEIH(state_prms, "R0_hsd", list(climob, sunob+permDST, dayob), census_pop)
  }

  ggplot() + geom_point(data = data.frame(date=seq(73, 396), sun=c(sunob[73:365], sunob[1:31])), aes(x=date, y=sun)) + 
    geom_rect(data=data.frame(x1=c(73), x2=c(396), y1=c(0.55), y2=c(0.6)), mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="red", alpha=0.2)+
    ylim(0.4, 0.7) + clean
  ggplot() + 
    geom_line(data = data.frame(date=seq(73, 365+31), hosp=DST_sim$hosp_traj[1:(396-72)]), aes(x = date, y=hosp), color = "#009E73", size=1.2) +
    #geom_line(data = data.frame(date=seq(73, 365+365), hosp=hd_sim$hosp_traj), aes(x = date, y=hosp), color = "#56B4E9", size=1.2) +
    #geom_line(data = data.frame(date=seq(73, 365+365), hosp=sansDST_sim$hosp_traj), aes(x = date, y=hosp), color = "#0072B2", linetype="dashed", size=1) +
    #geom_line(data = data.frame(date=seq(73, 365+365), hosp=permDST_sim$hosp_traj), aes(x = date, y=hosp), color = "#D55E00", linetype="dashed", size=1) +
    geom_line(data = state_hos, aes(x = date, y=hospitalizedCurrently), color = "#000000", size=1.2) + clean
  
}
