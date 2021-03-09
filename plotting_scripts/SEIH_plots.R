library(ggplot2)
library(gridExtra)

clean <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"),
               text = element_text(size=20))

load("covid_hosp_fit/eta0.05_var1_030817_sims.RDa") # load simulation results

CR_hosp <- readRDS("state_lv_data/censusReg_hospitalization.rds")
CR_hosp$date <- as.Date(CR_hosp$date, origin = "2020-01-01")
#CR_hosp <- CR_hosp[CR_hosp$date <= 365, ]
dates <- as.Date(seq(73, 365+31), origin = "2020-01-01")

plot_ls <- list()
for (census_reg in c("Northeast", "Midwest", "South", "West")){
  dst_df <- data.frame(date=dates, opt=dst_trajs[[census_reg]][,1],
                       min=apply(dst_trajs[[census_reg]], 1, FUN=min),
                       max=apply(dst_trajs[[census_reg]], 1, FUN=max))
  nodst_df <- data.frame(date=dates, opt=nodst_trajs[[census_reg]][,1],
                       min=apply(nodst_trajs[[census_reg]], 1, FUN=min),
                       max=apply(nodst_trajs[[census_reg]], 1, FUN=max))
  permdst_df <- data.frame(date=dates, opt=permdst_trajs[[census_reg]][,1],
                       min=apply(permdst_trajs[[census_reg]], 1, FUN=min),
                       max=apply(permdst_trajs[[census_reg]], 1, FUN=max))
  
  plot_ls[[census_reg]] <- 
    ggplot() + 
    geom_line(data = dst_df, aes(x=date, y=opt), color = "#D55E00", size=1.2) +
    geom_ribbon(data = dst_df, aes(x=date, ymin=min, ymax=max), fill = "#D55E00", alpha=0.5) +
    geom_line(data = nodst_df, aes(x = date, y=opt), color = "#7570B3", linetype="dashed", size=1) +
    geom_ribbon(data = nodst_df, aes(x=date, ymin=min, ymax=max), fill = "#7570B3", alpha=0.5) +
    geom_line(data = permdst_df, aes(x = date, y=opt), color = "#800020", linetype="dashed", size=1) +
    geom_ribbon(data = permdst_df, aes(x=date, ymin=min, ymax=max), fill = "#800020", alpha=0.55) +
    geom_line(data = subset(CR_hosp, Region == census_reg), aes(x = date, y=hospitalizedCurrently), color = "#000000", size=1.2) + clean
    #scale_x_date(labels = "%b") + clean
}

do.call(grid.arrange, c(plot_ls, list(ncol=2)))
#grid.arrange(p1, p2, p3, p4, ncol = 2)

# annual_cnt$wodst_reduc <- (annual_cnt$wodst-annual_cnt$wdst)/annual_cnt$wdst
# annual_cnt$permdst_reduc <- (annual_cnt$permdst-annual_cnt$wdst)/annual_cnt$wdst
no_states <- 50
permdst_plt <- data.frame(regime=rep("permdst", no_states+1))
permdst_plt$opt <- permdst_rc$X1
permdst_plt$low <- apply(permdst_rc, 1, FUN=min)
permdst_plt$high <- apply(permdst_rc, 1, FUN=max)
permdst_plt$juri <- row.names(permdst_rc)

nodst_plt <- data.frame(regime=rep("nodst", no_states+1))
nodst_plt$opt <- nodst_rc$X1
nodst_plt$low <- apply(nodst_rc, 1, FUN=min)
nodst_plt$high <- apply(nodst_rc, 1, FUN=max)
nodst_plt$juri <- row.names(nodst_rc)

plt_df <- rbind(permdst_plt, nodst_plt)

ordering <- c(intersect(order(nodst_plt$opt), which(nodst_plt$opt < permdst_plt$opt)),
              intersect(order(-permdst_plt$opt), which(nodst_plt$opt>= permdst_plt$opt)))
# ordering <- order(nodst_plt$opt-permdst_plt$opt)
plt_df$juri <- factor(plt_df$juri, levels = permdst_plt$juri[ordering], order = TRUE)

ggplot(plt_df, aes(x=juri, y=opt, fill=regime)) +
  geom_bar(stat="identity", color="black", position=position_dodge(), width = 0.8)+
  geom_errorbar(aes(ymin=low, ymax=high), width=.2,
                position=position_dodge(.8)) +
  scale_fill_manual(values=c('#7570B3','#800020')) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(-0.4, 0.8, 0.2)) +
  theme(axis.text.x = element_text(size=9), axis.text.y = element_text(size=15), axis.line.x = element_blank()) + clean


## state-level facet plots ##
covid_df <- readRDS("state_lv_data/state_hospitalization.rds")
covid_df$date <- as.numeric(as.Date(covid_df$date, format="%Y-%m-%d") - as.Date("2019-12-31", format="%Y-%m-%d"))

dat_mod_df <- data.frame(date=numeric(), obs=numeric(), pred=numeric(), state=character())

for (state_i in seq(no_states)){
  state_eg <- states[state_i]
  census_pop <- state_pop$pop[state_pop$code==state_eg]
  sunob <- all_state_sun[[state_eg]]/720
  dayob <- all_state_day[[state_eg]]/1440
  
  state_df <- subset(covid_df, state == state_eg)
  state_df <- state_df[order(state_df$date),]
  state_df <- state_df[!is.na(state_df$hospitalizedCurrently), ]
  state_df <- state_df[state_df$date <= 396, ]
  state_hos <- subset(state_df, select = c("date", "hospitalizedCurrently"))
  colnames(state_hos) <- c("date", "obs")
  
  s0_idx <- 2
  state_prms <- c(ss_fit[[state_eg]],
                  seed_sh[(no_cdv+1):(no_cdv+2+s0_idx-1)],
                  seed_sh[state2reg$CDV.Code[state2reg$State.Code == state_eg]],
                  seed_sh[-1:-(no_cdv+2+s0_idx-1)])
  
  cat(state_eg, state_prms, "\n")
  
  DST_sim <- sim_SEIH(state_prms, "R0_sd", list(sunob, dayob), census_pop)
  
  result <- merge(data.frame(date=seq(73, 365+31), pred=DST_sim$hosp_traj), state_hos, all = TRUE)
  result$state <- state_eg
  
  dat_mod_df <- rbind(dat_mod_df, result)
}

dat_mod <- ggplot(dat_mod_df) + 
  geom_line(aes(x=date, y=obs), color = "#000000", size=0.8) +
  geom_line(aes(x=date, y=pred), color = "#009E73", size=0.8) + clean

dat_mod + facet_wrap( ~ state, ncol=5, scales="free_y")

ggplot() + geom_point(data = data.frame(date=seq(73, 396), sun=c(sunob[73:365], sunob[1:31])), aes(x=date, y=sun)) + 
  geom_rect(data=data.frame(x1=c(73), x2=c(396), y1=c(0.55), y2=c(0.6)), mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="red", alpha=0.2)+
  ylim(0.4, 0.7) + clean
ggplot() + 
  geom_line(data = data.frame(date=seq(73, 365+31), hosp=DST_sim$hosp_traj[1:(396-72)]), aes(x = date, y=hosp), color = "#009E73", size=1.2) +
  #geom_line(data = data.frame(date=seq(73, 365+365), hosp=hd_sim$hosp_traj), aes(x = date, y=hosp), color = "#56B4E9", size=1.2) +
  #geom_line(data = data.frame(date=seq(73, 365+365), hosp=sansDST_sim$hosp_traj), aes(x = date, y=hosp), color = "#0072B2", linetype="dashed", size=1) +
  #geom_line(data = data.frame(date=seq(73, 365+365), hosp=permDST_sim$hosp_traj), aes(x = date, y=hosp), color = "#D55E00", linetype="dashed", size=1) +
  geom_line(data = state_hos, aes(x = date, y=hospitalizedCurrently), color = "#000000", size=1.2) + clean

