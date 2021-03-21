library(ggplot2)
library(gridExtra)

clean <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"),
               text = element_text(size=20))

load("covid_hosp_fit/eta0.05_bootsims_032114.RDa") # load simulation results

CR_hosp <- readRDS("state_lv_data/censusReg_hospitalization.rds")
CR_hosp$date <- as.Date(CR_hosp$date, origin = "2020-01-01")
#CR_hosp <- CR_hosp[CR_hosp$date <= 365, ]
dates <- as.Date(seq(73, 365+31), origin = "2020-01-01")
state2reg <- read.csv("state_lv_data/state2censusReg.csv")
state2reg <- state2reg[state2reg$State.Code %in% names(dst_trajs), ]

#err_quantile <- 0.95
plot_ls <- list()
for (census_reg in c("Northeast", "Midwest", "South", "West")){
  DSTtraj <- Reduce("+", dst_trajs[state2reg$State.Code[state2reg$Region == census_reg]])
  std_buf <- apply(DSTtraj, 1, FUN=sd)
  dst_df <- data.frame(date=dates, opt=DSTtraj[,1], min=DSTtraj[,1]-std_buf, max=DSTtraj[,1]+std_buf)
                       # min=apply(dst_trajs[[census_reg]], 1, FUN=quantile, probs = 1-err_quantile),
                       # max=apply(dst_trajs[[census_reg]], 1, FUN=quantile, probs = err_quantile))
  
  noDSTtraj <- Reduce("+", nodst_trajs[state2reg$State.Code[state2reg$Region == census_reg]])
  std_buf <- apply(noDSTtraj, 1, FUN=sd)
  nodst_df <- data.frame(date=dates, opt=noDSTtraj[,1], min=noDSTtraj[,1]-std_buf, max=noDSTtraj[,1]+std_buf)
  
  permDSTtraj <- Reduce("+", permdst_trajs[state2reg$State.Code[state2reg$Region == census_reg]])
  std_buf <- apply(permDSTtraj, 1, FUN=sd)
  permdst_df <- data.frame(date=dates, opt=permDSTtraj[,1], min=permDSTtraj[,1]-std_buf, max=permDSTtraj[,1]+std_buf)
  
  plot_ls[[census_reg]] <- 
    ggplot() + 
    geom_line(data = dst_df, aes(x=date, y=opt), color = "#D55E00", size=1.2) +
    geom_ribbon(data = dst_df, aes(x=date, ymin=min, ymax=max), fill = "#D55E00", alpha=0.3) +
    geom_line(data = nodst_df, aes(x = date, y=opt), color = "#7570B3", linetype="dashed", size=1) +
    geom_ribbon(data = nodst_df, aes(x=date, ymin=min, ymax=max), fill = "#7570B3", alpha=0.2) +
    geom_line(data = permdst_df, aes(x = date, y=opt), color = "#800020", linetype="dashed", size=1) +
    geom_ribbon(data = permdst_df, aes(x=date, ymin=min, ymax=max), fill = "#800020", alpha=0.2) +
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
permdst_plt$stdev <- apply(permdst_rc, 1, FUN=sd)
#permdst_plt$low <- apply(permdst_rc, 1, FUN=quantile, probs = 1-err_quantile)
permdst_plt$low <- permdst_plt$opt - permdst_plt$stdev #apply(permdst_rc, 1, FUN=min)
#permdst_plt$high <- apply(permdst_rc, 1, FUN=quantile, probs = err_quantile)
permdst_plt$high <- permdst_plt$opt + permdst_plt$stdev #apply(permdst_rc, 1, FUN=max)
permdst_plt$juri <- row.names(permdst_rc)

nodst_plt <- data.frame(regime=rep("nodst", no_states+1))
nodst_plt$opt <- nodst_rc$X1
nodst_plt$stdev <- apply(nodst_rc, 1, FUN=sd)
nodst_plt$low <- nodst_plt$opt - nodst_plt$stdev
nodst_plt$high <- nodst_plt$opt + nodst_plt$stdev
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

########### state-level facet plots ###########
states <- read.delim("states_49DC.tsv", header = FALSE)$V1
no_states <- length(states)

## Load COVID hospitalization data
covid_df <- readRDS("state_lv_data/state_hospitalization.rds")
covid_df$date <- as.numeric(as.Date(covid_df$date, format="%Y-%m-%d") - as.Date("2019-12-31", format="%Y-%m-%d"))

dat_mod_df <- NULL
for (state_i in seq(no_states)){
  state_eg <- states[state_i]
  
  state_df <- subset(covid_df, state == state_eg)
  state_df <- state_df[order(state_df$date),]
  state_df <- state_df[!is.na(state_df$hospitalizedCurrently), ]
  state_df <- state_df[state_df$date <= 396, ]
  state_hos <- subset(state_df, select = c("date", "hospitalizedCurrently"))
  colnames(state_hos) <- c("date", "obs")
  
  DST_sd <- apply(dst_trajs[[state_eg]], 1, FUN=sd)
  sansDST_sd <- apply(nodst_trajs[[state_eg]], 1, FUN=sd)
  permDST_sd <- apply(permdst_trajs[[state_eg]], 1, FUN=sd)
  
  result <- merge(data.frame(date=seq(73, 365+31),
                             pred=dst_trajs[[state_eg]][,1], predMin=pmax(dst_trajs[[state_eg]][,1]-DST_sd, 0), predMax=dst_trajs[[state_eg]][,1]+DST_sd,
                             sansDST_pred=nodst_trajs[[state_eg]][,1], sansDST_predMin=pmax(nodst_trajs[[state_eg]][,1]-sansDST_sd, 0), sansDST_predMax=nodst_trajs[[state_eg]][,1]+sansDST_sd,
                             permDST_pred=permdst_trajs[[state_eg]][,1], permDST_predMin=pmax(permdst_trajs[[state_eg]][,1]-permDST_sd, 0), permDST_predMax=permdst_trajs[[state_eg]][,1]+permDST_sd),
                  state_hos, all = TRUE)
  result$state <- state_eg
  
  dat_mod_df <- rbind(dat_mod_df, result)
}

dat_mod_df$date <- as.Date(dat_mod_df$date, origin = "2020-01-01")

dat_mod <- ggplot(dat_mod_df) + 
  geom_line(aes(x=date, y=obs), color = "#000000", size=0.8) +
  geom_line(aes(x=date, y=sansDST_pred), color = "#7570B3", size=0.5, linetype="dashed") +
  geom_ribbon(aes(x=date, ymin=sansDST_predMin, ymax=sansDST_predMax), fill = "#7570B3", alpha=0.2) +
  geom_line(aes(x=date, y=permDST_pred), color = "#800020", size=0.5, linetype="dashed") +
  geom_ribbon(aes(x=date, ymin=permDST_predMin, ymax=permDST_predMax), fill = "#800020", alpha=0.2) +
  geom_line(aes(x=date, y=pred), color = "#D55E00", size=0.8) +
  geom_ribbon(aes(x=date, ymin=predMin, ymax=predMax), fill = "#D55E00", alpha=0.3) + clean

dat_mod + facet_wrap( ~ state, ncol=9, scales="free") + theme(
  strip.text = element_text(margin = margin(0, 0, 0, 0), size=10, face="bold"),
  strip.background = element_rect(color="white", fill="#FFFFFF", size=1.5),
  axis.text.y = element_text(size=8),
  axis.text.x = element_text(size=8))

ggplot() + geom_point(data = data.frame(date=seq(73, 396), sun=c(sunob[73:365], sunob[1:31])), aes(x=date, y=sun)) + 
  geom_rect(data=data.frame(x1=c(73), x2=c(396), y1=c(0.55), y2=c(0.6)), mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="red", alpha=0.2)+
  ylim(0.4, 0.7) + clean
ggplot() + 
  geom_line(data = data.frame(date=seq(73, 365+31), hosp=DST_sim$hosp_traj[1:(396-72)]), aes(x = date, y=hosp), color = "#009E73", size=1.2) +
  #geom_line(data = data.frame(date=seq(73, 365+365), hosp=hd_sim$hosp_traj), aes(x = date, y=hosp), color = "#56B4E9", size=1.2) +
  #geom_line(data = data.frame(date=seq(73, 365+365), hosp=sansDST_sim$hosp_traj), aes(x = date, y=hosp), color = "#0072B2", linetype="dashed", size=1) +
  #geom_line(data = data.frame(date=seq(73, 365+365), hosp=permDST_sim$hosp_traj), aes(x = date, y=hosp), color = "#D55E00", linetype="dashed", size=1) +
  geom_line(data = state_hos, aes(x = date, y=hospitalizedCurrently), color = "#000000", size=1.2) + clean

## Plot year end Recovered ##

# sero_prvl <- read.csv("state_lv_data/Nationwide_Commercial_Laboratory_Seroprevalence_Survey.csv")
# y_end_sero <- subset(sero_prvl, Round == 12, select = c("Site", "Rate......Cumulative.Prevalence.", "Lower.CI..Cumulative.Prevalence.", "Upper.CI..Cumulative.Prevalence."))
# colnames(y_end_sero) <- c("State", "Rate", "Lower_CI", "Upper_CI")
# y_end_sero$Rate[y_end_sero$Rate == 777] <- NA 
# write.table(y_end_sero, file="state_lv_data/Sero_Jan2021.csv", quote = FALSE, sep = "\t", row.names = FALSE)
y_end_sero <- read.delim("state_lv_data/Sero_Jan2021.csv")
y_end_sero$Rate <- y_end_sero$Rate/100
y_end_sero$Lower_CI <- y_end_sero$Lower_CI/100
y_end_sero$Upper_CI <- y_end_sero$Upper_CI/100
ordering <- order(y_end_sero$Rate)
y_end_sero$Regime <- "sero"

yer_std <- apply(year_end_R, 1, sd)
y_end_Recv <- data.frame(State=row.names(year_end_R), Rate=year_end_R$X1,
                         Lower_CI=year_end_R$X1-yer_std, Upper_CI=year_end_R$X1+yer_std, Regime="recv")
sero_recv_df <- rbind(y_end_sero, y_end_Recv)
sero_recv_df$State <- factor(sero_recv_df$State, levels = y_end_sero$State[ordering], order = TRUE)

sero_recv_df <- subset(sero_recv_df, State %in% intersect(y_end_sero$State, row.names(year_end_R)))

ggplot(sero_recv_df, aes(x=State, y=Rate, fill=Regime)) +
  geom_bar(stat="identity", color="black", position=position_dodge(), width = 0.8)+
  geom_errorbar(aes(ymin=Lower_CI, ymax=Upper_CI), width=.2,
                position=position_dodge(.8)) +
  scale_fill_manual(values=c('#ff0000','#ffc100')) +
  #scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(-0.4, 0.8, 0.2)) +
  theme(axis.text.x = element_text(size=9), axis.text.y = element_text(size=15), axis.line.x = element_blank()) + clean

