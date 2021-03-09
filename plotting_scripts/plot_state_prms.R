library(ggplot2)
library(gridExtra)
library(ggrepel)

clean <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"),
               text = element_text(size=20))

states <- read.delim("states_49DC.tsv", header = FALSE)
no_states <- nrow(states)

state2reg <- read.csv("state_lv_data/state2censusReg.csv")
cdv <- unique(state2reg$Division)
no_cdv <- length(cdv)
state2reg$CDV.Code <- match(state2reg$Division, cdv)

#load("covid_hosp_fit/states_49DC.tsv_sd_020521_pso_iterative.rds") # ss_fit and seed_sh

ss_fit <- "0.004764315 0.003575784 0.004366139 0.003453841 0.005309574 0.009086865 0.004650126 0.008474114 0.007286027 0.003067197 0.01913177 0.001913116 0.006524789 0.003710029 0.003180887 0.003748295 0.002184048 0.007165196 0.00879047 0.006058895 0.01030153 0.003093634 0.002138411 0.005700072 0.003272083 0.002170342 0.09998269 0.004087742 0.008266404 0.01023782 0.004852979 0.01080531 0.004613433 0.002579945 0.002564202 0.002647211 0.001061218 0.007379454 0.007011503 0.004616346 0.002567075 0.003453841 0.003843991 0.002647731 0.009170044 0.005342982 0.001913764 0.0015825 0.005353082 0.002613343"
seed_sh <- "0.4549885 0.4712374 0.4852478 0.4779713 0.6246432 0.4864148 0.4750408 0.5055364 0.6053397 0.9998064 0.6811585 -35.16802 -147.2332 0.4378637"
ss_fit <- as.numeric(strsplit(ss_fit, " ", fixed=TRUE)[[1]])
names(ss_fit) <- states$V1
seed_sh <- as.numeric(strsplit(seed_sh, " ", fixed=TRUE)[[1]])

states$Init <- numeric(no_states)
states$s0 <- numeric(no_states)

for (state_i in seq(no_states)){
  state_code <- states$V1[state_i]
  states$Init[state_i] <- ss_fit[[state_code]]
  states$s0[state_i] <- seed_sh[state2reg$CDV.Code[state2reg$State.Code == state_code]]
}

#CR <- c("Northeast", "South", "Midwest", "West")

state_prms_df <- merge(states, state2reg, by.x = "V1", by.y = "State.Code")

ggplot(state_prms_df, aes(x=s0, y=Init, color=Region)) +
  geom_point(size=1) +
  scale_y_continuous(trans='log10') + 
  geom_text_repel(aes(label = V1), size = 3, max.overlaps = Inf) +
  scale_color_manual(values = c("#0066CC", "#00994C", "#FFD700", "#990000")) +
  xlim(0.4, 0.7) + clean
