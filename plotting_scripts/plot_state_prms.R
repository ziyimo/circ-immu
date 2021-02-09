library(ggplot2)
library(gridExtra)
library(ggrepel)

clean <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"),
               text = element_text(size=20))

states <- read.delim("states_49DC.tsv", header = FALSE)
no_states <- nrow(states)

param <- "0.00494648 0.004725592 0.002797547 0.006186709 0.00108099 0.02101321 0.009556145 0.08074675 0.005869293 0.004678179 0.03280633 0.0005922044 0.007022075 0.00294134 0.002480091 0.008191887 0.002571842 0.009529793 0.0003011255 0.009080876 0.02397457 0.002736044 0.0006965941 0.004541399 0.002859488 0.000799127 0.00101418 0.003802543 8.937205e-05 0.0280871 0.001989917 0.03000426 0.00312669 0.0003524225 0.000726298 0.001631182 0.0005225833 0.005743082 0.002155352 0.0055846 0.001112112 0.002967331 0.004022597 0.002340543 0.00492091 0.005791592 0.0005920849 0.0006699368 0.001011222 0.0001474013 0.4419011 0.7027351 0.4354307 0.4557757 0.3776401 0.5168677 0.6905534 0.8743672 0.4673612 0.4886965 0.0003147507 0.0207076 0.4622 0.4719123 0.4395842 0.4869069 0.4671157 0.4708707 5.732814e-05 0.4754487 0.5233149 0.4871234 0.4062399 0.436152 0.4508618 0.4509692 0.4453371 0.654378 0.09888624 0.5555276 0.4214297 0.5830996 0.439778 0.9935164 0.4414427 0.4576465 0.02298935 0.4782788 0.7353051 0.4715997 0.4480709 0.4465802 0.4742313 0.4230786 0.3008775 0.4584686 0.3567086 0.4144042 0.3719223 0.244507 0.01258636 1.161118 1.89173 -0.03399191 -54.28325 -255.3139 0.4132449"
param <- as.numeric(strsplit(param, " ", fixed=TRUE)[[1]])

states$Init <- param[1:no_states]
states$s0 <- param[(no_states+1):(2*no_states)]

state2reg <- read.csv("state_lv_data/state2censusReg.csv")
state2reg <- subset(state2reg, select=-c(State))
#CR <- c("Northeast", "South", "Midwest", "West")

state_prms_df <- merge(states, state2reg, by.x = "V1", by.y = "State.Code")

ggplot(state_prms_df, aes(x=s0, y=Init, color=Region)) +
  geom_point(size=1) +
  scale_y_continuous(trans='log10') + 
  geom_text_repel(aes(label = V1), size = 3, max.overlaps = Inf) +
  scale_color_manual(values = c("#0066CC", "#00994C", "#FFD700", "#990000")) + clean

## Plot state fit ##

#################