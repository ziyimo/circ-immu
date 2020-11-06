#!/usr/bin/env Rscript

install.packages("deSolve")
install.packages("reshape2")
install.packages("ggplot2")
install.packages("DEoptim")

library(deSolve)
library(reshape2)
library(ggplot2)
library(DEoptim)


### Functions ###

# SIRS loss function
sirsloss <- function(climvar,sunvar){
  xstart = c(S = 2999999, I = 1, R = 0)
  times = seq(1, 364 * 50, by = 1)
  paras = list(D = 5, 
               L = 40*7, 
               R0max = 2.2, 
               R0min = 1.2, 
               climob = climob, 
               sunob = sunob, 
               climvar = climvar, 
               sunvar = sunvar, 
               birthrate = 0,
               timestart = 0,
               timeend = 0, 
               R0fixed = 1)
  
  predictions <- as.data.frame(ode(xstart, times, SIRS_InterventionGlobal, paras))
  # Take only last year of predictions to speed up
  # Normally we use last 9 years
  years <- 1
  obsyearsmod <- subset(predictions, time > (max(predictions$time) - 364 * years))
  
  # Scale modelled to observed data
  obsyearsmod$I <- obsyearsmod$I / (max(obsyearsmod$I) - min(obsyearsmod$I))
  totinfected <- head(totinfected,364 * years)
  obsdata <- totinfected / (max(totinfected) - min(totinfected))
  scaler <- mean(obsdata) / mean(obsyearsmod$I)
  obsyearsmod$I <- obsyearsmod$I * scaler
  # Baker et al filtering method
  #good <- subset(out5y, I > 0.05 & observed_inf > 0.05)
  
  # calculate the sum of the squared error
  sum((obsyearsmod$I[-1] - obsdata[-1])^2)
  
  # Use binomial distribution
  #-sum(dbinom(as.integer(totinfected), size=as.integer(totpatients), prob=pmodel, log=TRUE))
  
}

# Wrapper function needed for optim()
sirsloss2 <- function(x) {
  sirsloss(climvar = x[1],sunvar = x[2])
}

# SIRS-climate-sunrise model
SIRS_InterventionGlobal <- function(time, state ,theta) {
  ## Parameters:
  climob <- theta[["climob"]]
  subob <- theta[["sunob"]]
  D <- theta[["D"]]
  L <- theta[["L"]]
  R0base <- theta[["R0max"]]
  R0min <- theta[["R0min"]]
  climvar <- theta[["climvar"]]
  sunvar <- theta[["sunvar"]]
  # the below four params are currently unused
  birthrate <- theta[["birthrate"]]
  timestart <- theta[["timestart"]]
  timeend <- theta[["timeend"]]
  R0fixed <- theta[["R0fixed"]]
  
  # get humidity and sunrise at time t
  climobt <- climob[time]
  sunobt <- sunob[time]
  
  ## States
  # susceptible people
  S <- state["S"]
  # Infectious people
  I <- state["I"]
  # Recovered
  R <- state["R"]
  # Total population
  N <- S + I + R
  
  # Ordinary differential equations
  # climvar is the climate dependence parameter
  R0 = exp(climvar*climobt + sunvar*sunobt + log(R0base - R0min)) + R0min
  beta = R0/D

  # birthrate is currently fixed to 0 and thus ignored
  dS <- birthrate*N + (R/L) -beta * S * I/N - S*birthrate
  dI <- beta * S * I/N - (I/D) - I*birthrate
  dR <- (I/D) - (R/L) - R*birthrate
  
  return(list(c(dS, dI, dR)))
}


### External inputs ###

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

fname <- tools::file_path_sans_ext(args[1])

# Load tsv
epi_mw <- read.csv(args[1], sep = "\t", header = FALSE,na.strings=c(""," ","NA"))

# get values as vectors
sunob<- as.numeric(epi_mw$V1)
climob <- as.numeric(epi_mw$V2)
totinfected <- as.numeric(epi_mw$V3)
totpatients <- as.numeric(epi_mw$V4)

# Repeat observations from a single year 50 times
climob = rep(head(climob,364), times = 50)
sunob = rep(head(sunob,364) / 10000, times = 50)
D
### Internal settings ###

# These settings are just for reference as they are currently hardcoded into the functions
# Starting population
pop = 3000000
# start with 1 infected or change here
xstart = c(S = 2999999, I = 1, R = 0)
# Run model for 50 years
times = seq(1, 364 * 50, by = 1)

### Run optimization ###

# initialisation parameters for climvar and sunvar 
starting_param_val <- c(-180,-10)

# Nelder-Mead optimization
outoptim <- paste(fname,".optim.nm.txt", sep = "")
sirs_optim <- optim(starting_param_val, sirsloss2)
sink(outoptim)
sirs_optim
sink()

# BGFS gradient optimization
#ss_optim <- optim(starting_param_val, ss2,method="L-BFGS-B")
#sink("bg1.txt")
#ss_optim
#sink()

### Plot optimization results ###

# Replace climvar and sunvar with fit params
paras = list(D = 5, 
             L = 40*7, 
             R0max = 2.2, 
             R0min = 1.2, 
             climob = climob, 
             sunob = sunob, 
             climvar = sirs_optim$par[1], 
             sunvar = sirs_optim$par[2], 
             birthrate = 0, 
             timestart = 0, 
             timeend = 0, 
             R0fixed= 0)

# Dataframe of S, I and R individuals at each time point
out = as.data.frame(ode(xstart, times, SIRS_InterventionGlobal, paras))

# Plot model output
outplot1 <- paste(fname,".model.infected.pdf", sep = "")
ggplot(out, aes(x= time,y = I)) + 
  geom_line() +
  xlab("Day") +
  ylab("Number of people")
ggsave(filename = outplot1)

# Plot modelled against observed data
out9y <- subset(out, time > (max(out$time) - 364 * 9))
out9y$I <- out9y$I / (max(out9y$I) - min(out9y$I))
infob <- infob / (max(infob) - min(infob))
scaler <- mean(infob) / mean(out9y$I)
out9y$I <- out9y$I * scaler
out9y$observed_inf <- infob
out9y <- melt(out9y, id.vars = c("time"))
names(out9y) <- c("Days", "Group", "Cases")
out9y$Group <- gsub("^S", "Susceptible", out9y$Group)
out9y$Group <- gsub("^I", "Infected", out9y$Group)
out9y$Group <- gsub("^R", "Recovered", out9y$Group)
out9y$Group <- gsub("^observed_inf", "Observed_infected", out9y$Group)
out9ys <- subset(out9y, Group == "Infected" | Group == "Observed_infected")
outplot2 <- paste(fname,".model.fit.pdf", sep = "")
ggplot(out9ys, aes(x= Days,y = Cases, color = Group)) + 
  geom_line() +
  xlab("Day") +
  ylab("Number of people")
ggsave(filename = outplot2)
