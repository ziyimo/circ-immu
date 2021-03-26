#!/usr/bin/env Rscript

# sinusoidal R0 model
R0_cos <- function(t, phi, R0base = 2, R0min = 1.2){ # a sinusoidal baseline model
  R0 <- (R0base - R0min)/2*cos(2*pi/364*(t - phi)) + (R0base + R0min)/2
  return(R0)
}

# single var R0 models

R0_hum <- function(h_t, alpha, R0base = 2, R0min = 1.2){ # Baker et al. exponential decay
  R0 <- exp(alpha*h_t + log(R0base - R0min)) + R0min
  return(R0)
}

R0_day <- function(d_t, alpha, d_0, R0base = 2, R0min = 1.2){
  R0 <- exp(alpha*(d_t-d_0)^2 + log(R0base - R0min)) + R0min
  return(R0)
}

# composite R0 models

R0_hd <- function(h_t, d_t, alpha_1, alpha_2, d_0, R0base = 2, R0min = 1.2){ # humidity + daytime
  R0 <- exp(alpha_1*h_t + alpha_2*(d_t-d_0)^2 + log(R0base - R0min)) + R0min
  return(R0)
}

R0_sd <- function(s_t, d_t, alpha_1, s_0, alpha_2, d_0, R0base = 2, R0min = 1.2){
  R0 <- exp(alpha_1*(s_t-s_0)^2 + alpha_2*(d_t-d_0)^2 + log(R0base - R0min)) + R0min
  return(R0)
}

R0_hsd <- function(h_t, s_t, d_t, alpha_1, alpha_2, s_0, alpha_3, d_0, R0base = 2, R0min = 1.2){
  R0 <- exp(alpha_1*h_t + alpha_2*(s_t-s_0)^2 + alpha_3*(d_t-d_0)^2 + log(R0base - R0min)) + R0min
  return(R0)
}