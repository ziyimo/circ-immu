#!/usr/bin/env Rscript

param_bounds <- list() # boundaries of parameters, first element is "c"

# sinusoidal R0 model
R0_cos <- function(t, phi, R0base = 2, R0min = 1.2){ # a sinusoidal baseline model
  R0 <- (R0base - R0min)/2*cos(2*pi/364*(t - phi)) + (R0base + R0min)/2
  return(R0)
}
#param_bounds[["cos"]] <- list(low = c(0, 1), high = c(1, 364))

# single var R0 models

R0_hum <- function(h_t, alpha, R0base = 2, R0min = 1.2){ # Baker et al. exponential decay
  R0 <- exp(alpha*h_t + log(R0base - R0min)) + R0min
  return(R0)
}
#param_bounds[["hum"]] <- list(low = c(0, -300), high = c(1, 0))
# > max(all_state_hum)
# [1] 0.02134766

R0_day <- function(d_t, alpha, d_0, R0base = 2, R0min = 1.2){
  R0 <- exp(alpha*(d_t-d_0)^2 + log(R0base - R0min)) + R0min
  return(R0)
}
#param_bounds[["day"]] <- list(low = c(0, -100, 0), high = c(1, 0, 1))
R0_sun <- R0_day

R0_sunAsym <- function(s_t, alpha_l, alpha_r, s_0, R0base = 2, R0min = 1.2){
  right <- as.numeric((s_t - s_0) > 0)
  left <- as.numeric((s_t - s_0) < 0)
  R0 <- exp(left*alpha_l*(s_t-s_0)^2 + right*alpha_r*(s_t-s_0)^2 + log(R0base - R0min)) + R0min
  return(R0)
}

# composite R0 models

R0_hd <- function(h_t, d_t, alpha_1, alpha_2, d_0, R0base = 2, R0min = 1.2){ # humidity + daytime
  R0 <- exp(alpha_1*h_t + alpha_2*(d_t-d_0)^2 + log(R0base - R0min)) + R0min
  return(R0)
}
#param_bounds[["hd"]] <- list(low = c(0, -300, -100, 0), high = c(1, 0, 0, 1))
R0_hs <- R0_hd

R0_sd <- function(s_t, d_t, alpha_1, s_0, alpha_2, d_0, R0base = 2, R0min = 1.2){
  R0 <- exp(alpha_1*(s_t-s_0)^2 + alpha_2*(d_t-d_0)^2 + log(R0base - R0min)) + R0min
  return(R0)
}
#param_bounds[["sd"]] <- list(low = c(0, -100, 0, -100, 0), high = c(1, 0, 1, 0, 1))

R0_hsd <- function(h_t, s_t, d_t, alpha_1, alpha_2, s_0, alpha_3, d_0, R0base = 2, R0min = 1.2){
  R0 <- exp(alpha_1*h_t + alpha_2*(s_t-s_0)^2 + alpha_3*(d_t-d_0)^2 + log(R0base - R0min)) + R0min
  return(R0)
}
#param_bounds[["hsd"]] <- list(low = c(0, -300, -100, 0, -100, 0), high = c(1, 0, 0, 1, 0, 1))

######## Plot R0 models #######
if(FALSE){
  library(ggplot2)
  library(gridExtra)
  
  plot(seq(364), R0_cos(seq(364), 120))
  
  clean <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                 text = element_text(size=20))
  
  sun_range = seq(0, 1, by=0.01)
  day_range = seq(0, 1, by=0.01)
  cli_range = seq(0, 0.03, by=0.0003)
  
  ggplot() + 
    geom_line(data = data.frame("var"= sun_range, "R0"= R0_sunAsym(sun_range, -50, -1, 0.6)),
              aes(x = var, y=R0), color = "#E69F00") +
    geom_line(data = data.frame("var"= sun_range, "R0"= R0_sunAsym(sun_range, -10, -10, 0.5)), 
              aes(x = var, y=R0), color = "#0072B2") +
    geom_line(data = data.frame("var"= sun_range, "R0"= R0_sunAsym(sun_range, -1, -1, 0.55)), 
              aes(x = var, y=R0), color = "#009E73") + clean
  
  ggplot() + 
    geom_line(data = data.frame("day"= day_range, "R0"= R0_day(day_range, -10, 0.4)), 
              aes(x = day, y=R0), color = "#E69F00") +
    geom_line(data = data.frame("day"= day_range, "R0"= R0_day(day_range, -100, 0.4)), 
              aes(x = day, y=R0), color = "#0072B2") +
    geom_line(data = data.frame("day"= day_range, "R0"= R0_day(day_range, -300, 0.4)), 
              aes(x = day, y=R0), color = "#009E73") + clean
  
  #grid.arrange(p1, p2, ncol = 1)
  
  hum_day_df <- expand.grid(cli_range, day_range)
  colnames(hum_day_df) <- c("hum", "day")
  
  sun_day_df <- expand.grid(sun_range, day_range)
  colnames(sun_day_df) <- c("sun", "day")
  
  hum_day_df$R0 <- R0_hd(hum_day_df$hum, hum_day_df$day, -100, -1)
  p1 <- ggplot(hum_day_df, aes(hum, day, z = R0)) + geom_contour_filled() + clean
  
  hum_day_df$R0 <- R0_hd(hum_day_df$hum, hum_day_df$day, -30, -5)
  p2 <- ggplot(hum_day_df, aes(hum, day, z = R0)) + geom_contour_filled() + clean
  
  sun_day_df$R0 <- R0_sd(sun_day_df$sun, sun_day_df$day, -100, -3)
  p3 <- ggplot(sun_day_df, aes(sun, day, z = R0)) + geom_contour_filled() + clean
  
  sun_day_df$R0 <- R0_sd(sun_day_df$sun, sun_day_df$day, -100, -100)
  p4 <- ggplot(sun_day_df, aes(sun, day, z = R0)) + geom_contour_filled() + clean
  
  grid.arrange(p1, p2, p3, p4, ncol = 2)
  grid.arrange(p2, p3, ncol = 1)
}
