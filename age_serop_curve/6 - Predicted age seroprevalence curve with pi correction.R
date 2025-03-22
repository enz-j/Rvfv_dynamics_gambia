# This version of the model simulates age-stratified seroprevalence dynamics in the cattle subpopulations, 
# incorporating a weekly-varying force of infection across distinct eco-regions. The model extends the 
# classical SIR framework by explicitly accounting for the effect of immunity on the force of infection.
# To reconcile discrepancies between predicted and observed seroprevalence levels, the model incorporates 
# seropositivity decay through a parameter (pi) that quantifies the decay of seropositivity in
# immune cattle over time.

rm(list = ls())

# Load libraries
library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)

# Load lambda vector
load("lambda_vec.RData")

# Constants and initial values
gamma = 0.0758493 * 7  # days
delta = 0.1161880 * 7  # days
pi = 0.0018155         # rate of removal of immune cattle

# Initial population setup
initial_N <- c(S.T = 35000, I.T = 1,
               S.M = 170000, I.M = 1,
               S.L = 85000, I.L = 1)

# Define SIR model function
SIR_model <- function(time, N, lambda_t) {
  with(as.list(c(N, lambda_t)), {
    dS.T <- -lambda_hn * S.T
    dI.T <- (delta / (delta + gamma)) * lambda_hn * S.T
    dS.M <- -lambda_m * S.M
    dI.M <- (delta / (delta + gamma)) * lambda_m * S.M
    dS.L <- -lambda_l * S.L
    dI.L <- (delta / (delta + gamma)) * lambda_l * S.L
    
    return(list(c(dS.T, dI.T, dS.M, dI.M, dS.L, dI.L)))
  })
}

# Extract population values at a specific time
get_N <- function(out_struc) {
  myN <- matrix(0, nrow = 6)
  myN <- c(out_struc[round(length(out_struc) / 7, 0), 2],
           out_struc[round(length(out_struc) / 7, 0), 3],
           out_struc[round(length(out_struc) / 7, 0), 4],
           out_struc[round(length(out_struc) / 7, 0), 5],
           out_struc[round(length(out_struc) / 7, 0), 6],
           out_struc[round(length(out_struc) / 7, 0), 7])
  return(myN)
}

# Set initial lambda and time parameters
lambda_hn <- lambda_vec[1, 1]
lambda_m <- lambda_vec[1, 2]
lambda_l <- lambda_vec[1, 3]
lambda_t <- c(lambda_hn, lambda_m, lambda_l)
param <- c(gamma, lambda_t)

# Initial time and ODE call
time <- seq(0, 1, by = 1)
N <- initial_N
out.serop <- ode(N, time, SIR_model, param)
all_out.serop <- out.serop

# Loop over time grid to run the ODE model for each week
time_grid <- seq(1, 555, by = 1)
for (weeks in 1:length(time_grid)) {
  lambda_hn <- lambda_vec[weeks, 1]
  lambda_m <- lambda_vec[weeks, 2]
  lambda_l <- lambda_vec[weeks, 3]
  
  lambda_t <- c(lambda_hn, lambda_m, lambda_l)
  param <- c(lambda_t, gamma)
  
  N <- get_N(out.serop)  # Initialize to the end of last integration period
  out.serop <- ode(y = N, times = c(weeks, weeks + 1), func = SIR_model, parms = param)
  all_out.serop <- rbind(all_out.serop, out.serop)  # Add the new output to the cumulative output
}

# Clean and prepare data for plotting
df <- as.data.frame(all_out.serop)
df <- distinct(df)
df$time <- 1:557

# Set plot theme
mytheme4 <- theme_bw() +
  theme(text = element_text(colour = "black")) +
  theme(panel.grid = element_line(colour = "#f6f6f6")) +
  theme(panel.background = element_rect(fill = "#FFFFFF"))
theme_set(mytheme4)

# Define data for comparison with observations
df.plot_ct.ag <- data.frame(
  x.M = c(50, 100, 154, 206, 258, 310, 362),
  Percentage.M = c(.207, .229, .211, .323, .334, .323, .457),
  lb.M = c(.126, .147, .141, .225, .243, .231, .344),
  ub.M = c(.316, .335, .302, .437, .439, .429, .573),
  x.L = c(52, 104, 156, 208, 260, 312, 364),
  Percentage.L = c(.289, .130, .361, .446, .239, .441, .492),
  lb.L = c(.094, .025, .203, .201, .075, .247, .257),
  ub.L = c(.597, .39, .551, .716, .525, .653, .728),
  x.T = c(54, 106, 158, 210, 262, 314, 366),
  Percentage.T = c(.360, .373, .378, .408, .58, .518, .438),
  lb.T = c(.257, .288, .297, .307, .481, .427, .344),
  ub.T = c(.477, .466, .466, .516, .667, .607, .536)
)

# Plot proportions
A <- ggplot() +
  geom_line(data = df, aes(x = time, y = I.T / (S.T + I.T), colour = "P.T")) +
  geom_line(data = df, aes(x = time, y = I.M / (S.M + I.M), colour = "P.M")) +
  geom_line(data = df, aes(x = time, y = I.L / (S.L + I.L), colour = "P.L")) +
  ylab(label = "Proportion of immune \nseropositive cattle") +
  xlab(label = "Age (weeks)") +
  theme(legend.justification = c(1, 0), legend.position = c(0.5, 0.8), legend.direction = "horizontal") +
  theme(legend.title = element_text(size = 11, face = "bold"),
        legend.background = element_rect(fill = '#FFFFFF', size = 0.4, linetype = "solid"),
        legend.text = element_text(size = 11),
        legend.key = element_rect(colour = "#FFFFFF", fill = '#C2C2C2', size = 0.35, linetype = "solid")) +
  theme(
    axis.title.x = element_text(size = 17),
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  ) +
  scale_colour_manual("Compartments",
                      breaks = c("P.M", "P.L", "P.T"),
                      values = c("darkgreen", "blue", "red")) +
  scale_y_continuous(limit = c(0, 1), breaks = c(0.2, 0.4, 0.6, 0.8, 1))

# Add observed data with error bars
B <- A + 
  geom_point(data = df.plot_ct.ag, aes(x = x.T, y = Percentage.T, colour = "P.T")) +
  geom_errorbar(data = df.plot_ct.ag, aes(x = x.T, ymin = lb.T, ymax = ub.T, colour = "P.T"), width = 3) +
  geom_point(data = df.plot_ct.ag, aes(x = x.M, y = Percentage.M, colour = "P.M")) +
  geom_errorbar(data = df.plot_ct.ag, aes(x = x.M, ymin = lb.M, ymax = ub.M, colour = "P.M"), width = 3) +
  geom_point(data = df.plot_ct.ag, aes(x = x.L, y = Percentage.L, colour = "P.L")) +
  geom_errorbar(data = df.plot_ct.ag, aes(x = x.L, ymin = lb.L, ymax = ub.L, colour = "P.L"), width = 3)

# Calculate the proportion of immune Cattle
dat <- as.Matrix(get_N(all_out.serop))
dat <- as.data.frame(dat)
P.T <- (dat[2, 1]) / (dat[1, 1] + dat[2, 1]) * 100
P.M <- (dat[4, 1]) / (dat[3, 1] + dat[4, 1]) * 100
P.L <- (dat[6, 1]) / (dat[5, 1] + dat[6, 1]) * 100

##########################################################################################
# Reconcilling the age-seroprevaelcne curves
##########################################################################################

# Force of infection model correction 
SIR_model_cr <- function(time, N, lambda_t) {
  with(as.list(c(N, lambda_t)), {
    dS.T <- -lambda_hn * S.T
    dI.T <- (delta / (delta + gamma)) * lambda_hn * S.T - pi * I.T
    dR.T <- pi * I.T
    dS.M <- -lambda_m * S.M
    dI.M <- (delta / (delta + gamma)) * lambda_m * S.M - pi * I.M
    dR.M <- pi * I.M
    dS.L <- -lambda_l * S.L
    dI.L <- (delta / (delta + gamma)) * lambda_l * S.L - pi * I.L
    dR.L <- pi * I.L
    
    return(list(c(dS.T, dI.T, dR.T, dS.M, dI.M, dR.M, dS.L, dI.L, dR.L)))
  })
}

# Correct the ODE model
get_N_cr <- function(out_struc) {
  myN <- matrix(0, nrow = 9)
  myN <- c(out_struc[round(length(out_struc) / 10, 0), 2],
           out_struc[round(length(out_struc) / 10, 0), 3],
           out_struc[round(length(out_struc) / 10, 0), 4],
           out_struc[round(length(out_struc) / 10, 0), 5],
           out_struc[round(length(out_struc) / 10, 0), 6],
           out_struc[round(length(out_struc) / 10, 0), 7],
           out_struc[round(length(out_struc) / 10, 0), 8],
           out_struc[round(length(out_struc) / 10, 0), 9],
           out_struc[round(length(out_struc) / 10, 0), 10])
  return(myN)
}

# Solve corrected ODE
param_cr <- c(gamma, pi, lambda_t)
out.serop_cr <- ode(N, time, SIR_model_cr, param_cr)
all_out.serop_cr <- out.serop_cr

# Time grid for correction
for (weeks in 1:length(time_grid)) {
  lambda_hn <- lambda_vec[weeks, 1]
  lambda_m <- lambda_vec[weeks, 2]
  lambda_l <- lambda_vec[weeks, 3]
  
  lambda_t <- c(lambda_hn, lambda_m, lambda_l)
  param_cr <- c(lambda_t, pi, gamma)
  
  N <- get_N_cr(out.serop_cr)
  out.serop_cr <- ode(y = N, times = c(weeks, weeks + 1), func = SIR_model_cr, parms = param_cr)
  all_out.serop_cr <- rbind(all_out.serop_cr, out.serop_cr)
}

# Proportions after correction (Plotting)
df_cr <- as.data.frame(all_out.serop_cr)
df_cr <- distinct(df_cr)
df_cr$time <- 1:557

C <- ggplot() +
  geom_line(data = df_cr, aes(x = time, y = I.T / (S.T + I.T), colour = "P.T")) +
  geom_line(data = df_cr, aes(x = time, y = I.M / (S.M + I.M), colour = "P.M")) +
  geom_line(data = df_cr, aes(x = time, y = I.L / (S.L + I.L), colour = "P.L")) +
  ylab(label = "Proportion of immune \nseropositive cattle") +
  xlab(label = "Age (weeks)") +
  theme(legend.justification = c(1, 0), legend.position = c(0.5, 0.8), legend.direction = "horizontal") +
  theme(legend.title = element_text(size = 11, face = "bold"),
        legend.background = element_rect(fill = '#FFFFFF', size = 0.4, linetype = "solid"),
        legend.text = element_text(size = 11),
        legend.key = element_rect(colour = "#FFFFFF", fill = '#C2C2C2', size = 0.35, linetype = "solid")) +
  theme(
    axis.title.x = element_text(size = 17),
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  ) +
  scale_colour_manual("Compartments",
                      breaks = c("P.M", "P.L", "P.T"),
                      values = c("darkgreen", "blue", "red")) +
  scale_y_continuous(limit = c(0, 1), breaks = c(0.2, 0.4, 0.6, 0.8, 1))

# Add observed data with error bars
D <- C + 
  geom_point(data = df.plot_ct.ag, aes(x = x.T, y = Percentage.T, colour = "P.T")) +
  geom_errorbar(data = df.plot_ct.ag, aes(x = x.T, ymin = lb.T, ymax = ub.T, colour = "P.T"), width = 3) + 
  geom_point(data = df.plot_ct.ag, aes(x = x.M, y = Percentage.M, colour = "P.M")) +
  geom_errorbar(data = df.plot_ct.ag, aes(x = x.M, ymin = lb.M, ymax = ub.M, colour = "P.M"), width = 3) + 
  geom_point(data = df.plot_ct.ag, aes(x = x.L, y = Percentage.L, colour = "P.L")) +
  geom_errorbar(data = df.plot_ct.ag, aes(x = x.L, ymin = lb.L, ymax = ub.L, colour = "P.L"), width = 3)

# Display results
ggarrange(B, D, nrow = 2, ncol = 1)
