# Hypothetical Movement Ban Simulation

# Clear the environment
rm(list = ls())

# Load necessary libraries
library(deSolve)      # For solving differential equations
library(ggplot2)      # For data visualization
library(reshape2)     # For data reshaping
library(patchwork)    # For combining multiple plots

### Initial conditions ###
initial_N <- c(S.T = 13165.62, I.T = 1.8447889, R.T = 10093.57,       # Transhumant cattle (T)
               S.M = 83395.99, I.M = 1.537987e-18, R.M = 55182.37,    # Resident cattle in Sahel (M)
               S.L = 67037.90, I.L = 9.404628, R.L = 42077.79)        # Resident cattle in river (L)

# Total population for each eco-region
N.T <- sum(initial_N[1:3])  
N.M <- sum(initial_N[4:6])  
N.L <- sum(initial_N[7:9])  

# Define model parameters (mean posterior values from ABC procedure)
mu <- 0.0016940        # Natural mortality rate
gamma <- 0.0758493     # RVFV-induced mortality rate
delta <- 0.1161880     # Recovery rate
b <- 0.0021754         # Birth rate
sf <- 2.0508335        # Scaling factor

# Calculate basic reproduction number for each eco-region (R0,i) following Anderson and May (1991)
# Using seroprevalence estimates
# Serop_subpop_T = 45%, Serop_subpop_M = 30%, Serop_subpop_L = 36%
# avg_serop_s.wet = ((30% x 117004) + (45%  x 35400))/(117004 + 35400) = 33% 
# avg_serop_s.dry = 0
# avg_serop_r.wet = 36%
# avg_serop_r.dry = ((36%  x 80433) + (45%  x 35400))/(80433 + 35400) = 39%

#R0.s.wet = 1/0.68 = 1.47
#R0.s.dry = 0
#R0.r.wet = 1/0.64 = 1.56
#R0.r.dry = 1/0.0.61 = 1.64

# Calculate basic reproduction number for each eco-region (R0,i) following Anderson and May (1991)
# Using seroprevalence estimates
beta_s.wet <- sf * ((1.47 * (mu + delta + gamma)) / (N.M + N.T))  # Wet season, susceptible
beta_r.wet <- (1.56 * (mu + delta + gamma)) / N.L                 # Wet season, recovered
beta_r.dry <- (1.64 * (mu + delta + gamma)) / (N.L + N.T)         # Dry season, recovered

# Define interaction matrices for wet and dry seasons
beta_w <- c(beta_s.wet, beta_s.wet, 0,
            beta_s.wet, beta_s.wet, 0,
            0, 0, beta_r.wet)
dim(beta_w) <- c(3, 3)

beta_d <- c(0, 0, 0,
            0, 0, 0,
            0, 0, beta_r.dry)
dim(beta_d) <- c(3, 3)

# Define parameter sets for wet and dry season models
param_w <- c(mu = mu,  
             b = b,
             gamma = gamma,                                                                    
             delta = delta,                                                           
             beta_s.w = beta_s.wet,
             beta_r.wet = beta_r.wet,
             sf = sf)

param_d <- c(mu = mu, 
             b = b,
             gamma = gamma,                                                                    
             delta = delta,                                                           
             beta_r.dry = beta_r.dry)

###################################################################################
# First simulate 20 years to quasi equilibrium
###################################################################################

# Define SIR model for wet season
SIR_model_w <- function(time, N, param_w) {
  with(as.list(c(N, param_w)), {
    
    dS.T <- b * N.T - mu * S.T - beta_w[1, 1] * S.T * I.T - beta_w[1, 2] * S.T * I.M - beta_w[1, 3] * S.T * I.L                              
    dI.T <- -mu * I.T - gamma * I.T - delta * I.T + beta_w[1, 1] * S.T * I.T + beta_w[1, 2] * S.T * I.M + beta_w[1, 3] * S.T * I.L
    dR.T <- delta * I.T - mu * R.T
    
    dS.M <- b * N.M - mu * S.M - beta_w[2, 1] * S.M * I.T - beta_w[2, 2] * S.M * I.M - beta_w[2, 3] * S.M * I.L                                 
    dI.M <- -mu * I.M - gamma * I.M - delta * I.M + beta_w[2, 1] * S.M * I.T + beta_w[2, 2] * S.M * I.M + beta_w[2, 3] * S.M * I.L 
    dR.M <- delta * I.M - mu * R.M
    
    dS.L <- b * N.L - mu * S.L - beta_w[3, 1] * S.L * I.T - beta_w[3, 2] * S.L * I.M - beta_w[3, 3] * S.L * I.L                            
    dI.L <- -mu * I.L - gamma * I.L - delta * I.L + beta_w[3, 1] * S.L * I.T + beta_w[3, 2] * S.L * I.M + beta_w[3, 3] * S.L * I.L  
    dR.L <- delta * I.L - mu * R.L
    
    return(list(c(dS.T, dI.T, dR.T, dS.M, dI.M, dR.M, dS.L, dI.L, dR.L)))
  })
}

# Define SIR model for dry season
SIR_model_d <- function(time, N, param_d) {
  with(as.list(c(N, param_d)), {
    
    dS.T <- b * N.T - mu * S.T - beta_d[1, 1] * S.T * I.T - beta_d[1, 2] * S.T * I.M - beta_d[1, 3] * S.T * I.L                              
    dI.T <- -mu * I.T - gamma * I.T - delta * I.T + beta_d[1, 1] * S.T * I.T + beta_d[1, 2] * S.T * I.M + beta_d[1, 3] * S.T * I.L
    dR.T <- delta * I.T - mu * R.T
    
    dS.M <- b * N.M - mu * S.M - beta_d[2, 1] * S.M * I.T - beta_d[2, 2] * S.M * I.M - beta_d[2, 3] * S.M * I.L                                 
    dI.M <- -mu * I.M - gamma * I.M - delta * I.M + beta_d[2, 1] * S.M * I.T + beta_d[2, 2] * S.M * I.M + beta_d[2, 3] * S.M * I.L 
    dR.M <- delta * I.M - mu * R.M
    
    dS.L <- b * N.L - mu * S.L - beta_d[3, 1] * S.L * I.T - beta_d[3, 2] * S.L * I.M - beta_d[3, 3] * S.L * I.L                            
    dI.L <- -mu * I.L - gamma * I.L - delta * I.L + beta_d[3, 1] * S.L * I.T + beta_d[3, 2] * S.L * I.M + beta_d[3, 3] * S.L * I.L  
    dR.L <- delta * I.L - mu * R.L
    
    return(list(c(dS.T, dI.T, dR.T, dS.M, dI.M, dR.M, dS.L, dI.L, dR.L)))
  })
}

# Helper function to extract the population values from the output
get_N <- function(out_struc) {
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

# Run the simulation for 20 years
time <- seq(0, 245, by = 1)  # First part of year (dry season)
N <- initial_N  # Initialize to almost entirely susceptible
out <- ode(N, time, SIR_model_d, param_d)

out.dat <- as.data.frame(out)
out.dat$I.hn[out.dat$I.hn < 1] <- 0
out.dat$I.m[out.dat$I.m < 1] <- 0
out.dat$I.l[out.dat$I.l < 1] <- 0

all_out <- out.dat

max_years <- 20

# Loop over each year and simulate wet and dry seasons
for (years in 1:max_years) {
  N <- get_N(out)  # Initialize to the end of last integration period
  
  # Wet season simulation
  day <- ((years - 1) * 365 + 1):((years - 1) * 365 + 120)  # Wet season
  out <- ode(N, day, SIR_model_w, param_w)
  
  out.dat <- as.data.frame(out)
  out.dat$I.hn[out.dat$I.hn < 1] <- 0
  out.dat$I.m[out.dat$I.m < 1] <- 0
  out.dat$I.l[out.dat$I.l < 1] <- 0
  all_out <- rbind(all_out, out.dat)  # Add the new output to the cumulative output
  
  # Update N for the dry season
  N <- get_N(out)
  
  # Dry season simulation
  day <- ((years - 1) * 365 + 121):((years - 1) * 365 + 120 + 245)  # Dry season
  out <- ode(N, day, SIR_model_d, param_d)
  
  out.dat <- as.data.frame(out)
  out.dat$I.hn[out.dat$I.hn < 1] <- 0
  out.dat$I.m[out.dat$I.m < 1] <- 0
  out.dat$I.l[out.dat$I.l < 1] <- 0
  all_out <- rbind(all_out, out.dat)  # Add the new output to the cumulative output
}

# The script ends here. You can proceed to visualize the results.