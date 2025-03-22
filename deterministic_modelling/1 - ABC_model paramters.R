# This script implements the Approximate Bayesian Computation (ABC) analysis to estimate parameter values of the 
# epidemiological model by fitting to observed seroprevalence data for the three cattle subpopulations

# Clear the environment

rm(list = ls())

# Load the necessary libraries
library(deSolve)
library(ggplot2)
library(reshape2)
library(abc)
library(EnvStats)
library(triangle)

# Simulate observed seroprevalences for 3 populations
observed_seroprevalence <- c(0.45, 0.3, 0.36)

### Initial conditions ###
simulated_serop <- function(parameters) {
  
  beta_s.w <- parameters[1]
  beta_r.w <- parameters[2]
  beta_r.d <- parameters[3]
  mu <- parameters[4]
  b <- parameters[5]
  gamma <- parameters[6]                                                                    
  delta <- parameters[7]                                                   
  sf <- parameters[8]     
  
  # Initial conditions
  initial_N <- c(S.T = 13165.62, I.T = 1.8447889, R.T = 10093.57,       
                 S.M = 83395.99, I.M = 1.537987e-18, R.M = 55182.37,  
                 S.L = 67037.90, I.L = 9.404628, R.L = 42077.79)
  
  N.T <- sum(initial_N[1:3])  
  N.M <- sum(initial_N[4:6])  
  N.L <- sum(initial_N[7:9])  
  
  param_w <- c(mu = mu, b = b, gamma = gamma, delta = delta, 
               beta_s.w = beta_s.w, beta_r.w = beta_r.w, sf = sf)
  
  param_d <- c(mu = mu, b = b, gamma = gamma, delta = delta, 
               beta_r.d = beta_r.d)
  
  # Define the ODE-based model
  SIR_model_w <- function(time, N, param_w) {
    with(as.list(c(N, param_w)), {
      # Rate of change for each compartment
      dS.T <- b * N.T - mu * S.T - sf * beta_s.w * S.T * I.T - sf * beta_s.w * S.T * I.M
      dI.T <- -mu * I.T - gamma * I.T - delta * I.T + sf * beta_s.w * S.T * I.T + sf * beta_s.w * S.T * I.M 
      dR.T <- delta * I.T - mu * R.T
      dS.M <- b * N.M - mu * S.M - sf * beta_s.w * S.M * I.T - sf * beta_s.w * S.M * I.M
      dI.M <- -mu * I.M - gamma * I.M - delta * I.M + sf * beta_s.w * S.M * I.T + sf * beta_s.w * S.M * I.M 
      dR.M <- delta * I.M - mu * R.M
      dS.L <- b * N.L - mu * S.L - beta_r.w * S.L * I.L
      dI.L <- -mu * I.L - gamma * I.L - delta * I.L + beta_r.w * S.L * I.L
      dR.L <- delta * I.L - mu * R.L
      
      output <- c(dS.T, dI.T, dR.T, dS.M, dI.M, dR.M, dS.L, dI.L, dR.L)
      list(output)  
    })
  }
  
  # Function to extract population sizes after simulation
  get_N <- function(out_struc) {
    myN <- matrix(0, nrow = 9)
    myN <- c(out_struc[round(length(out_struc)/10,0), 2], out_struc[round(length(out_struc)/10,0), 3], 
             out_struc[round(length(out_struc)/10,0), 4], out_struc[round(length(out_struc)/10,0), 5], 
             out_struc[round(length(out_struc)/10,0), 6], out_struc[round(length(out_struc)/10,0), 7], 
             out_struc[round(length(out_struc)/10,0), 8], out_struc[round(length(out_struc)/10,0), 9], 
             out_struc[round(length(out_struc)/10,0), 10])
    return(myN)
  }
  
  time <- seq(0, 245, by = 1)
  N <- initial_N
  
  suppressWarnings(out <- ode(N, time, SIR_model_w, param_w))
  all_out <- out
  
  max_years <- 20
  for (years in 1:max_years) {
    N <- get_N(out)
    day <- ((years-1)*365+1):((years-1)*365+120)
    suppressWarnings(try(out <- ode(N, day, SIR_model_w, param_w), TRUE))
    all_out <- rbind(all_out, out)
    N <- get_N(out)
    day <- ((years-1)*365+121):((years-1)*365+120+245)
    suppressWarnings(try(out <- ode(N, day, SIR_model_d, param_d), TRUE))
    all_out <- rbind(all_out, out)
  }
  
  df <- get_N(all_out)
  df <- as.data.frame(df)
  simulated_seroprevalence <- c(P.T = df[3,1] / (df[1,1] + df[2,1] + df[3,1]), 
                                P.M = df[6,1] / (df[4,1] + df[5,1] + df[6,1]),        
                                P.L = df[9,1] / (df[7,1] + df[8,1] + df[9,1]))
  return(simulated_seroprevalence)
}

# Set up prior parameters
n = 500000
prior_parameters = matrix(0, 8, 3)
fff <- 0.7

prior_parameters[1,3]<-1.674599e-06; prior_parameters[1,1]<-prior_parameters[1,3]-fff*prior_parameters[1,3]; prior_parameters[1,2]<-prior_parameters[1,3]+fff*prior_parameters[1,3]
prior_parameters[2,3]<-3.910813e-06; prior_parameters[2,1]<-prior_parameters[2,3]-fff*prior_parameters[2,3]; prior_parameters[2,2]<-prior_parameters[2,3]+fff*prior_parameters[2,3]
prior_parameters[3,3]<-2.854882e-06; prior_parameters[3,1]<-prior_parameters[3,3]-fff*prior_parameters[3,3]; prior_parameters[3,2]<-prior_parameters[3,3]+fff*prior_parameters[3,3]
prior_parameters[4,3]<-0.00164; prior_parameters[4,1]<-prior_parameters[4,3]-fff*prior_parameters[4,3]; prior_parameters[4,2]<-prior_parameters[4,3]+fff*prior_parameters[4,3]
prior_parameters[5,3]<-0.00215; prior_parameters[5,1]<-prior_parameters[5,3]-fff*prior_parameters[5,3]; prior_parameters[5,2]<-prior_parameters[5,3]+fff*prior_parameters[5,3]
prior_parameters[6,3]<-0.075; prior_parameters[6,1]<-prior_parameters[6,3]-fff*prior_parameters[6,3]; prior_parameters[6,2]<-prior_parameters[6,3]+fff*prior_parameters[6,3]
prior_parameters[7,3]<-0.125; prior_parameters[7,1]<-prior_parameters[7,3]-fff*prior_parameters[7,3]; prior_parameters[7,2]<-prior_parameters[7,3]+fff*prior_parameters[7,3]
prior_parameters[8,3]<-2; prior_parameters[8,1]<-prior_parameters[8,3]-fff*prior_parameters[8,3]; prior_parameters[8,2]<-prior_parameters[8,3]+fff*prior_parameters[8,3]

all_parameters <-NULL

# Progress bar initialization
pb <- txtProgressBar(min = 0, max = n, style = 3)

# Main simulation loop
incomplete <- 0
summary_stats <- matrix(0, ncol = 3, nrow = 0)
for (row in 1:n) {
  parameters <- c(
    beta_s.w = rtriangle(1, prior_parameters[1,1], prior_parameters[1,2], prior_parameters[1,3]),
    beta_r.w = rtriangle(1, prior_parameters[2,1], prior_parameters[2,2], prior_parameters[2,3]),
    beta_r.d = rtriangle(1, prior_parameters[3,1], prior_parameters[3,2], prior_parameters[3,3]),
    mu = rtriangle(1, prior_parameters[4,1], prior_parameters[4,2], prior_parameters[4,3]),  
    b = rtriangle(1, prior_parameters[5,1], prior_parameters[5,2], prior_parameters[5,3]),  
    gamma = rtriangle(1, prior_parameters[6,1], prior_parameters[6,2], prior_parameters[6,3]),                                                                  
    delta = rtriangle(1, prior_parameters[7,1], prior_parameters[7,2], prior_parameters[7,3]),                                                   
    sf = rtriangle(1, prior_parameters[8,1], prior_parameters[8,2], prior_parameters[8,3])
  )
  
  all_parameters <- rbind(all_parameters, parameters)
  
  outcomes <- simulated_serop(parameters)
  
  if (any(is.nan(outcomes))) {
    incomplete <- incomplete + 1
  } else {
    summary_stats <- rbind(summary_stats, outcomes)
  }
  
  setTxtProgressBar(pb, row)
}

# Perform ABC using a tolerance of 0.002
tolerance <- 0.002
abc_results <- abc(
  target = observed_seroprevalence,
  param = all_parameters,
  method = "rejection",
  tol = tolerance,
  n = n,
  sumstat = summary_stats
)

# Summary of the ABC results
options(digits = 10)
summary(abc_results)
