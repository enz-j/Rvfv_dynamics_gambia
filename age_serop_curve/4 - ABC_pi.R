# This script implements the Approximate Bayesian Computation (ABC) analysis to estimate the seropositivity decay
# rate parameter (pi) of the epidemiological model by fitting to observed seroprevalence data for the three cattle 
# subpopulations

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

# Load lambda values
load("lambda_vec.RData")

# Function to simulate seroprevalences based on given parameters
simulated_serop <- function(parameters) {
  pi <- parameters[1]
  gamma <- parameters[2]
  delta <- parameters[3]
  
  # Initial population values
  initial_N <- c(S.T = 35000, I.T = 1, R.T = 0,
                 S.M = 170000, I.M = 1, R.M = 0,
                 S.L = 85000, I.L = 1, R.L = 0) 
  
  # SIR model differential equations
  SIR_model.cr <- function(time, N, lambda_t) {
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
  
  # Function to extract population sizes from the output
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
  
  # Initial lambda values
  lambda_hn <- lambda_vec[1]
  lambda_m <- lambda_vec[1, 2]
  lambda_l <- lambda_vec[1, 3]
  
  lambda_t <- c(lambda_hn, lambda_m, lambda_l)
  param <- c(gamma, pi, lambda_t)
  
  time <- seq(0, 1, by = 1)
  N <- initial_N  # Initialize the populations to be almost entirely susceptible
  out.serop <- ode(N, time, SIR_model.cr, param)
  all_out.serop <- out.serop 
  
  time_grid <- seq(1, 555, by = 1)
  
  # Run the simulation for the specified number of weeks
  for (weeks in 1:length(time_grid)) {
    lambda_hn <- lambda_vec[weeks, 1]
    lambda_m <- lambda_vec[weeks, 2]
    lambda_l <- lambda_vec[weeks, 3]
    
    lambda_t <- c(lambda_hn, lambda_m, lambda_l)
    param <- c(lambda_t, pi, gamma)
    
    N <- get_N(all_out.serop)  # Initialize to the end of the last integration period
    out.serop <- ode(y = N, times = c(weeks, weeks + 1), func = SIR_model.cr, parms = param)
    all_out.serop <- rbind(all_out.serop, out.serop)  # Append the new output
  }
  
  # Extract the final populations
  df <- get_N(all_out.serop)
  df <- as.data.frame(df)
  
  # Calculate the simulated seroprevalence for each population
  simulated_seroprevalence <- c(P.T = (df[2, 1]) / (df[1, 1] + df[2, 1] + df[3, 1]),
                                P.M = (df[5, 1]) / (df[4, 1] + df[5, 1] + df[6, 1]),
                                P.L = (df[8, 1]) / (df[7, 1] + df[8, 1] + df[9, 1]))
  
  return(simulated_seroprevalence)
}

# Define the prior distribution for the parameters
prior_parameters <- matrix(0, 3, 3)
fff <- 0.7

# Set parameter bounds (± 50%)
prior_parameters[1, 3] <- 0.005
prior_parameters[1, 1] <- prior_parameters[1, 3] - fff * prior_parameters[1, 3]
prior_parameters[1, 2] <- prior_parameters[1, 3] + fff * prior_parameters[1, 3]

prior_parameters[2, 3] <- 0.0758493 * 7
prior_parameters[3, 3] <- 0.1161880 * 7

# Initialize variables for ABC
n <- 500000
summary_stats <- NULL
all_parameters <- NULL
incomplete <- 0

# Start the ABC sampling process
start_time <- Sys.time()
interval <- Sys.time()  # Interval in seconds

for (row in 1:n) {
  # Sample parameters from the prior distribution
  parameters <- c(
    pi = rtri(1, prior_parameters[1, 1], prior_parameters[1, 2], prior_parameters[1, 3]),
    gamma = prior_parameters[2, 3],
    delta = prior_parameters[3, 3]
  )
  
  all_parameters <- rbind(all_parameters, parameters)
  
  # Simulate the outcomes for the current parameters
  outcomes <- simulated_serop(parameters)
  
  # Check for incomplete or invalid outcomes (NaN values)
  if (any(is.nan(outcomes))) {
    incomplete <- incomplete + 1
  } else {
    summary_stats <- rbind(summary_stats, outcomes)
  }
  
  # Calculate elapsed time for progress tracking
  elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Update the tracker dynamically in the console
  cat(sprintf("\rIteration: %d/%d | Elapsed Time: %.2f seconds", row, n, elapsed_time))
  Sys.sleep(interval)  # Simulate work
}

# Set the tolerance for ABC
tolerance <- 0.002

# Perform the ABC analysis using rejection sampling
abc_results <- abc(
  target = observed_seroprevalence,
  param = all_parameters,
  method = "rejection",
  tol = tolerance,
  n = n,
  sumstat = summary_stats
)

# Print the summary of the ABC results
options(digits = 10)
summary(abc_results)
