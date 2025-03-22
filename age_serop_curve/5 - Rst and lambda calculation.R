# Estimate ecoclimatic region-specific Rst and R0 values across the live of a cattle (estimated as 10 years)

# Load necessary libraries
library(deSolve)      # For solving differential equations
library(ggplot2)      # For data visualization
library(reshape2)     # For data reshaping
library(patchwork)    # For combining multiple plots

# Clear the environment
rm(list = ls())

### Initial conditions ###
initial_N <- c(
  S.T = 13165.62, I.T = 1.8447889, R.T = 10093.57,       # Transhumant cattle (T)
  S.M = 83395.99, I.M = 1.537987e-18, R.M = 55182.37,    # Resident cattle in Sahel (M)
  S.L = 67037.90, I.L = 9.404628, R.L = 42077.79         # Resident cattle in river (L)
)

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

# Calculate the force of infection parameters for wet and dry seasons
beta_s.wet <- sf * ((R0.s.wet * (mu + delta + gamma)) / (N.M + N.T))  # Wet season, susceptible
beta_r.wet <- (R0.r.wet * (mu + delta + gamma)) / N.L                  # Wet season, recovered
beta_r.dry <- (R0.r.dry * (mu + delta + gamma)) / (N.L + N.T)          # Dry season, recovered

# Define matrices for beta parameters (force of infection) in wet and dry seasons
beta_w <- matrix(c(
  beta_s.wet, beta_s.wet, 0,
  beta_s.wet, beta_s.wet, 0,
  0, 0, beta_r.wet
), nrow = 3, byrow = TRUE)

beta_d <- matrix(c(
  beta_r.dry, 0, beta_r.dry,
  0, 0, 0,
  beta_r.dry, 0, beta_r.dry
), nrow = 3, byrow = TRUE)

### SIR Model for wet and dry seasons ###

# Wet season model function
SIR_model_w <- function(time, N, param_w) {
  with(as.list(c(N, param_w)), {
    # Force of infection calculations
    beta_w.S <- matrix(c(
      beta_s.wet * S.T, beta_s.wet * S.T, 0,
      beta_s.wet * S.M, beta_s.wet * S.M, 0,
      0, 0, beta_r.wet * S.L
    ), nrow = 3, byrow = TRUE)
    
    Trans <- matrix(c(-1 / denominator, 0, 0, 
                      0, -1 / denominator, 0, 
                      0, 0, -1 / denominator), nrow = 3, byrow = TRUE)
    
    vfinv.w <- -beta_w.S %*% Trans
    R0 <- max(eigen(vfinv.w)$values)
    
    # Infection pressures for each population
    lambda_hn <- beta_w[1, 1] * I.T + beta_w[1, 2] * I.M
    lambda_m <- beta_w[2, 1] * I.M + beta_w[2, 2] * I.T
    lambda_l <- beta_w[3, 3] * I.L
    
    # Differential equations for the model
    dS.T <- b * N.T - mu * S.T - beta_w[1, 1] * S.T * I.T - beta_w[1, 2] * S.T * I.M - beta_w[1, 3] * S.T * I.L
    dI.T <- -mu * I.T - gamma * I.T - delta * I.T + beta_w[1, 1] * S.T * I.T + beta_w[1, 2] * S.T * I.M + beta_w[1, 3] * S.T * I.L
    dR.T <- delta * I.T - mu * R.T
    dS.M <- b * N.M - mu * S.M - beta_w[2, 1] * S.M * I.T - beta_w[2, 2] * S.M * I.M - beta_w[2, 3] * S.M * I.L
    dI.M <- -mu * I.M - gamma * I.M - delta * I.M + beta_w[2, 1] * S.M * I.T + beta_w[2, 2] * S.M * I.M + beta_w[2, 3] * S.M * I.L
    dR.M <- delta * I.M - mu * R.M
    dS.L <- b * N.L - mu * S.L - beta_w[3, 1] * S.L * I.T - beta_w[3, 2] * S.L * I.M - beta_w[3, 3] * S.L * I.L
    dI.L <- -mu * I.L - gamma * I.L - delta * I.L + beta_w[3, 1] * S.L * I.T + beta_w[3, 2] * S.L * I.M + beta_w[3, 3] * S.L * I.L
    dR.L <- delta * I.L - mu * R.L
    
    output <- c(dS.T, dI.T, dR.T, dS.M, dI.M, dR.M, dS.L, dI.L, dR.L)
    
    list(output, lambda_hn = lambda_hn, lambda_m = lambda_m, lambda_l = lambda_l, R0 = R0)
  })
}

# Dry season model function (similar to wet season with adjusted parameters)
SIR_model_d <- function(time, N, param_d) {
  with(as.list(c(N, param_d)), {
    # Force of infection calculations
    beta_d.S <- matrix(c(
      beta_r.dry * S.T, 0, beta_r.dry * S.T,
      0, 0, 0,
      beta_r.dry * S.L, 0, beta_r.dry * S.L
    ), nrow = 3, byrow = TRUE)
    
    Trans <- matrix(c(-1 / denominator, 0, 0, 
                      0, -1 / denominator, 0, 
                      0, 0, -1 / denominator), nrow = 3, byrow = TRUE)
    
    vfinv.d <- -beta_d.S %*% Trans
    R0 <- max(eigen(vfinv.d)$values)
    
    # Infection pressures for each population
    lambda_hn <- beta_d[1, 1] * I.T + beta_d[1, 3] * I.L
    lambda_m <- 0 * I.M
    lambda_l <- beta_d[3, 1] * I.L + beta_d[3, 3] * I.T
    
    # Differential equations for the model
    dS.T <- b * N.T - mu * S.T - beta_d[1, 1] * S.T * I.T - beta_d[1, 2] * S.T * I.M - beta_d[1, 3] * S.T * I.L
    dI.T <- -mu * I.T - gamma * I.T - delta * I.T + beta_d[1, 1] * S.T * I.T + beta_d[1, 2] * S.T * I.M + beta_d[1, 3] * S.T * I.L
    dR.T <- delta * I.T - mu * R.T
    dS.M <- b * N.M - mu * S.M - beta_d[2, 1] * S.M * I.T - beta_d[2, 2] * S.M * I.M - beta_d[2, 3] * S.M * I.L
    dI.M <- -mu * I.M - gamma * I.M - delta * I.M + beta_d[2, 1] * S.M * I.T + beta_d[2, 2] * S.M * I.M + beta_d[2, 3] * S.M * I.L
    dR.M <- delta * I.M - mu * R.M
    dS.L <- b * N.L - mu * S.L - beta_d[3, 1] * S.L * I.T - beta_d[3, 2] * S.L * I.M - beta_d[3, 3] * S.L * I.L
    dI.L <- -mu * I.L - gamma * I.L - delta * I.L + beta_d[3, 1] * S.L * I.T + beta_d[3, 2] * S.L * I.M + beta_d[3, 3] * S.L * I.L
    dR.L <- delta * I.L - mu * R.L
    
    output <- c(dS.T, dI.T, dR.T, dS.M, dI.M, dR.M, dS.L, dI.L, dR.L)
    
    list(output, lambda_hn = lambda_hn, lambda_m = lambda_m, lambda_l = lambda_l, R0 = R0)
  })
}

# Function to extract the current state of populations
get_N <- function(out_struc) {
  myN <- c(
    out_struc[round(length(out_struc) / 14, 0), 2],
    out_struc[round(length(out_struc) / 14, 0), 3],
    out_struc[round(length(out_struc) / 14, 0), 4],
    out_struc[round(length(out_struc) / 14, 0), 5],
    out_struc[round(length(out_struc) / 14, 0), 6],
    out_struc[round(length(out_struc) / 14, 0), 7],
    out_struc[round(length(out_struc) / 14, 0), 8],
    out_struc[round(length(out_struc) / 14, 0), 9],
    out_struc[round(length(out_struc) / 14, 0), 10]
  )
  return(myN)
}

# Time vector from 0 to 245 days, with a time step of 1 day
time <- seq(0, 245, by = 1)

# Initial population values (ensure 'intial_N' is correctly defined in your script)
N <- intial_N

# Running the ODE model for the dry season (using SIR_model_d)
out_wk <- ode(N, time, SIR_model_d, param_d)

# Extract data at 7-day intervals for the dry season output
all_out_wk <- out_wk[seq(7, nrow(out_wk), by = 7), ]

# Define the number of years for the simulation
max_years <- 10

# Loop through the years, simulating both wet and dry seasons
for (years in 1:max_years) {
  # Initialize populations at the end of the previous simulation period
  N <- get_N(out_wk)
  
  # Update population totals for each region (T, M, L)
  N.T <- sum(N[1:3])
  N.M <- sum(N[4:6])
  N.L <- sum(N[7:9])
  
  # Wet season (day 1 to day 120 of the year)
  wet_start <- (years - 1) * 365 + 1
  wet_end <- (years - 1) * 365 + 120
  week <- seq(wet_start, wet_end, by = 1) # Wet season time period
  
  # Run the ODE model for the wet season (using SIR_model_w)
  out_wk <- ode(N, week, SIR_model_w, param_w)
  
  # Extract data at 7-day intervals for the wet season output
  all_out_wk1 <- out_wk[seq(7, nrow(out_wk), by = 7), ]
  all_out_wk <- rbind(all_out_wk, all_out_wk1)
  
  # Update populations at the end of the wet season
  N <- get_N(out_wk)
  N.T <- sum(N[1:3])
  N.M <- sum(N[4:6])
  N.L <- sum(N[7:9])
  
  # Dry season (day 121 to day 245 of the year)
  dry_start <- (years - 1) * 365 + 121
  dry_end <- (years - 1) * 365 + 120 + 245
  week <- seq(dry_start, dry_end, by = 1) # Dry season time period
  
  # Run the ODE model for the dry season (using SIR_model_d)
  out_wk <- ode(N, week, SIR_model_d, param_d)
  
  # Extract data at 7-day intervals for the dry season output
  all_out_wk1 <- out_wk[seq(7, nrow(out_wk), by = 7), ]
  all_out_wk <- rbind(all_out_wk, all_out_wk1)
}

########################## Plotting the results for 10 years ##########################

# Load necessary libraries
library(ggplot2)
library(ggpubr)

# Convert the output into a data frame
df <- as.data.frame(all_out_wk)

# Add a time variable to represent weeks
df$time <- 1:555

# Set a custom theme for the plot
mytheme4 <- theme_bw() +
  theme(text = element_text(colour = "black")) +
  theme(panel.grid = element_line(colour = "#F6F6F6")) +
  theme(panel.background = element_rect(fill = "#FFFFFF"))
theme_set(mytheme4)

# Plotting R0, lambda_s, and lambda_r with ggplot2
title <- bquote("Rst")

lamR <- ggplot(df, aes(x = time)) +
  ggtitle(bquote(atop(bold(.(title))))) +
  geom_line(aes(y = R0, colour = "R.st")) +
  geom_line(aes(y = lambda_m, colour = "lambda.s")) +
  geom_line(aes(y = lambda_l, colour = "lambda.r")) +
  
  # Highlight the wet and dry seasons
  annotate("rect", xmin = 451, xmax = 468, ymin = 0, ymax = 2, alpha = .5, fill = "antiquewhite") +
  annotate("rect", xmin = 468, xmax = 503, ymin = 0, ymax = 2, alpha = .5, fill = "azure2") +
  annotate("text", x = 459, y = 1.9, label = "wet", size = 4) +
  annotate("text", x = 485, y = 1.9, label = "dry", size = 4) +
  
  annotate("rect", xmin = 503, xmax = 520, ymin = 0, ymax = 2, alpha = .5, fill = "antiquewhite") +
  annotate("rect", xmin = 520, xmax = 555, ymin = 0, ymax = 2, alpha = .5, fill = "azure2") +
  annotate("text", x = 511, y = 1.9, label = "wet", size = 4) +
  annotate("text", x = 537, y = 1.9, label = "dry", size = 4) +
  
  ylab(expression(R["0,st"] ~ " and " ~ lambda["i"])) +
  xlab(label = "Time (weeks)") +
  
  # Customize the legend and appearance
  theme(legend.justification = c(1, 0), legend.position = c(0.4, 0.851), legend.direction = "horizontal") +
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.background = element_rect(fill = '#FFFFFF', size = 0.2, linetype = "solid"),
        legend.text = element_text(size = 11),
        legend.key = element_rect(colour = "#FFFFFF", fill = '#C2C2C2', size = 0.4, linetype = "solid")) +
  theme(
    axis.title.x = element_text(size = 17),  # Increase x-axis title size
    axis.title.y = element_text(size = 17),  # Increase y-axis title size
    axis.text.x = element_text(size = 16),   # Increase x-axis tick labels
    axis.text.y = element_text(size = 16)    # Increase y-axis tick labels
  ) +
  
  # Add horizontal line for R0 = 1
  geom_hline(yintercept = 1, colour = "black", linetype = "longdash", size = 0.4) +
  annotate("text", x = 572, y = 1.09, label = expression(R["st"] == 1), parse = TRUE) +
  
  # Customize color scale
  scale_colour_manual("Key Metrics", breaks = c("R.st", "lambda.s", "lambda.r"),
                      values = c("darkgreen", "blue", "red"),
                      labels = c("R.st" = expression(R["st"]), "lambda.s" = expression(lambda["s"]), "lambda.r" = expression(lambda["r"]))) +
  
  # Customize y and x axis limits and breaks
  scale_y_continuous(limit = c(0, 2), breaks = c(0.02, 0.1, 1, 1.5, 2)) +
  scale_x_continuous(limit = c(0, 572), breaks = c(0, 200, 400, 555)) +
  
  # Apply a logarithmic transformation to y-axis
  coord_trans(y = "sqrt")
