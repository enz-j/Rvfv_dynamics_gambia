# Function to compute next configuration for dry season
get_next_configuration <- function(out_row) {
  # Define dry season events
  events <- c(
    b * sum(out_row[c(2, 3, 4)]),  # Host birth
    mu * out_row[2],  # Death of S.T
    mu * out_row[3],  # Death of I.T (unrelated to infection)
    mu * out_row[4],  # Death of R.T
    beta_r.d * out_row[2] * out_row[3],  # Infection from I.T
    beta_s.d * out_row[2] * out_row[6],  # Infection from I.M
    beta_r.d * out_row[2] * out_row[9],  # Infection from I.L
    gamma * out_row[3],  # Death of I.T due to infection
    delta * out_row[3],  # Recovery of I.T
    b * sum(out_row[c(5, 6, 7)]),  # Host birth
    mu * out_row[5],  # Death of S.M
    mu * out_row[6],  # Death of I.M (unrelated to infection)
    mu * out_row[7],  # Death of R.M
    beta_s.d * out_row[5] * out_row[3],  # Infection from I.T
    beta_s.d * out_row[5] * out_row[6],  # Infection from I.M
    beta_s.d * out_row[5] * out_row[9],  # Infection from I.L
    gamma * out_row[6],  # Death of I.M due to infection
    delta * out_row[6],  # Recovery of I.M
    b * sum(out_row[c(8, 9, 10)]),  # Host birth
    mu * out_row[8],  # Death of S.L
    mu * out_row[9],  # Death of I.L (unrelated to infection)
    mu * out_row[10],  # Death of R.L
    beta_r.d * out_row[8] * out_row[3],  # Infection from I.T
    beta_s.d * out_row[8] * out_row[6],  # Infection from I.M
    beta_r.d * out_row[8] * out_row[9],  # Infection from I.L
    gamma * out_row[9],  # Death of I.L due to infection
    delta * out_row[9]   # Recovery of I.L
  )
  
  # Apply Poisson distribution
  K <- rpois(length(events), events * tau)
  
  # Update state based on event outcomes
  out_row[c(2, 3, 4)] <- out_row[c(2, 3, 4)] + c(K[1] - K[2] - K[5] - K[7], -K[3] + K[5] + K[7] - K[8] - K[9], -K[4] + K[9])
  out_row[c(5, 6, 7)] <- out_row[c(5, 6, 7)] + c(K[10] - K[11], -K[12] - K[17] - K[18], -K[13] + K[18])
  out_row[c(8, 9, 10)] <- out_row[c(8, 9, 10)] + c(K[19] - K[20] - K[23] - K[25], -K[21] + K[23] + K[25] - K[26] - K[27], -K[22] + K[27])
  
  # Ensure populations do not go negative
  out_row[2:10] <- pmax(out_row[2:10], 0)
  
  # Update time
  out_row[1] <- out_row[1] + tau
  
  return(list(out_row = out_row))
}

# Initialize parameters
mu <- 0.0016940    
gamma <- 0.0758493  
delta <- 0.1161880  
b <- 0.0021754      
sf <- 2.0508335

# Infection rates
beta_s.w <- sf * 0.0000014
beta_r.d <- 0.0000032
beta_s.d <- 0
beta_r.w <- 0.0000037

# Initial conditions
initial_conditions <- c(0, 13165.62, 1.8447889, 10093.57, 83395.99, 1.537987e-18, 55182.37, 67037.90, 9.404628, 42077.79)

# Run a single step of simulation
out_row <- initial_conditions

# Simulation function with extinction tracking
run_simulation_with_extinction <- function(nsim) {
  all_results <- list()
  extinction_times <- numeric(nsim)
  extinction_count <- 0
  
  for (sim in 1:nsim) {
    out <- matrix(0, nrow = 1, ncol = 10)
    out[1, ] <- c(0, 13165.62, 1.8447889, 10093.57,
                  83395.99, 1.537987e-18, 55182.37,
                  67037.90, 9.404628, 42077.79)
    
    out_row <- out[1, ]
    extinction_time <- NA
    
    for (year in 1:20) {
      # Dry season simulation (245 days)
      for (day in ((year-1)*365+1):((year-1)*365+245)) {
        while (out_row[1] <= day) {
          mylist <- get_next_configuration(out_row)
          out_row <- mylist$out_row
          out <- rbind(out, out_row)
          
          if (out_row[3] == 0 && out_row[6] == 0 && out_row[9] == 0) {
            extinction_time <- out_row[1]
          }
        }
      }
      
      # Wet season simulation (120 days)
      for (day in ((year-1)*365+245.01):((year-1)*365+245+120)) {
        while (out_row[1] <= day) {
          mylist <- get_next_configuration(out_row)
          out_row <- mylist$out_row
          out <- rbind(out, out_row)
          
          if (out_row[3] == 0 && out_row[6] == 0 && out_row[9] == 0) {
            extinction_time <- out_row[1]
          }
        }
      }
    }
    
    # Track extinction occurrences and times
    if (!is.na(extinction_time)) {
      extinction_count <- extinction_count + 1
      extinction_times[sim] <- extinction_time
    }
    
    # Convert to data frame and add a column for the simulation number
    out_df <- as.data.frame(out)
    out_df$Simulation <- sim
    all_results[[sim]] <- out_df
  }
  
  extinction_prob <- extinction_count / nsim
  list(results = do.call(rbind, all_results), extinction_prob = extinction_prob, extinction_times = extinction_times)
}

# Run the simulations with extinction tracking
nsim <- 1000

simulation_result <- run_simulation_with_extinction(nsim)
result.rep <- simulation_result$results

# Function to count the unique simulations where I.T = 0 and I.T, I.M, and I.L = 0
count_extinctions <- function(simulation_result) {
  
  # Initialize counters
  count_I.T_zero <- 0
  count_all_zero <- 0
  
  # Iterate over all simulations
  for (sim in unique(result.rep$Simulation)) {
    sim_data <- subset(result.rep, Simulation == sim)
    
    # Check if I.T = 0 in any time step for this simulation
    if (any(sim_data$I.T == 0)) {
      count_I.T_zero <- count_I.T_zero + 1
    }
    
    # Check if I.T, I.M, and I.L = 0 in any time step for this simulation
    if (any(sim_data$I.T == 0 & sim_data$I.M == 0 & sim_data$I.L == 0)) {
      count_all_zero <- count_all_zero + 1
    }
  }
  
  return(list(count_I.T_zero = count_I.T_zero, count_all_zero = count_all_zero))
}

# Count occurrences of I.T = 0 and I.T, I.M, and I.L = 0
extinction_counts <- count_extinctions(simulation_result)


library(ggplot2)
library(tidyr)
library(dplyr)

# Create the ggplot for simulation results
ggplot(result.rep, aes(x = time)) +
  geom_line(aes(y=(I.T/(S.T+I.T+R.T)),colour="I.T"), size=0.5) +
  geom_line(aes(y=(I.L/(S.L+I.L+R.L)),colour="I.L"), size=0.5) +
  geom_line(aes(y=(I.M/(S.M+I.M+R.M)),colour="I.M"), size=0.5) +
  
  facet_wrap(~ Simulation) +
  labs(title = 'Simulation Results', x = 'Time (Days)', y = 'Proportion infected') +
  theme_minimal() +
  theme(legend.position = "bottom") + 
  scale_colour_manual("Compartments",
                      breaks=c("I.T","I.L", "I.M"),
                      values=c("red", "blue", "darkgreen"))+
  coord_trans(y='sqrt')