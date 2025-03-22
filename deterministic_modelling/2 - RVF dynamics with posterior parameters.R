# This script simulates Rift Valley fever (RVF) dynamics in three cattle subpopulations: T, M and L, which is modelled across two 
# distinct eco-regions: Sahelian (s) and Gambia river (r). It utilizes a Susceptible-Infected-Recovered (SIR) model with separate 
# dynamics for the wet and dry seasons.

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

# Wet season and dry season beta values following standard expression from the model
beta_s.wet <- sf * ((1.47 * (mu + delta + gamma)) / (N.M + N.T))
beta_r.wet <- (1.56 * (mu + delta + gamma)) / N.L
beta_r.dry <- (1.64 * (mu + delta + gamma)) / (N.L + N.T)

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
# Deterministic simulation: 20 years to quasi-equilibrium
###################################################################################

# Wet season model 

SIR_model_w <- function(time, N, param_w) {
  with(as.list(c(N, param_w)), {
    # Rate of change for each compartment
    dS.T <- b * N.T - mu * S.T - beta_s.wet * S.T * I.T - beta_s.wet * S.T * I.M                            
    dI.T <- -mu * I.T - gamma * I.T - delta * I.T + beta_s.wet * S.T * I.T + beta_s.wet * S.T * I.M 
    dR.T <- delta * I.T - mu * R.T
    dS.M <- b * N.M - mu * S.M - beta_s.wet * S.M * I.T - beta_s.wet * S.M * I.M                              
    dI.M <- -mu * I.M - gamma * I.M - delta * I.M + beta_s.wet * S.M * I.T + beta_s.wet * S.M * I.M 
    dR.M <- delta * I.M - mu * R.M
    dS.L <- b * N.L - mu * S.L - beta_r.wet * S.L * I.L                            
    dI.L <- -mu * I.L - gamma * I.L - delta * I.L + beta_r.wet * S.L * I.L  
    dR.L <- delta * I.L - mu * R.L
    
    return(list(c(dS.T, dI.T, dR.T, dS.M, dI.M, dR.M, dS.L, dI.L, dR.L)))
  })
}

# Dry season model 

SIR_model_d <- function(time, N, param_d) {
  with(as.list(c(N, param_d)), {
    # Rate of change for each compartment
    dS.T <- b * N.T - mu * S.T - beta_r.dry * S.T * I.T - beta_r.dry * S.T * I.L                              
    dI.T <- -mu * I.T - gamma * I.T - delta * I.T + beta_r.dry * S.T * I.T + beta_r.dry * S.T * I.L
    dR.T <- delta * I.T - mu * R.T
    dS.M <- b * N.M - mu * S.M  # No transmission in dry season for migrant cattle
    dI.M <- -mu * I.M - gamma * I.M - delta * I.M
    dR.M <- delta * I.M - mu * R.M
    dS.L <- b * N.L - mu * S.L - beta_r.dry * S.L * I.T - beta_r.dry * S.L * I.L                            
    dI.L <- -mu * I.L - gamma * I.L - delta * I.L + beta_r.dry * S.L * I.T + beta_r.dry * S.L * I.L  
    dR.L <- delta * I.L - mu * R.L
    
    return(list(c(dS.T, dI.T, dR.T, dS.M, dI.M, dR.M, dS.L, dI.L, dR.L)))
  })
}

# Function to extract compartment values from the output structure

get_N <- function(out_struc) {
  N <- matrix(0, nrow = 9)
  N <- c(out_struc[round(length(out_struc)/10, 0), 2],
         out_struc[round(length(out_struc)/10, 0), 3], 
         out_struc[round(length(out_struc)/10, 0), 4],
         out_struc[round(length(out_struc)/10, 0), 5],
         out_struc[round(length(out_struc)/10, 0), 6],
         out_struc[round(length(out_struc)/10, 0), 7],
         out_struc[round(length(out_struc)/10, 0), 8],
         out_struc[round(length(out_struc)/10, 0), 9],
         out_struc[round(length(out_struc)/10, 0), 10])
  return(N)
}

# Define simulation time for dry season (245 days)
time <- seq(0, 245, by = 1)

N <- initial_N   

# Simulate for the first year (dry season)
out <- ode(N, time, SIR_model_d, param_d)

# Combine outputs for further years

all_out <- out

max_years <- 20  # Simulation period (in years)

# Loop through 20 years, alternating between dry and wet seasons

for (years in 1:max_years) {
  N <- get_N(out)  # Initialize with the end of last integration period
  
  # Total population for each eco-region
  N.T <- sum(N[1:3])
  N.M <- sum(N[4:6])
  N.L <- sum(N[7:9])
  
  # Simulate wet season
  day <- ((years - 1) * 365 + 1):((years - 1) * 365 + 120)
  out <- ode(N, day, SIR_model_w, param_w)
  all_out <- rbind(all_out, out)  # Append the new output
  
  # Initialize for dry season
  N <- get_N(out)
  
  # Simulate dry season
  day <- ((years - 1) * 365 + 121):((years - 1) * 365 + 120 + 245)
  out <- ode(N, day, SIR_model_d, param_d)
  all_out <- rbind(all_out, out)  # Append the new output
}

# Visualization: Plotting the infection dynamics (Figure 3)

title <- bquote("RVF infection dynamics")

# Create main plot (Figure 3)
Aqu <- ggplot(all_out, aes(x = time)) +
  ggtitle(bquote(atop(bold(.(title))))) +
  geom_line(aes(y = I.T / (S.T + I.T + R.T), colour = "I.T")) +
  geom_line(aes(y = I.L / (S.L + I.L + R.L), colour = "I.L")) +
  geom_line(aes(y = I.M / (S.M + I.M + R.M), colour = "I.M")) +
  ylab(label = "Proportion of infectious cattle") +
  xlab(label = "Time (days)") +
  theme(legend.position = c(0.95, 0.85),
        legend.direction = "horizontal",
        legend.title = element_text(size = 7, face = "bold"),
        legend.background = element_rect(fill = '#FFFFFF', size = 0.4, linetype = "solid"),
        legend.text = element_text(size = 7)) +
  scale_colour_manual("Compartments", 
                      breaks = c("I.T", "I.L", "I.M"),
                      values = c("red", "blue", "darkgreen")) +
  scale_y_continuous(limit = c(0, 0.2), breaks = c(0.02, 0.04, 0.08, 0.2)) +
  coord_trans(y = "sqrt")

# Inset plot (for Figure 3 inset)
df1 <- all_out[6817:7546, , drop = FALSE]
df1 <- as.data.frame(df1)
df1$time <- 6571:7300

Aqu2 <- ggplot(df1, aes(x = time)) +
  geom_line(aes(y = I.T / (S.T + I.T + R.T), colour = "I.H")) +
  geom_line(aes(y = I.L / (S.L + I.L + R.L), colour = "I.L")) +
  geom_line(aes(y = I.M / (S.M + I.M + R.M), colour = "I.M")) +
  annotate("rect", xmin = 6571, xmax = 6691, ymin = 0, ymax = 0.2, alpha = 0.5, fill = "antiquewhite") +
  annotate("rect", xmin = 6691, xmax = 6936, ymin = 0, ymax = 0.2, alpha = 0.5, fill = "azure2") +
  annotate("text", x = 6630, y = 0.19, label = "Wet Season") +
  annotate("text", x = 6820, y = 0.19, label = "Dry Season") +
  annotate("rect", xmin = 6936, xmax = 7056, ymin = 0, ymax = 0.2, alpha = 0.5, fill = "antiquewhite") +
  annotate("rect", xmin = 7056, xmax = 7300, ymin = 0, ymax = 0.2, alpha = 0.5, fill = "azure2") +
  annotate("text", x = 7000, y = 0.19, label = "Wet Season") +
  annotate("text", x = 7180, y = 0.19, label = "Dry Season") +
  ylab(label = "Proportion of infectious cattle") +
  xlab(label = "Time (days)") +
  scale_colour_manual("Compartments", 
                      breaks = c("I.H", "I.L", "I.M"),
                      values = c("red", "blue", "darkgreen")) +
  scale_y_continuous(limit = c(0, 0.2), breaks = c(0.01, 0.03, 0.05, 0.1, 0.2)) +
  scale_x_continuous(limit = c(6571, 7300)) +
  coord_trans(y = "sqrt")


# Combine main and inset plots of Figure 3
library(patchwork)
 Scho + 
  plot_layout(tag_level = "new") &
  inset_element(Aqu2, left = 0.5, bottom = 0.55, right = 0.99, top = 0.99)

