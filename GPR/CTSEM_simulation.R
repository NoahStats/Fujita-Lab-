# Load libraries
library(ctsem)
library(Matrix)

set.seed(42)

# --- 1. Define the True Parameters ---
# Negative values on diagonal (a11, a22) = stability
# Off-diagonal (a21) = Cross-lag (L1 influences L2)
drift_matrix <- matrix(c(-0.5,  0.0, 
                          0.3, -0.8), nrow=2, byrow=TRUE)

diffusion_matrix <- diag(0.2, 2) # Random noise (innovation)
intercepts <- c(10, 20)          # The "mean" levels for the person

# --- 2. Generate Irregular Time Points ---
# 50 observations over 100 hours with random gaps
n_obs <- 50
time_points <- cumsum(runif(n_obs, 0.5, 4.0)) 

# --- 3. Simulate the Continuous Trajectory ---
# We use a simple Euler discretization for the simulation
dt <- 0.01
full_time <- seq(0, max(time_points), by=dt)
latent_data <- matrix(0, nrow=length(full_time), ncol=2)
latent_data[1, ] <- intercepts # Starting point

for(i in 2:length(full_time)) {
  # Change = (Drift * (Current - Mean)) * dt + Noise
  diffusion_noise <- rnorm(2, 0, sqrt(dt)) %*% diffusion_matrix
  change <- drift_matrix %*% (latent_data[i-1, ] - intercepts) * dt + diffusion_noise
  latent_data[i, ] <- latent_data[i-1, ] + change
}

# --- 4. Subsample at our "Irregular" Wearable Pings ---
obs_indices <- sapply(time_points, function(t) which.min(abs(full_time - t)))
observed_data <- latent_data[obs_indices, ]

# Final Dataframe in Long Format
df <- data.frame(
  id = 1,
  time = time_points,
  V1 = observed_data[, 1] + rnorm(n_obs, 0, 0.1), # Add measurement error
  V2 = observed_data[, 2] + rnorm(n_obs, 0, 0.1)
)

print(head(df))

library(ggplot2)
library(tidyr)

# Reshape data for plotting
df_plot <- df %>%
  pivot_longer(cols = c(V1, V2), names_to = "Variable", values_to = "Value")

ggplot(df_plot, aes(x = time, y = Value, color = Variable)) +
  geom_line(alpha = 0.5, linetype = "dashed") +  # Connect the irregular points
  geom_point(size = 2) +                         # Show the actual measurement points
  theme_minimal() +
  labs(title = "Irregularly Sampled Time Series (N=1)",
       subtitle = "Simulated Wearable Data: V1 influences V2 over time",
       x = "Time (Hours)",
       y = "Measured Value") +
  scale_color_manual(values = c("V1" = "#E41A1C", "V2" = "#377EB8"))



# 1. Define the model (Continuous-Time bivariate)
# 'stanct' tells ctsem to treat time as continuous
model_spec <- ctModel(
  type = 'ct', 
  n.latent = 2, 
  latentNames = c("L1", "L2"),
  manifestNames = c("V1", "V2"),
  LAMBDA = diag(2) # Map V1 to L1 and V2 to L2
)

# 2. Fit the model to our simulated data
# optimize = TRUE uses Maximum Likelihood (fast for N=1)
# stationary = FALSE is safer for short time series
fit <- ctStanFit(
  ctstanmodel = model_spec,
  datalong = df, 
  optimize = TRUE
)

# 3. View the results
results <- summary(fit)
print(results$parmatrices)


# Calculate discrete-time parameters for a 1-hour lag
discrete_pars <- ctStanDiscretePars(fit, times = 1)
print(discrete_pars)


# 1. Generate the Kalman Filter estimates
# This calculates the expected 'path' based on the fitted parameters
k_estimates <- ctKalman(fit, timerange = c(0, max(df$time)), timestep = 0.1)

# 2. Plot the results
# ctsem has a built-in plot function for Kalman objects
plot(k_estimates, 
     plotcontrol = list(
       xlab = "Time (Hours)", 
       ylab = "Latent Level",
       main = "Model Estimates (Lines) vs. Observed Data (Points)"
     ))
