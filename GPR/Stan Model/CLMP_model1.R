setwd('/home/moosehunter/R/Fujita Lab/GPR/')
library(openesm)
library(rlang)
library(dplyr)
library(ggplot2)
library(tidyr)

#simulate the data 
set.seed(123)

N <- 30
T <- 40

beta_x <- 0.5
beta_y <- 0.6
delta  <- -0.3

sigma_x <- 0.5
sigma_y <- 0.5

# Create smooth time-varying cross-lag
time <- 2:T
gamma_true <- 0.4 * sin( (time-1)/T * pi )  # smooth nonlinear pattern

# Storage
x <- matrix(0, N, T)
y <- matrix(0, N, T)

# Initial values
x[,1] <- rnorm(N, 0, 1)
y[,1] <- rnorm(N, 0, 1)

for (t in 2:T) {
  x[,t] <- beta_x * x[,t-1] +
           delta  * y[,t-1] +
           rnorm(N, 0, sigma_x)

  y[,t] <- beta_y * y[,t-1] +
           gamma_true[t-1] * x[,t-1] +
           rnorm(N, 0, sigma_y)
}


#estimation 
library(rstan)

stan_data <- list(
  N = N,
  T = T,
  x = x,
  y = y
)

fit <- stan(
  file = "Stan Model/CLPM_model1.stan",
  data = stan_data,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4
)

print(fit, pars = c("beta_x","beta_y","delta","eta","rho","sigma_x","sigma_y"))

#recovering parameters
posterior_gamma <- extract(fit)$gamma
gamma_est <- apply(posterior_gamma, 2, mean)

plot(time, gamma_true, type="l", lwd=2, col="black",
     ylim=range(c(gamma_true, gamma_est)),
     ylab="gamma(t)", xlab="time")

lines(time, gamma_est, col="red", lwd=2)
legend("topright", legend=c("True","Estimated"),
       col=c("black","red"), lwd=2)


time_obs <- 2:T
time_grid <- seq(2, T, length.out = 200)
rbf_kernel <- function(x1, x2, eta, rho) {
  outer(x1, x2, function(a, b)
    eta^2 * exp(-(a - b)^2 / (2 * rho^2)))
}

library(MASS)

eta_post <- extract(fit)$eta
rho_post <- extract(fit)$rho


n_draws <- nrow(posterior_gamma)

gamma_smooth <- matrix(NA, n_draws, length(time_grid))

for (s in 1:n_draws) {

  K     <- rbf_kernel(time_obs, time_obs, eta_post[s], rho_post[s])
  K_inv <- solve(K + 1e-6 * diag(nrow(K)))  # jitter for stability

  K_star <- rbf_kernel(time_grid, time_obs,
                       eta_post[s], rho_post[s])

  gamma_smooth[s, ] <- K_star %*% K_inv %*% posterior_gamma[s, ]
}

gamma_mean <- apply(gamma_smooth, 2, mean)
gamma_low  <- apply(gamma_smooth, 2, quantile, 0.025)
gamma_high <- apply(gamma_smooth, 2, quantile, 0.975)
plot(time_grid, gamma_mean,
     type = "l",
     lwd = 2,
     ylim = range(c(gamma_low, gamma_high)),
     ylab = "gamma(t)",
     xlab = "time")

lines(time_grid, gamma_low, lty = 2)
lines(time_grid, gamma_high, lty = 2)

points(time_obs,
       colMeans(posterior_gamma),
       pch = 16)
