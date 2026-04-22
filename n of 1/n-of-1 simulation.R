setwd("/home/moosehunter/R/Fujita Lab")


library(rstan)
library(foreach)
library(doParallel)
library(dplyr)

# n-of-1  -----------------------------------------------------------------


set.seed(123)
library(rstan)

T <- 10
t <- 1:T
treat_time <- 6

# Structural parameters
beta_0 <- 0
beta_1 <- 0.1
beta_2 <- 2
beta_3 <- 0.2

# AR parameters
phi   <- 0.3
sigma <- 0.4

# Indicators
I_treat <- ifelse(t >= treat_time, 1, 0)
t_post  <- pmax(0, t - treat_time)

# Structural mean
mu <- beta_0 +
  beta_1 * t +
  beta_2 * I_treat +
  beta_3 * t_post

# Generate y with AR term
y <- numeric(T)

# Stationary initialization
y[1] <- rnorm(
  1,
  mean = mu[1] / (1 - phi),
  sd   = sigma / sqrt(1 - phi^2)
)

for (tt in 2:T) {
  y[tt] <- mu[tt] + phi * y[tt - 1] + rnorm(1, 0, sigma)
}

dat <- data.frame(
  time = t,
  y = y,
  mu = mu,
  treatment = I_treat
)

dat

ggplot(dat, aes(x = time)) +
  geom_line(aes(y = y), alpha = 0.4) +
  geom_point(aes(y = y), size = 2) +
  
  
  geom_vline(xintercept = treat_time, linetype = "dashed") +
  
  labs(
    title = "N-of-1 Data with Level Change, Trend Change, and AR(1) term and Normal errors",
    subtitle = "Solid line = true mean; correlated errors",
    x = "Time",
    y = "Outcome"
  ) +
  theme_minimal()



stan_data <- list(
  T = nrow(dat),
  y = dat$y,
  time = dat$time,
  I_treat = dat$treatment,
  t_post = ifelse(dat$time >= 10, dat$time - 10, 0)
)

fit <- stan(
  file = "n-of-1.stan",
  data = stan_data,
  chains = 4,
  iter = 2000,
  cores = 4
)

print(fit, pars = c(
  "beta_0", "beta_1", "beta_2", "beta_3", "phi", "sigma"
))






# hierarchical n-of-1 T = 10, treat = 5  ----------------------------------------------------

set.seed(123)
library(tidyverse)

# -----------------------
# Design
# -----------------------
N  <- 10
T1 <- 10
time <- 1:T1
treat_time <- 6

I_treat <- ifelse(time >= treat_time, 1, 0)
t_post  <- pmax(0, time - treat_time)

# -----------------------
# Population-level parameters
# -----------------------
beta_mu <- c(0, 0.1, 2, 0.2)
beta_sd <- c(1, 0.5, 1, 0.5)

phi_mu <- 0.3
phi_sd <- 0.2
sigma <- 0.4


# -----------------------
# Individual-level parameters
# -----------------------
beta_0 <- rnorm(N, beta_mu[1], beta_sd[1])
beta_1 <- rnorm(N, beta_mu[2], beta_sd[2])
beta_2 <- rnorm(N, beta_mu[3], beta_sd[3])
beta_3 <- rnorm(N, beta_mu[4], beta_sd[4])

phi <- tanh(rnorm(N, phi_mu, phi_sd))  # constrain to (-1,1)

# -----------------------
# Generate data
# -----------------------
y  <- matrix(NA, N, T1)
mu <- matrix(NA, N, T1)

for (i in 1:N) {
  
  mu[i, ] <-
    beta_0[i] +
    beta_1[i] * time +
    beta_2[i] * I_treat +
    beta_3[i] * t_post
  
  y[i, 1] <- rnorm(
    1,
    mean = mu[i, 1] / (1 - phi[i]),
    sd   = sigma / sqrt(1 - phi[i]^2)
  )
  
  for (t in 2:T1) {
    y[i, t] <- mu[i, t] + phi[i] * y[i, t - 1] + rnorm(1, 0, sigma)
  }
}



mu_pop <-
  beta_mu["beta_0"] +
  beta_mu["beta_1"] * time +
  beta_mu["beta_2"] * I_treat +
  beta_mu["beta_3"] * t_post

y_pop <- numeric(T1)

y_pop[1] <- rnorm(
  1,
  mean = mu_pop[1] / (1 - phi_mu),
  sd   = sigma / sqrt(1 - phi_mu^2)
)

for (t in 2:T1) {
  y_pop[t] <- mu_pop[t] + phi_mu * y_pop[t - 1] + rnorm(1, 0, sigma)
}



df <- expand_grid(
  id = factor(1:N),
  time = time
)

df$y <- as.vector(t(y))

df_pop <- tibble(
  time = time,
  y_pop = y_pop
)


ggplot() +
  # Individual observations (POINTS ONLY)
  geom_point(
    data = df,
    aes(x = time, y = y),
    alpha = 0.6
  ) +
  
  # Population-level AR trajectory (LINE)
  geom_line(
    data = df_pop,
    aes(x = time, y = y_pop),
    linewidth = 1.6
  ) +
  
  geom_vline(xintercept = treat_time, linetype = "dashed") +
  
  labs(
    title = "Population-Level AR Process with Individual Observations",
    subtitle = "Line = population AR trajectory, points = individual data",
    x = "Time",
    y = "y"
  ) +
  theme_minimal()




stan_data1 <- list(
  N = N,
  T = T1, #T1 = 10 
  y = y,
  time = time,
  I_treat = I_treat,
  t_post = t_post
)

#T1 = 10 and one estimation 
fit0 <- stan( 
  file = "model/hierarchical_n-of-1.stan",
  data = stan_data1,
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 11
)
fit0


# true group level data vs estimated group level data 
mu_group_true <- beta_0_mu +
  beta_1_mu * time +
  beta_2_mu * I_treat +
  beta_3_mu * t_post

post <- extract(fit0)

beta_0_mu_draws <- post$beta_0_mu
beta_1_mu_draws <- post$beta_1_mu
beta_2_mu_draws <- post$beta_2_mu
beta_3_mu_draws <- post$beta_3_mu

mu_group_draws <-
  outer(beta_0_mu_draws, rep(1, length(time))) +
  outer(beta_1_mu_draws, time) +
  outer(beta_2_mu_draws, I_treat) +
  outer(beta_3_mu_draws, t_post)

mu_group_ci <- apply(
  mu_group_draws,
  2,
  quantile,
  probs = c(0.025, 0.5, 0.975)
)

df_plot <- data.frame(
  time = time,
  true = mu_group_true,
  est  = mu_group_ci[2, ],
  lo   = mu_group_ci[1, ],
  hi   = mu_group_ci[3, ]
)

ggplot(df_plot, aes(x = time)) +
  
  # 95% credible interval (estimated)
  geom_ribbon(
    aes(ymin = lo, ymax = hi),
    fill = "steelblue",
    alpha = 0.25
  ) +
  
  # Estimated group-level mean (posterior median) — DASHED
  geom_line(
    aes(y = est),
    color = "steelblue4",
    linetype = "dashed",
    linewidth = 1.4
  ) +
  
  # True group-level mean — SOLID
  geom_line(
    aes(y = true),
    color = "gray20",
    linewidth = 1.2
  ) +
  
  # Intervention time
  geom_vline(
    xintercept = treat_time -1,
    color = "firebrick",
    linetype = "dotted",
    linewidth = 1
  ) +
  
  labs(
    title = "Group-Level Time Series",
    subtitle = "Dashed = estimated mean (posterior median); solid = true mean; shaded = 95% CI",
    x = "Time",
    y = "Outcome"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )






# model recovery n = (3,5,10), T1 = c(7,11,15) ----------------------------
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

R <- 100                      # number of simulation replications
N_grid <- c(3, 5, 10)
T_grid <- c(7, 11, 15)

beta_true <- c(
  beta_0_mu = 5,
  beta_1_mu = 0.1,
  beta_2_mu = 2,
  beta_3_mu = 0.2
)

phi   <- 0.3
sigma <- 0.5
treat_time <- 6


#define function
run_one_sim <- function(N, T1) {
  
  time <- 0:(T1 - 1)
  
  I_treat <- ifelse(seq_along(time) > treat_time, 1, 0)
  t_post  <- pmax(0, seq_along(time) - treat_time)
  
  # individual-level parameters
  beta_0 <- rnorm(N, beta_true["beta_0_mu"], 1)
  beta_1 <- rnorm(N, beta_true["beta_1_mu"], 0.5)
  beta_2 <- rnorm(N, beta_true["beta_2_mu"], 1)
  beta_3 <- rnorm(N, beta_true["beta_3_mu"], 0.5)
  
  y <- matrix(NA, N, T1)
  
  for (i in 1:N) {
    mu <- beta_0[i] +
      beta_1[i] * time +
      beta_2[i] * I_treat +
      beta_3[i] * t_post
    
    eps <- numeric(T1)
    eps[1] <- rnorm(1, 0, sigma / sqrt(1 - phi^2))
    
    for (t in 2:T1) {
      eps[t] <- phi * eps[t - 1] + rnorm(1, 0, sigma)
    }
    
    y[i, ] <- mu + eps
  }
  
  stan_data <- list(
    N = N,
    T = T1,
    y = y,
    time = time,
    I_treat = I_treat,
    t_post = t_post
  )
  
  fit <- stan(
    file = "model/hierarchical_n-of-1.stan",
    data = stan_data,
    chains = 2,
    iter = 2000,
    warmup = 1000,
    refresh = 0
  )
  
  post <- rstan::extract(fit)
  
  # posterior means
  est <- c(
    beta_0_mu = mean(post$beta_0_mu),
    beta_1_mu = mean(post$beta_1_mu),
    beta_2_mu = mean(post$beta_2_mu),
    beta_3_mu = mean(post$beta_3_mu)
  )
  
  # CI for coverage
  ci <- apply(
    cbind(
      post$beta_0_mu,
      post$beta_1_mu,
      post$beta_2_mu,
      post$beta_3_mu
    ),
    2,
    quantile,
    probs = c(0.025, 0.975)
  )
  
  list(est = est, ci = ci)
}

#run 
results <- foreach(
  N = N_grid,
  .combine = rbind,
  .packages = c("rstan")
) %:%
  foreach(
    T1 = T_grid,
    .combine = rbind
  ) %:%
  foreach(
    r = 1:R,
    .combine = rbind
  ) %dopar% {
    
    out <- run_one_sim(N, T1)
    
    data.frame(
      N = N,
      T = T1,
      param = names(beta_true),
      true = beta_true,
      est  = out$est,
      lo   = out$ci[1, ],
      hi   = out$ci[2, ]
    )
  }

#compute bias, RMSE and covergae
summary_stats <- results %>%
  group_by(N, T, param) %>%
  summarise(
    Bias = mean(est - true),
    RMSE = sqrt(mean((est - true)^2)),
    Coverage = mean(true >= lo & true <= hi),
    .groups = "drop"
  )

#stop cluster 
stopCluster(cl)




# new ---------------------------------------------------------------------

library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = 1)

stan_model_compiled <- stan_model("model/hierarchical_n-of-1.stan")

run_one_sim <- function(N, T1) {
  
  time <- 0:(T1 - 1)
  
  I_treat <- ifelse(time >= treat_time, 1, 0)
  t_post  <- pmax(0, time - treat_time)
  
  beta_0 <- rnorm(N, beta_true["beta_0_mu"], 1)
  beta_1 <- rnorm(N, beta_true["beta_1_mu"], 0.5)
  beta_2 <- rnorm(N, beta_true["beta_2_mu"], 1)
  beta_3 <- rnorm(N, beta_true["beta_3_mu"], 0.5)
  
  y <- matrix(NA, N, T1)
  
  for (i in 1:N) {
    mu <- beta_0[i] +
      beta_1[i] * time +
      beta_2[i] * I_treat +
      beta_3[i] * t_post
    
    eps <- numeric(T1)
    eps[1] <- rnorm(1, 0, sigma / sqrt(1 - phi^2))
    
    for (t in 2:T1) {
      eps[t] <- phi * eps[t - 1] + rnorm(1, 0, sigma)
    }
    
    y[i, ] <- mu + eps
  }
  
  stan_data <- list(
    N = N,
    T = T1,
    y = y,
    time = time,
    I_treat = I_treat,
    t_post = t_post
  )
  
  fit <- sampling(
    stan_model_compiled,
    data = stan_data,
    chains = 1,
    iter = 2000,
    warmup = 1000,
    refresh = 0
  )
  
  post <- rstan::extract(fit)
  
  est <- c(
    beta_0_mu = mean(post$beta_0_mu),
    beta_1_mu = mean(post$beta_1_mu),
    beta_2_mu = mean(post$beta_2_mu),
    beta_3_mu = mean(post$beta_3_mu)
  )
  
  ci <- apply(
    cbind(
      post$beta_0_mu,
      post$beta_1_mu,
      post$beta_2_mu,
      post$beta_3_mu
    ),
    2,
    quantile,
    probs = c(0.025, 0.975)
  )
  
  list(est = est, ci = ci)
}


library(doParallel)
library(foreach)

n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

results <- foreach(
  N = N_grid,
  .combine = rbind,
  .packages = c("rstan")
) %:%
  foreach(
    T1 = T_grid,
    .combine = rbind
  ) %:%
  foreach(
    r = 1:R,
    .combine = rbind
  ) %dopar% {
    
    out <- run_one_sim(N, T1)
    
    data.frame(
      N = N,
      T = T1,
      param = names(beta_true),
      true = beta_true,
      est  = out$est,
      lo   = out$ci[1, ],
      hi   = out$ci[2, ]
    )
  }

stopCluster(cl)


# cmdstan RMSE and coverage---------------------------------------------------------------------

library(cmdstanr)
library(posterior)
library(doParallel)
library(foreach)
library(dplyr)

set.seed(123)

mod <- cmdstan_model("model/hierarchical_n-of-1.stan")

run_one_sim <- function(N, T1) {
  
  time <- 0:(T1 - 1)
  
  I_treat <- ifelse(time >= treat_time, 1, 0)
  t_post  <- pmax(0, time - treat_time)
  
  beta_0 <- rnorm(N, beta_true["beta_0_mu"], 1)
  beta_1 <- rnorm(N, beta_true["beta_1_mu"], 0.5)
  beta_2 <- rnorm(N, beta_true["beta_2_mu"], 1)
  beta_3 <- rnorm(N, beta_true["beta_3_mu"], 0.5)
  
  y <- matrix(NA, N, T1)
  
  for (i in 1:N) {
    mu <- beta_0[i] +
      beta_1[i] * time +
      beta_2[i] * I_treat +
      beta_3[i] * t_post
    
    eps <- numeric(T1)
    eps[1] <- rnorm(1, 0, sigma / sqrt(1 - phi^2))
    
    for (t in 2:T1) {
      eps[t] <- phi * eps[t - 1] + rnorm(1, 0, sigma)
    }
    
    y[i, ] <- mu + eps
  }
  
  stan_data <- list(
    N = N,
    T = T1,
    y = y,
    time = time,
    I_treat = I_treat,
    t_post = t_post
  )
  
  fit <- mod$sample(
    data = stan_data,
    chains = 1,
    iter_sampling = 1000,
    iter_warmup = 1000,
    refresh = 0
  )
  
  draws <- posterior::as_draws_df(fit$draws())
  
  est <- c(
    beta_0_mu = mean(draws$beta_0_mu),
    beta_1_mu = mean(draws$beta_1_mu),
    beta_2_mu = mean(draws$beta_2_mu),
    beta_3_mu = mean(draws$beta_3_mu)
  )
  
  ci <- apply(
    cbind(
      draws$beta_0_mu,
      draws$beta_1_mu,
      draws$beta_2_mu,
      draws$beta_3_mu
    ),
    2,
    quantile,
    probs = c(0.025, 0.975)
  )
  
  ci <- apply(
    cbind(
      draws$beta_0_mu,
      draws$beta_1_mu,
      draws$beta_2_mu,
      draws$beta_3_mu
    ),
    2,
    quantile,
    probs = c(0.025, 0.975)
  )
  
  list(est = est, ci = ci)
}

R <- 100
N_grid <- c(3, 5, 10)
T_grid <- c(7, 11, 15)

beta_true <- c(
  beta_0_mu = 5,
  beta_1_mu = 0.1,
  beta_2_mu = 2,
  beta_3_mu = 0.2
)

phi   <- 0.3
sigma <- 0.5
treat_time <- 6

n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

results <- foreach(
  N = N_grid,
  .combine = rbind
) %:%
  foreach(
    T1 = T_grid,
    .combine = rbind
  ) %:%
  foreach(
    r = 1:R,
    .combine = rbind
  ) %dopar% {
    
    out <- run_one_sim(N, T1)
    
    data.frame(
      N = N,
      T = T1,
      param = names(beta_true),
      true = beta_true,
      est  = out$est,
      lo   = out$ci[1, ],
      hi   = out$ci[2, ]
    )
  }



#compute bias, RMSE and covergae
summary_stats <- results %>%
  group_by(N, T, param) %>%
  summarise(
    Bias = mean(est - true),
    RMSE = sqrt(mean((est - true)^2)),
    Coverage = mean(true >= lo & true <= hi),
    .groups = "drop"
  )

write.csv(summary_stats, "summary_stats.csv", row.names = FALSE)
saveRDS(summary_stats, "summary_stats.rds")


summary_stats_RMSE <- summary_stats %>%
  group_by(param) %>%      # group by parameter
  arrange(RMSE, .by_group = TRUE) %>%  # sort by RMSE within each parameter
  ungroup()                # remove grouping

View(summary_stats_RMSE)

write.csv(summary_stats_RMSE, "summary_stats_RMSE.csv", row.names = FALSE)
saveRDS(summary_stats_RMSE, "summary_stats_RMSE.rds")


summary_stats_CI <- summary_stats %>%
  group_by(param) %>%      # group by parameter
  arrange(Coverage, .by_group = TRUE) %>%  # sort by RMSE within each parameter
  ungroup()                # remove grouping

View(summary_stats_CI)

write.csv(summary_stats_CI, "summary_stats_CI.csv", row.names = FALSE)
saveRDS(summary_stats_CI, "summary_stats_CI.rds")

RMSEstopCluster(cl)
stopCluster(cl)







# final plz ---------------------------------------------------------------

set.seed(123)
library(tidyverse)


# Design
N  <- 10
T  <- 10
time <- 1:T
treat_time <- 6

I_treat <- ifelse(time >= treat_time, 1, 0)
t_post  <- pmax(0, time - treat_time)

# Population-level parameters

beta_mu <- c(0, 0.1, 2, 0.2)
beta_sd <- c(1, 0.5, 1, 0.5)

phi_mu <- 0.3
phi_sd <- 0.2

sigma <- 0.4

# Individual-level parameters
beta_0 <- rnorm(N, beta_mu[1], beta_sd[1])
beta_1 <- rnorm(N, beta_mu[2], beta_sd[2])
beta_2 <- rnorm(N, beta_mu[3], beta_sd[3])
beta_3 <- rnorm(N, beta_mu[4], beta_sd[4])

phi <- tanh(rnorm(N, phi_mu, phi_sd))  # constrain to (-1,1)

# Generate data

y <- matrix(NA, N, T)
mu <- matrix(NA, N, T)

for (i in 1:N) {
  
  mu[i, ] <-
    beta_0[i] +
    beta_1[i] * time +
    beta_2[i] * I_treat +
    beta_3[i] * t_post
  
  # stationary initialization
  y[i, 1] <- rnorm(
    1,
    mean = mu[i, 1] / (1 - phi[i]),
    sd   = sigma / sqrt(1 - phi[i]^2)
  )
  
  for (t in 2:T) {
    y[i, t] <- mu[i, t] + phi[i] * y[i, t - 1] + rnorm(1, 0, sigma)
  }
}

# long format
df <- expand_grid(
  id = factor(1:N),
  time = time
) |>
  mutate(y = as.vector(t(y)))


mu_pop <-
  beta_mu[1] +
  beta_mu[2] * time +
  beta_mu[3] * I_treat +
  beta_mu[4] * t_post

y_pop <- numeric(T)
y_pop[1] <- mu_pop[1] / (1 - phi_mu)

for (t in 2:T) {
  y_pop[t] <- mu_pop[t] + phi_mu * y_pop[t - 1]
}

df_pop <- tibble(time = time, y_pop = y_pop)



ggplot() +
  # individual trajectories (thin)
  geom_line(
    data = df,
    aes(x = time, y = y, group = id),
    alpha = 0.4,
    linewidth = 0.4
  ) +
  
  # population trajectory (thick)
  geom_line(
    data = df_pop,
    aes(x = time, y = y_pop),
    linewidth = 1.5
  ) +
  
  geom_vline(
    xintercept = treat_time,
    linetype = "dashed"
  ) +
  
  labs(
    title = "Hierarchical Dynamic Regression",
    subtitle = "Thin lines = individuals, thick line = population",
    x = "Time",
    y = "y"
  ) +
  
  theme_minimal()



library(rstan)

y_mat <- matrix(df$y, nrow = Tt, ncol = N)

stan_data <- list(
  N = N,
  T = Tt,
  y = y_mat
)


library(rstan)

# reshape y: [N, T]
y_mat <- matrix(df$y, nrow = N, ncol = T, byrow = TRUE)

stan_data <- list(
  N = N,
  T = T,
  y = y_mat,
  time = time,
  I_treat = I_treat,
  t_post = t_post
)

library(rstan)
fit <- stan(
  file = "model/hierarchical_n-of-1.stan",
  data = stan_data,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4
)

fit

print(
  fit,
  pars = c("beta_mu", "beta_sd", "phi_mu", "phi_sd", "sigma")
)







# simulation1 fixed parameters and idfferent n and T -------------------------------------------------------------
library(cmdstanr)
library(posterior)
library(doParallel)
library(foreach)
library(dplyr)

stopCluster(cl)

cl <- makeCluster(n_cores)
registerDoParallel(cl)


set.seed(123)

mod <- cmdstan_model("model/hierarchical_n-of-1.stan")

#parameters
beta_mu <- c(
  beta_0 = 0,
  beta_1 = 0.1,
  beta_2 = 2,
  beta_3 = 0.2
)

beta_sd <- c(
  beta_0 = 1,
  beta_1 = 0.5,
  beta_2 = 1,
  beta_3 = 0.5
)

phi_mu <- 0.3
phi_sd <- 0.05
sigma <- 0.4

# single simulation 
run_one_sim <- function(N, T) {
  
  time <- 1:T
  treat_time <- ceiling(T / 2)
  
  I_treat <- ifelse(time >= treat_time, 1, 0)
  t_post  <- pmax(0, time - treat_time)
  
  beta_0 <- rnorm(N, beta_mu["beta_0"], beta_sd["beta_0"])
  beta_1 <- rnorm(N, beta_mu["beta_1"], beta_sd["beta_1"])
  beta_2 <- rnorm(N, beta_mu["beta_2"], beta_sd["beta_2"])
  beta_3 <- rnorm(N, beta_mu["beta_3"], beta_sd["beta_3"])
  
  phi <- rnorm(N, phi_mu, phi_sd)
  
  y <- matrix(NA, N, T)
  
  for (i in 1:N) {
    mu <- beta_0[i] +
      beta_1[i] * time +
      beta_2[i] * I_treat +
      beta_3[i] * t_post
    
    y[i, 1] <- rnorm(
      1,
      mean = mu[1] / (1 - phi[i]),
      sd   = sigma / sqrt(1 - phi[i]^2)
    )
    
    for (t in 2:T) {
      y[i, t] <- mu[t] + phi[i] * y[i, t - 1] + rnorm(1, 0, sigma)
    }
  }
  
  stan_data <- list(
    N = N,
    T = T,
    y = y,
    time = time,
    I_treat = I_treat,
    t_post = t_post
  )
  
  fit <- mod$sample(
    data = stan_data,
    chains = 1,
    iter_warmup = 1000,
    iter_sampling = 1000,
    refresh = 0
  )
  
  draws <- as_draws_df(fit$draws())
  
  param_names <- c(
    "beta_0_mu", "beta_1_mu", "beta_2_mu", "beta_3_mu", "phi_mu"
  )
  
  est <- draws |>
    summarise(across(all_of(param_names), mean))
  
  ci <- draws |>
    summarise(
      across(
        all_of(param_names),
        list(
          lo = ~quantile(.x, 0.025),
          hi = ~quantile(.x, 0.975)
        )
      )
    )
  
  list(est = est, ci = ci)
}


#simulation setting 
R <- 100

N_grid <- c(3, 5, 10)
T_grid <- c(8, 14, 20)

#foreach 
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

#main loop 
results <- foreach(
  N = N_grid,
  .combine = rbind,
  .packages = c("cmdstanr", "posterior", "dplyr")
) %:%
  foreach(T = T_grid, .combine = rbind) %:%
  foreach(r = 1:R, .combine = rbind) %dopar% {
    
    out <- run_one_sim(N, T)
    
    param_names <- c(
      "beta_0_mu", "beta_1_mu", "beta_2_mu", "beta_3_mu", "phi_mu"
    )
    
    data.frame(
      N = N,
      T = T,
      param = param_names,
      true = c(beta_mu, phi_mu),
      est  = as.numeric(out$est[1, param_names]),
      lo   = as.numeric(out$ci[1, paste0(param_names, "_lo")]),
      hi   = as.numeric(out$ci[1, paste0(param_names, "_hi")])
    )
  }



#summary 
summary_stats <- results |>
  group_by(N, T, param) |>
  summarise(
    Bias = mean(est - true),
    RMSE = sqrt(mean((est - true)^2)),
    Coverage = mean(lo <= true & true <= hi),
    .groups = "drop"
  )

summary_stats

summary_stats_RMSE <- summary_stats %>%
  group_by(param) %>%      # group by parameter
  arrange(RMSE, .by_group = TRUE) %>%  # sort by RMSE within each parameter
  ungroup()                # remove grouping

View(summary_stats_RMSE)

View(summary_stats)

write.csv(
  summary_stats,
  file = "results/fina_bias_RMSE.csv",
  row.names = FALSE
)

write.csv(
  summary_stats_RMSE,
  file = "results/final_bias_RMSE_byParameters.csv",
  row.names = FALSE
)

saveRDS(
  summary_stats,
  file = "results_bias_RMSE.rds"
)

stopCluster(cl)



# verying parameters 
# 1. Update the function definition to handle the randomized inputs
run_one_sim_recovery <- function(N, T, true_beta_mu, true_phi_mu) {
  
  time <- 1:T
  treat_time <- ceiling(T / 2)
  I_treat <- ifelse(time >= treat_time, 1, 0)
  t_post  <- pmax(0, time - treat_time)
  
  # Individual parameters drawn from the RANDOMIZED population means
  # Note: beta_sd, phi_sd, and sigma need to be defined or passed in
  beta_0 <- rnorm(N, true_beta_mu["beta_0"], 1)   # Using sd=1
  beta_1 <- rnorm(N, true_beta_mu["beta_1"], 0.5) 
  beta_2 <- rnorm(N, true_beta_mu["beta_2"], 1)
  beta_3 <- rnorm(N, true_beta_mu["beta_3"], 0.5)
  phi    <- rnorm(N, true_phi_mu, 0.05)
  sigma  <- 0.4
  
  y <- matrix(NA, N, T)
  
  for (i in 1:N) {
    mu <- beta_0[i] + beta_1[i]*time + beta_2[i]*I_treat + beta_3[i]*t_post
    
    y[i, 1] <- rnorm(1, mu[1] / (1 - phi[i]), sigma / sqrt(1 - phi[i]^2))
    
    for (t in 2:T) {
      y[i, t] <- mu[t] + phi[i] * y[i, t-1] + rnorm(1, 0, sigma)
    }
  }
  
  stan_data <- list(N=N, T=T, y=y, time=time, I_treat=I_treat, t_post=t_post)
  
  # Run model
  fit <- mod$sample(data = stan_data, chains = 1, iter_warmup = 500, 
                    iter_sampling = 500, refresh = 0, show_messages = FALSE)
  
  # Extract population means
  draws <- fit$summary(c("beta_0_mu", "beta_1_mu", "beta_2_mu", "beta_3_mu", "phi_mu"))
  return(draws$mean)
}

# 2. Export the new function and objects to the cluster
parallel::clusterExport(cl, c("run_one_sim_recovery", "mod"))

# 3. Run the loop
results_recovery <- foreach(
  N = N_grid,
  .combine = rbind,
  .packages = c("cmdstanr", "posterior", "dplyr")
) %:%
  foreach(T = T_grid, .combine = rbind) %:%
  foreach(r = 1:R, .combine = rbind) %dopar% {
    
    # Generate random truths for this iteration
    current_true_beta_mu <- c(beta_0 = runif(1, -1, 1), beta_1 = runif(1, 0, 0.5),
                              beta_2 = runif(1, 1, 3), beta_3 = runif(1, 0, 0.5))
    current_true_phi_mu <- runif(1, 0.1, 0.7)
    
    # CALL THE UPDATED FUNCTION
    est_mu <- run_one_sim_recovery(N, T, current_true_beta_mu, current_true_phi_mu)
    
    data.frame(
      N = N, T = T, Iteration = r,
      param = c("beta_0_mu", "beta_1_mu", "beta_2_mu", "beta_3_mu", "phi_mu"),
      true_val = c(current_true_beta_mu, current_true_phi_mu),
      est_val = est_mu
    )
  }

results_recovery
View(results_recovery)

library(dplyr)

recovery_correlations <- results_recovery %>%
  group_by(param, N, T) %>%
  summarise(
    # Pearson correlation between truth and estimation
    recovery_rate = cor(true_val, est_val, use = "complete.obs"),
    .groups = "drop"
  )

# View the results
print(recovery_correlations)
View(recovery_correlations)

write.csv(
  recovery_correlations,
  file = "results/parameter_recovery.csv",
  row.names = FALSE
)
# varying parameters: MC bias, RMSE, coverage and correlation  ------------

library(cmdstanr)
library(posterior)
library(doParallel)
library(foreach)
library(dplyr)


cl <- makeCluster(n_cores)
registerDoParallel(cl)


set.seed(123)

mod <- cmdstan_model("model/hierarchical_n-of-1.stan")


run_one_sim <- function(N, T) {
  
  # ---- time structure ----
  time <- 1:T
  treat_time <- ceiling(T / 2)
  I_treat <- ifelse(time >= treat_time, 1, 0)
  t_post  <- pmax(0, time - treat_time)
  
  # ---- draw TRUE population-level parameters ----
  true <- list(
    beta_0_mu = rnorm(1, 0, 1),
    beta_1_mu = rnorm(1, 0.1, 0.5),
    beta_2_mu = rnorm(1, 2, 1),
    beta_3_mu = rnorm(1, 0.2, 0.5),
    phi_mu    = rnorm(1, 0.3, 0.05)
  )
  
  beta_sd <- c(1, 0.5, 1, 0.5)
  sigma   <- 0.4
  
  # ---- individual-level parameters ----
  beta_0 <- rnorm(N, true$beta_0_mu, beta_sd[1])
  beta_1 <- rnorm(N, true$beta_1_mu, beta_sd[2])
  beta_2 <- rnorm(N, true$beta_2_mu, beta_sd[3])
  beta_3 <- rnorm(N, true$beta_3_mu, beta_sd[4])
  
  phi <- rnorm(N, true$phi_mu, 0.05)
  
  # ---- simulate data ----
  y <- matrix(NA, N, T)
  
  for (i in 1:N) {
    mu <- beta_0[i] +
      beta_1[i] * time +
      beta_2[i] * I_treat +
      beta_3[i] * t_post
    
    y[i, 1] <- rnorm(
      1,
      mu[1] / (1 - phi[i]),
      sigma / sqrt(1 - phi[i]^2)
    )
    
    for (t in 2:T) {
      y[i, t] <- mu[t] + phi[i] * y[i, t - 1] + rnorm(1, 0, sigma)
    }
  }
  
  # ---- fit Stan ----
  fit <- mod$sample(
    data = list(
      N = N,
      T = T,
      y = y,
      time = time,
      I_treat = I_treat,
      t_post = t_post
    ),
    chains = 1,
    iter_warmup = 1000,
    iter_sampling = 1000,
    refresh = 0
  )
  
  draws <- posterior::as_draws_df(fit$draws())
  
  params <- names(true)
  
  est <- draws |>
    dplyr::summarise(dplyr::across(all_of(params), mean))
  
  ci <- draws |>
    dplyr::summarise(
      dplyr::across(
        all_of(params),
        list(lo = ~quantile(.x, 0.025),
             hi = ~quantile(.x, 0.975))
      )
    )
  
  data.frame(
    param = params,
    true  = unlist(true),
    est   = as.numeric(est[1, params]),
    lo    = as.numeric(ci[1, paste0(params, "_lo")]),
    hi    = as.numeric(ci[1, paste0(params, "_hi")])
  )
}



R <- 100
N_grid <- c(3, 5, 10)
T_grid <- c(8, 14, 20)

results <- foreach(
  N = N_grid,
  .combine = rbind,
  .packages = c("cmdstanr", "posterior", "dplyr")
) %:%
  foreach(T = T_grid, .combine = rbind) %dopar% {
    
    # Run the simulation R times
    sim_res <- lapply(1:R, function(r) {
      run_one_sim(N, T) |> mutate(iter = r) # Track which iteration it was
    })
    
    # Return the raw data frame instead of summarizing it
    bind_rows(sim_res) |> mutate(N = N, T = T)

    
    # ---- Bias / RMSE / Coverage ----
    summary_stats <- sim_res |>
      group_by(param) |>
      summarise(
        Bias = mean(est - true),
        RMSE = sqrt(mean((est - true)^2)),
        Coverage = mean(lo <= true & true <= hi),
        .groups = "drop"
      )
    
    # ---- Parameter recovery (correlation) ----
    recovery <- sim_res |>
      group_by(param) |>
      summarise(
        Recovery = cor(true, est),
        .groups = "drop"
      )
    
    left_join(summary_stats, recovery, by = "param") |>
      mutate(N = N, T = T)
  }

summary_stats
stopCluster


# bias plot 
plot_data <- data.frame(
  param = rep(c("beta_0_mu", "beta_1_mu", "beta_2_mu", "beta_3_mu", "phi_mu"), 2),
  Scenario = rep(c("Worst Case (N=3, T=8)", "Best Case (N=10, T=20)"), each = 5),
  Bias = c(-0.029, 0.026, 0.096, 0.035, -0.147,  # Worst Case
           -0.033, -0.019, -0.064, 0.001, -0.038) # Best Case
)

ggplot(plot_data, aes(x = param, y = Bias, fill = Scenario)) +
  geom_bar(stat = "identity", position = "dodge", color = "white") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.5) +
  geom_text(aes(label = round(Bias, 3)), 
            position = position_dodge(width = 0.9), 
            vjust = ifelse(plot_data$Bias >= 0, -0.5, 1.5), 
            size = 3.5) +
  scale_fill_manual(values = c("Best Case (N=10, T=20)" = "#00BFC4", 
                               "Worst Case (N=3, T=8)" = "#F8766D")) +
  labs(
    title = "Bias Comparison",
    x = "Population Parameter",
    y = "Bias (Estimate - True)",
    fill = "Simulation Scenario"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# correlation plot 
ggplot(results, aes(x = true, y = est)) +
  geom_point(alpha = 0.4, color = "midnightblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  facet_wrap(~param, scales = "free") +
  labs(
    title = "Parameter Recovery: True vs. Estimated",
    subtitle = "Red dashed line indicates perfect recovery",
    x = "True Value (Simulated)",
    y = "Posterior Mean (Estimated)"
  ) +
  theme_minimal()

# bias of phi plot 
phi_results <- results %>%
  filter(param == "phi_mu") %>%
  mutate(bias = est - true)

ggplot(phi_results, aes(x = true, y = bias, color = as.factor(T))) +
  geom_smooth(method = "loess", se = FALSE) + # This creates the "curves"
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Bias of Phi",
    x = "True Phi Value",
    y = "Bias (Estimated - True)",
    color = "Time Points (T)"
  ) +
  theme_minimal()


# Load necessary libraries
library(ggplot2)
library(dplyr)

# 1. Create the dataset
df <- data.frame(
  param = rep(c("beta_0", "beta_1", "beta_2", "beta_3", "phi"), each = 9),
  N = rep(c(3, 3, 3, 5, 5, 5, 10, 10, 10), 5),
  T_val = rep(c(8, 14, 20), 15),
  recovery_r = c(
    0.644548, 0.703339, 0.686669, 0.747429, 0.796795, 0.713623, 0.865372, 0.855740, 0.877437,
    0.403611, 0.373734, 0.464697, 0.390244, 0.583387, 0.495871, 0.584623, 0.663579, 0.610352,
    0.705306, 0.629383, 0.676521, 0.737269, 0.767881, 0.814471, 0.841320, 0.874615, 0.889190,
    0.092331, 0.292757, 0.420269, 0.347843, 0.414120, 0.456205, 0.567879, 0.583562, 0.680659,
    0.678714, 0.792473, 0.796886, 0.786716, 0.911130, 0.921408, 0.889768, 0.935992, 0.961770
  )
)

# 2. Create the plot
ggplot(df, aes(x = factor(T_val), y = recovery_r, color = factor(N), group = N)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  facet_wrap(~param, scales = "free_y") +
  scale_color_viridis_d(name = "Sample Size (N)") +
  labs(
    title = "Parameter Recovery Rate",
    x = "Time Points (T)",
    y = " Correlation "
  ) +
  theme_minimal(base_size = 14) + # Larger text for slides
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1, "lines")
  )









# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)
# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Create the dataset
df <- data.frame(
  param = rep(c("beta_0_mu", "beta_1_mu", "beta_2_mu", "beta_3_mu", "phi_mu"), each = 9),
  recovery_rate = c(
    0.644548, 0.703339, 0.686669, 0.747429, 0.796795, 0.713623, 0.865372, 0.855740, 0.877437,
    0.403611, 0.373734, 0.464697, 0.390244, 0.583387, 0.495871, 0.584623, 0.663579, 0.610352,
    0.705306, 0.629383, 0.676521, 0.737269, 0.767881, 0.814471, 0.841320, 0.874615, 0.889190,
    0.092331, 0.292757, 0.420269, 0.347843, 0.414120, 0.456205, 0.567879, 0.583562, 0.680659,
    0.678714, 0.792473, 0.796886, 0.786716, 0.911130, 0.921408, 0.889768, 0.935992, 0.961770
  )
)

# 2. Extract Min and Max per parameter
summary_df <- df %>%
  group_by(param) %>%
  summarise(
    Worst = min(recovery_rate),
    Best = max(recovery_rate)
  ) %>%
  pivot_longer(cols = c(Worst, Best), names_to = "Scenario", values_to = "Value")

# 3. Create the Bar Plot
ggplot(summary_df, aes(x = param, y = Value, fill = Scenario)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("Worst" = "#e41a1c", "Best" = "#4daf4a")) +
  labs(
    title = "Worst vs Best Recovery Performance",
    subtitle = "Comparing simulation extremes for each parameter",
    x = "Parameter",
    y = "Recovery Rate (r)"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 1))




# varying parameters: MC bias, RMSE, coverage and correlation  ------------

library(cmdstanr)
library(posterior)
library(doParallel)
library(foreach)
library(dplyr)


cl <- makeCluster(n_cores)
registerDoParallel(cl)


set.seed(123)

mod <- cmdstan_model("model/hierarchical_n-of-1.stan")


run_one_sim <- function(N, T) {
  
  # ---- time structure ----
  time <- 1:T
  treat_time <- ceiling(T / 2)
  I_treat <- ifelse(time >= treat_time, 1, 0)
  t_post  <- pmax(0, time - treat_time)
  
  # ---- draw TRUE population-level parameters ----
  true <- list(
    beta_0_mu = rnorm(1, 0, 1),
    beta_1_mu = rnorm(1, 0.1, 0.5),
    beta_2_mu = rnorm(1, 2, 1),
    beta_3_mu = rnorm(1, 0.2, 0.5),
    phi_mu    = rnorm(1, 0.3, 0.05)
  )
  
  beta_sd <- c(1, 0.5, 1, 0.5)
  sigma   <- 0.4
  
  # ---- individual-level parameters ----
  beta_0 <- rnorm(N, true$beta_0_mu, beta_sd[1])
  beta_1 <- rnorm(N, true$beta_1_mu, beta_sd[2])
  beta_2 <- rnorm(N, true$beta_2_mu, beta_sd[3])
  beta_3 <- rnorm(N, true$beta_3_mu, beta_sd[4])
  
  phi <- rnorm(N, true$phi_mu, 0.05)
  
  # ---- simulate data ----
  y <- matrix(NA, N, T)
  
  for (i in 1:N) {
    mu <- beta_0[i] +
      beta_1[i] * time +
      beta_2[i] * I_treat +
      beta_3[i] * t_post
    
    y[i, 1] <- rnorm(
      1,
      mu[1] / (1 - phi[i]),
      sigma / sqrt(1 - phi[i]^2)
    )
    
    for (t in 2:T) {
      y[i, t] <- mu[t] + phi[i] * y[i, t - 1] + rnorm(1, 0, sigma)
    }
  }
  
  # ---- fit Stan ----
  fit <- mod$sample(
    data = list(
      N = N,
      T = T,
      y = y,
      time = time,
      I_treat = I_treat,
      t_post = t_post
    ),
    chains = 1,
    iter_warmup = 1000,
    iter_sampling = 1000,
    refresh = 0
  )
  
  draws <- posterior::as_draws_df(fit$draws())
  
  params <- names(true)
  
  est <- draws |>
    dplyr::summarise(dplyr::across(all_of(params), mean))
  
  ci <- draws |>
    dplyr::summarise(
      dplyr::across(
        all_of(params),
        list(lo = ~quantile(.x, 0.025),
             hi = ~quantile(.x, 0.975))
      )
    )
  
  data.frame(
    param = params,
    true  = unlist(true),
    est   = as.numeric(est[1, params]),
    lo    = as.numeric(ci[1, paste0(params, "_lo")]),
    hi    = as.numeric(ci[1, paste0(params, "_hi")])
  )
}



R <- 100
N_grid <- c(3, 5, 10)
T_grid <- c(8, 14, 20)

results <- foreach(
  N = N_grid,
  .combine = rbind,
  .packages = c("cmdstanr", "posterior", "dplyr")
) %:%
  foreach(T = T_grid, .combine = rbind) %dopar% {
    
    # Run the simulation R times
    sim_res <- lapply(1:R, function(r) {
      run_one_sim(N, T) |> mutate(iter = r) # Track which iteration it was
    })
    
    # Return the raw data frame instead of summarizing it
    bind_rows(sim_res) |> mutate(N = N, T = T)
    
    
    # ---- Bias / RMSE / Coverage ----
    summary_stats <- sim_res |>
      group_by(param) |>
      summarise(
        Bias = mean(est - true),
        RMSE = sqrt(mean((est - true)^2)),
        Coverage = mean(lo <= true & true <= hi),
        .groups = "drop"
      )
    
    # ---- Parameter recovery (correlation) ----
    recovery <- sim_res |>
      group_by(param) |>
      summarise(
        Recovery = cor(true, est),
        .groups = "drop"
      )
    
    left_join(summary_stats, recovery, by = "param") |>
      mutate(N = N, T = T)
  }

summary_stats
stopCluster







# 2-2 ---------------------------------------------------------------------


library(cmdstanr)
library(posterior)
library(doParallel)
library(foreach)
library(dplyr)


cl <- makeCluster(n_cores)
registerDoParallel(cl)


set.seed(123)

mod <- cmdstan_model("model/hierarchical_n-of-1.stan")


run_one_sim <- function(N, T) {
  
  # ---- time structure ----
  time <- 1:T
  treat_time <- ceiling(T / 2)
  I_treat <- ifelse(time >= treat_time, 1, 0)
  t_post  <- pmax(0, time - treat_time)
  
  # ---- draw TRUE population-level parameters ----
  true <- list(
    beta_0_mu = runif(1, -1, 1),
    beta_1_mu = runif(1, 0, 0.5),
    beta_2_mu = runif(1, 1, 3),
    beta_3_mu = runif(1, 0, 0.5),
    phi_mu    = runif(1, 0.1, 0.7)
  )
  
  beta_sd <- c(1, 0.5, 1, 0.5)
  sigma   <- 0.4
  
  # ---- individual-level parameters ----
  beta_0 <- rnorm(N, true$beta_0_mu, beta_sd[1])
  beta_1 <- rnorm(N, true$beta_1_mu, beta_sd[2])
  beta_2 <- rnorm(N, true$beta_2_mu, beta_sd[3])
  beta_3 <- rnorm(N, true$beta_3_mu, beta_sd[4])
  phi <- pmax(pmin(rnorm(N, true$phi_mu, 0.05), 0.95), -0.95)

  # ---- simulate data ----
  y <- matrix(NA, N, T)
  
  for (i in 1:N) {
    mu <- beta_0[i] +
      beta_1[i] * time +
      beta_2[i] * I_treat +
      beta_3[i] * t_post
    
    y[i, 1] <- rnorm(
      1,
      mu[1] / (1 - phi[i]),
      sigma / sqrt(1 - phi[i]^2)
    )
    
    for (t in 2:T) {
      y[i, t] <- mu[t] + phi[i] * y[i, t - 1] + rnorm(1, 0, sigma)
    }
  }
  
  # ---- fit Stan ----
  fit <- mod$sample(
    data = list(N=N, T=T, y=y, time=time, I_treat=I_treat, t_post=t_post),
    chains = 1, iter_warmup = 300, iter_sampling = 300, # Dropped slightly for speed
    max_treedepth = 8, refresh = 0, show_messages = FALSE
  )
  
  draws <- posterior::as_draws_df(fit$draws())
  
  params <- names(true)
  
  est <- draws |>
    dplyr::summarise(dplyr::across(all_of(params), mean))
  
  ci <- draws |>
    dplyr::summarise(
      dplyr::across(
        all_of(params),
        list(lo = ~quantile(.x, 0.025),
             hi = ~quantile(.x, 0.975))
      )
    )
  
  data.frame(
    param = params,
    true  = unlist(true),
    est   = as.numeric(est[1, params]),
    lo    = as.numeric(ci[1, paste0(params, "_lo")]),
    hi    = as.numeric(ci[1, paste0(params, "_hi")])
  )
}



R <- 50
N_grid <- c(3, 5, 10)
T_grid <- c(8, 14, 20)



# We only return the RAW data frames here so we can use them for plots later.
results_raw <- foreach(
  N = N_grid,
  .combine = rbind,
  .packages = c("cmdstanr", "posterior", "dplyr")
) %:%
  foreach(T = T_grid, .combine = rbind) %dopar% {
    
    # Run the simulation R times
    sim_list <- lapply(1:R, function(r) {
      run_one_sim(N, T) |> mutate(iter = r)
    })
    
    # Combine the list into one data frame and return it
    # This is the "raw" data with 100 iterations per scenario
    bind_rows(sim_list) |> mutate(N = N, T = T)
  }
# Save the object
write.csv(results_raw, "simulation_results_raw.csv", row.names = FALSE)
saveRDS(results_raw, "simulation_results_raw.rds")

library(dplyr)
results_raw <- readRDS("simulation_results_raw.rds")
# To load it back later, you would use:
# results_raw <- readRDS("simulation_results_raw.rds")
# ---- 1. Calculate the Summary Table (Bias, RMSE, Correlation) ----
# We do this calculation outside the loop now
summary_stats_table <- results_raw %>%
  group_by(param, N, T) %>%
  summarise(
    Bias = mean(est - true),
    RMSE = sqrt(mean((est - true)^2)),
    Coverage = mean(lo <= true & true <= hi),
    Recovery = cor(true, est),
    .groups = "drop"
  )

print(summary_stats_table)
View(summary_stats_table)
write.csv(summary_stats_table, "results/recovery.csv", row.names = FALSE)


# parameter recovery plot 
library(ggplot2)

# Ensure N and T are treated as factors for clean plotting
plot_data <- summary_stats_table %>%
  mutate(
    N_label = paste0("N = ", N),
    T_factor = as.factor(T)
  )

ggplot(plot_data, aes(x = T_factor, y = Recovery, group = as.factor(N), color = as.factor(N))) +
  # Lines and points to match your screenshot style
  geom_line(size = 1) +
  geom_point(size = 3) +
  # Separate plot for each parameter
  facet_wrap(~param, scales = "free_y") +
  # Professional styling
  theme_minimal(base_size = 14) +
  scale_color_viridis_d(option = "viridis", end = 0.9) + # Matches the yellow/green/purple theme
  labs(
    title = "Parameter Recovery Rate (Correlation)",
    x = "Time Points (T)",
    y = "Correlation (Recovery Rate)",
    color = "Sample Size (N)"
  ) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )
# ---- 2. Parameter Recovery Plot (Correlation) ----
# Filter for the specific scenario
cor_labels <- data.frame(
  param = c("beta_0_mu", "beta_1_mu", "beta_2_mu", "beta_3_mu", "phi_mu"),
  label = c("r = 0.87", "r = 0.71", "r = 0.86", "r = 0.68", "r = 0.95"),
  # These coordinates place the text in the top-left corner of each facet
  # Adjust x and y based on your parameter ranges if needed
  x = -Inf, 
  y = Inf
)

# 2. Generate the plot
results_raw %>%
  filter(T == 14, N == 10) %>%
  ggplot(aes(x = true, y = est)) +
  geom_point(alpha = 0.5, color = "midnightblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  # Add the correlation labels
  geom_text(data = cor_labels, aes(x = x, y = y, label = label), 
            hjust = -0.2, vjust = 1.5, fontface = "italic", size = 4.5, inherit.aes = FALSE) +
  facet_wrap(~param, scales = "free") +
  labs(
    title = "Parameter Recovery (T = 14, N = 10)",
    x = "True Value", 
    y = "Estimated Value"
  ) +
  theme_minimal()

ggplot(results_raw %>% filter(param == "phi_mu"), aes(x = true, y = est)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_point(alpha = 0.4, color = "midnightblue") +
  # This splits the plot: Rows are N, Columns are T
  facet_grid(N ~ T, labeller = label_both) + 
  labs(
    title = "Phi Recovery by Scenario",
    x = "True Phi",
    y = "Estimated Phi"
  ) +
  theme_minimal()

# ---- 3. Phi Bias Curves (Professor's Request) ----
# 横軸: phiの真値, 縦軸: バイアス (推定値 - 真値)
phi_data <- results_raw %>%
  filter(param == "phi_mu") %>%
  mutate(iter_bias = est - true)

ggplot(phi_data, aes(x = true, y = iter_bias, color = as.factor(T))) +
  geom_smooth(method = "loess", se = TRUE) + # 曲線を描画
  geom_hline(yintercept = 0, linetype = "dotted") +
  facet_wrap(~N, labeller = label_both) + # サンプルサイズ(N)ごとに分割
  labs(
    title = "Phi Estimation Bias by Sample Size",
    x = "True Phi Value ",
    y = "Bias ",
    color = "Time points (T)"
  ) +
  theme_minimal()
