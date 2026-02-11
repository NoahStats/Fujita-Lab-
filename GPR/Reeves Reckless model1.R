setwd('/home/moosehunter/R/Fujita Lab/GPR/')
library(openesm)
library(rlang)
library(dplyr)
library(ggplot2)
library(tidyr)

ls = list_datasets()
View(ls)
d1 = get_dataset('0034')
d2 = d1$data
View(d2)

# all id
d_all_id_clean <- d2 %>%
  dplyr::select(completion_date, reckless, id) %>%
  filter(!is.na(completion_date),
         !is.na(reckless)) %>%
  group_by(id) %>%
  mutate(
    time_day = as.numeric(
      difftime(
        completion_date,
        as.POSIXct(
          as.Date(min(completion_date)),
          tz = attr(completion_date, "tzone")
        ),
        units = "days"
      )
    )
  ) %>%
  ungroup()

# uncentered reckless for each id 
ggplot(d_all_id_clean,
       aes(x = time_day, y = reckless)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ id, scales = "free_x") +
  theme_minimal()

# C002 making the midnight of the 1st day 0, and making 1 day 1
d_C002 <- d2 %>%
  filter(id == 'C002') %>%
  dplyr::select(reckless, completion_date)%>%
  filter(!is.na(completion_date)) %>%
  mutate(
    time_day = as.numeric(
      difftime(
        completion_date,
        as.POSIXct(
          as.Date(min(completion_date)),
          tz = attr(completion_date, "tzone")
        ),
        units = "days"
      )
    )
  ) %>%
  mutate(centered_reckless = reckless - mean(reckless))


# uncentered C002 
ggplot(d_C002, aes(x = time_day, y = reckless)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(
    x = "day",
    y = "Reckless",
    title = "Time Series of Reckless of C002"
  )

# centered C002 
ggplot(d_C002, aes(x = time_day, y = centered_reckless)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(
    x = "day",
    y = "Reckless",
    title = "Time Series of Reckless of C002"
  )


# model 1  with train = 84, test = 21 20%
# prior specified by ChatGPT 
# matern 3/2 for trend, periodic kernel  p = 7 and Gaussian noise 
head(d_C002)
d_C002_train = d_C002 %>%
  slice(1:84)
d_C002_test = d_C002 %>%
  slice(85:104)

N_train <- 84

x_train <- d_C002$time_day[1:N_train]
y_train <- d_C002$centered_reckless[1:N_train]

x_new <- d_C002_test$time_day


stan_data1<- list(
  N = nrow(d_C002_train),
  x = d_C002_train$time_day,
  y = d_C002_train$centered_reckless, 
  p = 7)



library(cmdstanr)
mod1 <- cmdstan_model(
  'Stan Model/Reeves_model_1.stan',
  force_recompile = TRUE
)

fit1 <- mod1$sample(
  data = stan_data1,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  refresh = 100
)

library(posterior)
draws <- fit1$draws()


draws <- fit1$draws()
theta <- as_draws_df(draws)
d <- theta[1, ]

alpha_trend <- d$alpha_trend
rho_trend   <- d$rho_trend
alpha_per   <- d$alpha_per
rho_per     <- d$rho_per
sigma       <- d$sigma


x <- d_C002_train$time_day
y <- d_C002_train$centered_reckless

k_matern32 <- function(x1, x2, alpha, rho) {
  d <- abs(outer(x1, x2, "-"))
  alpha^2 * (1 + sqrt(3) * d / rho) * exp(-sqrt(3) * d / rho)
}

k_periodic <- function(x1, x2, alpha, rho, p) {
  d <- abs(outer(x1, x2, "-"))
  alpha^2 * exp(-2 * sin(pi * d / p)^2 / rho^2)
}

draws <- fit1$draws() |> posterior::as_draws_df()

theta <- draws |>
  summarise(
    alpha_trend = mean(alpha_trend),
    rho_trend   = mean(rho_trend),
    alpha_per   = mean(alpha_per),
    rho_per     = mean(rho_per),
    sigma       = mean(sigma)
  )

x_train <- d_C002_train$time_day
y_train <- d_C002_train$centered_reckless

x_grid <- seq(min(x_train), max(x_train), length.out = 300)

# Covariance matrices
K <- k_matern32(x_train, x_train,
                theta$alpha_trend, theta$rho_trend) +
     k_periodic(x_train, x_train,
                theta$alpha_per, theta$rho_per, p = 7)

diag(K) <- diag(K) + theta$sigma^2 + 1e-6

K_star <- k_matern32(x_grid, x_train,
                     theta$alpha_trend, theta$rho_trend) +
          k_periodic(x_grid, x_train,
                     theta$alpha_per, theta$rho_per, p = 7)

K_starstar <- k_matern32(x_grid, x_grid,
                         theta$alpha_trend, theta$rho_trend) +
              k_periodic(x_grid, x_grid,
                         theta$alpha_per, theta$rho_per, p = 7)

# Cholesky
L <- chol(K)

# Posterior mean
alpha <- backsolve(L, forwardsolve(t(L), y_train))
mu_star <- K_star %*% alpha

# Posterior variance
v <- forwardsolve(t(L), t(K_star))
Sigma_star <- K_starstar - t(v) %*% v
sd_star <- sqrt(pmax(diag(Sigma_star), 0))

plot_df <- data.frame(
  time = x_grid,
  mean = as.vector(mu_star),
  lower = mu_star - 1.96 * sd_star,
  upper = mu_star + 1.96 * sd_star
)


p = ggplot() +
  geom_ribbon(data = plot_df,
              aes(time, ymin = lower, ymax = upper),
              fill = "steelblue", alpha = 0.25) +
  geom_line(data = plot_df,
            aes(time, mean),
            color = "steelblue", linewidth = 1) +
  geom_point(data = d_C002_train,
             aes(time_day, centered_reckless),
             color = "black", size = 1.8) +
  labs(
    title = "GP posterior with between-points uncertainty",
    x = "Time (days)",
    y = "Centered reckless score"
  ) +
  theme_minimal()

ggsave('figures/model1.png',
  plot = p)
