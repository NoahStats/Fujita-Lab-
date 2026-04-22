# model 2 
# Matern 1/2 for trend, weekly periodic kernel and white noise 
setwd("//wsl.localhost/Ubuntu/home/moosehunter/R/Fujita Lab/GPR")
library(openesm)
library(rlang)
library(dplyr)
library(ggplot2)
library(rstan)
library(tidyr)



d1 = get_dataset('0034')
d2 = d1$data

# all id
View(d2)
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
  mutate(centered_reckless = reckless - mean(reckless))%>%
  mutate(normalized_reckless = centered_reckless/sd(reckless))


# model 1  with train = 84, test = 21 20%
# prior specified by ChatGPT 
# matern 3/2 for trend, periodic kernel  p = 7 and Gaussian noise 
d_C002_train = d_C002 %>%
  slice(1:84)
d_C002_test = d_C002 %>%
  slice(85:104)

N_train <- 84

x_train <- d_C002$time_day[1:N_train] 
y_train <- d_C002$normalized_reckless[1:N_train]

x_new <- d_C002_test$time_day


stan_data1<- list(
  N = nrow(d_C002_train),
  x = d_C002_train$time_day,
  y = d_C002_train$normalized_reckless, # use normalized reckless score
  p = 7)

fit1 = stan(file = 'Stan Model/Reeves reckless model2.stan', data = stan_data1)
save.image(file = 'output/Reeves reckless model2.Rdata')

library(ggmcmc)
load('output/Reeves reckless model2.Rdata')

write.table(data.frame(summary(fit1)$summary),
            file = 'output/Reeves reckless model2 fit-summary.txt',
            sep = '\t', quote = FALSE, col.names = NA)

ggmcmc(ggs(fit1, inc_warmup = TRUE, stan_include_auxiliar = TRUE), 
       file = 'output/Reeves reckless model2 traceplot.pdf', plot = 'traceplot')


# Extract the summary table
posterior_summary <- summary(fit1)$summary

# Extract the EAP (the 'mean' column)
eap_estimates <- posterior_summary[, "mean"]

# Extract individual values correctly
alpha_trend <- eap_estimates["alpha_trend"]
rho_trend   <- eap_estimates["rho_trend"]
alpha_per   <- eap_estimates["alpha_per"]
rho_per     <- eap_estimates["rho_per"]
sigma       <- eap_estimates["sigma"]

# Now you can use them:
print(alpha_trend)

x <- d_C002_train$time_day
y <- d_C002_train$normalized_reckless


k_matern12 <- function(x1, x2, alpha, rho) {
  d <- abs(outer(x1, x2, "-"))
  alpha^2  * exp(- d / rho)
}

k_periodic <- function(x1, x2, alpha, rho, p) {
  d <- abs(outer(x1, x2, "-"))
  alpha^2 * exp(-2 * sin(pi * d / p)^2 / rho^2)
}


x_train <- d_C002_train$time_day
y_train <- d_C002_train$normalized_reckless

x_grid <- seq(min(x_train), max(x_train), length.out = 300)

# Covariance matrices
K <- k_matern12(x_train, x_train,
                alpha_trend, rho_trend) +
     k_periodic(x_train, x_train,
                alpha_per, rho_per, p = 7)

diag(K) <- diag(K) + sigma^2 + 1e-6

K_star <- k_matern12(x_grid, x_train,
                     alpha_trend, rho_trend) +
          k_periodic(x_grid, x_train,
                     alpha_per, rho_per, p = 7)

K_starstar <- k_matern12(x_grid, x_grid,
                         alpha_trend, rho_trend) +
              k_periodic(x_grid, x_grid,
                         alpha_per, rho_per, p = 7)

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
             aes(time_day, normalized_reckless),
             color = "black", size = 1.8) +
  labs(
    title = "GP posterior with between-points uncertainty",
    x = "Time (days)",
    y = "Centered reckless score"
  ) +
  theme_minimal()
p
ggsave('figures/model2_training.png',
  plot = p)





# prediction 

x_test <- x_new
y_test <- d_C002_test$centered_reckless

# Cross-covariance between test and train
K_star_test <- k_matern12(x_test, x_train,
                          theta$alpha_trend, theta$rho_trend) +
               k_periodic(x_test, x_train,
                          theta$alpha_per, theta$rho_per, p = 7)

# Covariance among test points
K_starstar_test <- k_matern12(x_test, x_test,
                              theta$alpha_trend, theta$rho_trend) +
                   k_periodic(x_test, x_test,
                              theta$alpha_per, theta$rho_per, p = 7)

# Posterior mean
mu_test <- K_star_test %*% alpha

# Posterior variance
v_test <- forwardsolve(t(L), t(K_star_test))
Sigma_test <- K_starstar_test - t(v_test) %*% v_test
sd_test <- sqrt(pmax(diag(Sigma_test), 0))

# Put into dataframe
test_df <- data.frame(
  time = x_test,
  mean = as.vector(mu_test),
  lower = mu_test - 1.96 * sd_test,
  upper = mu_test + 1.96 * sd_test,
  observed = y_test
)
library(ggplot2)

p_test <- ggplot(test_df, aes(x = time)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "darkorange", alpha = 0.25) +
  geom_line(aes(y = mean),
            color = "darkorange", linewidth = 1) +
  geom_point(aes(y = observed),
             color = "black", size = 2) +
  labs(
    title = "GP Test Predictions",
    x = "Time (days)",
    y = "Centered reckless score"
  ) +
  theme_minimal()

ggplot(test_df, aes(x = time)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "darkorange", alpha = 0.25) +
  geom_line(aes(y = mean),
            color = "darkorange", linewidth = 1) +
  geom_point(aes(y = observed),
             color = "black", size = 2) +
  labs(
    title = "GP Test Predictions",
    x = "Time (days)",
    y = "Centered reckless score"
  ) +
  theme_minimal()
p_test
ggsave("figures/model2_test_predictions.png",
       plot = p_test,
       width = 7, height = 5)




