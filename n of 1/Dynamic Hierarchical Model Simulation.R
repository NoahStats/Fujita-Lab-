setwd("//wsl.localhost/Ubuntu/home/moosehunter/R/Fujita Lab")
library(ggplot2)
library(rstan)
library(ggmcmc)
library(dplyr)

set.seed(1)

# true model  -----------------------------------------------------------
#based on Hierarchical Dynamic Models  Marina Silva Paez and Dani Gamerman


n = 10 # sample size
T = 300 # observations (time)

# parameters for observation equations
sd_y = 0.01

# parameters for structural equations 
level = 10
sd_level = 0.8
sd_mu = 0.5
A = 0.5
P = 50 
phi = 5
sd_A = 0.01
sd_P = 0.5
sd_phi = 0.8
sd_season = 0.8

# parameters for system equations
sd_mu_t = 0.3
 # seasonal period 
p = 0.01
mean_amp = 1
sd_amp = 0.01

# initial values for group level parameters
season0 = 0.3
amp0 = 0.8 
mu0= 10


data = list()

for(i in 1:n){
  # time dependent parameters (mean parameters)
  mu_t = numeric(T) ; mu_t[1] =  mu0

  for(t in 2:T){
    # System Equations
    mu_t[t] = mu_t[t-1] + rnorm(1,0,sd_mu_t)
  }
  
  # structural equations
  level_i = level + rnorm(1,0,sd_level)
  mu_it = mu_t + rnorm(t,0,sd_mu)
  A_i = A + rnorm(1,0,sd_A)
  P_i = P + rnorm(1,0, sd_P)
  phi_i = phi + rnorm(1,0, sd_phi)
  prob_it = rbinom(300,1,p)
  amp_it = abs(rnorm(T,mean_amp, sd_amp))
  
  season_it = A_i * sin(2 * pi * (1:T) /P_i + phi_i) + rnorm(t,0,sd_season)
  
  # observation equations
  y = level_i + mu_it + season_it + prob_it*amp_it + rnorm(t,0,sd_y) 
  data[[i]] = data.frame(mouse = i, time = 1:T, y = y)
}
sim_data = do.call(rbind,data)


ggplot(data = sim_data, aes(x = time, y = y, color = factor(mouse))) + 
  geom_line(alpha = 0.8) + 
  labs(x = 'Time',
       y = '',
       color = 'Mouse')



# Dynamic Hierarchical Model 1 ChatGPT  -----------------------------------------------------

library(rstan)

y_array <- array(NA, dim = c(n, T))
for (i in 1:n) {
  y_array[i, ] <- sim_data$y[sim_data$mouse == i]
}

stan_data <- list(
  n = n,
  T = T,
  y = y_array,
  p = p
)

fit <- stan(
  file = "Mice_SSmodel.stan",
  data = stan_data,
  chains = 4,
  iter = 2000,
  cores = 4
)



print(fit, pars = c("level", "A", "P", "phi", "sd_y", "sd_mu_t"))
library(ggmcmc)

get_stanmodel(fit)
get_elapsed_time(fit)
get_inits(fit)


save.image(file = 'result_DHM1.RData')
saveRDS(fit, file = "fit_DHM1.rds")

# Dynamic Hierarchical Models 2 with no spikes (me) ---------------------------
y_array2 <- array(NA, dim = c(n, T))
for (i in 1:n) {
  y_array[i, ] <- sim_data$y[sim_data$mouse == i]
}

stan_data <- list(
  n = n,
  T = T,
  y = y_array
)

fit2 <- stan(
  file = "Mice_SSmodel2.stan",
  data = stan_data,
  chains = 4,
  iter = 4000,
  cores = 11
)
print(fit2)

saveRDS(fit2, file = "fit_DHM2.rds")

# Dynamic Hierarchical Models 3 with no levels(centered) and spikes (me) ---------------------------


# GP regression 
# https://mc-stan.org/docs/stan-users-guide/gaussian-processes.html#gaussian-process-regression

#make centered data 
centered_sim_data = sim_data %>%
  group_by(mouse) %>% 
  mutate(
    y_mean = mean(y),
    y_centered = y - y_mean
  ) %>% 
  ungroup()

y_array3 <- array(NA, dim = c(n, T))
for (i in 1:n) {
  y_array3[i, ] <- centered_sim_data$y_centered[centered_sim_data$mouse == i]
}

stan_data3 <- list(
  n = n,
  T = T,
  y = y_array3
)

fit3 <- stan(
  file = "Mice_SSmodel3.stan",
  data = stan_data3,
  chains = 4,
  iter = 2000,
  cores = 11
)

View(summary(fit3)$summary)
saveRDS(fit3, file = "fit_DHM3.rds")


# 4 Dynamic hierarchical Model with no level, spikes,
# no hierarchy in Period P 
# non-centered A_i and phi_i 
# Dynamic Hierarchical Model 4 marginalize mu_it (no more individual trends)--------------------

centered_sim_data = sim_data %>%
  group_by(mouse) %>% 
  mutate(
    y_mean = mean(y),
    y_centered = y - y_mean
  ) %>% 
  ungroup()

y_array4 <- array(NA, dim = c(n, T))
for (i in 1:n) {
  y_array4[i, ] <- centered_sim_data$y_centered[centered_sim_data$mouse == i]
}

stan_data4 <- list(
  n = n,
  T = T,
  y = y_array4
)

fit4 <- stan(
  file = "Mice_SSmodel4.stan",
  data = stan_data4,
  chains = 4,
  iter = 2000,
  cores = 11
)
saveRDS(fit4, file = "fit_DHM4.rds")

View(summary(fit4)$summary)

# 5 Dynamic hierarchical Model with no level, spikes,
# no hierarchy in A, P and Phi 
# marginalize mu_it (no more individual trends)--------------------

centered_sim_data = sim_data %>%
  group_by(mouse) %>% 
  mutate(
    y_mean = mean(y),
    y_centered = y - y_mean
  ) %>% 
  ungroup()

y_array5 <- array(NA, dim = c(n, T))
for (i in 1:n) {
  y_array5[i, ] <- centered_sim_data$y_centered[centered_sim_data$mouse == i]
}

stan_data5 <- list(
  n = n,
  T = T,
  y = y_array5
)

fit5 <- stan(
  file = "Mice_SSmodel5.stan",
  data = stan_data5,
  chains = 4,
  iter = 2000,
)  cores = 11

save.image(file = 'output/fit_DHM5.RData')
View(summary(fit5)$summary)
ggmcmc(ggs(fit5, inc_warmup = TRUE, stan_include_auxiliar = TRUE),
       file = 'output.fit-traceplot.pdf', plot = 'traceplot')


# 6 simple AR(1) model with no hierarchy (no more individual trends)--------------------

centered_sim_data = sim_data %>%
  group_by(mouse) %>% 
  mutate(
    y_mean = mean(y),
    y_centered = y - y_mean
  ) %>% 
  ungroup()

y_array6 <- array(NA, dim = c(n, T))
for (i in 1:n) {
  y_array6[i, ] <- centered_sim_data$y_centered[centered_sim_data$mouse == i]
}

stan_data6 <- list(
  n = n,
  T = T,
  y = y_array6
)

fit6 <- stan(
  file = "Mice_SSmodel6.stan",
  data = stan_data6,
  chains = 4,
  iter = 2000,
  cores = 11
  )

#saving outputs
save.image(file = 'output/fit_DHM6.RData')
View(summary(fit6)$summary)
ggmcmc(ggs(fit6, inc_warmup = TRUE, stan_include_auxiliar = TRUE),
       file = 'output.fit-traceplot.pdf', plot = 'traceplot')
summary(fit6)
class(fit6)


#extracting parameters 
ms6 = rstan::extract(fit6)
fit6_mu_t = rstan::extract(fit6, pars = 'mu_t')$mu_t
dim(fit6_mu_t)
fit6_mu_t_EAP = colMeans(fit6_mu_t)

# mean trend
fit6_pred = rstan::extract(fit6, pars = 'y_rep')$y_rep
dim(fit6_pred)
pred_EAP <- apply(fit6_pred, c(2,3), mean)
View(pred_EAP)
fit6_pred = as.data.frame(t(pred_EAP))
fit6_pred

fit6_mean_trend <- fit6_pred %>%
  mutate(time = row_number()) %>%
  pivot_longer(
    cols = starts_with("V"),
    names_to = "mouse",
    values_to = "mean_trend"
  ) %>%
  mutate(mouse = as.integer(sub("V", "", mouse))) %>%
  select(mouse, time, mean_trend)

#long data format as in sim_data
fit6_long = fit6_pred %>% 
  mutate(time = row_number()) %>%
  pivot_longer(
    cols = starts_with('V'),
    names_to = 'mouse',
    values_to = 'y'
  ) %>% 
  mutate(mouse = as.integer(sub('V','',mouse)))%>% 
  select(mouse, time,y)

sim_and_fit6_long = centered_sim_data%>% 
  left_join(
    fit6_long, 
    by = c('mouse', 'time'),
    suffix = c('_true', '_pred')
  )%>%
  left_join(
    fit6_mean_trend,
    by = c("mouse", "time")
  )

#plot 

sim_and_fit6_long %>%
  filter(mouse == 1) %>%
  ggplot(aes(x = time)) +
  geom_line(aes(y = y_centered), alpha = 0.8, color = "black") +
  geom_line(aes(y = mean_trend), color = "red", linewidth = 1) +
  labs(
    x = "Time",
    y = "",
    title = "Mouse 1: Observed vs Mean Trend"
  )

# comparing true trend vs EAP trend 

fit6_res =fit6_res =mu_t
  ggplot(sim_and_fit6_long, aes(x = time)) +
  geom_line(aes(y = y_centered), alpha = 0.4) +
  geom_line(aes(y = mean_trend), color = "red", linewidth = 1) +
  facet_wrap(~ mouse, scales = "free_y") +
  labs(
    x = "Time",
    y = "",
    title ="EPA trend and observed data"
  ) +
  theme_minimal()
ggsave(file = 'output/fit6_res.png', plot = fit6_res)

# comparing true trend vs EAP trend vs observed data
centered_mu_t = mu_t - mean(mu_t)

true_trend_df <- data.frame(
  time = 1:300,
  true_centered_trend = centered_mu_t
)

sim_and_fit6_long <- sim_and_fit6_long %>%
  left_join(true_trend_df, by = "time")

fit6_res = ggplot(sim_and_fit6_long, aes(x = time)) +
  geom_line(aes(y = y_centered), alpha = 0.35) +
  geom_line(aes(y = mean_trend), color = "red", linewidth = 1) +
  geom_line(
    aes(y = true_centered_trend),
    color = "blue",
    linetype = "dashed",
    linewidth = 1
  ) +
  facet_wrap(~ mouse, scales = "free_y") +
  labs(
    x = "Time",
    y = "",
    title = "Observed Data, Posterior Mean Trend, and True Centered Trend"
  ) +
  theme_minimal()
ggsave(file = 'figures/fit6_res.png', plot = fit6_res)

