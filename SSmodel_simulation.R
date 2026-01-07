
library(ggplot2)

# true model  -------------------------------------------------------------

set.seed(123)

#parameters 
n_mice = 10
T = 300
P = 50 # seasonal period

sigma_trend_global = 0.02
sigma_obs = 0.5

data = list()

for(i in 1:n_mice){
  
  #mouse-specific parameters 
  sigma_trend = abs(rnorm(1,0,sigma_trend_global))
  A = rnorm(1,0.5, 0.1)
  phi = runif(1,0,2*pi)
  p_spike = runif(1,0.1, 0.8)
  
  # latent states
  mu = numeric(T)
  v = numeric(T)
  mu[1] = 0 
  v[1] = 0 
  
  for(t in 2:T){
    v[t] = v[t-1] + rnorm(1,0,sigma_trend)
    mu[t] = mu[t-1] + v[t]
  }
  
  # seasonality 
  s = A * sin(2*pi*(1:T)/P + phi)
  
  # spikes 
  spikes = rbinom(T,1,p_spike) * rlnorm(T,0,0.8)
  
  # observed data
  y = mu + s + spikes + rnorm(T,0,sigma_obs)
   
  data[[i]] = data.frame(
    mouse = i, 
    time = 1:T, 
    y = y
  )
}
print(data)

sim_data = do.call(rbind, data)


ggplot(sim_data, aes(x = time, y = y, color = factor(mouse))) +
  geom_line(alpha = 0.8) +
  labs(
    x = "Time",
    y = "Serotonin concentration",
    color = "Mouse"
  ) +
  theme_minimal()



# state space models  -----------------------------------------------------

# GP regression 
