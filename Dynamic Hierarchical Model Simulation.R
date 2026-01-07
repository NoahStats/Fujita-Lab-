
library(ggplot2)
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

# state space models  -----------------------------------------------------

# GP regression 

