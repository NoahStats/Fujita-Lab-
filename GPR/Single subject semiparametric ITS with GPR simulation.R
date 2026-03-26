setwd("//wsl.localhost/Ubuntu/home/moosehunter/R/Fujita Lab/GPR")
library(rstan)
library(tidyverse)
library(MASS)

set.seed(123)

# generating a time series data  ------------------------------------------
# n = 1, T = 140(70 pre and 70 post treatments)
# randomly sample a time points from the 2hour windows 5 times a day 
# 10:12 , 12:14, 14:16, 16:18, 18:20 
# a linear model for the mean function 
# GP for the circadian rhythm and the auto-regressive structure 

#randomly draw time points from the fixed windows 5 days a day for 4 weeks
time_points = numeric(140) #unit is 'day' 
for(day in 1:28){
  print(1*day)
  time_points[5*(day-1)+1] = runif(1,10,12)/24 + (day-1)
  time_points[5*(day-1)+2] = runif(1,12,14)/24 + (day-1)
  time_points[5*(day-1)+3] = runif(1,14,16)/24 + (day-1)
  time_points[5*(day-1)+4] = runif(1,16,18)/24 + (day-1)
  time_points[5*(day-1)+5] = runif(1,18,20)/24 + (day-1)
}

# define the mean function 
linear_term = numeric(140) #linear term for the tratement effects

intercept = 0; slope = 0 # basline intercept and slope 
intercept_treat  = 0.15; slope_treat = 0.005 # treatment effects

for(t in 1:140){ 
  print(t)
  if(t < 70)
  linear_term[t] = intercept + slope*time_points[t]
  
  else{
  linear_term[t] = intercept + slope*time_points[t]+
    intercept_treat + slope_treat*(time_points[t] - 14)
  }
}
plot(time_points, linear_term)

# define the kernel 
OU_term = numeric(140)
Periodic_term = numeric(140)

OU_amp = 0.5; OU_scale = 20 # hyparparameters for OU
Per_amp = 0.05; Per_length_scale = 20 #hyperparameters for Periodic kernel
#OU_scale = 0.05; Per_length_scale = 0.005

OU_kernel = function(t_1,t_2){ 
  OU_amp^2 * exp(-abs(t_1 - t_2)/OU_scale)
} 

Periodic_kernel = function(t_1, t_2){ #circadian rhythm
  Per_amp^2 * exp(-1/(2*Per_length_scale^2) * sin(pi*abs((t_1-t_2)/1))^2) 
}

# define covariance matrices using the kernels and the inputs
OU_cov_matrix = matrix(0,nrow = 140, ncol  = 140) 
Per_cov_matrix = matrix(0, nrow = 140, ncol =140)
for(s in 1:140){ 
  for(t in 1:140){
    OU_cov_matrix[s,t] = OU_kernel(time_points[s], time_points[t])
    Per_cov_matrix[s,t] = Periodic_kernel(time_points[s], time_points[t])
  }
}

#sample OU_term and Periodic_term from MVN distributions
diag(Per_cov_matrix) <- diag(Per_cov_matrix) + 1e-12 # ensure positive definitess
diag(OU_cov_matrix) <- diag(OU_cov_matrix) + 1e-12 # ensure positive definitess

OU_Cholesky = chol(OU_cov_matrix)
Per_Cholesky = chol(Per_cov_matrix)

OU_term = OU_Cholesky %*% 
  mvrnorm(n = 1, mu = rep(0,140), Sigma = diag(140))
Periodic_term = Per_Cholesky %*%
  mvrnorm(n = 1, mu = rep(0,140), Sigma = diag(140))

# experiment with different scale parameters in Periodic kernel
ggplot(data = data.frame(time_points, Periodic_term),
       aes(x= time_points, y = Periodic_term))+
  geom_line()

# Add the linear term, OU and the periodic term
y = linear_term + OU_term + Periodic_term +
  mvrnorm(1,mu = rep(0,140), Sigma = 0.1*diag(140))

df = data.frame(
  time = time_points, 
  y = y, 
  linear = linear_term,
  OU = OU_term, 
  Periodic = Periodic_term
)

# Visualization 
df_long = pivot_longer(df, cols = c(y,linear,OU, Periodic), 
                       names_to = 'variable',
                       values_to = 'value')
ggplot(df_long, aes(x = time, y = value, color = variable)) +
  geom_line() +
  geom_point(data = subset(df_long, variable == "y")) +
  scale_x_continuous(breaks = seq(1, max(df_long$time), by = 1)) +
  theme_minimal() #OU takes over in the baseline 
  # linear treatment effect is a decnet APPROXIMATION, 
  #esp when the OU effect is dominant. 



# Model estimation and recovery rate --------------------------------------

data = list(T = nrow(y), T = time_points, Y = y)
fit1 = stan(file = "Stan Model/Single subject semiparametric ITS1.stan")

