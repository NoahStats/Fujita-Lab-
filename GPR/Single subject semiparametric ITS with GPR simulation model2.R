#model2; no periodicity
#OU_amp = 0.8; OU_scale = 0.5 # hyparparameters for OU
#intercept = 0; slope = 0 # basline intercept and slope 
#intercept_treat  = 0.15; slope_treat = 0.005 # treatment effects
#sigma = 0.1

# Windows Rstuido
setwd("//wsl.localhost/Ubuntu/home/moosehunter/R/Fujita Lab/GPR")
#vscode
setwd('/home/moosehunter/R/Fujita Lab/GPR') 
library(rstan)
library(tidyverse)
library(MASS)
library(ggmcmc)
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

OU_amp = 0.8; OU_scale = 0.5 # hyparparameters for OU

OU_kernel = function(t_1,t_2){ 
  OU_amp^2 * exp(-abs(t_1 - t_2)/OU_scale)
} 

# define covariance matrices using the kernels and the inputs
OU_cov_matrix = matrix(0,nrow = 140, ncol  = 140) 
for(s in 1:140){ 
  for(t in 1:140){
    OU_cov_matrix[s,t] = OU_kernel(time_points[s], time_points[t])
  }
}

#sample OU_term and Periodic_term from MVN distributions
diag(OU_cov_matrix) <- diag(OU_cov_matrix) + 1e-12 # ensure positive definitess

OU_Cholesky = chol(OU_cov_matrix)

OU_term = t(U_Cholesky) %*% 
  mvrnorm(n = 1, mu = rep(0,140), Sigma = diag(140))

# experiment with different scale parameters 

ggplot(data = data.frame(time_points, OU_term),
       aes(x= time_points, y = OU_term))+
  geom_line()

# Add the linear term, OU and the periodic term
y = linear_term + OU_term  +
  mvrnorm(1,mu = rep(0,140), Sigma = 0.1*diag(140))
y = as.numeric(y)

df = data.frame(
  time = time_points, 
  y = y, 
  linear = linear_term,
  OU = OU_term
  )

# Visualization 
df_long = pivot_longer(df, cols = c(y,linear,OU), 
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

stan_data = list(T = length(y), time_points = time_points, Y = y)
fit1 = stan(file = "Stan Model/Single subject semiparametric ITS2.stan",data = stan_data )

save.image(file = 'output/single subject Semiparametric ITS model2.Rdata')
write.table(data.frame(summary(fit1)$summary),
            file = 'output/Single subject semiparametric ITS model2-summary',
            sep = '\t', quote = FALSE, col.names = NA)
ggmcmc(ggs(fit1, inc_warmup = TRUE, stan_include_auxiliar = TRUE),
       file = 'output/Single subject semiparametric ITS model2-traceplot.pdf', 
       plot = 'traceplot')
ggmcmc(ggs(fit1), file = 'output/Single subject semiparametric ITS model2-ggmcmc.pdf')
