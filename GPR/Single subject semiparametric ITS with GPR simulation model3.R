#model3; no periodicity; same as model2 with different parameters
#OU_amp = 1; OU_scale = 2 # hyparparameters for OU
#intercept = 0; slope = 0 # basline intercept and slope 
#intercept_treat  = 0.4; slope_treat = 0.1 # treatment effects
#sigma = 0.5

setwd("//wsl.localhost/Ubuntu/home/moosehunter/R/Fujita Lab/GPR")
setwd('/home/moosehunter/R/Fujita Lab/GPR') # VS code
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
intercept_treat  = 0.4; slope_treat = 0.1 # treatment effects

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

OU_amp = 1; OU_scale = 2 # hyparparameters for OU

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

OU_term = t(OU_Cholesky) %*% 
  mvrnorm(n = 1, mu = rep(0,140), Sigma = diag(140))

# experiment with different scale parameters 

ggplot(data = data.frame(time_points, OU_term),
       aes(x= time_points, y = OU_term))+
  geom_line()

# Add the linear term, OU and the periodic term
y = linear_term + OU_term  +
  mvrnorm(1,mu = rep(0,140), Sigma = 0.5*diag(140))
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



# Baysian Model estimation and recovery rate --------------------------------------

stan_data = list(T = length(y), time_points = time_points, Y = y)
#the model is the same as stan.model2; only the parameters changed
stan_file <- normalizePath("Stan Model/Single subject semiparametric ITS2.stan")
fit1 = stan(file = stan_file,data = stan_data )

save.image(file = 'output/single subject Semiparametric ITS model3.Rdata')
write.table(data.frame(summary(fit1)$summary),
            file = 'output/Single subject semiparametric ITS model3-summary',
            sep = '\t', quote = FALSE, col.names = NA)
ggmcmc(ggs(fit1, inc_warmup = TRUE, stan_include_auxiliar = TRUE),
       file = 'output/Single subject semiparametric ITS model3-traceplot.pdf', 
       plot = 'traceplot')
ggmcmc(ggs(fit1), file = 'output/Single subject semiparametric ITS model3-ggmcmc.pdf')

# REML via mgcv -----------------------------------------------------------
library(mgcv)
library(tidyverse)

# ── Prep the data ──────────────────────────────────────────────────────────────
df <- df %>%
  mutate(
    treat       = as.numeric(row_number() >= 70),          # 0 = pre, 1 = post
    time_post   = ifelse(treat == 1, time - time_points[70], 0)  # time since treatment
  )

# ── Fit with REML ──────────────────────────────────────────────────────────────
# m = c(1, OU_scale_init) : Matern 1/2 ≡ OU kernel
# The second element is the range parameter (starting value / fixed value)
# Leave it free by letting mgcv estimate it, or fix it if you prefer.

fit_reml <- gam(
  y ~ 
    # --- linear ITS terms (parametric) ---
    time +                   # baseline slope
    treat +                  # treatment intercept shift
    time_post +              # treatment slope change
    # --- OU / autocorrelation structure (nonparametric GP) ---
    s(time, bs = "gp", m = c(1, 0.5), k = 40),
  data   = df,
  method = "REML",           # <-- REML smoothing parameter selection
  family = gaussian()
)

summary(fit_reml)

# Fixed effects (intercept, slope, treat intercept/slope)
coef(fit_reml)

# Smoothing/variance parameters
gam.object <- fit_reml
fit_reml$sp          # smoothing parameter (~ 1/OU_amp^2 * sigma^2)
fit_reml$sig2        # residual variance (sigma^2)

sigma2 <- fit_reml$sig2        # residual variance (sigma^2)
lambda <- fit_reml$sp          # smoothing parameter
implied_amp2 <- sigma2 / lambda
implied_amp  <- sqrt(implied_amp2)
implied_amp
