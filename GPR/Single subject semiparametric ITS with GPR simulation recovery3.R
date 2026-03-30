#Recovery3: checking parameter recovery rate by correlation
# assume positive treatment effects in the stan code

# randomly generate 50 datasets, compute the correlation for each parameter
# standardize the data by the baseline data
# use the same model as model 2~4

#Y_t ~ GP(mu_t, k(t,t') + sigma^2 I)
# mu_t = b0 + b1*t + b2*1(t>14) + b3*1(t>14)*(t-14)
# k(t,t') = alpha^2 * exp(-|t-t'|/l)
# 28 days total,5 times a day, 14 days pre and 14 days post treatment

#generate standardized(with respect to basleine)datasets by the following method
# 1. sample b1 ~ Unif(-1/14, 1/14) ; slope taking max of 1sd on day 14
# 2. set b0 = -7*b1 ; to keep the basline mean 0 
# 3. sample b2 ~ Unif(0, 0.6); basline change
# 4. b3 ~ Unif(0,1/14); max 1 * SD change for slope                                                                                       
# 5. sample alpha  ~ Unif(0.3, 0.8 )
# 6. set sigma = sqrt(1-alpha^2); basline variance alpha^2 + sigma^2 = 1
# 7. sample l ~ Unif(1,2)


setwd("//wsl.localhost/Ubuntu/home/moosehunter/R/Fujita Lab/GPR") #Rstudio
setwd('/home/moosehunter/R/Fujita Lab/GPR') # VScode
library(rstan)
library(ggplot2)
library(MASS)
library(foreach)
library(doParallel)
library(ggmcmc)
library(purrr)
set.seed(123)

# generating a time series data  ------------------------------------------
# n = 1, T = 140(70 pre and 70 post treatments)
# randomly sample a time points from the 2hour windows 5 times a day 
# 10:12 , 12:14, 14:16, 16:18, 18:20 
# a linear model for the mean function 
# GP for the circadian rhythm and the auto-regressive structure 

n_sim = 50 # number of datasets generated
records = vector('list', n_sim) #making a list with 50, this is efficient conpared to list()

#define the exponential kernel
OU_kernel = function(t_1,t_2, alpha, l){ 
  alpha^2 * exp(-abs(t_1 - t_2)/l)
} 

### simulate the parameters and the datasets 50 times
for(k in 1:n_sim){
  # 1. sample b1 ~ Unif(-1/14, 1/14) ; slope taking max of 1sd on day 14
  # 2. set b0 = -7*b1 ; to keep the basline mean 0      
  # 3. sample b2 ~ Unif(0, 0.6); basline change
  # 4. b3 ~ Unif(0,1.5/14); max 1.5 * SD change for slope                                                                                       
  # 5. sample alpha  ~ Unif(0.3, 0.8 )
  # 6. set sigma = sqrt(1-alpha^2); basline variance alpha^2 + sigma^2 = 1
  # 7. sample l ~ Unif(1,2)
  
  # sample the parameters
  b1 = runif(1, -1/14,1/14) ; b0 = -7 * b1 #base line
  b2 = runif(1, 0,0.6); b3 = runif(1, 0, 1/14) #treatment effects
  alpha = runif(1, 0.3, 0.8)
  sigma = sqrt(1-alpha^2) #this ensures that the var in the baseline alpha2 + sigma2 = 1
  l = runif(1, 1,2)
  
  #randomly draw time points from the fixed windows 5 days a day for 4 weeks
  
  time_points = numeric(140) #unit is 'day' 
  for(day in 1:28){
    time_points[5*(day-1)+1] = runif(1,10,12)/24 + (day-1)
    time_points[5*(day-1)+2] = runif(1,12,14)/24 + (day-1)
    time_points[5*(day-1)+3] = runif(1,14,16)/24 + (day-1)
    time_points[5*(day-1)+4] = runif(1,16,18)/24 + (day-1)
    time_points[5*(day-1)+5] = runif(1,18,20)/24 + (day-1)
  }
  
  
  # define the mean function mu
  mu = numeric(140) #linear term for the tratement effects
  for(t in 1:140){ 
    if(t < 70)
      mu[t] = b0 + b1*time_points[t]
    
    else{
      mu[t] = b0 + b1*time_points[t]+
        b2 + b3*(time_points[t] - 14)
    }
  }
  
  # define the GP(OU) term 
  GP = numeric(140)
  
  # define covariance matrices using the kernels and the inputs
  OU_cov_matrix = matrix(0,nrow = 140, ncol  = 140) 
  for(s in 1:140){ 
    for(t in 1:140){
      OU_cov_matrix[s,t] = OU_kernel(time_points[s], time_points[t], alpha, l)
    }
  }
  
  #sample OU_term and Periodic_term from MVN distributions
  diag(OU_cov_matrix) <- diag(OU_cov_matrix) + 1e-12 # ensure positive definitess
  
  OU_Cholesky = chol(OU_cov_matrix)
  
  GP = t(OU_Cholesky) %*% rnorm(140)
  
  # define the error 
  error = rnorm(140,0, sigma)
  
  # total data
  Y = mu + GP + error 
  
  # recoding the parameter and the data
  records[[k]] = list(
    true = list(b0 = b0, b1 = b1, b2 = b2, b3 =  b3,
                alpha = alpha,l = l,sigma = sigma),
    time_points = time_points, 
    Y = as.vector(Y)
  )
}


# Visualization 
df = data.frame(
  time = records[[5]]$time_points, 
  y = records[[5]]$Y
)

ggplot(df, aes(x = time, y = y)) +
  geom_line() +
  geom_point() +
  theme_minimal() 


#Model estimation and recovery rates --------------------------------------
n_cores = detectCores()-1 #the usual three lines
c1 = makeCluster(n_cores)
registerDoParallel(c1)

#compile the file beforehand
stan_model_obj <- stan_model(file = "Stan Model/Single subject semiparametric ITS with GPR simulation recovery3.stan")

clusterExport(c1, varlist = c('stan_model_obj', 'records')) #sends objects to each worker?
clusterEvalQ(c1, library(rstan)) #runs code on each worker? 

#iterate over k datasets using rstan
results = foreach( 
  k = 1:n_sim, 
  .packages = 'rstan', 
  .errorhandling = 'pass'
) %dopar% {
  
  rec = records[[k]]#extract the k-th dataset
  
  stan_data = list( #prepare the stan df for the k-th data
    T = length(rec$Y), 
    time_points = rec$time_points, 
    Y = rec$Y
  ) 
  
  fit = sampling(
    object = stan_model_obj, 
    data = stan_data, 
    chains = 4, 
    iter = 2000, 
    warmup = 1000, 
    cores = 1,# no sub-parallelization
    refresh = 0, 
    seed = k
  )
  
  post = as.data.frame(fit) #extract the posterior the k-th fit
  
  #RETURN OBJECT1
  list(
    sim  = k,
    true = rec$true,
    Posterior_summary = summary(fit)$summary)
}

stopCluster(c1) 

#calculate the recovery rates 
params <- c("b0", "b1", "b2", "b3", "alpha", "l", "sigma")

ok <- sapply(results, function(r) {
  !inherits(r, "error") &&
    !is.null(r$Posterior_summary)
})

cat(sprintf("Successful fits: %d / %d\n", sum(ok), n_sim))

recovery <- sapply(params, function(p) {
  true_vals <- as.numeric(sapply(results[ok], function(r) r$true[[p]]))
  est_vals  <- as.numeric(sapply(results[ok], function(r) r$Posterior_summary[p, "mean"]))
  cor(true_vals, est_vals)
})

saveRDS(results,  "output/Single subject semiparametric ITS recovery3_results.rds")
saveRDS(recovery, "output/Single subject semiparametric ITS recovery3_correlations.rds")



library(ggplot2)
library(tidyverse)

params <- c("b0", "b1", "b2", "b3", "alpha", "l", "sigma")

# ── Extract true and estimated values into a long data frame ──────────────────
ok <- sapply(results, function(r) !inherits(r, "error") && !is.null(r$Posterior_summary))

recovery_df <- map_dfr(params, function(p) {
  tibble(
    parameter = p,
    true      = sapply(results[ok], function(r) r$true[[p]]),
    estimated = sapply(results[ok], function(r) r$Posterior_summary[p, "mean"])
  )
})

# ── Compute correlations per parameter ───────────────────────────────────────
cor_df <- recovery_df %>%
  group_by(parameter) %>%
  summarise(r = cor(true, estimated), .groups = "drop") %>%
  mutate(label = sprintf("r = %.3f", r))

# ── Scatter plots ─────────────────────────────────────────────────────────────
p = ggplot(recovery_df, aes(x = true, y = estimated)) +
  geom_point(alpha = 0.7, size = 2, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # identity line
  geom_text(
    data    = cor_df,
    aes(label = label),
    x       = -Inf, y = Inf,
    hjust   = -0.1, vjust = 1.5,
    size    = 4, color = "black", inherit.aes = FALSE
  ) +
  facet_wrap(~ parameter, scales = "free") +
  labs(
    title    = "Parameter Recovery: True vs. Posterior Mean",
    subtitle = "Dashed red line = perfect recovery (slope 1, intercept 0)",
    x        = "True value",
    y        = "Posterior mean"
  ) +
  theme_bw()
ggsave(file = 'figures/Single subjec semiparametric ITS with GPR recovery3.png',p)
