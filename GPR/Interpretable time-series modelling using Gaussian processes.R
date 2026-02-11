setwd('/home/moosehunter/R/Fujita Lab/GPR/')
# https://hendersontrent.github.io/posts/2024/05/gaussian-process-time-series/
library(Rcpp)
library(dplyr)
library(tidyr)
library(ggplot2)
library(MASS)
library(TSA)

# 2. Generating some data 
set.seed(123)
y = 3 * sin(2 * seq(0, 4 * pi, length.out = 100)) + runif(100)*2
trend = 0.08 * seq(from = 1, to = 100, by = 1)
y = y + trend
tmp = data.frame(timepoint = 1:length(y), y = y)

# Draw plot 
ggplot(data = tmp ) + 
    geom_point(aes(x = timepoint, y = y), color = 'black') + 
    geom_line(aes(x = timepoint, y = y), color = 'black') + 
    labs(x = 'Timepoint', 
    y = 'values')  + 
    theme_bw()

# 3. Decomposing a time series for modeling 
# we are interested in constructing a model that can accurately 
# capture temporal properties and generate sensible predictions on which to do inference on. 


#4. Gausian Processes 
x1 = seq(from = -2, to = 2, length.out = 100)

#4.1 SE Kernel or RBF kernel 

Rcpp::cppFunction("
NumericMatrix cov_exp_quad(NumericVector xa, NumericVector xb, double sigma, double l) {
  int n1 = xa.size();
  int n2 = xb.size();
  NumericMatrix K(n1, n2);
  
  for (int i = 0; i < n1; ++i) {
    for (int j = 0; j < n2; ++j) {
      double diff = xa[i] - xb[j];
      double exp_term = exp(-0.5 * pow(diff / l, 2));
      K(i, j) = pow(sigma, 2) * exp_term;
    }
  }
  
  return K;
}
")

mat_exp_quad = cov_exp_quad(x1, x1, 1,1)
mat_exp_quad[1:5, 1:5]


sample_kernel <- function(n, k = 5, mu = NULL, X, Sigma){
  
  # Prepare mean vector if NULL
  
  if(is.null(mu)){
    mu <- integer(length = n)
  }
  
  # Produce dynamic column names, draw from multivariate normal, and convert to dataframe
  
  thenames <- list()
  
  for(i in 1:k){
    thenames[[i]] <- paste0("Sample ", i)
  }
  
  thenames <- unlist(thenames)
  y <- as.data.frame(t(MASS::mvrnorm(n = k, mu = mu, Sigma = Sigma)))
  colnames(y) <- thenames
  draws <- data.frame(x = X)
  draws <- cbind(draws, y)
  
  # Wrangle into long format for ggplot
  
  draws <- reshape(draws, direction = "long",
                   v.names = "y",
                   varying = 2:ncol(draws),
                   times = names(draws)[2:ncol(draws)],
                   timevar = "sample_ind")
  
  # Draw plot
  
  p <- ggplot(data = draws, aes(x = x, y = y, colour = sample_ind)) +
      geom_line(linewidth = 0.7) +
      labs(x = "x",
           y = "y = f(x)",
           colour = NULL) +
      theme_bw() +
      theme(legend.position = "bottom")
  
  return(p)
}

p <- sample_kernel(n = nrow(mat_exp_quad), k = 5, mu = NULL, X = x1, Sigma = mat_exp_quad)
print(p)

mat_exp_quad2 <- cov_exp_quad(x1, x1, 1, 0.3)
p_b <- sample_kernel(n = nrow(mat_exp_quad2), k = 5, mu = NULL, X = x1, Sigma = mat_exp_quad2)
print(p_b)


plot_covariance <- function(matrix){
  
  mat_df <- as.data.frame(matrix)
  mat_df$x <- seq_len(nrow(mat_df))
  
  mat_df <- reshape(mat_df, 
                    direction = "long",
                    v.names = "values",
                    varying = 1:(ncol(mat_df) - 1),
                    times = 1:nrow(mat_df),
                    idvar = "x")
  
  p <- ggplot(data = mat_df) +
    geom_tile(aes(x = x, y = time, fill = values)) +
    labs(title = "Heatmap of covariance matrix",
         x = "x",
         y = "x",
         fill = "k(x,x')") +
    scale_fill_viridis_c(limits = c(0, 1),
                         breaks = seq(from = 0, to = 1, by = 0.2)) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  return(p)
}

plot_covariance(matrix = mat_exp_quad)

plot_covariance(matrix = mat_exp_quad2)

# 4.2  Rational Quadratic Kernel 
cppFunction("
NumericMatrix cov_rat_quad(NumericVector xa, NumericVector xb, double sigma, double alpha, double l) {
  int n1 = xa.size();
  int n2 = xb.size();
  NumericMatrix K(n1, n2);
  
  for (int i = 0; i < n1; ++i) {
    for (int j = 0; j < n2; ++j) {
      double diff = xa[i] - xb[j];
      double norm_sq = pow(diff, 2);
      double exp_term = pow(1 + norm_sq / (2 * alpha * pow(l, 2)), -alpha);
      K(i, j) = pow(sigma, 2) * exp_term;
    }
  }
  
  return K;
}
")

mat_rat_quad = cov_rat_quad(x1, x1, 1,1,1)
mat_rat_quad[1:5, 1:5]

p1 = sample_kernel(n = nrow(mat_rat_quad), k =5, mu = NULL, X = x1, Sigma = mat_rat_quad)
print(p1)

mat_rat_quad2 <- cov_rat_quad(x1, x1, 1, 1, 0.3)
p1_b <- sample_kernel(n = nrow(mat_rat_quad2), k = 5, mu = NULL, X = x1, Sigma = mat_rat_quad2)
print(p1_b)

mat_rat_quad3 <- cov_rat_quad(x1, x1, 1, 0.1, 1)
p1_c <- sample_kernel(n = nrow(mat_rat_quad3), k = 5, mu = NULL, X = x1, Sigma = mat_rat_quad3)
print(p1_c)

# 4.3 Periodic Kernel # 
cppFunction("
NumericMatrix cov_periodic(NumericVector xa, NumericVector xb, double sigma, double l, double p) {
  int n1 = xa.size();
  int n2 = xb.size();
  NumericMatrix K(n1, n2);
  
  for (int i = 0; i < n1; ++i) {
    for (int j = 0; j < n2; ++j) {
      double diff = std::abs(xa[i] - xb[j]);
      double sin_term = std::sin(M_PI * diff / p);
      double exp_term = std::exp(-2 / pow(l, 2) * std::pow(sin_term, 2));
      K(i, j) = pow(sigma, 2) * exp_term;
    }
  }
  
  return K;
}
")

mat_periodic <- cov_periodic(x1, x1, 1, 2, 1)
mat_periodic[1:5, 1:5]

p2 <- sample_kernel(n = nrow(mat_periodic), k = 5, mu = NULL, X = x1, Sigma = mat_periodic)
print(p2)

# 4.4 whitenoise

#5 GPR for Time series analysis 

Sigma_exp_quad = cov_exp_quad(y, y, 1, length(y))# trend kernel
Sigma_periodic = cov_periodic(y, y, 1,1,25) # periodic kernel
Sigma_white_noise = cov_noise(x = 1:nrow(Sigma_exp_quad), sigma = 0.01) # white noise kernel 
Sigma_comb = Sigma_exp_quad + Sigma_periodic + Sigma_white_noise # combined kernel 

peri = TSA::periodogram(y) # Do Fourier Transform to decompose signal
dd = data.frame(freq = peri$freq, spec = peri$spec) # Extract frequency and spctreal values 
order  = dd[order(-dd$spec), ] # rearrange to descending order 
time = 1/ order$f # convert frequency to origiinal scale
time = time[!time %in% length(y)] 
time[1]

# 5.2 Sampling from the prior 
p4 = sample_kernel(n = nrow(Sigma_comb), k = 5, mu = NULL, X = y, Sigma = Sigma_comb)
print(p4)

plot_covariance(matrix = Sigma_comb) + 
  scale_fill_viridis_c()

# 5.3 Predicting from the posterior 

cov_summed = function(xa, xb, sigma_1 = 1, sigma_2 = 1, 
l_1 = 1, l_2 = 1, p = 1){
  Sigma_exp_quad = cov_exp_quad(xa, xb, sigma_1, l_1)
  Sigma_periodic = cov_periodic(xa, xb, sigma_2, l_2, p)
  Sigma = Sigma_exp_quad + Sigma_periodic
  return(Sigma)
}

GP <- function(x, xprime, y, sigma = 0.05, cov_function, ...){
  Sigma_11 <- cov_function(x, x, ...) + diag(sigma ^ 2, length(x))
  Sigma_12 <- cov_function(x, xprime, ...)
  Sigma_inv <- t(solve(Sigma_11, Sigma_12))
  Sigma_22 <- cov_function(xprime, xprime, ...)
  Sigma_2 <- Sigma_22 - (Sigma_inv %*% Sigma_12) # Posterior covariance matrix
  mu_2 <- Sigma_inv %*% y # Posterior mean vector
  gp <- list(x, xprime, y, mu_2, Sigma_2)
  names(gp) <- c("x", "xprime", "y", "mu", "Sigma")
  gp <- structure(gp, class = c("GaussProc", "list"))
  return(gp)
}

mod <- GP(x = 1:length(y), xprime = seq(from = 1, to = 100, by = 10), 
         y = y, sigma = 0.8, cov_function = cov_summed, 
         sigma_1 = 5, sigma_2 = 1, l_1 = 75, l_2 = 1, p = 25)

plot.GaussProc <- function(x, data, draws = 100, ...){
  
  stopifnot(inherits(x, "GaussProc") == TRUE)
  
  # Generate dataframe of draws from posterior
  
  thenames <- list()
    
  for(i in 1:draws){
    thenames[[i]] <- paste0("Sample ", i)
  }
  
  thenames <- unlist(thenames)
  posterior <- as.data.frame(t(MASS::mvrnorm(n = draws, mu = x$mu, Sigma = x$Sigma)))
  colnames(posterior) <- thenames
  
  posterior <- posterior %>%
    mutate(timepoint = x$xprime) %>%
    pivot_longer(cols = 1:ncol(posterior), names_to = "sample", values_to = "y")
  
  posterior_mean <- posterior %>%
    reframe(mu = mean(y), .by = "timepoint")

  posterior_sd <- data.frame(timepoint = x$xprime,
                             sigma = sqrt(diag(x$Sigma))) # SD is square root of variance

  posterior_summary <- merge(x = posterior_mean, y = posterior_sd, by = "timepoint")
  multiplier <- 1.96 # Multiplier for 95% interval over normal distribution
  posterior_summary$lower <- posterior_summary$mu - (multiplier * posterior_summary$sigma)
  posterior_summary$upper <- posterior_summary$mu + (multiplier * posterior_summary$sigma)
  
  # Draw plot
  
  p <- ggplot(data = posterior_summary) +
    geom_ribbon(aes(x = timepoint, ymin = lower, ymax = upper), fill = "steelblue2", alpha = 0.5) +
    geom_line(data = data, aes(x = timepoint, y = y), colour = "black") +
    geom_point(data = data, aes(x = timepoint, y = y), colour = "black") +
    geom_line(aes(x = timepoint, y = mu), colour = "steelblue2", size = 1) +
    labs(title = "Posterior mean and 95% prediction interval",
         subtitle = paste0("Interval calculated over ", draws, " samples"),
         x = "Timepoint",
         y = "Value") +
    theme_bw()
  
  return(p)
}

p6 = plot(mod, data = tmp, 100)
print(p6)

mod2 <- GP(x = 1:length(y), xprime = seq(from = 1, to = 100, by = 5),
           y = y, sigma = 0.8, cov_function = cov_summed, 
           sigma_1 = 5, sigma_2 = 1, l_1 = 75, l_2 = 1, p = 25)

p7 <- plot(mod2, tmp, 100)
print(p7)

mod3 <- GP(x = 1:length(y), xprime = 1:length(y),
           y = y, sigma = 0.8, cov_function = cov_summed, 
           sigma_1 = 5, sigma_2 = 1, l_1 = 75, l_2 = 1, p = 25)

p8 <- plot(mod3, tmp, 100)
print(p8)

mod4 <- GP(x = 1:length(y), xprime = 1:length(y),
           y = y, sigma = 1.2, cov_function = cov_summed, 
           sigma_1 = 5, sigma_2 = 1, l_1 = 75, l_2 = 1, p = 25)

p9 <- plot(mod4, tmp, 100)
print(p9)

p10 <- sample_kernel(n = nrow(mod4$Sigma), k = 5, X = mod4$x, mu = mod4$mu, Sigma = mod4$Sigma) +
  labs(title = "5 realisations of the posterior distribution")

print(p10)
