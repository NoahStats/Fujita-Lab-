data {
  int<lower=1> T;
  vector[T] y;
  vector[T] time;
  vector[T] I_treat;
  vector[T] t_post;
}

parameters {
  real beta_0;
  real beta_1;
  real beta_2;
  real beta_3;
  real<lower=-1, upper=1> phi;
  real<lower=0> sigma;
}

model {
  vector[T] mu;

  mu = beta_0
      + beta_1 * time
      + beta_2 * I_treat
      + beta_3 * t_post;
      
   // Priors
  beta_0 ~ normal(0, 5);
  beta_1 ~ normal(0, 1);
  beta_2 ~ normal(0, 2);
  beta_3 ~ normal(0, 1);
  phi    ~ uniform(-1, 1);
  sigma  ~ exponential(1);

  // Initial observation (stationary)
  y[1] ~ normal(mu[1] / (1 - phi), sigma / sqrt(1 - phi^2));

  // AR term
  for (t in 2:T) {
    y[t] ~ normal(mu[t] + phi * y[t - 1], sigma);
  }
}
