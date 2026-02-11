data {
  int<lower=1> n;        // number of individuals
  int<lower=1> T;        // number of time points
  vector[T] y[n];        // observed data
}

parameters {
  // System-level trend
  vector[T] mu_t;
  real<lower=0> sd_mu_t;

  // Population-level parameters
  real A;              // amplitude
  real P;              // period
  real phi;            // phase


  // Observation noise
  real<lower=0> sd_y;
}


model {

  // System trend evolution
  for (t in 2:T) {
    mu_t[t] ~ normal(mu_t[t-1], sd_mu_t);
  }


  // Observation model
  for (i in 1:n) {
    for (t in 1:T) {
      real season = A * sin(2 * pi() / P * t + phi);
      y[i][t] ~ normal(mu_t[t] + season, sd_y);
    }
  }
}

generated quantities {
  matrix[n, T] y_rep;

  for (i in 1:n) {
    for (t in 1:T) {
      real season = A * sin(2 * pi() / P * t + phi);
      y_rep[i, t] = normal_rng(mu_t[t] + season, sd_y);
    }
  }
}
