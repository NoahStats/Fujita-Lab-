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
  real level;          // population mean
  real A;              // amplitude
  real P;              // period
  real phi;            // phase
  real<lower=0> sd_level;
  real<lower=0> sd_trend_it;
  real<lower=0> sd_A;
  real<lower=0> sd_P;
  real<lower=0> sd_phi;

  // Individual-level parameters
  vector[n] level_i;
  vector[n] A_i;
  vector[n] P_i;
  vector[n] phi_i;

  // Individual deviations from system trend
  matrix[n, T] mu_it;

  // Observation noise
  real<lower=0> sd_y;
}

model {
  // System trend evolution
  for (t in 2:T) {
    mu_t[t] ~ normal(mu_t[t-1], sd_mu_t);
  }

  // Hierarchical priors
  level_i ~ normal(level, sd_level);
  A_i ~ normal(A, sd_A);
  P_i ~ normal(P, sd_P);
  phi_i ~ normal(phi, sd_phi);

  for (i in 1:n) {
    for (t in 1:T) {
      mu_it[i, t] ~ normal(mu_t[t], sd_trend_it);
    }
  }

  // Observation model
  for (i in 1:n) {
    for (t in 1:T) {
      real season = A_i[i] * sin(2 * pi() / P_i[i] * t + phi_i[i]);
      real mu0 = level_i[i] + mu_it[i, t] + season;
      y[i][t] ~ normal(mu0, sd_y);
    }
  }
}

generated quantities {
  matrix[n, T] y_rep;

  for (i in 1:n) {
    for (t in 1:T) {
      real season = A_i[i] * sin(2 * pi() / P_i[i] * t + phi_i[i]);
      real mu0 = level_i[i] + mu_it[i, t] + season;
      y_rep[i, t] = normal_rng(mu0, sd_y);
    }
  }
}
