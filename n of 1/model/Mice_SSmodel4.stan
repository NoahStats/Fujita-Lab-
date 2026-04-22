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
  real<lower=0> sd_A;
  real<lower=0> sd_phi;

  // Individual-level parameters
  vector[n] A_raw;
  vector[n] phi_raw;

  // Observation noise
  real<lower=0> sd_y;
}

transformed parameters{
    vector[n] A_i   = A   + sd_A   * A_raw;
    vector[n] phi_i = phi + sd_phi * phi_raw;
}

model {
    A_raw ~ std_normal();
    phi_raw ~ std_normal();

  // System trend evolution
  for (t in 2:T) {
    mu_t[t] ~ normal(mu_t[t-1], sd_mu_t);
  }


  // Observation model
  for (i in 1:n) {
    for (t in 1:T) {
      real season = A_i[i] * sin(2 * pi() / P * t + phi_i[i]);
      y[i][t] ~ normal(mu_t[t] + season, sd_y);
    }
  }
}

generated quantities {
  matrix[n, T] y_rep;

  for (i in 1:n) {
    for (t in 1:T) {
      real season = A_i[i] * sin(2 * pi() / P * t + phi_i[i]);
      y_rep[i, t] = normal_rng(mu_t[t] + season, sd_y);
    }
  }
}
