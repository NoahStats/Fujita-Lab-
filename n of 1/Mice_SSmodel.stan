data {
  int<lower=1> n;           // number of mice
  int<lower=1> T;           // number of time points
  vector[T] y[n];           // observations
  real<lower=0,upper=1> p;  // spike probability
}

parameters {
  // ----- Latent state -----
  vector[T] mu_t;
  real<lower=0> sd_mu_t;

  // ----- Population-level parameters -----
  real level;
  real A;
  real<lower=1> P;
  real phi;

  real<lower=0> sd_level;
  real<lower=0> sd_A;
  real<lower=0> sd_P;
  real<lower=0> sd_phi;
  real<lower=0> sd_mu;

  // ----- Individual-level parameters -----
  vector[n] level_i;
  vector[n] A_i;
  vector<lower=1>[n] P_i;
  vector[n] phi_i;

  // ----- Noise parameters -----
  real<lower=0> sd_y;
  real<lower=0> sd_amp;
  real<lower=0> mean_amp;
}

transformed parameters {
  matrix[n, T] mu_it;

  for (i in 1:n)
    for (t in 1:T)
      mu_it[i, t] = mu_t[t];
}

model {
  // ===== Priors =====
  level ~ normal(10, 5);
  A ~ normal(0, 2);
  P ~ normal(50, 10);
  phi ~ normal(0, 5);

  sd_mu_t ~ exponential(1);
  sd_level ~ exponential(1);
  sd_A ~ exponential(1);
  sd_P ~ exponential(1);
  sd_phi ~ exponential(1);
  sd_mu ~ exponential(1);
  sd_y ~ exponential(1);
  sd_amp ~ exponential(1);
  mean_amp ~ normal(1, 0.5);

  // ===== State equation =====
  mu_t[1] ~ normal(0, 5);
  for (t in 2:T)
    mu_t[t] ~ normal(mu_t[t-1], sd_mu_t);

  // ===== Hierarchical structure =====
  level_i ~ normal(level, sd_level);
  A_i ~ normal(A, sd_A);
  P_i ~ normal(P, sd_P);
  phi_i ~ normal(phi, sd_phi);

  // ===== Observation model (mixture) =====
  for (i in 1:n) {
    for (t in 1:T) {
      real season = A_i[i] * sin(2 * pi() * t / P_i[i] + phi_i[i]);
      real mu0 = level_i[i] + mu_it[i, t] + season;

      target += log_mix(
        p,
        normal_lpdf(y[i][t] | mu0 + mean_amp, sd_y),
        normal_lpdf(y[i][t] | mu0, sd_y)
      );
    }
  }
}

generated quantities {
  matrix[n, T] y_rep;

  for (i in 1:n) {
    for (t in 1:T) {
      real season = A_i[i] * sin(2 * pi() * t / P_i[i] + phi_i[i]);
      real mu0 = level_i[i] + mu_it[i, t] + season;
      y_rep[i, t] = normal_rng(mu0, sd_y);
    }
  }
}
 