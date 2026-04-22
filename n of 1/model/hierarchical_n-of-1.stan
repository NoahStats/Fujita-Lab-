data {
  int<lower=1> N;          // number of individuals
  int<lower=1> T;          // time points
  matrix[N, T] y;

  vector[T] time;
  vector[T] I_treat;
  vector[T] t_post;
}

parameters {
  // ---- Population-level means ----
  real beta_0_mu;
  real beta_1_mu;
  real beta_2_mu;
  real beta_3_mu;

  // ---- Population-level SDs ----
  real<lower=0> beta_0_sd;
  real<lower=0> beta_1_sd;
  real<lower=0> beta_2_sd;
  real<lower=0> beta_3_sd;

  // ---- Individual-level (non-centered) ----
  vector[N] beta_0_raw;
  vector[N] beta_1_raw;
  vector[N] beta_2_raw;
  vector[N] beta_3_raw;

  // ---- AR(1) hierarchy (NORMAL, NOT TANH) ----
  real phi_mu;                 
  real<lower=0> phi_sd;        
  vector[N] phi_raw;           

  // ---- Noise ----
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] beta_0 = beta_0_mu + beta_0_sd * beta_0_raw;
  vector[N] beta_1 = beta_1_mu + beta_1_sd * beta_1_raw;
  vector[N] beta_2 = beta_2_mu + beta_2_sd * beta_2_raw;
  vector[N] beta_3 = beta_3_mu + beta_3_sd * beta_3_raw;

  // DIRECT normal, no transformation
  vector[N] phi = phi_mu + phi_sd * phi_raw;
}

model {
  // ---- Priors ----
  beta_0_mu ~ normal(0, 5);
  beta_1_mu ~ normal(0, 1);
  beta_2_mu ~ normal(0, 5);
  beta_3_mu ~ normal(0, 1);

  beta_0_sd ~ exponential(1);
  beta_1_sd ~ exponential(1);
  beta_2_sd ~ exponential(1);
  beta_3_sd ~ exponential(1);

  beta_0_raw ~ std_normal();
  beta_1_raw ~ std_normal();
  beta_2_raw ~ std_normal();
  beta_3_raw ~ std_normal();

  phi_mu  ~ normal(0, 0.5);
  phi_sd  ~ exponential(1);
  phi_raw ~ std_normal();

  sigma ~ exponential(1);

  // ---- Likelihood ----
  for (i in 1:N) {
    vector[T] mu_i;

    mu_i =
      beta_0[i]
      + beta_1[i] * time
      + beta_2[i] * I_treat
      + beta_3[i] * t_post;

    // ---- Stationary initialization (matches R) ----
    y[i, 1] ~ normal(
      mu_i[1] / (1 - phi[i]),
      sigma / sqrt(1 - square(phi[i]))
    );

    // ---- UNCENTERED AR(1) recursion ----
    for (t in 2:T) {
      y[i, t] ~ normal(
        mu_i[t] + phi[i] * y[i, t - 1],
        sigma
      );
    }
  }
}
