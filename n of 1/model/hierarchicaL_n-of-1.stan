data {
  int<lower=1> N;
  int<lower=1> T;
  matrix[N, T] y;
  vector[T] time;
  vector[T] I_treat;
  vector[T] t_post;
}

parameters {
  // population means
  real beta_0_mu;
  real beta_1_mu;
  real beta_2_mu;
  real beta_3_mu;
  real phi_mu;

  // population SDs
  real<lower=0> beta_0_sd;
  real<lower=0> beta_1_sd;
  real<lower=0> beta_2_sd;
  real<lower=0> beta_3_sd;
  real<lower=0> phi_sd;

  // individual parameters
  vector[N] beta_0;
  vector[N] beta_1;
  vector[N] beta_2;
  vector[N] beta_3;
  vector[N] phi_raw;

  real<lower=0> sigma;
}

transformed parameters {
  vector[N] phi;
  phi = tanh(phi_raw);
}

model {
  // priors
  beta_0_mu ~ normal(0, 5);
  beta_1_mu ~ normal(0, 5);
  beta_2_mu ~ normal(0, 5);
  beta_3_mu ~ normal(0, 5);
  phi_mu    ~ normal(0, 1);

  beta_0_sd ~ exponential(1);
  beta_1_sd ~ exponential(1);
  beta_2_sd ~ exponential(1);
  beta_3_sd ~ exponential(1);
  phi_sd    ~ exponential(1);

  beta_0 ~ normal(beta_0_mu, beta_0_sd);
  beta_1 ~ normal(beta_1_mu, beta_1_sd);
  beta_2 ~ normal(beta_2_mu, beta_2_sd);
  beta_3 ~ normal(beta_3_mu, beta_3_sd);

  phi_raw ~ normal(phi_mu, phi_sd);
  sigma   ~ exponential(1);

  // likelihood
  for (i in 1:N) {
    real mu_prev;

    mu_prev =
      beta_0[i] +
      beta_1[i] * time[1] +
      beta_2[i] * I_treat[1] +
      beta_3[i] * t_post[1];

    y[i,1] ~ normal(mu_prev / (1 - phi[i]),
                    sigma / sqrt(1 - square(phi[i])));

    for (t in 2:T) {
      real mu_t;
      mu_t =
        beta_0[i] +
        beta_1[i] * time[t] +
        beta_2[i] * I_treat[t] +
        beta_3[i] * t_post[t];

      y[i,t] ~ normal(mu_t + phi[i] * y[i,t-1], sigma);
    }
  }
}
