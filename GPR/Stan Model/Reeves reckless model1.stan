data {
  int<lower=1> N;
  array[N] real x;
  vector[N] y;
  real<lower=0> p;
}

transformed data {
  vector[N] mu = rep_vector(0, N);
}

parameters {
  real<lower=0> alpha_trend;
  real<lower=0> rho_trend;

  real<lower=0> alpha_per;
  real<lower=0> rho_per;

  real<lower=0> sigma;
}

model {
  matrix[N, N] K;
  matrix[N, N] L_K;

  K =
    gp_matern32_cov(x, alpha_trend, rho_trend)
    + gp_periodic_cov(x, alpha_per, rho_per, p);

  for (n in 1:N)
    K[n, n] += square(sigma) + 1e-6;

  L_K = cholesky_decompose(K);

  alpha_trend ~ normal(0, 1);
  rho_trend   ~ lognormal(0, 1);
  alpha_per   ~ normal(0, 0.5);
  rho_per     ~ lognormal(0, 0.5);
  sigma       ~ normal(0, 0.3);

  y ~ multi_normal_cholesky(mu, L_K);
}
