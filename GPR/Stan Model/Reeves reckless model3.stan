
// model2 : matern 1/2 for trend, periodic kernel and gaussian noise
functions{
  matrix daily_periodic_kernel(real[] x, real alpha, real rho, real p) {
    int N = size(x);
    matrix[N, N] K;
    real sq_rho = square(rho);
    real inv_2sq_rho = -1.0 / (2.0 * sq_rho); // Your specific -1/2 multiplier
    
    for (i in 1:N) {
      for (j in i:N) {
        real delta = fabs(x[i] - x[j]);
        // Your specific formula
        K[i, j] = square(alpha) * exp(inv_2sq_rho * pow(sin(pi() * delta / p), 2));
        K[j, i] = K[i, j]; 
      }
    }
    return K;
  }
  }
  
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
    gp_exponential_cov(x, alpha_trend, rho_trend)
    + daily_periodic_kernel(x, alpha_per, rho_per, p);
//gp_periodic_cov uses different parametrizations from 
// Rasumussen! 
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
