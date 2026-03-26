data {
  int<lower=1> N;
  int<lower=2> T;
  matrix[N, T] x;
  matrix[N, T] y;
}

parameters {
  real beta_x;
  real beta_y;
  real delta;

  vector[T-1] gamma;

  real<lower=0> eta;
  real<lower=0> rho;

  real<lower=0> sigma_x;
  real<lower=0> sigma_y;
}

transformed parameters {
  matrix[T-1, T-1] K;

  for (i in 1:(T-1)) {
    for (j in 1:(T-1)) {
      real d = i - j;
      K[i,j] = eta^2 * exp( - square(d) / (2 * square(rho)) );
    }
  }
}

model {
  // Priors
  beta_x ~ normal(0,1);
  beta_y ~ normal(0,1);
  delta  ~ normal(0,1);

  eta ~ exponential(1);
  rho ~ exponential(1);

  sigma_x ~ exponential(1);
  sigma_y ~ exponential(1);

  gamma ~ multi_normal(rep_vector(0, T-1), K);

  for (i in 1:N) {
    for (t in 2:T) {

      x[i,t] ~ normal(
        beta_x * x[i,t-1] +
        delta  * y[i,t-1],
        sigma_x
      );

      y[i,t] ~ normal(
        beta_y * y[i,t-1] +
        gamma[t-1] * x[i,t-1],
        sigma_y
      );
    }
  }
}
