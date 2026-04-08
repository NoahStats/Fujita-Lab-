//recovery3;  assume positive treatment effects b2 and b3

//MODEL
//Y_t ~ GP(mu_t, k(t,t') + sigma^2 I) 0<=t<=28
// mu_t = b0 + b1*t + b2*1(t>14) + b3*1(t>14)*(t-14)
// k(t,t') = alpha^2 * exp(-|t-t'|/l)
// 28 days total,5 times a day, 14 days pre and post treatment respectively

// PRIORS
//  intercept ~ std_normal(); given the standardization by the baseline
//  intercept_treat ~ std_normal();  given the standardization by the baseline
//  slope ~ normal(0, 1/14);  given the standardization by the baseline
//  slope_treat ~ normal(0, 1/14);given the standardization by the baseline
// alpha ~ std_normal(); recommendation from the Stan webpage
// l ~ generalized_inverse_gaussian(2, 4,3); Prior predictive check 
// sigma ~ std_normal() recommendation from the Stan webpage


functions {
  //define functions above the data block  
  real generalized_inverse_gaussian_lpdf(real x, int p,
                                        real a, real b) {
    return p * 0.5 * log(a / b)
      - log(2 * modified_bessel_second_kind(p, sqrt(a * b)))
      + (p - 1) * log(x)
      - (a * x + b / x) * 0.5;
 }
}


data{
  int T;
  array[T] real time_points;
  vector[T] Y;
}


transformed data{
  vector[T] D;  // intervention indicator
  vector[T] T_post; // time since intervention
  vector[T] time_vec;
  
  real cutpoint = 14;
  
  for(t in 1:T ){
    time_vec[t] = time_points[t];
    D[t] = time_points[t] > cutpoint ? 1.0 : 0.0;
    T_post[t] = fmax(0, time_points[t] - cutpoint);
  }
}

parameters{
//linear terms; b0, b1, b2 and b3
//OU terms; alpha, l
  real b0; 
  real b1; 
  real<lower = 0> b2; 
  real<lower = 0> b3;
  real<lower = 0> alpha ; 
  real<lower = 0> l ; 
  real<lower = 0> sigma ; 
}


model{
//define the kernel here!
  matrix[T,T] K;
  matrix[T,T] L_K;
  vector[T] mu;
  real sq_sigma = square(sigma);

  K = gp_exponential_cov(time_points,alpha,l);

  for(t in 1:T){
    K[t,t] = K[t,t] + sq_sigma;
  }

  L_K = cholesky_decompose(K);

  mu = b0 + b1 * time_vec + 
       b2 * D + b3 * T_post .* D;
       
  Y ~ multi_normal_cholesky(mu, L_K );

//priors for the slope and the intercept
//The data is standaardized by the mean and sd of the baseline data
// 1SD = 1/14
  b0 ~ std_normal();
  b2 ~ std_normal();
  b1 ~ normal(0, 1.0/14); // given the standardization
  b3 ~ normal(0, 1.0/14);
  
  alpha ~ std_normal();// prior for amplitude;scale ~ normal(0,1) or t! 
  l ~ generalized_inverse_gaussian(2, 4,3); // prior for length-scale 
  sigma ~ std_normal();// prior for the error sigma ; sigma ~ normal(0,1)

// large length-scale captures low frequnecy, almost linear, variations
// and results in non-identifiability
// It is recommended to assign a prior with mass on the left tail (near zero) 
// to capture only high frequency variations
// let the mean capture the low frequency
  
}
