//model4; same as model 2 and 3 but with priors on the 
// intercept and the slope
// parameters are different as well 
//  intercept ~ std_normal();
//  intercept_treat ~ std_normal();
//  slope ~ normal(0, 0.07); // given the standardization
//  slope_treat ~ std_normal(0, 0.07);
//#OU_amp = 0.8; OU_scale = 0.5 # hyparparameters for OU
//#intercept = 0; slope = 0 # basline intercept and slope 
//#intercept_treat  = 0.15; slope_treat = 0.005 # treatment effects
//#sigma = 0.1
// OU_length_scale ~ generalized_inverse_gaussian(2, 4,3);

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
//Define mean functions and kernels in this block
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
//linear terms; intercept, slope, intercept_treat, slope_treat
//OU terms; OU_amp, OU_scale
//Periodic terms; Per_amp, Per_scale
  real intercept; 
  real slope; 
  real intercept_treat; 
  real slope_treat;
  real<lower = 0> OU_amp ; 
  real<lower = 0> OU_length_scale ; 
  real<lower = 0> sigma ; 
}


model{
//define the kernel here
  matrix[T,T] K;
  matrix[T,T] L_K;
  vector[T] mu;
  real sq_sigma = square(sigma);

  K = gp_exponential_cov(time_points,OU_amp,OU_length_scale);

  for(t in 1:T){
    K[t,t] = K[t,t] + sq_sigma;
  }

  L_K = cholesky_decompose(K);

// Y = a * 1(t<70) + b * t * 1(t<70) 
//    + c * 1(t>=70) + d * t * 1(t>=70)
 //    + f1 + f2 + epsilon
// f1 ~ GP(0,OU kernel), f2 ~ GP(0, periodic kernel)
// epsilon ~ N(0,\sigma^2)
  mu = intercept + slope * time_vec + 
       intercept_treat * D + slope_treat * T_post .* D;
       
  Y ~ multi_normal_cholesky(mu, L_K );

//prior for the slope and the intercept
//since the data is standaardized, 
// slope of 2SD is given by (2.0/14 days) = 0.14 
// 1SD = 0.07
  intercept ~ std_normal();
  intercept_treat ~ std_normal();
  slope ~ normal(0, 0.07); // given the standardization
  slope_treat ~ normal(0, 0.07);
// prior for amplitude;scale ~ normal(0,1) or t! 
// prior for the error sigma ; sigma ~ normal(0,1)
  sigma ~ std_normal();
  OU_amp ~ std_normal();
// prior for length-scale 
//amp ~ invGamma(5,5),soft right tailed,works when is it pure GP 

//library(invgamma)
//x <- seq(0.01, 10, length.out = 1000)  
//plot(x, dinvgamma(x, shape = 5, rate = 5), type = "l")

// but doesn't work when GP is a component because 
// large length-scale captures low frequnecy, almost linear, variations
// and results in non-identifiability
// It is recommended to assign a prior with mass on the left tail (near zero) 
// to capture high frequency variation
// let the mean capture the low frequency
  OU_length_scale ~ generalized_inverse_gaussian(2, 4,3);
}
