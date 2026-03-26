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
  vector[T] time_points;
  vector[T] Y;
}

transformed data{
//Define mean functions and kernels in this block
  vector[T] D;  // intervention indicator
  vector[T] T_post // time since intervention
  
  real cutpoint = 14; 
  
  for(t in T ){
    D[t] = time_points[t] > cutpoint; 
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
  real OU_amp <lower = 0>; 
  real OU_length_scale<lower = 0>; 
  real Per_amp<lower = 0>; 
  real Per_length_scale<lower = 0> ; 
  real sigma<lower = 0>; 
}


model{
//define the kernel here
matrix[T,T] L_K;
matrix[T,T] K = gp_exponental_cov(time_points,OU_amp,OU_length_scale)+ 
  gp_periodic_cov(time_points,Per_amp, Per_length_scale,1);
real sq_sigma = square(sigma);

//diagonal elements 
for(t in 1:T){
  K[t,t] = K[t,t] + sq_sigma;
}

L_K = cholesky_decompose(K);

// Y = a * 1(t<70) + b * t * 1(t<70) 
//    + c * 1(t>=70) + d * t * 1(t>=70)
 //    + f1 + f2 + epsilon
// f1 ~ GP(0,OU kernel), f2 ~ GP(0, periodic kernel)
// epsilon ~ N(0,\sigma^2)

Y ~ multi_normal_cholesky(intercept + slope * time_points + 
  intercept_treat * D + slope_treat * T_post * D, L_K );

// prior for amplitude;scale ~ normal(0,1) or t! 
// prior for the error sigma ; sigma ~ normal(0,1)
sigma ~ std_normal();
OU_amp ~ std_normal();
Per_amp ~ std_normal();
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
OU_length_scale ~ generalized_inverse_gaussian(p, a,b);
Per_length_scale ~ geralized_inverse_gaussan(p,a,b);

}
