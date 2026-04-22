setwd("//wsl.localhost/Ubuntu/home/moosehunter/R/Fujita Lab/GPR")
library(GPBayes) #BesselK function
library(GIGrvg) #random number generator for GIG
library(GeneralizedHyperbolic) # quantile function for GIG
library(MASS)
set.seed(123)

# The shape parameter denoted by lambda or p takes integers. 
# p > 0 for light right tail; p< 0 for heavy right tail; p = 0 for symmetry
#a,b are real numbers 
# large a shrinks the distributions toward zero
# large b avoids near zero 
# since we want to limit the length scale to high frequency (small length scales),
# we need p>0, large a and small b 


# OU_length_scale <= 3 based on prior predictive checking
# I tried OU_length_scale \in [0,10] and below 3 works best
# p =2, a = 4, b = 3

# rep_length_scale must be smaller than 0.1 at least! 
# I tried rep-lenth_scla from 10 to 0.0002

# It is not exactly prior predictive check because I didn't 
# sample from the prior everytime, 
# but I set the hyperparameters manually and generated 
# sampled from prior GP and checked it visually 
# 


library(GeneralizedHyperbolic)
p_vals = 1:5
a_vals = seq(0.1, 5, length.out = 20)
b_vals = seq(0.1, 5, length.out = 20)

GIG = function(x, p, a, b){
  
  (a/b)^(p/2) * (1/(2*besselK(p, sqrt(a*b)))) * x^(p-1) * exp(-(a*x + b/x)/2)
  
} 

x = seq(0, 5, length.out = 1000)

plot(x,GIG(x,p = 2, a = 4, b=3    ))
#p <=4 is required to ensure light right tail




