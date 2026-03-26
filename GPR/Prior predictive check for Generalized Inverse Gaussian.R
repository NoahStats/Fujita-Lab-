setwd("//wsl.localhost/Ubuntu/home/moosehunter/R/Fujita Lab/GPR")
library(GPBayes) #BesselK function
library(GIGrvg) #random number generator for GIG
library(MASS)
set.seed(123)

# The shape parameter denoted by lambda or p takes integers. 
# p > 0 for light right tail; p< 0 for heavy right tail; p = 0 for symmetry
#a,b are real numbers 
# large a shrinks the distributions toward zero
# large b avoids near zero 
# since we want to limit the length scale to high frequency (small length scales),
# we need p>0, large a and small b 

#list of hyperparameters ro run grid search on
p = 1:5
a = seq(0,5,length.out = 10)
b = seq(0,5, length.out= 10)

OU_amp = mvrnorm(1,mu = rep(0,5*10*10),Sigma = diag(5*10*10)) #empty vector to store amp hyperparameters from normal
Per_amp = mvrnorm(1,mu = rep(0,5*10*10),Sigma = diag(5*10*10))

#use mean and mode of GIG to find appropriate sets of hyperparameters 
GIG_mode = function(p,a,b){
  
}
for(i in p){
  for(j in a){
    for(k in b){
  
    }
  }
}


GIG = function(x,p,a,b){
  (a/b)^(p/2) * (1/(2*BesselK(p,sqrt(a*b)))) * x^(p-1) * exp(-(a*x + b/x)/2)
}

x = seq(0,3,length.out = 1000)
plot(x,GIG(x,5,0.5,0.5))



