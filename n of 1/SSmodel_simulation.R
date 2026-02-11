####true model 
t = seq(1,100) # observation points 

n = 10 # sample size 
x0 = 0.1#initial value for the state model
a0 = 0.8 # the group level slope
b0 = 0.1 # group level intercept

data = data.frame(X = (1:n),x1 = x0 +rnorm(n),
 a = a0 + rnorm(n), b0 = b0 + rnorm(n))
res = data.frame(X = )
x = numeric(100)
x[1] = x0 # initial value
for(t in 1:(length(x)-1)){
    x[t+1] = a* x[t] + b 
} # state equation 
x

# observation equation 
s = 0.1 
m = 4
for(t in 1:length(t)){
    y[t] = x[t] + rnorm(1, m, s)
}

# data simulation 

# model for estimation 

#sample size n = 5 
# n = 10 
# n = 15 
# n = 20 
# n = 50 
# n = 100 