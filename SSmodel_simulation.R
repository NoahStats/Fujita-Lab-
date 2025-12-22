####true model 
t = seq(1,100) # observations points 

x0 = 0.1#initial value for the latent variable
a = 0.5 # inital values for the state equation parameters 
b = 0.1 

x = numeric(100)
x[1] = x0 # initial value
for(t in 1:(length(x)-1)){
    x[t+1] = a* x[t] + b 
} # state equation 

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