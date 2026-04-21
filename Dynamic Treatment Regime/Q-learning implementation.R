# simulation 1, setting 6, of Murray et.al. 2018

# true parameters; first half for main effect and second half for treatment effects
beta1 = c(0.644, 0.006, 0.369, 0.019) #stage 1 Q-function 
beta2 = c(0.00, 0.00, -0.50, 0.00, 0.00, 0.25, 0.00, 0.50, 0.50) # stage 2 Q-function
delta = c(0.1, 0.1) # o2 outcome model p(o2 | o1, a1) 

# sample 
set.seed(67);n = 200 # {-1,1} coding for covariates and assignment
o1 = 2* rbinom(n, 1, 0.5) -1 ; a1 = 2* rbinom(n, 1, 0.5) -1 
o2 = 2* rbinom(n, 1, 1/(1 + exp(-cbind(o1, a1) %*% delta))) -1 
a2 = 2 * rbinom(n, 1, 0.5) -1
y = rnorm(n, cbind(1, o1, a1, o1*a1,o2, a2,a2*o1,a2*a1, a2*o2 ) %*% beta2, 1)

my_dataset = list(data = cbind(o1, a1, o2, a2,y), beta1 = beta1, beta2 = beta2, delta = delta)

# check
View(my_dataset$data)
cor(my_dataset$data[,'a2'], my_dataset$data[,'y'])
boxplot(my_dataset$data[,'y'] ~my_dataset$data[,'a1'], 
        xlab = 'a1', 
        ylab = 'final payoff') 
boxplot(my_dataset$data[,'y'] ~my_dataset$data[,'a2'], 
        xlab = 'a2', 
        ylab = 'final payoff') 
head(my_dataset$data)


