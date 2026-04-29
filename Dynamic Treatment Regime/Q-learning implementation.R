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

attach(my_dataset$data)
# stage 2 regression 
stage2 = lm(y ~ 1 + o1 + a1 + o1*a1 + o2 + a2*(1  + o1 + a1 + o2))
(coef2 = coefficients(stage2))
beta2_hat_lm = c(b20= coef2[1],b21 = coef2['o1'], b22 = coef2['a1'], b23 = coef2['o1:a1'], b24 = coef2['o2'] ) 
phi2_hat_lm = c(phi20 = coef2['a2'], phi21 = coef2['o1:a2'], 
                phi22 = coef2['a1:a2'], phi23 = coef2['o2:a2'])
# compute stage 1 psuedo-outcome 
Q_2_max = abs(cbind(1,o1, a1, o2) %*% phi2_hat_lm) #a_2 \in {-1,1}
y_1_hat = 0 + cbind(1,o1, a1, o1*a1, o2)%*%beta2_hat_lm + Q_2_max

# stage 1 regression 
stage1 = lm(y_1_hat ~ 1 + o1 + a1*(1 + o1))
(coef1 = coefficients(stage1))
beta1_hat_lm = c(b10 = coef1[1], coef1['o1'])
phi1_hat_lm = c(phi10 = coef1['a1'], phi11 = coef1['o1:a1'])

#identify optimal regimes for stage 1 
d1_opt_lm = ifelse(cbind(1,o1)%*%phi1_hat_lm >=0 ,1,-1) #sign
#identify optimal regimes for stage 2
d2_opt_lm = ifelse(cbind(1,o1, a1, o2) %*% phi2_hat_lm >=0, 1,-1)

#true optimal regime for stage 1 
true_d1_opt = ifelse(cbind(1, o1) %*% phi1 >= 0, 1, -1)
true_d2_opt = ifelse(cbind(1,o1, a1, o2) %*% phi2 >=0, 1,-1)

#Proportion of optimal action 
(mean(d1_opt_lm == true_d1_opt))
(mean(d2_opt_lm == true_d2_opt))

