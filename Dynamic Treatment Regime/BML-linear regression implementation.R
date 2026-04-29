# Based on simulation 1, setting 6, of Murray et.al. 2018
library(MASS)
library(mvtnorm)
# true parameters; first half for main effect and second half for treatment effects
beta1 = c(0.644, 0.006, -0.369, 0.019) #stage 1 Q-function 
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
boxplot(my_dataset$data[,'y'] ~my_dataset$data[,'a2'], 
        xlab = 'a2', 
        ylab = 'final payoff') 
# BIG sampler with g prior-------------------------------------------------------------
p1 = 4; p2 = 9 #dims of the parameters 
S = 2000
BETA1  = matrix(NA, nrow = S, ncol = p1)
SIGMA1 = rep(NA, S)
BETA2  = matrix(NA, nrow = S, ncol = p2)
SIGMA2 = rep(NA, S)
X1 = model.matrix(~ o1 + a1 + o1*a1)
X2 = model.matrix(~ o1 + a1 + o1*a1 + o2 + a2*(1 + o1 + a1 + o2))
#use model.matrix when constructing design matrices! 
treat_cols2 = grep('a2', colnames(X2)) # store the col names of treatment effects
treat_cols1 = grep('a1', colnames(X1))
main_cols2 = grep('a2', colnames(X2), invert = TRUE)
main_cols1 = grep('a1', colnames(X1), invert = TRUE)
fit2 = lm(y ~ o1 + a1 + o1:a1 + o2 + a2 + a2:o1 + a2:a1 + a2:o2)
Beta_2_ols  = coef(fit2)
sigma2_2_ols = summary(fit2)$sigma^2
SSRg_2 = as.numeric(t(y) %*% (diag(n) - n/(n+1) * (X2 %*% solve(t(X2)%*%X2)) 
                              %*% t(X2)) %*%y)

for(s in 1:S){
  #1. sample THETA.2 from stage 2 posterior
  sigma2_2 = 1/rgamma(1,shape = (n+1)/2,
                      rate = (sigma2_2_ols + SSRg_2)/2) 
  Beta_2 = mvrnorm(1, mu = n/(n+1) * Beta_2_ols, 
                   Sigma = n/(n+1) * sigma2_2 * solve(t(X2)%*%X2))
  BETA2[s, ] = Beta_2
  SIGMA2[s] = sigma2_2
  
  interaction = cbind(1, o1, a1, o2) %*% 
    Beta_2[treat_cols2]
  a2_opt = ifelse(interaction >= 0, 1, -1)
  
  #2. Step2: generate a counterfactual dataset 
  X2_new =  cbind(1,o1,a1,o1*a1,o2,
                  a2_opt, a2_opt*o1, a2_opt*a1, a2_opt*o2)
  y_mis = rmvt(1, df = n+1, delta = (n)/(n+1) *X2_new %*% Beta_2_ols, 
               sigma = (sigma2_2_ols + SSRg_2)/(n+1)*
                 (diag(n) + X2_new %*% solve(t(X2)%*%X2) %*% t(X2_new)),
               type = 'shifted'
  )
  y_counterfactual = y 
  for(i in 1:n){
    if(a2[i] != a2_opt[i]){
      y_counterfactual[i] = y_mis[i]
    }
  }
  
  #3 Step3: sample from stage 1 posterior using y2_counterfactual
  fit1 = lm(y_counterfactual ~ o1 + a1 + o1:a1)
  Beta_1_ols  = coef(fit1)
  sigma2_1_ols = summary(fit1)$sigma^2
  SSRg_1 = t(y_counterfactual) %*% 
    (diag(n) - n/(n+1) * (X1 %*% solve(t(X1)%*%X1)) %*% t(X1)) %*%
    y_counterfactual
  
  sigma2_1 = 1/rgamma(1,shape = (n+1)/2,
                      rate = (sigma2_1_ols + SSRg_1)/2) 
  
  Beta_1 = mvrnorm(1, mu = n/(n+1) * Beta_1_ols, 
                   Sigma = n/(n+1) * sigma2_1 * solve(t(X1)%*%X1))
  BETA1[s, ] = Beta_1
  SIGMA1[s] = sigma2_1
}


BETA2_EAP = colMeans(BETA2);BETA1_EAP = colMeans(BETA1)
a2_opt_hat = ifelse(cbind(1, o1, a1, o2)%*%
                      BETA2_EAP[treat_cols2]>=0, 
                    1,-1)
a1_opt_hat = ifelse(cbind(1, o1)%*%
                      BETA1_EAP[treat_cols1]>=0, 
                    1,-1)
#true optimal regime for stage 1 
true_d1_opt = ifelse(cbind(1, o1) %*% beta1[3:4] >= 0, 1, -1)
true_d2_opt = ifelse(cbind(1,o1, a1, o2) %*% beta2[6:9] >=0, 1,-1) 
mean(a2_opt_hat == true_d2_opt)
mean(a1_opt_hat == true_d1_opt)

# probability of optimal treatments for individual i
i = 30
input2_i = X2[i,treat_cols2]/X2[i,'a2'] #input
input1_i = X1[i,treat_cols1]/X1[i,'a1']
treat2_i = BETA2[,treat_cols2]%*%input2_i 
treat1_i = BETA1[,treat_cols1]%*%input1_i 
(prob_one_is_better_for_a2 = mean(treat2_i>=0))
(prob_one_is_better_for_a1 = mean(treat1_i >=0))
 