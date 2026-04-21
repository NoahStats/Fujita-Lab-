setwd("//wsl.localhost/Ubuntu/home/moosehunter/R/Fujita Lab/Dynamic Treatment Regime/") #Rstudio
setwd('/home/moosehunter/R/Fujita Lab/Dynamic Treatment Regime/') # VScode

library(DTRreg) #Wallace, Erica Moodie and David Stephens! 
# https://cran.r-project.org/web/packages/DTRreg/refman/DTRreg.html#DTRreg

library(DTRlearn2)
#works for OWL (Outcome Weighted Learning)
# https://cran.r-project.org/web/packages/DTRlearn2/refman/DTRlearn2.html

library(DynTxRegime) #bmiData
# https://cran.r-project.org/web/packages/DynTxRegime/refman/DynTxRegime.html#.newClassificationObj


data(bmiData) #DynTxRegime
data(twoStageCont) #DTRreg
data(adhd) #DTRlearn2

attach(twoStageCont)

# DTRreg for Q-learning, G methods and Dynamic WOLS  by Moodie and Stephens-------------------------------------------------------------------------

# THis estimates the parameters of the BLIP FUNCTION
# using G-estimation, Dynamic Weighted OLS and Q-learning
# 2-stage DTRreg

# A list specifying covariates of the blip functions in order
# X1 for the first point nd X2 for the second point
blip.mod = list( ~X1, ~ X2)  #iinteraction between patients and treatments

# A list specifying the treatment models for each stage in order
treat.mod = list(A1 ~ X1, A2 ~ 1) # Propensity score

# A list specifying the covariates of the treatment-free model
tf.mod = list( ~ X1, ~ X2) # baseline outcome (no treatments)

#G-estimation
mod1 = DTRreg(twoStageCont$Y, blip.mod, treat.mod, tf.mod,
              data = twoStageCont, method = 'gest',var.estim = 'bootstrap') #method = gest

# Q-learning 
mod2 =　DTRreg(twoStageCont$Y, blip.mod, treat.mod, tf.mod, 
               data = twoStageCont, method = 'qlearn') #method = qlearn

# Dynamic WOLS 
mod3 = DTRreg(twoStageCont$Y, blip.mod, treat.mod, tf.mod, 
              data = twoStageCont, method = 'dwols') #method = dwols



# DTRlearn2 ---------------------------------------------------------------

data(adhd) 
n = length(adhd$a1) # sample size

View(cbind(adhd$o11, adhd$o12, adhd$o13, adhd$o14))
H1 = scale(cbind(adhd$o11, adhd$o12, adhd$o13, adhd$o14))# standardizing  0-1 data??
# history at stage 1 or initial covariates 
View(H1)


H2 = as.data.frame(scale(cbind(H1, adhd$a1, H1 * adhd$a1, adhd$r, adhd$o22, adhd$r * adhd$a1,
                 adhd$o22 * adhd$a1)))#hstory at stage2
View(H2)
colnames(H2)[12] = "r*a1"
colnames(H2)[13] = "o22*a1"
# r for response to the treatment: 0 or 1 (Interim payoff)

# q learning
fit_ql = ql(H = list(H1, H2), AA = list(adhd$a1, adhd$a2), 
            R = list(rep(0,n),y), pi = list(rep(0.5, n), rep(0.5, n)), 
            k = 2, m = 3, lasso = TRUE) ##????  Lasso for Q learning 
