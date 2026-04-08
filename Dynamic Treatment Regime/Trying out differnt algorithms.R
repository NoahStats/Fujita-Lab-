setwd("//wsl.localhost/Ubuntu/home/moosehunter/R/Fujita Lab/Dynamic Treatment Regime/") #Rstudio
setwd('/home/moosehunter/R/Fujita Lab/Dynamic Treatment Regime/') # VScode

library(DTRreg) #Wallace, Erica Moodie and David Stephens! 
# https://cran.r-project.org/web/packages/DTRreg/refman/DTRreg.html#DTRreg

library(DynTxRegime) #bmiData
# https://cran.r-project.org/web/packages/DynTxRegime/refman/DynTxRegime.html#.newClassificationObj


data(bmiData) #DynTxRegime
data(twoStageCont)


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
