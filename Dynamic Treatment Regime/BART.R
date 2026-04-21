setwd("//wsl.localhost/Ubuntu/home/moosehunter/R/Fujita Lab/Dynamic Treatment Regime/") #Rstudio
setwd('/home/moosehunter/R/Fujita Lab/Dynamic Treatment Regime/') # VScode
library(BayesTree)
set.seed(1)


# Bart from BayesTree -----------------------------------------------------

f = function(x){
  10*sin(pi*x[,1]*x[,2] + 20*(x[,3] -.5)^2 +
           10*x[,4] + 5*x[,5])
}
sigma = 1.0 # y = f(x) + sigma*z, z ~ N(0,1)
n = 100
x = matrix(runif(n*10), n, 10)

Ey = f(x) # f(x) without error
y = Ey + sigma*rnorm(n) #observed data

lmFit = lm(y ~., data.frame(x,y))

## run BART
bartFit = bart(x,y,ndpost = 200) 
#ntree = 200 by default
#ndpost = number of posterior draws. It is 1000 by default

plot(bartFit)

fitmat = cbind(y, Ey, lmFit$fitted.values, 
               bartFit$yhat.train.mean)
colnames(fitmat) = c('y','Ey','lm','bart')
print(cor(fitmat))


# partial dependence plot -------------------------------------------------
f = function(x) { return(.5 * x[,1] + 2 * x[,2] * x[,3])}

sigma = .2 # y = f(x) + sigma*x
n = 100 #number of observations

x = matrix(2*runif(n*3)-1,ncol=3) ; colnames(x) = c('rob','hugh','ed')
Ey = f(x)
y = Ey +  sigma*rnorm(n)
lmFit = lm(y~.,data.frame(x,y))
par(mfrow=c(1,3))

##pdbart: one dimensional partial dependence plot
pdb1 = pdbart(x,y,xind=c(1,2),
              levs=list(seq(-1,1,.2),seq(-1,1,.2)),pl=FALSE,
              keepevery=10,ntree=100,nskip=100,ndpost=200) #should run longer!
plot(pdb1,ylim=c(-.6,.6))
##pd2bart: two dimensional partial dependence plot
set.seed(99)

pdb2 = pd2bart(x,y,xind=c(2,3),
               levquants=c(.05,.1,.25,.5,.75,.9,.95),pl=FALSE,
               ntree=100,keepevery=10,verbose=FALSE,nskip=100,ndpost=200) #should run longer!
plot(pdb2)
##compare BART fit to linear model and truth = Ey
fitmat = cbind(y,Ey,lmFit$fitted,pdb1$yhat.train.mean)
colnames(fitmat) = c('y','Ey','lm','bart')
print(cor(fitmat))
## plot.bart(pdb1) displays the BART run used to get the plot.
