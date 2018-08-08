
################################################################
################### Lambda Sq EM        ########################
################### this cannot work for p > n #################
################################################################

BlassoemR <- function(X=X,Y=Y,niter,nburn,nrep,nthin){


library(Rcpp) ; library(RcppArmadillo);
sourceCpp("C:\\Users\\chetkar\\Workspace\\BioInformatics\\BLassoem\\Blassoem.cpp");

#library(bayesm)
n <- length(Y)
p <- dim(X)[2]

#### data scaling

J <- matrix (1, nrow=n, ncol=1)
X  <- (X - J%*%colMeans(X))/(J%*%apply(X,2,sd))
XX <- t(X)%*%X    
Y.til <- Y - mean(Y)

#### starting values 

ini <- lm(Y.til ~ X) # use ridge for better
res <-  ini$residuals
b.est <- ini$coefficients

#### starting values 

lambda <- p*sqrt(sum(res^2)/(n-p))/sum(abs(b.est))
lambda.sq <- lambda^2
sig.em <- runif(1,0.1,10)
tau.em <- rexp(p, rate=lambda.sq/2)
tau.sq  <- rexp(p, rate=lambda.sq/2)
mean.be <- rep(0,p)
cov.be <- sig.em*diag(tau.em)
beta.em.p <- rmvnorm(1,mean=mean.be,sigma=cov.be)
sigsc_hyper=0.0;
sigsh_hyper =0.0;


emo <-BLassoem( X,Y,lambda.sq,sig.em,sigsc_hyper,sigsh_hyper,tau.em,tau.sq,beta.em.p, niter,nrep);
 
return(emo);
}