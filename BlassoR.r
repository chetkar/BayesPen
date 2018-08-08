
################################################################
################### Lambda Sq EM        ########################
################### this cannot work for p > n #################
################################################################

library(Rcpp) ; library(RcppArmadillo);


	
BlassoemR <- function(X=X,Y=Y,sigsh_hyper,sigsc_hyper,niter,nburn,nrep,nthin){


sourceCpp("Blassoem.cpp");


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
#sig.em <- runif(1,0.1,10)
sig.em  <- sum(res^2)/n
#tau.em <- rexp(p, rate=lambda.sq/2)
#tau.sq  <- rexp(p, rate=lambda.sq/2)
#mean.be <- rep(0,p)
#cov.be <- sig.em*diag(tau.em)
#beta.em.p <- rmvnorm(1,mean=mean.be,sigma=cov.be)
beta.em.p <-solve(t(X)%*%X,t(X)%*%Y)
tau.em  <- tau.sq   <-beta.em.p^2

Blassoem <-BLassoem( X,Y,lambda.sq,sig.em,sigsh_hyper,sigsc_hyper,tau.em,tau.sq,beta.em.p, niter,nrep);
 
beta.em     <- Blassoem$Betahist
sigsq.em    <- Blassoem$sigma2hist
tausq.em  <- Blassoem$Tausqhist
lambda.em   <- sqrt(Blassoem$Lambda2hist)
 
burnin <- nburn
 
### Inference 

beta.post.mean <- apply(beta.em[-(1:burnin),],2,mean)
beta.post.median <- apply(beta.em[-(1:burnin),],2,median)
beta.post.sd <- apply(beta.em[-(1:burnin),],2,sd)

sigsq.post.mean <- mean(sigsq.em[-(1:burnin)])
sigsq.post.median <- median(sigsq.em[-(1:burnin)])
sigsq.post.sd <- sd(sigsq.em[-(1:burnin)])

tausq.post.mean <- apply(tausq.em[-(1:burnin),],2,mean)
tausq.post.median <- apply(tausq.em[-(1:burnin),],2,median)
tausq.post.sd <- apply(tausq.em[-(1:burnin),],2,sd)

lambda.post.mean <- mean(sqrt(lambda.em[-(1:burnin)]))
lambda.post.median <- median(sqrt(lambda.em[-(1:burnin)]))
lambda.post.sd <- sd(sqrt(lambda.em[-(1:burnin)]))

return.list <- list(beta.post.mean=beta.post.mean,beta.post.median=beta.post.median,beta.post.sd=beta.post.sd,sigsq.post.mean=sigsq.post.mean,
sigsq.post.median=sigsq.post.median,sigsq.post.sd=sigsq.post.sd,tausq.post.mean=tausq.post.mean,tausq.post.median=tausq.post.median,tausq.post.sd=tausq.post.sd,lambda.post.mean=lambda.post.mean,lambda.post.median=lambda.post.median,lambda.post.sd=lambda.post.sd,beta.MCMC.store = beta.em, sigsq.MCMC.store = sigsq.em, tausq.MCMC.store = tausq.em,lambda.MCMC.store = sqrt(lambda.em))

return(return.list)

}
	

################################################################
################### Lambda Sq         ##########################
################################################################
################################################################

	
Blassofixedlambda <- function(X=X,Y=Y,lambda,sigsc_hyper,sigsh_hyper,niter,nburn,nthin){


sourceCpp("Blassofixedlambda.cpp");

#library(bayesm)
n <- length(Y)
p <- dim(X)[2]

#### data scaling

J <- matrix (1, nrow=n, ncol=1)
X  <- (X - J%*%colMeans(X))/(J%*%apply(X,2,sd))
XX <- t(X)%*%X    
Y.til <- Y - mean(Y)

BLassofixedlambda <- BLassofixedlambda(X,Y,lambda,sigsc_hyper,sigsh_hyper,niter,nburn,nthin)

beta.p      <-BLassofixedlambda$Betahist
sigsq.post  <-Blassofixedlambda$sigma2hist
tausq.post  <-Blassofixedlambda$Tausqhist

burnin =nburn

### Inference 

beta.post.mean <- apply(beta.p[-(1:burnin),],2,mean)
beta.post.median <- apply(beta.p[-(1:burnin),],2,median)
beta.post.sd <- apply(beta.p[-(1:burnin),],2,sd)

sigsq.post.mean <- mean(sigsq.post[-(1:burnin)])
sigsq.post.median <- median(sigsq.post[-(1:burnin)])
sigsq.post.sd <- sd(sigsq.post[-(1:burnin)])

tausq.post.mean <- apply(tausq.post[-(1:burnin),],2,mean)
tausq.post.median <- apply(tausq.post[-(1:burnin),],2,median)
tausq.post.sd <- apply(tausq.post[-(1:burnin),],2,sd)

return.list <- list(beta.post.mean=beta.post.mean,beta.post.median=beta.post.median,beta.post.sd=beta.post.sd,sigsq.post.mean=sigsq.post.mean,
sigsq.post.median=sigsq.post.median,sigsq.post.sd=sigsq.post.sd,tausq.post.mean=tausq.post.mean,tausq.post.median=tausq.post.median,tausq.post.sd=tausq.post.sd, beta.MCMC.store = beta.p, sigsq.MCMC.store = sigsq.post, tausq.MCMC.store = tausq.post)

return(return.list)

}


Blassopriorlambda <- function(X=X,Y=Y,a,b,sigsc_hyper,sigsh_hyper,niter,nburn,nthin){


sourceCpp("BLassopriorlambda.cpp");

#library(bayesm)
n <- length(Y)
p <- dim(X)[2]

#### data scaling

J <- matrix (1, nrow=n, ncol=1)
X  <- (X - J%*%colMeans(X))/(J%*%apply(X,2,sd))
XX <- t(X)%*%X    
Y.til <- Y - mean(Y)

BLassopriorlambda <- BLassopriorlambda(X,Y,a,b,sigsh_hyper,sigsc_hyper,niter)

beta.post      <-BLassopriorlambda$Betahist
sigsq.post  <-BLassopriorlambda$sigma2hist
tausq.post  <-BLassopriorlambda$Tausqhist
lambda.post   <- sqrt(BLassopriorlambda$Lambda2hist)
 
burnin =nburn

### Inference 

beta.post.mean <- apply(beta.post[-(1:burnin),],2,mean)
beta.post.median <- apply(beta.post[-(1:burnin),],2,median)
beta.post.sd <- apply(beta.post[-(1:burnin),],2,sd)

sigsq.post.mean <- mean(sigsq.post[-(1:burnin)])
sigsq.post.median <- median(sigsq.post[-(1:burnin)])
sigsq.post.sd <- sd(sigsq.post[-(1:burnin)])

tausq.post.mean <- apply(tausq.post[-(1:burnin),],2,mean)
tausq.post.median <- apply(tausq.post[-(1:burnin),],2,median)
tausq.post.sd <- apply(tausq.post[-(1:burnin),],2,sd)

lambda.post.mean <- mean(sqrt(lambda.post[-(1:burnin)]))
lambda.post.median <- median(sqrt(lambda.post[-(1:burnin)]))
lambda.post.sd <- sd(sqrt(lambda.post[-(1:burnin)]))

return.list <- list(beta.post.mean=beta.post.mean,beta.post.median=beta.post.median,beta.post.sd=beta.post.sd,sigsq.post.mean=sigsq.post.mean,
sigsq.post.median=sigsq.post.median,sigsq.post.sd=sigsq.post.sd,tausq.post.mean=tausq.post.mean,tausq.post.median=tausq.post.median,tausq.post.sd=tausq.post.sd,lambda.post.mean=lambda.post.mean,lambda.post.median=lambda.post.median,lambda.post.sd=lambda.post.sd,beta.MCMC.store = beta.post, sigsq.MCMC.store = sigsq.post, tausq.MCMC.store = tausq.post,lambda.MCMC.store = sqrt(lambda.post))

return(return.list)

 
}