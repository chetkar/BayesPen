	wkdir = "C:/Users/chetkar/Workspace/BioInformatics"
	setwd(wkdir);
	library(mvtnorm);
	library(SuppDists);
	library(MCMCpack);
	library(MASS);
	#FIX PARAMETERS
 
    n           = 100
    beta        = c(rep(0,10),rep(2,10),rep(0,10),rep(2,10))
    p           = length(beta)
    sig.ep      = 1
    rho         = 0.5
    Cor.X       = matrix(rho, nrow=p, ncol=p)
    diag(Cor.X) = 1             
                               
   
 
    # EXPLANATORY MATRIX
    X     = rmvnorm(n, sigma=Cor.X, method="chol")
    J      = matrix (1, nrow=n, ncol=1)
    X1      = (X - J%*%colMeans(X))/(J%*%apply(X,2,sd)) #STANDARDIZED X
   
    # RESPONSE VECTOR
    Y      = X1 %*% beta + sig.ep*rnorm(n)
	

	sourceCpp("C:\\Users\\chetkar\\Workspace\\BioInformatics\\BLasso_fixed_lambda\\BLasso_Fixed_Lambda.cpp")
	
	
	J <- matrix (1, nrow=n, ncol=1)
	X  <- (X - J%*%colMeans(X))/(J%*%apply(X,2,sd))
	
	sigsc_hyper =0.1;
	sigsh_hyper =0.1;
	
	c1 <- BLassofl(X,Y,1,sigsc_hyper,sigsh_hyper,1000,500,1)
	c2 <- BLasso_Prior_Lambda(X,Y,100,10,sigsc_hyper,sigsh_hyper,1000,500,1)
	c3 <- BlassoemR(X,Y,20,10,sigsc_hyper,sigsh_hyper,500,1)
	