#include <iostream>
#include<RcppArmadillo.h>
#include<math.h>

//Sampling Scheme from Inverse Gaussian Distribution
//[[Rcpp::export]]
Rcpp::NumericVector inversegauss(int k,Rcpp::NumericVector mu_ig_v,Rcpp::NumericVector lambda_ig_v){

//variable declaration
Rcpp::NumericVector x_ig_v(k);
Rcpp::NumericVector UniRand   = Rcpp::runif(k); //sampling a vector from U(0,1)
Rcpp::NumericVector NormRand  = Rcpp::rnorm(k); //sampling a vector from N(0,1)

for(int i =0; i<k ; i++){

    //variable declaration
    double mu_ig = mu_ig_v(i);
    double lambda_ig = lambda_ig_v(i);
    double y_ig;
    double x_ig;

        // returns 0, when we sample from inverse gauss(0.0,lambda)
        if(mu_ig!=0.0) {

            y_ig      = pow(NormRand[i],2) ;
            x_ig      = mu_ig + (pow(mu_ig,2)*y_ig)/(2*lambda_ig) - (mu_ig/(2*lambda_ig))*sqrt(4*mu_ig*lambda_ig*y_ig+ pow(mu_ig,2)*pow(y_ig,2));

                if (UniRand[i] >= mu_ig/(mu_ig+x_ig)) {
                        x_ig_v[i] = (pow(mu_ig,2))/x_ig;
                } else {
                            x_ig_v[i] = x_ig;
                        }

        } else {
            x_ig_v[i] = 0.0;
        }

}
return x_ig_v;
}


//Sampling Beta from Multivariate Normal
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
Rcpp::NumericVector sampleBetainC(double sigma2, Rcpp::NumericMatrix Lambda,Rcpp::NumericMatrix XtX,Rcpp::NumericVector XtYminusYbar){

int n=XtX.nrow(), p=XtX.ncol();
arma::mat xtx(XtX.begin(),n,p,false) ;
arma::mat lambda(Lambda.begin(),p,p,false) ;
arma::vec xtyminusybar(XtYminusYbar.begin(),XtYminusYbar.size(),false);

arma::mat cov_beta  =  sigma2*arma::inv(xtx + lambda);
arma::vec mean_beta =  (1/sigma2)*cov_beta*xtyminusybar ;

arma::mat U(p,p) ;
arma::vec s(p) ;
arma::mat V(p,p);
arma::svd(U,s,V,cov_beta);
arma::mat d(p,p);
s = arma::pow(s,.5) ;
d = arma::diagmat(s);

Rcpp::NumericVector mean = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(mean_beta));
Rcpp::NumericMatrix cov  = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(cov_beta));

arma::vec beta(p);

arma::vec normrand = Rcpp::as<arma::vec>(Rcpp::rnorm(p)) ;
//arma::vec normrand(p, arma::fill::randn) ;

//beta = mean_beta + std::sqrt(sigma2)*arma::chol(cov_beta)*normrand;
beta = mean_beta +U*d*normrand;

return Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(beta));

}


//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
Rcpp::List BElasticNetH(Rcpp::NumericMatrix X,Rcpp::NumericVector Y, double a1,double b1, double a2,double b2, double asig, double bsig, int niter){

// Variable Declaration
int ncol = X.ncol();
int nrow = X.nrow();

arma::mat XX(nrow,ncol);
arma::mat XtX(ncol,ncol);
arma::colvec XtY(ncol);
arma::colvec YY(nrow);
arma::colvec XX1(ncol,ncol);
arma::colvec tausqe(ncol);
arma::colvec tauinvsqe(ncol);
arma::colvec onecol(ncol);
Rcpp::NumericMatrix dtauinv(ncol,ncol);
Rcpp::NumericMatrix xtx(ncol,ncol);
Rcpp::NumericVector xty(ncol);
Rcpp::NumericVector betae(ncol);
arma::colvec beta(ncol);
arma::colvec mu1(ncol);
arma::colvec mu2(ncol);
Rcpp::NumericVector Mu1(ncol);
Rcpp::NumericVector Mu2(ncol);
Rcpp::NumericVector Tau(ncol);

int N=niter;

arma::mat Betahist(N,ncol);
arma::mat Tausqhist(N,ncol);
Rcpp::NumericVector sigma2hist(N);
Rcpp::NumericVector lambda1sqhist(N);
Rcpp::NumericVector lambda2hist(N);

// Data Scaling
        XX      = Rcpp::as<arma::mat>(X);
        XX      = XX - arma::cumsum(XX,0)/nrow;
        XX1     = XX*XX ;
        XX1     = arma::cumsum(XX1,0)/(nrow-1);
        XX1     = arma::pow(XX,-1/2);
        XX      = XX*XX1;


// Basic Variables

        XtX     = XX.t()*XX;
        xtx     = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(XtX));
        YY      = Rcpp::as<arma::colvec>(Y);
        YY      = YY - arma::cumsum(YY)/nrow;
        XtY     = XX.t()*YY;
        xty     = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(XtY));

// Starting Values
double lambda1sq = R::rgamma(1,a1,1/b1);
double lambda2   = R::rgamma(1,a2,1/b2);

double sigsqe    = R::runif(.1,10);

    for (int i=0;i <ncol;i++) {
        tausqe[i] = R::rgamma(1,2/lambda1sq);
     tauinvsqe[i] = 1/tausqe[i];
    }


// MCMC Looping

    for (int i =0; i <niter; i ++) {

    // #beta

        dtauinv = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(arma::diagmat(tauinvsqe + lambda2*onecol)));
        betae   = sampleBetainC(sigsqe,dtauinv,xtx,xty);

    // #tausq
        beta    = Rcpp::as<arma::colvec>(betae);
        mu1     = lambda1sq*sigsqe*arma::pow(beta,-2)
        mu1     = arma::pow(mu1,.5);
        mu2     = lambda2*onecol;
        Mu1     = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(mu1));
        Mu2     = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(mu2));
        Tau     = inversegauss(ncol,Mu1,Mu2);
    tauinvsqe   = Rcpp::as<arma::colvec>(Tau);
        tausqe  = arma::pow(tauinvsqe,-1) ;

    // #sigsq

double  shsig   = (ncol+nrow -1)/2 + asig;
double  scsig   = .5*( YY - XX*beta).t()*( YY - XX*beta) + 0.5*beta.t()*arma::diagmat(tauinvsqe +lambda2*onecol)*beta + bsig;
        sigsqe  = 1/R::rgamma(shsig,1/scsig);

    // #lambda 1
double  shlam1  = ncol + a1;
double  sclam1  = .5*arma::accu(tausqe) + b1;
    lambda1sq   = R::rgamma(1,shlam1,1/sclam1);

    // #lambda 2
double  shlam2  = .5*ncol + a2;
double  sclam2  = .5*arma::accu(arma::pow(beta,2)) + b2;
    lambda2     = R::rgamma(1,shlam2,1/sclam2);


        Betahist.row(i)     =  beta.t() ;
        Tausqhist.row(i)    =  arma::pow(tauinvsqe.t(),-1) ;
        sigma2hist[i]       =  sisqe;
        lambda1sqhist[i]    =  lambda1sq ;
        lambda2hist[i]      =  lambda2;

    }

 return Rcpp::List::create(Rcpp::Named("Betahist")=Betahist,Rcpp::Named("Tausqhist")=Tausqhist,Rcpp::Named("sigma2hist")=sigma2hist,Rcpp::Named("Lambda1sqhist")=lambda1sqhist,Rcpp::Named("Lambdahist")=lambda2);
}
