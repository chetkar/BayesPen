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
Rcpp::List BLassopl(Rcpp::NumericMatrix X,Rcpp::NumericVector Y,double a_shape,double b_rate, double sigsc_hyper, double sigsh_hyper,int niter, int nburn,int nthin)
{

// 0. Variable Declaration
   //     RNGScope scope;

        int ncol = X.ncol() ;
        int nrow = X.nrow() ;
        double sigma2 ;
        double sigma2_shape;
        double sigma2_scale;
        int N = niter/nthin;

        Rcpp::NumericVector tau_inv(ncol) ;
        Rcpp::NumericMatrix XtX(ncol,ncol) ;
        Rcpp::NumericVector Beta(nrow);
        Rcpp::NumericVector XtY(ncol) ;
        Rcpp::NumericVector Mu(ncol);
        Rcpp::NumericMatrix TauInverse(ncol,ncol);
        Rcpp::NumericVector Lambdasquare(ncol) ;


        arma::mat       XX(nrow,ncol) ;
        arma::colvec    one(nrow,arma::fill::ones);
        arma::colvec    onecol(ncol,arma::fill::ones);
        arma::rowvec    XX_mean(ncol) ;
        arma::mat       XttX(ncol,ncol) ;
        arma::mat       dtau_inv(ncol,ncol);
        arma::colvec    XttY(ncol) ;
        arma::colvec    Beta2(ncol) ;
        arma::colvec    YY(nrow) ;
        arma::colvec    Lambda2(ncol);


        arma::colvec    mu(ncol) ;
        arma::colvec    Beta1(ncol);
        arma::colvec    tau_inverse(ncol);

        arma::mat Betahist(N,ncol);
        arma::mat Tausqhist(N,ncol);
        Rcpp::NumericVector sigma2hist(N);
        Rcpp::NumericVector Lambda2hist(N);

        arma::colvec Betahat(ncol);
        arma::colvec Tausqhat(ncol);
        double sigma2hat;
        double Lambda2hat;
        int count;



// 1. Data Scaling
        XX      =    Rcpp::as<arma::mat>(X) ;
        XX_mean =    arma::sum(XX,0)/nrow;

        YY      =    Rcpp::as<arma::colvec>(Y) ;
        YY      =    YY - arma::accu(YY)*one/nrow ;

        XX      =    XX - one*XX_mean;
        X       =    Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(XX)) ;

        XttX    =    XX.t()*XX ;
        XtX     =    Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(XttX)) ;

        XttY    =    XX.t()*YY ;
        XtY     =    Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(XttY)) ;

// 2. Starting Values

        sigma2 = R::runif(0,1)*9.9 + 0.1;

        double lambda_sq   = R::rgamma(a_shape,1/b_rate);
//        double lambda_sq = std::pow(lambda,2);

        for ( int i =0 ; i < ncol;i++)
        {
            tau_inv[i] = 1/R::rgamma(1,2/lambda_sq) ;
            Beta[i]= 0 ;
            Lambda2[i] = lambda_sq;


        }

        tau_inverse =    Rcpp::as<arma::colvec>(tau_inv) ;
        dtau_inv   =    arma::diagmat(tau_inverse) ;



// 3. For Posterior

// 4a Displaying in the R window when the job started

// 4. MCMC Sampling

        tau_inverse =    Rcpp::as<arma::colvec>(tau_inv) ;
        dtau_inv   =    arma::diagmat(tau_inverse) ;



// 3. For Posterior

// 4a Displaying in the R window when the job started

// 4. MCMC Sampling

///*
for (int i=0 ;i < niter ; i++)
    {
//*/
//int i=1;

    TauInverse =    Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(dtau_inv)) ;

    Beta       =    sampleBetainC(sigma2,TauInverse,XtX,XtY) ;
    Beta1      =    Rcpp::as<arma::colvec>(Beta) ;

    Beta2      =    arma::square(Beta1) ;



    double a    =   ncol + a_shape;
    double b    =   Rcpp::as<double>(Rcpp::wrap(.5*arma::accu(arma::pow(tau_inverse,-1)))) +b_rate;


    lambda_sq   =   R::rgamma(a_shape,1/b_rate);

    Lambda2     =   onecol*lambda_sq;
//*/
    mu         =    arma::sqrt(sigma2*Lambda2%arma::pow(Beta2,-1)) ;
    Mu         =    Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(mu)) ;



    Lambdasquare =  Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(Lambda2));
    tau_inv    =    inversegauss(ncol,Mu,Lambdasquare) ;

    tau_inverse =    Rcpp::as<arma::colvec>(tau_inv) ;
    dtau_inv   =    arma::diagmat(tau_inverse) ;


    sigma2_shape    =(ncol + nrow -1)/2 +sigsh_hyper;
    sigma2_scale    =Rcpp::as<double>(Rcpp::wrap(.5*(YY.t() -(XX*Beta1).t())*(YY -(XX*Beta1)) + .5*(Beta1.t()*(dtau_inv*Beta1)))) +sigsc_hyper;


    sigma2          = 1/(R::rgamma(sigma2_shape,1/sigma2_scale) );

    if( (1+i)%nthin == 0)
        {
            Betahist.row(i)     = Beta1.t() ;
            Tausqhist.row(i)    = arma::pow(tau_inverse.t(),-1) ;
            sigma2hist[i]       = sigma2;
            Lambda2hist[i]      = lambda_sq ;

            if( (1+i) > nburn ) {
                    Betahat         = Betahat + Beta1 ;
                    Tausqhat        = Tausqhat + arma::pow(tau_inverse,-1);
                    sigma2hat       = sigma2hat + sigma2;
                    count           = count + 1;
                    Lambda2hat      = Lambda2hat + lambda_sq;
                                }
                }
  //  /*
    }
    //*/

// 5. Full Conditional Posterior
    Betahat     = Betahat/count;
    Tausqhat    = Tausqhat/count;
    sigma2hat   = sigma2hat/count;
    Lambda2hat  = Lambda2hat/count;


// 6. Inference

// 7. Return

///*
return Rcpp::List::create(Rcpp::Named("Betahat")=Betahat,Rcpp::Named("Tausqhat")=Tausqhat,Rcpp::Named("sigma2hat")=sigma2hat,Rcpp::Named("Betahist")=Betahist,Rcpp::Named("Tausqhist")=Tausqhist,Rcpp::Named("sigma2hist")=sigma2hist);
//*/
//return 0;
}
