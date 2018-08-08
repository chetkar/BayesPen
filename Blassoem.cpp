#include <iostream>
#include <RcppArmadillo.h>


//Sampling Scheme from Inverse Gaussian Distribution
//[[Rcpp::export]]
double rinvgauss(double mu, double lambda){

    //variable declaration
    double mu_ig = mu;
    double lambda_ig = lambda;
    double y_ig;
    double x_ig;



            y_ig      = pow(R::rnorm(0,1),2) ;
            x_ig      = mu_ig + (pow(mu_ig,2)*y_ig)/(2*lambda_ig) - (mu_ig/(2*lambda_ig))*sqrt(4*mu_ig*lambda_ig*y_ig+ pow(mu_ig,2)*pow(y_ig,2));

                if (R::runif(0,1) >= mu_ig/(mu_ig+x_ig)) {
                        x_ig = (pow(mu_ig,2))/x_ig;
                } else {
                            x_ig = x_ig;
                        }

return x_ig;
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
Rcpp::List BLassoem(Rcpp::NumericMatrix X,Rcpp::NumericVector Y,double lambdasq,double sigem, double sigsc_hyper,double sigsh_hyper,Rcpp::NumericVector tauem,Rcpp::NumericVector tausq,Rcpp::NumericVector betaem,int niter,int nrep,int nrep0)
{

// Variable Declaration
int ncol = X.ncol() ;
int nrow = X.nrow() ;

arma::colvec J(nrow, arma::fill::ones) ;
arma::colvec one(ncol,arma::fill::ones);
arma::mat xx(nrow,ncol) ;
arma::rowvec xx_mean(ncol) ;
arma::rowvec xx_sd(ncol);
arma::colvec YY(ncol);
arma::mat XttX(ncol,ncol);
arma::colvec XttY(ncol);
Rcpp::NumericMatrix XtX;
Rcpp::NumericVector XtY;

// Data Scaling
        xx          = Rcpp::as<arma::mat>(Rcpp::wrap(X)) ;

        YY      =    Rcpp::as<arma::colvec>(Y) ;

        XttX    =    xx.t()*xx ;
        XtX     =    Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(XttX)) ;

        XttY    =    xx.t()*YY ;
        XtY     =    Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(XttY)) ;


    double sigsq = R::runif(0.1,10);
    double sum_em =0.0;

    arma::colvec betaest(ncol);

    Rcpp::NumericVector Beta(ncol);

    int ntot = nrep + 1 ;
    int ntotby2 =nrep0;

    arma::mat Tausqemp(ntot,ncol, arma::fill::zeros) ;
    arma::mat  Betaemp(ntot,ncol,arma::fill::zeros);
    Rcpp::NumericVector sigsqemp(ntot);

    int ntotal = niter + 1 ;

    arma::mat Tausqem(ntotal,ncol,arma::fill::zeros) ;
    arma::mat  Betaem(ntotal,ncol,arma::fill::zeros);
    Rcpp::NumericVector sigsqem(ntotal);
    Rcpp::NumericVector lambdaem(ntotal);

    Rcpp::NumericMatrix dtauinv(ncol,ncol);
 //   Rcpp::NumericVector Beta(ncol);
    double z, sigshape, sigscale;

     Rcpp::NumericMatrix d(ncol,ncol);

    Tausqem.row(0)          = Rcpp::as<arma::rowvec>(tauem) ;
    Betaem.row(0)           = Rcpp::as<arma::rowvec>(betaem) ;
    sigsqem[0]              = sigem;
    lambdaem[0]             = lambdasq;

    Tausqemp.row(0)          = Rcpp::as<arma::rowvec>(tauem) ;
    Betaemp.row(0)           = Rcpp::as<arma::rowvec>(betaem) ;
    sigsqemp[0]              = sigem;

    arma::colvec Betahat(ncol,arma::fill::zeros);
    arma::colvec Tausqhat(ncol,arma::fill::zeros);
    double sigma2hat;
    double Lambda2hat;
    int count =0, ind=0;
    arma::mat dtaui(ncol,ncol,arma::fill::zeros);
 Rprintf("One is %f\n",2.0);
 for (int i = 0 ; i < niter ; i ++ )
 {

        for (int j =0; j < nrep ; j++ )
        {

                for (int l=0 ; l <ncol; l++) {
                    for(int m=0; m< ncol; m++){
                        if( l ==m ){
                        dtauinv(l,l) = 1/tauem[l] ;
                        } else
                        {
                            dtauinv(l,m) = 0 ;
                        }
                    }
                }


        Beta            = sampleBetainC(sigem,dtauinv,XtX,XtY) ;
        ind         = j +1;
        Betaemp.row(ind) = Rcpp::as<arma::rowvec>(Beta) ;
        Rprintf("iteration is %f\n",sum(Betaemp.row(j)));

                for ( int k=0 ; k < ncol ;k++) {

                    double temp = (lambdasq*sigem) /(Betaemp(ind,k)*Betaemp(ind,k)) ;
                    z =  rinvgauss(std::sqrt(temp),lambdasq);
                    tauem[k] = 1/z;
                    //Rprintf("temp : %f\n",std::sqrt(temp)) ;
                    //Rprintf("1/tau[i] : %f\n",1/z) ;

                }

        Tausqemp.row(ind) = Rcpp::as<arma::rowvec>(tauem) ;

        }


int row = i +1;

        sum_em=0.0 ;
        for (int l =ntotby2;l<ntot;l++) {
            for (int m=0;m<ncol;m++){
                sum_em = sum_em + Tausqemp(l,m);
            }

        }
        Rprintf("Tau0 sum : %f\n",sum_em) ;
        lambdasq    = ( 2*ncol*(ntot-ntotby2 ) )/sum_em;
        lambdaem[row] = lambdasq;


       for (int l=0 ; l <ncol; l++) {
                for(int m=0; m< ncol; m++){
                    if( l ==m ){
                    dtauinv(l,l) = 1/tausq[l] ;
                    } else
                    {
                    dtauinv(l,m) = 0.0 ;
                    }
                }
            }


        Beta            = sampleBetainC(sigsq,dtauinv,XtX,XtY) ;
        Betaem.row(row) = Rcpp::as<arma::rowvec>(Beta) ;
        Rprintf("Beta em sum : %f\n",arma::accu(Betaem.row(row))) ;

        for ( int k=0 ; k < ncol ;k++) {

        double t = lambdasq*sigsq/(Betaem(row,k)*Betaem(row,k)) ;
        //Rprintf("t : %f\n",std::sqrt(t)) ;
        //Rprintf("1/tau[i] : %f\n",1/z) ;

        z =    rinvgauss(std::sqrt(t),lambdasq);

        tausq[k] = 1/z;
        tauem[k] = 1/z;
       }

       if (i==0){
        Lambda2hat = lambdasq ;
        Tausqhat = Rcpp::as<arma::colvec>(tausq) ;
        sigma2hat = sigsq;
        Betahat  = Rcpp::as<arma::colvec>(Beta) ;
        d = dtauinv;
        }

        Tausqem.row(row) = Rcpp::as<arma::rowvec>(tausq) ;
        Rprintf("Tau sq sum : %f\n",arma::accu(Tausqem.row(row))) ;

            for (int l=0 ; l <ncol; l++) {
                for(int m=0; m< ncol; m++){
                    if( l ==m ){
                    dtaui(l,l) = 1/tausq[l] ;
                    } else
                    {
                    dtaui(l,m) = 0 ;
                    }
                }
            }



      sigshape        = (nrow + ncol -1)/2 + sigsh_hyper;
      betaest         =  Rcpp::as<arma::colvec>(Beta);

      sigscale = .5*Rcpp::as<double>(Rcpp::wrap((YY.t() -(xx*betaest).t())*(YY -(xx*betaest)) )) + sigsc_hyper;
      Rprintf("Sigscale : %f\n",sigscale) ;

      sigscale = sigscale + .5*Rcpp::as<double>(Rcpp::wrap( (betaest.t()*(dtaui)*betaest))) ;
      Rprintf("Sum dtau %f\n",arma::accu(dtaui));
      Rprintf("Sigscale : %f\n",sigscale) ;

      sigsq           = 1/R::rgamma(sigshape,1/sigscale) ;
      Rprintf("Sigma : %f\n",sigsq) ;
      sigsqem[row]    = sigsq;
      sigem           = sigsq;


}



return Rcpp::List::create(Rcpp::Named("Betahist")=Betaem,Rcpp::Named("Tausqhist")=Tausqemp,Rcpp::Named("sigma2hist")=sigsqem,Rcpp::Named("Lambda2hist")=lambdaem,Rcpp::Named("Lambdahat")=Lambda2hat,Rcpp::Named("d")=d);
//return sigma2;
}
