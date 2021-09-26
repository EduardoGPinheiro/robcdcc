//' @importFrom Rcpp evalCpp

#include <stdio.h>
#include <stdlib.h>
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// set seed
void set_seed(double seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(std::floor(std::fabs(seed)));
}

// [[Rcpp::depends(RcppArmadillo)]]
arma::mat mvrnormArma(int n, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return Y * arma::chol(sigma);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat simCDCC_C(arma::vec phi, arma::mat S, int nobs, double seed){
  double a = phi[0];
  double b = phi[1];
  double intcp = (1-a-b);
  
  int ndim = S.n_cols;
  
  arma::mat rt(nobs, ndim);
  arma::mat Q = S;
  arma::mat Qs = eye(ndim, ndim);
  arma::mat IQs = eye(ndim, ndim);
  arma::mat R(ndim, ndim);  
  arma::mat SS = S; 
  arma::mat A(ndim, ndim);
  arma::mat B(ndim, ndim);

  for(int i=0; i<ndim; i++){
    Qs(i,i) = sqrt(Q(i,i));
    IQs(i,i) = 1.0 / Qs(i,i);
  }
  
  set_seed(seed);                                      
  for(int t=0; t<nobs;t++){
    A = intcp * S;
    B = a * diagmat(Qs);
    B = B * SS; 
    B = B * diagmat(Qs);
    
    Q = b * Q;
    Q = Q + B; 
    Q = Q + A;
    
    for(int i=0; i<ndim; i++){
      Qs(i,i) = sqrt(Q(i,i));
      IQs(i,i) = 1.0 / Qs(i,i);
    }
    
    R = diagmat(IQs) * Q;
    R = R * diagmat(IQs);
    
    rt.row(t) = mvrnormArma(1, R);
    SS = rt.row(t).t() * rt.row(t); 
    
//    Rcout << "Q" << std::endl << Q << std::endl;
//    Rcout << "rt" << std::endl << rt.row(t) << std::endl;
    
  }
  
  return rt;
}
