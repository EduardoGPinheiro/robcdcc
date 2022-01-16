//' @importFrom Rcpp evalCpp

#include <stdio.h>
#include <stdlib.h>
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;
  
  // [[Rcpp::depends(RcppArmadillo)]]
  // [[Rcpp::export]]
  arma::mat unconditional_correlation(arma::mat tilde_epsilon, int nobs){
    mat S = zeros(2, 2);
    
    double S11 = 0;
    double S22 = 0;
    double S12 = 0;
    
    for(int t = 0; t<nobs; t++){
      S11 += tilde_epsilon(t,0) * tilde_epsilon(t,0);
      S22 += tilde_epsilon(t,1) * tilde_epsilon(t,1);
      S12 += tilde_epsilon(t,0) * tilde_epsilon(t,1);
    }    
    
    S(0,0) = 1.0;
    S(1,1) = 1.0; 
    S(1,0) = S(0,1) = S12 / (sqrt(S11 * S22));
    
    return S;
  }
  
  
  // [[Rcpp::depends("RcppArmadillo")]]
  double loglikelihoodCDCC_C(double alpha, double beta, 
                             arma::mat rt, int nobs){
    mat rts = zeros(nobs, 2);
    mat S = zeros(2, 2);
    mat S0 = S;
    
    double dt = 0;
    double lkh = 0;
    double r11 = 0;
    double r22 = 0;
    double intercepto = 1-alpha-beta;
    
    NumericVector Q11(nobs);
    NumericVector Q22(nobs);
    
    Q11[0] = 1.0;
    Q22[0] = 1.0;
    
    arma::mat Qt = eye(2, 2);
    arma::mat Rt = eye(2, 2); 
    arma::mat IRt = eye(2, 2);
    arma::mat Qs = eye(2, 2);
    arma::mat IQs = eye(2, 2); 
    arma::mat rr = eye(2, 2); 
    
    r11 = rt(0,0) * rt(0,0);
    r22 = rt(0,1) * rt(0,1);
    rts.rows(0,0) = rt.rows(0,0);
    
    for(int t=1; t < nobs; t++){
      Q11[t] = intercepto + alpha * r11 * Q11[t-1] +
        beta * Q11[t-1];
      Q22[t] = intercepto + alpha * r22 * Q22[t-1] + 
        beta * Q22[t-1];
      
      rts(t,0) = sqrt(Q11[t]) * rt(t,0);
      rts(t,1) = sqrt(Q22[t]) * rt(t,1);
      r11 = rt(t,0) * rt(t,0);
      r22 = rt(t,1) * rt(t,1);
    }
    
    S = unconditional_correlation(rts, nobs);
    S0 = (1-alpha-beta) * S;
    
    Qt = S;
    Rt = S;
    IRt = arma::inv(S);
    
    for(int t=0; t < nobs; t++){
      rr = rt.row(t).t() * rt.row(t);
      dt = (rt.row(t) * IRt * rt.row(t).t()).eval()(0,0);
      lkh += log(arma::det(Rt)) + dt; 

      Qt = S0 + 
        alpha * diagmat(Qs) * rr * diagmat(Qs) +
        beta * Qt; 
      
      for(int i=0; i < 2; i++){
        Qs(i,i) = sqrt(Qt(i,i)); 
        IQs(i,i) = 1.0 / Qs(i,i);
      }
      
      Rt = diagmat(IQs) * Qt * diagmat(IQs);
      IRt = arma::inv(Rt);
    }
    
    return(lkh / nobs);
  }
  
  // [[Rcpp::depends("RcppArmadillo")]]
  // [[Rcpp::export]]
  double compositeCDCC_C(double alpha, double beta, arma::mat St, 
                         int nobs, int ndim){
    double lkh = 0;
    mat Stb = zeros(nobs, 2);
    
    for (int i = 0; i < ndim - 1; i++){
      Stb = St.cols(i,i+1);

      lkh += loglikelihoodCDCC_C(alpha, beta, Stb, nobs);
    }
    
    return(lkh / (1.0 * ndim));
  }
  
  
  // [[Rcpp::depends("RcppArmadillo")]]
  // [[Rcpp::export]]
  arma::mat calc_Qs(arma::vec phi, arma::mat rt){
    int nobs = rt.n_rows;
    int ndim = rt.n_cols;
    
    mat rts = zeros(nobs, ndim);

    double a = phi[0];
    double b = phi[1];
    
    double rii = 0;
    double intercepto = 1-a-b;
    
    arma::mat Q(ndim, ndim);
    arma::mat Qs = eye(ndim, ndim);
    arma::mat Rt(ndim, ndim);
    
    for(int t=0; t < nobs; t++){
      for(int i=0; i < ndim; i++){
        rii = rt(t,i) * rt(t,i);
        rts(t, i) = sqrt(Qs(i,i)) * rt(t,i);
        
        Qs(i,i) = intercepto + 
          a * rii * Qs(i,i) +
          b * Qs(i,i);
      }
    }
    
    return rts;
  } 
  
  // [[Rcpp::depends("RcppArmadillo")]]
  // [[Rcpp::export]]
  arma::mat calc_Rt_C(arma::vec phi, 
                      arma::mat rt,
                      arma::mat S){
    int nobs = rt.n_rows;
    int ndim = rt.n_cols;
    
    double a = phi[0];
    double b = phi[1];
    double intercepto = 1-a-b;

    arma::mat Q = S;
    arma::mat Qs = eye(ndim, ndim);
    arma::mat Rt = S;
    arma::mat IQs = eye(ndim, ndim);
    
    arma::mat rr(ndim, ndim);

    for(int t=0; t < nobs; t++){
      rr = rt.row(t).t() * rt.row(t); 

      Q = intercepto * S + 
        a * diagmat(Qs) * rr * diagmat(Qs) + 
        b * Q;
      
      for(int i=0; i < ndim; i++){
        Qs(i,i) = sqrt(Q(i,i));
        IQs(i,i) = 1.0 / Qs(i,i);
      }
    }
    
    Rt = diagmat(IQs) * Q * diagmat(IQs);
    return Rt;
  } 
  