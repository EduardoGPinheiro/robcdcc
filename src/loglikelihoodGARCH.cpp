#include <stdio.h>
#include <stdlib.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double loglikelihoodGARCH_C(double alpha1, double beta1, 
                            NumericVector rt, int nobs, double h){
  double ht = h;
  double alpha0 = h * (1-alpha1-beta1);
  double lkh = 0.0; 
  
  for(int t=0; t < nobs; t++){
    lkh += pow(rt[t], 2.0) / ht + log(ht);
    
    ht = alpha0 +  alpha1 * pow(rt[t], 2.0) + beta1 * ht;
  }
  
  return(lkh / nobs); 
}

// [[Rcpp::export]]
NumericVector calc_ht_C(double omega, double alpha1, double beta1, 
                        NumericVector rt, int nobs){
  NumericVector ht(nobs+1);
  ht[0] = omega / (1-alpha1-beta1);
  
  for(int t=0; t < nobs; t++){
    ht[t+1] = omega + alpha1 * pow(rt[t], 2.0) + beta1 * ht[t];
  }
  
  return ht;
}