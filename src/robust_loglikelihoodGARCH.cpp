#include <stdio.h>
#include <stdlib.h>
#include <Rcpp.h>

using namespace Rcpp;

/* Auxiliary functions --------------------------------------------------------*/
double rho_g(double u){
  double res = -u + 0.8260 * ((5.0) * log(1.0 + exp(u) / 2.0));
  
  return res;
}

double rc_g(double x, double chisq){
  double res = 0;
  
  if(x <= chisq){
    res = 1.0;
  }else{
    res = chisq / x;
  }
  
  return res;
}

// [[Rcpp::export]]
double robust_loglikelihoodGARCH_C(double alpha1, 
                                   double beta1, 
                                   NumericVector rt, 
                                   int nobs, 
                                   double h, 
                                   double cy, 
                                   double chisq){
  double ht = h;
  double alpha0 = h * (1-alpha1-beta1);
  double lkh = 0.0;
  
  for(int t=0; t < nobs; t++){
    lkh += rho_g(log(pow(rt[t], 2.0) / ht));
    
    ht = alpha0 + 
      alpha1 * cy * rc_g(pow(rt[t], 2.0) / ht, chisq) * pow(rt[t], 2.0) +
      beta1 * ht;
  }
  
  return(lkh / nobs); 
}

// [[Rcpp::export]]
NumericVector robust_calc_ht_C(double omega, 
                               double alpha1, 
                               double beta1, 
                               NumericVector rt, 
                               int nobs,
                               double cy, 
                               double chisq){
  NumericVector ht(nobs);

  ht[0] = omega / (1-alpha1-beta1);
  
  for(int t=0; t < (nobs-1); t++){
    ht[t+1] = omega + alpha1 * cy * rc_g(pow(rt[t], 2.0) / ht[t], chisq) * pow(rt[t], 2.0) +
      beta1 * ht[t];
  }
  
  return ht;
}
