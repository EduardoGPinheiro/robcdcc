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
  double loglikelihoodCDCC_C(double alpha, double beta, arma::mat epsilon, 
                             arma::mat Dt, int nobs){
    mat tilde_epsilon = zeros(nobs, 2);
    mat S = zeros(2, 2);
    mat S0 = S;
    
    double dt = 0;
    double detHt = 0;
    double lkh = 0;
    double ee12 = 0;
    double ee22 = 0;
    double intercepto = 1-alpha-beta;
    
    NumericVector Q11(nobs);
    NumericVector Q22(nobs);
    NumericVector Q12(nobs);
    
    NumericVector H12(nobs);
    NumericVector IH11(nobs);
    NumericVector IH22(nobs);
    NumericVector IH12(nobs);
    
    Q11[0] = 1.0;
    Q22[0] = 1.0;
    Q12[0] = 0;
    H12[0] = 0;
    
    IH11[0] = 1.0;
    IH22[0] = 1.0;
    IH12[0] = 0;
    
    ee12 = epsilon(0,0) * epsilon(0,0);
    ee22 = epsilon(0,1) * epsilon(0,1);
    tilde_epsilon.rows(0,0) = epsilon.rows(0,0);
    
    for(int t=1; t < nobs; t++){
      Q11[t] = intercepto + alpha * ee12 * Q11[t-1] + beta * Q11[t-1];
      Q22[t] = intercepto + alpha * ee22 * Q22[t-1] + beta * Q22[t-1];
      
      tilde_epsilon(t,0) = sqrt(Q11[t]) * epsilon(t,0);
      tilde_epsilon(t,1) = sqrt(Q22[t]) * epsilon(t,1);
      ee12 = epsilon(t,0) * epsilon(t,0);
      ee22 = epsilon(t,1) * epsilon(t,1);
    }
    
    S = unconditional_correlation(tilde_epsilon, nobs);
    S0 = (1-alpha-beta) * S;
    
    Q12[0] = S(0,1);
    H12[0] = sqrt(Dt(0,0) * Dt(0,1)) * S(0,1);
    
    for(int t=0; t < (nobs-1); t++){
      detHt = Dt(t,0) * Dt(t,1) - H12[t] * H12[t];
      
      IH11[t] = Dt(t,1) / detHt;
      IH22[t] = Dt(t,0) / detHt;
      IH12[t] = -H12[t] / detHt;
      
      dt = Dt(t,0) * epsilon(t,0) * epsilon(t,0) * IH11[t] + 
        2 * sqrt(Dt(t,0) * Dt(t,1)) * epsilon(t,0) * epsilon(t,1) * IH12[t] + 
        Dt(t,1) * epsilon(t,1) * epsilon(t,1) * IH22[t];
      lkh += log(detHt) + dt;
      
      Q12[t+1] = S0(0,1) + 
        alpha * (sqrt(Q11[t] * Q22[t]) * epsilon(t,0) * epsilon(t,1)) + 
        beta * Q12[t];
      
      H12[t+1] = sqrt(Dt(t+1,0) * Dt(t+1, 1)) * Q12[t+1] / sqrt(Q11[t+1] * Q22[t+1]); 
    }
    
    return(lkh / nobs);
  }
  
  // [[Rcpp::depends("RcppArmadillo")]]
  // [[Rcpp::export]]
  double compositeCDCC_C(double alpha, double beta, arma::mat epsilon, 
                         arma::mat Dt, int nobs, int ndim){
    double lkh = 0;
    mat epsilon_biv = zeros(nobs, 2);
    mat Dt_biv = zeros(nobs, 2);
    
    for (int i = 0; i < ndim - 1; i++){
      epsilon_biv = epsilon.cols(i,i+1);
      Dt_biv = Dt.cols(i, i+1);
      
      lkh += loglikelihoodCDCC_C(alpha, beta, epsilon_biv, Dt_biv, nobs);
    }
    
    return(lkh / (1.0 * ndim));
  }
  
  // [[Rcpp::depends("RcppArmadillo")]]
  // [[Rcpp::export]]
  arma::mat calc_Rt_C(double alpha, double beta, arma::mat epsilon, int nobs){
    mat tilde_epsilon = zeros(nobs, 2);
    mat S = zeros(2, 2);
    mat S0 = S;
    
    double detRt = 0;
    double lkh = 0;
    double ee12 = 0;
    double ee22 = 0;
    double intercepto = 1-alpha-beta;
    
    NumericVector Q11(nobs);
    NumericVector Q22(nobs);
    NumericVector Q12(nobs);
    
    NumericVector R12(nobs);
    NumericVector IR11(nobs);
    NumericVector IR22(nobs);
    NumericVector IR12(nobs);
    
    Q11[0] = 1.0;
    Q22[0] = 1.0;
    Q12[0] = 0;
    R12[0] = 0;
    
    IR11[0] = 1.0;
    IR22[0] = 1.0;
    IR12[0] = 0;
    
    ee12 = epsilon(0,0) * epsilon(0,0);
    ee22 = epsilon(0,1) * epsilon(0,1);
    tilde_epsilon.rows(0,0) = epsilon.rows(0,0);
    
    for(int t=1; t < nobs; t++){
      Q11[t] = intercepto + alpha * ee12 * Q11[t-1] + beta * Q11[t-1];
      Q22[t] = intercepto + alpha * ee22 * Q22[t-1] + beta * Q22[t-1];
      
      tilde_epsilon(t,0) = sqrt(Q11[t]) * epsilon(t,0);
      tilde_epsilon(t,1) = sqrt(Q22[t]) * epsilon(t,1);
      ee12 = epsilon(t,0) * epsilon(t,0);
      ee22 = epsilon(t,1) * epsilon(t,1);
    }
    
    S = unconditional_correlation(tilde_epsilon, nobs);
    S0 = (1-alpha-beta) * S;
    
    Q12[0] = S(0,1);
    R12[0] = S(0,1);
    
    for(int t=0; t < (nobs-1); t++){
      detRt = 1.0 - R12[t] * R12[t];
      
      IR11[t] = 1.0 / detRt;
      IR22[t] = 1.0 / detRt;
      IR12[t] = -R12[t] / detRt;
      
      Q12[t+1] = S0(0,1) + 
        alpha * (sqrt(Q11[t] * Q22[t]) * epsilon(t,0) * epsilon(t,1)) + 
        beta * Q12[t];
      
      R12[t+1] = Q12[t+1] / sqrt(Q11[t+1] * Q22[t+1]); 
    }
    
    return R12;
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
  double calc_portfolio_variance(arma::vec phi, arma::vec hat_phi, 
                                 arma::mat rt, arma::mat cont_rt,
                                 arma::mat S, arma::mat hat_S){
    int nobs = rt.n_rows;
    int ndim = rt.n_cols;
    
    double a = phi[0];
    double b = phi[1];
    
    double hat_a = hat_phi[0];
    double hat_b = hat_phi[1];
    
    double dt = 0;
    double detRt = 0;
    double intercepto = 1-a-b;
    double hat_intercepto = 1-hat_a-hat_b;
    double port_num = 0;
    double port_denom = 0; 
    double port_var = 0;
    
    arma::mat Q = S;
    arma::mat Qs = eye(ndim, ndim);
    arma::mat Rt = S;
    arma::mat IRt = S;
    arma::mat IQs = eye(ndim, ndim);
    
    arma::mat hat_Q = hat_S;
    arma::mat hat_Qs = eye(ndim, ndim);
    arma::mat hat_Rt = hat_S;
    arma::mat hat_IRt = arma::inv_sympd(hat_Rt); 
    arma::mat hat_IQs = eye(ndim, ndim);
    
    arma::mat rr(ndim, ndim);
    arma::mat cont_rr(ndim, ndim);
    
    for(int t=0; t < nobs; t++){
      // Known Ht
      Q = intercepto * S + 
        a * diagmat(Qs) * rr * diagmat(Qs) + 
        b * Q;
      
      for(int i=0; i < ndim; i++){
        Qs(i,i) = sqrt(Q(i,i));
        IQs(i,i) = 1.0 / Qs(i,i);
      }
      
      Rt = diagmat(IQs) * Q * diagmat(IQs);
      rr = rt.row(t).t() * rt.row(t); 
      
      // estimated correlation matrix 
      hat_Q = hat_intercepto * hat_S + 
        hat_a * diagmat(hat_Qs) * cont_rr * diagmat(hat_Qs) + 
        hat_b * hat_Q;
      
      for(int i=0; i < ndim; i++){
        hat_Qs(i,i) = sqrt(hat_Q(i,i));
        hat_IQs(i,i) = 1.0 / hat_Qs(i,i);
      }
      
      hat_Rt = diagmat(hat_IQs) * hat_Q * diagmat(hat_IQs);
      cont_rr= cont_rt.row(t).t() * cont_rt.row(t);
      
      // calculating excess portfolio variance
      hat_IRt = arma::inv(hat_Rt); 
      
      port_num = trace(hat_IRt * Rt * hat_IRt)/ ndim;
      port_denom = (trace(hat_IRt) / ndim) * (trace(hat_IRt) / ndim);
      
      port_var+= port_num / port_denom;  
    }
    
    return port_var / nobs;
  } 
  
  
