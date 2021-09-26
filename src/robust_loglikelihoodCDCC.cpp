#include <stdio.h>
#include <stdlib.h>
#include <RcppArmadillo.h>


using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends("RcppArmadillo")]]
arma::vec my_rank(arma::vec v){
  uvec indices = sort_index(v);
  uvec indices2 = sort_index(indices);
  
  vec res = conv_to< vec >::from(indices2);
  res +=1;
  
  return res;
}

// [[Rcpp::depends("RcppArmadillo")]]
arma::mat spearman_corr(arma::mat St, int nobs){
  
  mat mrho = zeros(2, 2);
  mat rank_windows = zeros(nobs, 2);
  
  double mu1 = 0;
  double mu2 = 0;
  double den = 1 / (1.0 * (nobs-1));
  
  rank_windows.col(0) = my_rank(St.col(0));
  rank_windows.col(1) = my_rank(St.col(1));
  
  // calculating mean
  for(int n = 0; n < nobs; n++){
    mu1 += rank_windows(n, 0);
    mu2 += rank_windows(n, 1);
  }
  
  mu1 = mu1 / nobs;
  mu2 = mu2 / nobs;
  
  // subtracting from the mean
  rank_windows.col(0) += -mu1;
  rank_windows.col(1) += -mu2;
  
  // cross product
  mrho =   den * rank_windows.t() * rank_windows;
  
  // correlation matrix
  mrho(0,1) = mrho(0,1) / (sqrt(mrho(1,1) * mrho(0,0))); 
  mrho(1,0) = mrho(0,1);
  mrho(0,0) = 1.0;
  mrho(1,1) = 1.0;
  
  return(mrho);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat corr_reweighted_C(arma::mat St, int nobs){
  mat windows = zeros(nobs, 2);
  mat Ct = zeros(2, 2);
  mat SCt = zeros(2, 2);
  mat Lt = zeros(nobs, 1);
  mat RC = zeros(2, 2);
  mat diagRC = zeros(2, 2);
  mat S = zeros(2, 2);
  
  double cond = 0;
  double den = 0;
  
  for(int t=0; t < nobs; t++){
    if(t <= 30 / 2){
      windows(t, 0) = 0;
      windows(t, 1) = 30;
    }else{if(t > 30/2 && t < (nobs - 30 / 2)){
      windows(t,0) = t - (30 / 2);
      windows(t,1) = t + (30 / 2);
    }else{
      windows(t,0) = nobs - 30 - 1;
      windows(t,1) = nobs-1;
    }
    }
    
    Ct =  spearman_corr(St.rows(windows(t,0), windows(t,1)), 30 + 1);
    SCt = 2 * sin((3.141593 / 6) * Ct);
    cond = (St.rows(t,t) * inv(SCt) * St.rows(t,t).t()).eval()(0,0); 
    
    if(cond <= 5.991465){
      Lt(t,0) = 1;
    }else{
      Lt(t,0) = 0;
    }
    
    den += Lt(t,0);
  }
  
  for(int t=0; t < nobs; t++){
    RC += (St.rows(t,t).t() * St.rows(t,t) * Lt(t,0));
  }
  
  RC = (1.0257 / den) * RC;
  
  for(int i=0; i < 2; i++){
    diagRC(i,i) = sqrt(1/RC(i,i)); 
  }
  
  S = diagRC * RC * diagRC;
  
  return S;
}

/* Auxiliary functions -------------------------------------------------------*/
double rho(double u){
  double res = ((6.0) * log(1.0 + u / 2.0));
  
  return res;
}

double rc(double x, double k){
  double res = 0;
  
  if(x <= k){
    res = 1.0;
  }else{
    res =  k / x;
  }
  
  return res;
}


/* Bivariate robust loglikelihood --------------------------------------------*/
double robust_loglikelihoodCDCC_C(double alpha, double beta, arma::mat rt, 
                                  int nobs, double cy1, double chisq1,
                                  double cy2, double chisq2){
  mat rts = zeros(nobs, 2);
  mat S = zeros(2, 2);
  mat S0 = S;
  
  double dt = 0;
  double detRt = 0;
  double lkh = 0;
  double r11 = 0;
  double r22 = 0;
  double intercepto = 1-alpha-beta;
  
  NumericVector Q11(nobs);
  NumericVector Q22(nobs);
  NumericVector Q12(nobs);
  
  NumericVector R11(nobs);
  NumericVector R22(nobs);
  NumericVector R12(nobs);
  
  NumericVector IR11(nobs);
  NumericVector IR22(nobs);
  NumericVector IR12(nobs);
  
  Q11[0] = 1.0;
  Q22[0] = 1.0;
  Q12[0] = 0;
  
  R11[0] = 1.0;
  R22[0] = 1.0;
  R12[0] = 0;
  
  IR11[0] = 1.0;
  IR22[0] = 1.0;
  IR12[0] = 0;
  
  r11 = rt(0,0) * rt(0,0);
  r22 = rt(0,1) * rt(0,1);
  rts.rows(0,0) = rt.rows(0,0);
  
  for(int t=1; t < nobs; t++){
    Q11[t] = intercepto + alpha * cy1 * rc(r11, chisq1) * r11 * Q11[t-1] +
      beta * Q11[t-1];
    Q22[t] = intercepto + alpha * cy1 * rc(r22, chisq1) * r22 * Q22[t-1] + 
      beta * Q22[t-1];
    
    rts(t,0) = sqrt(Q11[t]) * rt(t,0);
    rts(t,1) = sqrt(Q22[t]) * rt(t,1);
    r11 = rt(t,0) * rt(t,0);
    r22 = rt(t,1) * rt(t,1);
  }
  
  S = corr_reweighted_C(rts, nobs);
  S0 = (1-alpha-beta) * S;
  
  Q12[0] = S(0,1);
  R12[0] = S(0,1);
  
  for(int t=0; t < (nobs-1); t++){
    detRt = 1.0 - R12[t] * R12[t];
    
    IR11[t] = 1.0 / detRt;
    IR22[t] = 1.0 / detRt;
    IR12[t] = -R12[t] / detRt;
    
    dt = rt(t,0) * rt(t,0) * IR11[t] + 
      2.0 * rt(t,0) * rt(t,1) * IR12[t] + 
      rt(t,1) * rt(t,1) * IR22[t];
    lkh += log(detRt) + .8258 * rho(dt);
    
    Q12[t+1] = S0(0,1) + 
      alpha * cy2 * rc(dt, chisq2) * (sqrt(Q11[t] * Q22[t]) * rt(t,0) * rt(t,1)) + 
      beta * Q12[t];
    
    Q11[t+1] = intercepto + 
      alpha * cy2 * rc(dt, chisq2) * (rt(t,0) * rt(t,0)) * Q11[t] + 
      beta * Q11[t]; 

    Q22[t+1] = intercepto + 
      alpha * cy2 * rc(dt, chisq2) * (rt(t,1) * rt(t,1)) * Q22[t] + 
      beta * Q22[t]; 
    
    R12[t+1] = Q12[t+1] / sqrt(Q11[t+1] * Q22[t+1]); 
    R11[t+1] = 1.0;
    R22[t+1] = 1.0;
    
  }
  
  return(lkh / nobs);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double robust_compositeCDCC_C(double alpha, double beta, arma::mat St,
                              int nobs, int ndim, double cy1, double chisq1,
                              double cy2, double chisq2){
  double lkh = 0;
  mat Stb = zeros(nobs, 2);
  
  for (int i = 0; i < ndim - 1; i++){
    Stb = St.cols(i,i+1);
    
    lkh += robust_loglikelihoodCDCC_C(alpha, beta, Stb, nobs, 
                                      cy1, chisq1, cy2, chisq2);
  }
  
  return(lkh / (1.0 * ndim));
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat robust_calc_Qs(arma::vec phi, arma::mat rt, 
                         double cy1, double chisq1){
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
        a * cy1 * rc(rii, chisq1) * rii * Qs(i,i) +
        b * Qs(i,i);
    }
  }
  
  return rts;
} 


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double robust_calc_portfolio_variance(arma::vec phi, arma::vec hat_phi, 
                                      arma::mat rt, arma::mat cont_rt,
                                      arma::mat S, arma::mat hat_S,
                                      double cy2, double chisq2){
  int nobs = rt.n_rows;
  int ndim = rt.n_cols;
  
  double a = phi[0];
  double b = phi[1];
  
  double hat_a = hat_phi[0];
  double hat_b = hat_phi[1];
  
  double dt = 0;
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
    
    // Estimated Ht
    dt = (cont_rt.row(t) * hat_IRt * cont_rt.row(t).t()).eval()(0,0);
    
    hat_Q = hat_intercepto * hat_S + 
      hat_a * cy2 * rc(dt, chisq2) * diagmat(hat_Qs) * cont_rr * diagmat(hat_Qs) + 
      hat_b * hat_Q;
    
    for(int i=0; i < ndim; i++){
      hat_Qs(i,i) = sqrt(hat_Q(i,i));
      hat_IQs(i,i) = 1.0 / hat_Qs(i,i);
    }
    
    hat_Rt = diagmat(hat_IQs) * hat_Q * diagmat(hat_IQs);
    cont_rr= cont_rt.row(t).t() * cont_rt.row(t);
    
    // calculating excess portfolio variance
    hat_IRt = arma::inv(hat_Rt); 
    
    port_num = trace(hat_IRt * Rt * hat_IRt) / ndim;
    port_denom = (trace(hat_IRt) / ndim) * (trace(hat_IRt) / ndim);
    
    port_var+= port_num / port_denom;    
    }
  
  return port_var / nobs;;
} 
