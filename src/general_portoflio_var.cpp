#include <stdio.h>
#include <stdlib.h>
#include <RcppArmadillo.h>


using namespace Rcpp;
using namespace arma;

double rc2(double x, double k){
  double res = 0;
  
  if(x <= k){
    res = 1.0;
  }else{
    res =  k / x;
  }
  
  return res;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec geral_calc_portfolio_variance(arma::vec phi, 
                                        arma::vec q_phi, arma::vec r_phi,  
                                        arma::mat rt, arma::mat burn_rt,
                                        arma::mat cont_rt,
                                        arma::mat S, 
                                        arma::mat q_S, arma::mat r_S,
                                        double cy2, double chisq2){
  int nobs = rt.n_rows;
  int ndim = rt.n_cols;
  int bnobs = burn_rt.n_rows - 1;
  
  double a = phi[0];
  double b = phi[1];
  
  double q_a = q_phi[0];
  double q_b = q_phi[1];
  
  double r_a = r_phi[0];
  double r_b = r_phi[1];
  
  double dt = 0;
  double intercepto = 1-a-b;
  double q_intercepto = 1-q_a-q_b;
  double r_intercepto = 1-r_a-r_b;
  
  double q_port_num = 0;
  double q_port_denom = 0; 
  double q_port_var = 0;
  
  double r_port_num = 0;
  double r_port_denom = 0; 
  double r_port_var = 0;
  
  arma::mat Q = S;
  arma::mat Qs = eye(ndim, ndim);
  arma::mat Rt = S;
  arma::mat IRt = S;
  arma::mat IQs = eye(ndim, ndim);
  
  arma::mat q_Q = q_S;
  arma::mat q_Qs = eye(ndim, ndim);
  arma::mat q_Rt = q_S;
  arma::mat q_IRt = arma::inv(q_Rt); 
  arma::mat q_IQs = eye(ndim, ndim);
  
  arma::mat r_Q = r_S;
  arma::mat r_Qs = eye(ndim, ndim);
  arma::mat r_Rt = r_S;
  arma::mat r_IRt = arma::inv(r_Rt); 
  arma::mat r_IQs = eye(ndim, ndim);
  
  arma::mat rr(ndim, ndim);
  arma::mat cont_rr(ndim, ndim);
  arma::mat burn_rr(ndim, ndim);
  
  double q_fro=0.0;
  double r_fro=0.0;
  double gmv = 0.0; 
  double q_gmv = 0.0;
  double r_gmv = 0.0;
  
  arma::vec w(ndim); 
  arma::vec q_hat_w(ndim); 
  arma::vec r_hat_w(ndim); 
  
  arma::vec ones_vec(ndim); 
  arma::vec results(7);
  
  ones_vec += 1.0;
  // burn in 
  for(int t=0; t < bnobs; t++){
    burn_rr = burn_rt.row(t).t() * burn_rt.row(t); 
    
    // Known Ht
    Q = intercepto * S + 
      a * diagmat(Qs) * burn_rr * diagmat(Qs) + 
      b * Q;
    
    for(int i=0; i < ndim; i++){
      Qs(i,i) = sqrt(Q(i,i));
      IQs(i,i) = 1.0 / Qs(i,i);
    }
    
    Rt = diagmat(IQs) * Q * diagmat(IQs);
  }
  
  for(int t=0; t < nobs; t++){
    rr = rt.row(t).t() * rt.row(t); 
    cont_rr= cont_rt.row(t).t() * cont_rt.row(t);
    dt = (cont_rt.row(t) * r_IRt * cont_rt.row(t).t()).eval()(0,0);
    
    // Known Ht
    Q = intercepto * S + 
      a * diagmat(Qs) * rr * diagmat(Qs) + 
      b * Q;
    
    for(int i=0; i < ndim; i++){
      Qs(i,i) = sqrt(Q(i,i));
      IQs(i,i) = 1.0 / Qs(i,i);
    }
    
    Rt = diagmat(IQs) * Q * diagmat(IQs);
    
    // QMVn 
    q_Q = q_intercepto * q_S + 
      q_a * diagmat(q_Qs) * cont_rr * diagmat(q_Qs) + 
      q_b * q_Q;
    
    for(int i=0; i < ndim; i++){
      q_Qs(i,i) = sqrt(q_Q(i,i));
      q_IQs(i,i) = 1.0 / q_Qs(i,i);
    }
    
    q_Rt = diagmat(q_IQs) * q_Q * diagmat(q_IQs);
    
    // Robust
    r_Q = r_intercepto * r_S + 
      r_a * cy2 * rc2(dt, chisq2) * diagmat(r_Qs) * cont_rr * diagmat(r_Qs) + 
      r_b * r_Q;
    
    for(int i=0; i < ndim; i++){
      r_Qs(i,i) = sqrt(r_Q(i,i));
      r_IQs(i,i) = 1.0 / r_Qs(i,i);
    }
    
    r_Rt = diagmat(r_IQs) * r_Q * diagmat(r_IQs);
    
    //    Rcout << "q_diag : " << q_test_diag / ndim << "in" << t<<"\n";
    //    Rcout << "r_diag : " << q_test_diag / ndim << "in" << t<<"\n";
  }
  
  // calculating excess portfolio variance
  q_IRt = arma::inv(q_Rt); 
  r_IRt = arma::inv(r_Rt);
  IRt = arma::inv(Rt); 
  
  q_port_num = trace(q_IRt * Rt * q_IRt) / ndim;
  q_port_denom = (trace(q_IRt) / ndim) * (trace(q_IRt) / ndim);
  
  r_port_num = trace(r_IRt * Rt * r_IRt) / ndim;
  r_port_denom = (trace(r_IRt) / ndim) * (trace(r_IRt) / ndim);
  
  q_port_var = q_port_num / q_port_denom;
  r_port_var = r_port_num / r_port_denom;
  
  // calculate Frobenius loss function
  q_fro = norm(q_Rt - Rt, "fro");
  r_fro = norm(r_Rt - Rt, "fro"); 
  
  // calculate hat weights
  q_hat_w = (q_IRt * ones_vec) / (ones_vec.t() * q_IRt * ones_vec).eval()(0,0);
  r_hat_w = (r_IRt * ones_vec) / (ones_vec.t() * r_IRt * ones_vec).eval()(0,0);
  w = (IRt * ones_vec) / (ones_vec.t() * IRt * ones_vec).eval()(0,0);
  
  // Calculate variance in GMV
  gmv = (w.t() * Rt * w).eval()(0,0);
  q_gmv = (q_hat_w.t() * Rt * q_hat_w).eval()(0,0);
  r_gmv = (r_hat_w.t() * Rt * r_hat_w).eval()(0,0);
  
  results(0) = q_port_var;
  results(1) = r_port_var;
  results(2) = q_fro; 
  results(3) = r_fro;
  results(4) = gmv;
  results(5) = q_gmv;
  results(6) = r_gmv; 
  
  return results;
} 

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double real_portfolio_var(arma::vec phi, arma::mat S, arma::mat rt, int n){
  int nobs = rt.n_rows;
  int ndim = rt.n_cols;
  
  double a = phi[0];
  double b = phi[1];
  
  double dt = 0;
  double intercepto = 1-a-b;
  
  arma::mat Q = S;
  arma::mat Qs = eye(ndim, ndim);
  arma::mat Rt = S;
  arma::mat IRt = S;
  arma::mat IQs = eye(ndim, ndim);
  
  arma::mat rr(ndim, ndim);
  
  double real_var = 0;
  
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
    
    if(t >= n){
      //    Rcout << "q_diag : " << q_test_diag / ndim << "in" << t<<"\n";
      //    Rcout << "r_diag : " << q_test_diag / ndim << "in" << t<<"\n";
      
      //    Rcout << "q_port_var : " << q_port_var << "in" << t<<"\n";
      //    Rcout << "r_port_var : " << r_port_var << "in" << t<<"\n";
      
      IRt = arma::inv(Rt);
      real_var += 1.0 / (trace(IRt) / ndim); 
    }
  }
  
  return real_var / (nobs-n);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec qmvn_calc_portfolio_variance(arma::vec phi, 
                                        arma::vec q_phi, 
                                        arma::mat rt, arma::mat burn_rt,
                                        arma::mat cont_rt,
                                        arma::mat S, 
                                        arma::mat q_S){
  int nobs = rt.n_rows;
  int ndim = rt.n_cols;
  int bnobs = burn_rt.n_rows - 1;
  
  double a = phi[0];
  double b = phi[1];
  
  double q_a = q_phi[0];
  double q_b = q_phi[1];
  
  double dt = 0;
  double intercepto = 1-a-b;
  double q_intercepto = 1-q_a-q_b;

  double q_port_num = 0;
  double q_port_denom = 0; 
  double q_port_var = 0;
 
  arma::mat Q = S;
  arma::mat Qs = eye(ndim, ndim);
  arma::mat Rt = S;
  arma::mat IRt = S;
  arma::mat IQs = eye(ndim, ndim);
  
  arma::mat q_Q = q_S;
  arma::mat q_Qs = eye(ndim, ndim);
  arma::mat q_Rt = q_S;
  arma::mat q_IRt = arma::inv(q_Rt); 
  arma::mat q_IQs = eye(ndim, ndim);
  
  arma::mat rr(ndim, ndim);
  arma::mat cont_rr(ndim, ndim);
  arma::mat burn_rr(ndim, ndim);
  
  double q_fro=0.0;
  double gmv = 0.0; 
  double q_gmv = 0.0;

  arma::vec w(ndim); 
  arma::vec q_hat_w(ndim); 
  
  arma::vec ones_vec(ndim); 
  arma::vec results(7);
  
  ones_vec += 1.0;
  // burn in 
  for(int t=0; t < bnobs; t++){
    burn_rr = burn_rt.row(t).t() * burn_rt.row(t); 
    
    // Known Ht
    Q = intercepto * S + 
      a * diagmat(Qs) * burn_rr * diagmat(Qs) + 
      b * Q;
    
    for(int i=0; i < ndim; i++){
      Qs(i,i) = sqrt(Q(i,i));
      IQs(i,i) = 1.0 / Qs(i,i);
    }
    
    Rt = diagmat(IQs) * Q * diagmat(IQs);
  }
  
  for(int t=0; t < nobs; t++){
    rr = rt.row(t).t() * rt.row(t); 
    cont_rr= cont_rt.row(t).t() * cont_rt.row(t);

    // Known Ht
    Q = intercepto * S + 
      a * diagmat(Qs) * rr * diagmat(Qs) + 
      b * Q;
    
    for(int i=0; i < ndim; i++){
      Qs(i,i) = sqrt(Q(i,i));
      IQs(i,i) = 1.0 / Qs(i,i);
    }
    
    Rt = diagmat(IQs) * Q * diagmat(IQs);
    
    // QMVn 
    q_Q = q_intercepto * q_S + 
      q_a * diagmat(q_Qs) * cont_rr * diagmat(q_Qs) + 
      q_b * q_Q;
    
    for(int i=0; i < ndim; i++){
      q_Qs(i,i) = sqrt(q_Q(i,i));
      q_IQs(i,i) = 1.0 / q_Qs(i,i);
    }
    
    q_Rt = diagmat(q_IQs) * q_Q * diagmat(q_IQs);
    
    //    Rcout << "q_diag : " << q_test_diag / ndim << "in" << t<<"\n";
    //    Rcout << "r_diag : " << q_test_diag / ndim << "in" << t<<"\n";
  }
  
  // calculating excess portfolio variance
  q_IRt = arma::inv(q_Rt); 
  IRt = arma::inv(Rt); 
  
  q_port_num = trace(q_IRt * Rt * q_IRt) / ndim;
  q_port_denom = (trace(q_IRt) / ndim) * (trace(q_IRt) / ndim);
  
  q_port_var = q_port_num / q_port_denom;

  // calculate Frobenius loss function
  q_fro = norm(q_Rt - Rt, "fro");

  // calculate hat weights
  q_hat_w = (q_IRt * ones_vec) / (ones_vec.t() * q_IRt * ones_vec).eval()(0,0);

  // Calculate variance in GMV
  gmv = (w.t() * Rt * w).eval()(0,0);
  q_gmv = (q_hat_w.t() * Rt * q_hat_w).eval()(0,0);

  results(0) = q_port_var;
  results(2) = q_fro; 
  results(4) = gmv;
  results(5) = q_gmv;

  return results;
} 
