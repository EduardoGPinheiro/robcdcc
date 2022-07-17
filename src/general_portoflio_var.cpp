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
List geral_calc_portfolio_variance_C(arma::vec phi, 
                                    arma::vec q_phi, arma::vec r_phi,  
                                    arma::mat rt, arma::mat burn_rt,
                                    arma::mat q_rt,
                                    arma::mat r_rt,
                                    arma::mat S, 
                                    arma::mat Dt, arma::mat q_Dt, 
                                    arma::mat r_Dt,
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
  
  double r_dt = 0;
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
  arma::mat Ht = S;
  arma::mat IHt = S;
  arma::mat IQs = eye(ndim, ndim);
  
  arma::mat q_Q = q_S;
  arma::mat q_Qs = eye(ndim, ndim);
  arma::mat q_Rt = q_S;
  arma::mat q_Ht = q_S;
  arma::mat q_IRt = arma::inv(q_Rt); 
  arma::mat q_IHt = arma::inv(q_Ht); 
  arma::mat q_IQs = eye(ndim, ndim);
  
  arma::mat r_Q = r_S;
  arma::mat r_Qs = eye(ndim, ndim);
  arma::mat r_Rt = r_S;
  arma::mat r_Ht = r_S;
  arma::mat r_IRt = arma::inv(r_Rt);
  arma::mat r_IHt = arma::inv(r_Ht); 
  arma::mat r_IQs = eye(ndim, ndim);
  
  arma::mat rr(ndim, ndim);
  arma::mat q_rr(ndim, ndim);
  arma::mat r_rr(ndim, ndim);
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
  arma::vec wyc(nobs);
  
  ones_vec += 1.0;
  
  arma::mat sqrt_Dt(ndim, ndim);
  arma::mat sqrt_qDt(ndim, ndim);
  arma::mat sqrt_rDt(ndim, ndim);
  
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
    q_rr= q_rt.row(t).t() * q_rt.row(t);
    r_rr= r_rt.row(t).t() * r_rt.row(t);
    
    r_dt = (r_rt.row(t) * r_IRt * r_rt.row(t).t()).eval()(0,0);
    wyc[t] = rc2(r_dt, chisq2);
      
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
      q_a * diagmat(q_Qs) * q_rr * diagmat(q_Qs) + 
      q_b * q_Q;
    
    for(int i=0; i < ndim; i++){
      q_Qs(i,i) = sqrt(q_Q(i,i));
      q_IQs(i,i) = 1.0 / q_Qs(i,i);
    }
    
    q_Rt = diagmat(q_IQs) * q_Q * diagmat(q_IQs);
    
    // Robust
    r_Q = r_intercepto * r_S + 
      r_a * cy2 * rc2(r_dt, chisq2) * diagmat(r_Qs) * r_rr * diagmat(r_Qs) + 
      r_b * r_Q;
    
    for(int i=0; i < ndim; i++){
      r_Qs(i,i) = sqrt(r_Q(i,i));
      r_IQs(i,i) = 1.0 / r_Qs(i,i);
    }
    
    r_Rt = diagmat(r_IQs) * r_Q * diagmat(r_IQs);
    r_IRt = arma::inv(r_Rt);
  }
  
  for(int i=0; i < ndim; i++){
    sqrt_Dt(i,i) = sqrt(Dt(i,i));
    sqrt_qDt(i,i) = sqrt(q_Dt(i,i));
    sqrt_rDt(i,i) = sqrt(r_Dt(i,i));
  }
  
  // calculating excess portfolio variance
  q_Ht = sqrt_qDt * q_Rt * sqrt_qDt;
  r_Ht = sqrt_rDt * r_Rt * sqrt_rDt;
  Ht = sqrt_Dt * Rt * sqrt_Dt;
  
  q_IHt = arma::inv(q_Ht); 
  r_IHt = arma::inv(r_Ht);
  IHt = arma::inv(Ht); 
  
  q_port_num = trace(q_IHt * Ht * q_IHt) / ndim;
  q_port_denom = (trace(q_IHt) / ndim) * (trace(q_IHt) / ndim);
  
  r_port_num = trace(r_IHt * Ht * r_IHt) / ndim;
  r_port_denom = (trace(r_IHt) / ndim) * (trace(r_IHt) / ndim);
  
  q_port_var = q_port_num / q_port_denom - 1 / (trace(IHt) / ndim);
  r_port_var = r_port_num / r_port_denom - 1 / (trace(IHt) / ndim);
  
  // calculate Frobenius loss function
  q_fro = norm(q_Ht - Ht, "fro");
  r_fro = norm(r_Ht - Ht, "fro"); 
  
  // calculate hat weights
  q_hat_w = (q_IHt * ones_vec) / (ones_vec.t() * q_IHt * ones_vec).eval()(0,0);
  r_hat_w = (r_IHt * ones_vec) / (ones_vec.t() * r_IHt * ones_vec).eval()(0,0);
  w = (IHt * ones_vec) / (ones_vec.t() * IHt * ones_vec).eval()(0,0);
  
  // Calculate variance in GMV
  gmv = (w.t() * Ht * w).eval()(0,0);
  q_gmv = (q_hat_w.t() * Ht * q_hat_w).eval()(0,0);
  r_gmv = (r_hat_w.t() * Ht * r_hat_w).eval()(0,0);
  
  results(0) = q_port_var;
  results(1) = r_port_var;
  results(2) = q_fro; 
  results(3) = r_fro;
  results(4) = gmv;
  results(5) = q_gmv;
  results(6) = r_gmv; 
  
  return Rcpp::List::create(Rcpp::Named("portfolio_metrics") = results,
                            Rcpp::Named("w") = wyc,
                            Rcpp::Named("Ht") = Ht,
                            Rcpp::Named("q_Ht") = q_Ht,
                            Rcpp::Named("r_Ht") = r_Ht);
} 

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double calc_real_portfolio_variance_C(arma::vec phi, 
                                      arma::mat burn_rt,
                                      arma::mat rt, 
                                      arma::mat S, 
                                      arma::mat Dt){
  int nobs = rt.n_rows;
  int ndim = rt.n_cols;
  int bnobs = burn_rt.n_rows - 1;
  
  double a = phi[0];
  double b = phi[1];
  
  double intercepto = 1-a-b;
  double port_var=0;
  
  arma::mat Q = S;
  arma::mat Qs = eye(ndim, ndim);
  arma::mat Rt = S;
  arma::mat Ht = S;
  arma::mat IHt = S;
  arma::mat IQs = eye(ndim, ndim);

  arma::mat burn_rr(ndim, ndim);
  arma::mat rr(ndim, ndim);

  arma::mat sqrt_Dt(ndim, ndim);

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

    // Known Ht
    Q = intercepto * S + 
      a * diagmat(Qs) * rr * diagmat(Qs) + 
      b * Q;
    
    for(int i=0; i < ndim; i++){
      Qs(i,i) = sqrt(Q(i,i));
      IQs(i,i) = 1.0 / Qs(i,i);
    }
    
    Rt = diagmat(IQs) * Q * diagmat(IQs);
  }
  
  for(int i=0; i < ndim; i++){
    sqrt_Dt(i,i) = sqrt(Dt(i,i));
  }
  
  // calculating excess portfolio variance
  Ht = sqrt_Dt * Rt * sqrt_Dt;
  IHt = arma::inv(Ht); 
  port_var = 1 / (trace(IHt) / ndim);

  return port_var;
} 
