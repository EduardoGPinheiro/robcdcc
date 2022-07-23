fc = function(delta, N){
  set.seed(1)
  chisq = rchisq((10^5), N)
  w = as.list(chisq) %>% 
    lapply(function(x){x * min(1.0, qchisq(delta, N)/x)}) %>% 
    unlist
  
  res = mean(chisq) / mean(w)
  
  return(res)
}

biv_robust_control = function(delta){
  rc = robust_control = case_when(
    delta == .99 ~ c(1.0185, 1.0101, .99), 
    delta == .975 ~ c(1.0465, 1.0257, .975), 
    delta == .95 ~ c(1.0953, 1.0526, .95), 
    delta == .90 ~ c(1.2030, 1.111, .90), 
    TRUE ~ c(fc(delta, 1), fc(delta, 2), delta))
  return(rc)
}

corr_reweighted = function(rt, delta=.975){
  # Get dimentions and correction factors
  N = ncol(rt)
  
  if(N==2){
    robust_control = biv_robust_control(delta)  
    cy2 = robust_control[2]
  } else {
    cy2 = fc(delta, N)  
  }
  
  chisq2 = qchisq(delta, N)
  corr_reweighted_C(St=rt, cy2=cy2, chisq2=chisq2)
}

#' Estimate GARCH(1, 1) parameters
#'
#' @param rt return matrix
#'
#' @return A list containing the estimated parameters, standardized residuals 
#' and conditional variance for each series
#' 
#' @importFrom dplyr %>%
#' @importFrom purrr map

estimateGARCH = function(rt){
  res = apply(rt, 2, optimGARCH)
  par_df = map(res,1) %>% rlist::list.rbind(.)
  std_residuals = map(res, 2) %>% rlist::list.cbind(.)
  ri = seq(1, nrow(par_df))
  par_df$ri = ri
  Dt = map(res, 3) %>% rlist::list.cbind(.)
  
  return(list(par_df=par_df, std_residuals=std_residuals, Dt=Dt))
}

#' Calculate ht for GARCH(1,1)
#'
#' @param rt returns
#' @param eta GARCH parameters
#'
#' @return volatility (ht)

calc_ht = function(rt, eta) {
  rt = rt - eta[1]
  nobs = length(rt)

  if(is.vector(rt)==TRUE){
    return(calc_ht_C(eta[2], eta[3], eta[4], rt, nobs))
  }else{
    print('rt in calc_ht is not a vector')
  }
}

#' Calculate Dt for GARCH(1,1), a diagonal matrix containing estimated 
#' volatilities for each series in each time.
#'
#' @param rt returns
#' @param eta GARCH parameters
#'
#' @return DataFrame of volatilities (Dt)

calc_Dt = function(rt, eta_df){
  eta_df = eta_df %>% select(mu, omega, alpha, beta)
  eta_lst = eta_df %>% apply(FUN=list, MARGIN=1) %>% map(1) %>% lapply(unlist)
  rt_lst = rt %>% apply(FUN=list, MARGIN=2) %>% map(1)
  
  Dt = mapply(
    calc_ht,
    rt=rt_lst,
    eta=eta_lst,
    SIMPLIFY = FALSE
  ) %>%
    rlist::list.cbind(.) 
  
  return(Dt)
}

#' Calculate devolatilized returns for GARCH(1,1).
#'
#' @param rt returns
#' @param eta GARCH parameters
#'
#' @return devolatilized returns

calc_devolatilized_returns = function(rt, eta_df){
  nobs = nrow(rt)
  Dt = calc_Dt(rt=rt, eta_df=eta_df)
  
  for(i in 1:ncol(rt)){
    rt[, i] = rt[, i] - eta_df$mu[i]
    rt[, i] =  rt[, i] / sqrt(Dt[1:nobs, i])
  }
  
  return(rt)
}

#' Estimate GARCH(1,1)-cDCC(1, 1) parameters
#'
#' @param rt return matrix
#'
#' @return A list containing the estimated parameters, standardized residuals 
#' and conditional variance for each series

estimateCDCC = function(rt, ini_par=c(.05, .93)){
  # GARCH estimation
  estimated_GARCH_result = estimateGARCH(rt)
  
  eta = estimated_GARCH_result$par_df
  epsilon = estimated_GARCH_result$std_residuals
  Dt = estimated_GARCH_result$Dt
  
  # cDCC estimation
  phi = optimCDCC(epsilon, ini_par=ini_par)
  
  return(list(eta=eta, phi=phi, epsilon=epsilon, Dt=Dt))
}

#' Estimate cDCC(1, 1) parameters
#'
#' @param rt returns matrix
#' 
#' @return The estimated parameters for cDCC(1,1)

high_dimension_estimateCDCC = function(rt){
  phi = optimCDCC(rt)
  
  return(phi)
}

calc_Qs = function(rt, phi) {
  return(calc_Qs_C(phi=phi, rt=rt))
}

calc_Rt = function(rt, phi, S) {
  phi = as.numeric(phi)
  return(calc_Rt_C(phi, rt, S))
}

#' Estimate GARCH(1,1) parameters using BIP-GARCH specification and 
#' M-estimation
#'
#' @param rt returns
#' @param cy correction factor
#' @param chisq Quantile for Chi-Squared 
#'
#' @return devolatilized returns

robust_estimateGARCH = function(rt, cy, chisq, k=30) {
  res = apply(rt, 2, robust_optimGARCH, cy, chisq, k)
  par_df = map(res, 1) %>% rlist::list.rbind(.)
  ri = seq(1, nrow(par_df))
  par_df$ri = ri
  std_residuals = map(res, 2) %>% rlist::list.cbind(.)
  Dt = map(res, 3) %>% rlist::list.cbind(.)
  
  return(list(par_df = par_df, std_residuals = std_residuals, Dt=Dt))
}

#' Calculate devolatilized returns using BIP-GARCH specification
#'
#' @param rt returns
#' @param eta_df GARCH parameters
#' @param delta Quantile for Chi-Squared 
#'
#' @return devolatilized returns

calc_robust_devolatilized_returns = function(rt, eta_df, delta=.975){
  robust_control = biv_robust_control(delta)  
  cy1 = robust_control[1]
  chisq1 = qchisq(delta, 1) # df = 1
  
  nobs = nrow(rt)
  Dt = robust_calc_Dt(rt=rt, eta_df=eta_df, delta=delta)$Dt
  
  for(i in 1:ncol(rt)){
    rt[, i] = rt[, i] - eta_df$mu[i]
    rt[, i] =  rt[, i] / sqrt(Dt[1:nobs, i])
  }
  
  return(rt)
}

#' Calculate ht using BIP-GARCH
#'
#' @param rt returns
#' @param eta GARCH parameters
#' @param delta Quantile for Chi-Squared 
#'
#' @return volatility (ht)

robust_calc_ht = function(rt, eta, delta) {
  robust_control = biv_robust_control(delta)  
  nobs = length(rt)
  rt = rt - eta[1]

  # robust parameters
  cy = robust_control[1]
  chisq = qchisq(delta, 1)
  
  if(is.vector(rt)==TRUE){
    return(robust_calc_ht_C(eta[2], eta[3], eta[4], rt, nobs, cy, chisq))
  }else{
    print('rt in calc_ht is not a vector')
  }
}

#' Calculate Dt for BIP-GARCH(1,1), a diagonal matrix containing estimated 
#' volatilities for each series in each time.
#'
#' @param rt returns
#' @param eta_df BIP-GARCH parameters dataframe
#' @param delta BIP-GARCH parameters
#'
#' @return DataFrame of volatilities (Dt) and volatility filter values for 
#' each time 

robust_calc_Dt = function(rt, eta_df, delta=.975){
  eta_df = eta_df %>% select(mu, omega, alpha, beta)
  eta_lst = eta_df %>% apply(FUN=list, MARGIN=1) %>% map(1) %>% lapply(unlist)
  rt_lst = rt %>% apply(FUN=list, MARGIN=2) %>% map(1)
  
  Dt_results = mapply(
    robust_calc_ht,
    rt=rt_lst,
    eta=eta_lst,
    MoreArgs = list(delta = delta),
    SIMPLIFY = FALSE
  ) 
  
  Dt = Dt_results %>% map(1) %>% rlist::list.cbind(.)
  w = Dt_results %>% map(2) %>% rlist::list.cbind(.)
  
  return(list(Dt=Dt, w=w))
}

#' Estimate BIP-cDCC(1, 1) parameters
#'
#' @param rt returns matrix
#' @param delta chi-squared quantile for the definition of the weighting scheme 
#' @param k size of the local window for the reweighted mean vector and 
#' covariance matrix estimators
#' @param ini_par initial values for the parameters to be optimized over.
#' 
#' @return A list containing the estimated parameters for 
#' BIP_GARCH(1,1)-BIP-cDCC(1,1), standardized residuals and conditional 
#' variance for each series

robust_estimateCDCC = function(rt, delta=.975, k=30, ini_par=c(0.05, .93)){
  robust_control = biv_robust_control(delta)  
  cy1 = robust_control[1]
  cy2 = robust_control[2]
  delta = robust_control[3]
  
  # quantile for chi square distribution
  chisq1 = qchisq(delta, 1) # df = 1
  chisq2 = qchisq(delta, 2) # df = 2
  
  # GARCH robust estimation
  robust_estimated_GARCH_result = robust_estimateGARCH(rt, cy1, chisq1, k=k)
  
  eta = robust_estimated_GARCH_result$par_df
  epsilon = robust_estimated_GARCH_result$std_residuals
  Dt = robust_estimated_GARCH_result$Dt
  
  # cDCC robust estimation
  phi = robust_optimCDCC(rt=epsilon, cy1=cy1, chisq1=chisq1, 
                         cy2=cy2, chisq2=chisq2, ini_par = ini_par)
  
  return(list(eta=eta, phi=phi, epsilon=epsilon, Dt=Dt))
}

#' Estimate BIP-cDCC(1, 1) parameters
#'
#' @param rt returns matrix
#' @param delta chi-squared quantile for the definition of the weighting scheme 
#' 
#' @return The estimated parameters for BIP_GARCH(1,1)-BIP-cDCC(1,1)

high_dimension_robust_estimateCDCC = function(rt, delta){
  robust_control = biv_robust_control(delta)  
  cy1 = robust_control[1]
  cy2 = robust_control[2]
  delta = robust_control[3]
  
  # quantile for chi square distribution
  chisq1 = qchisq(delta, 1) # df = 1
  chisq2 = qchisq(delta, 2) # df = 2
  
  # cDCC robust estimation
  phi = robust_optimCDCC(rt, cy1, chisq1, cy2, chisq2)
  
  return(phi)  
}

robust_calc_Qs = function(rt, phi, delta=0.975) {
  robust_control = biv_robust_control(delta)  
  cy1 = robust_control[1]
  chisq1 = qchisq(delta, 1) 
  
  return(robust_calc_Qs_C(phi=phi, rt=rt, cy1=cy1, chisq1=chisq1))
}

robust_calc_Rt = function(rt, phi, S, delta=0.975) {
  N = ncol(rt)
  if(N==2){
    robust_control = biv_robust_control(delta)  
    cy2 = robust_control[2]
  } else {
    cy2 = fc(delta, N)  
  }
  
  chisq2 = qchisq(delta, N)
  
  return(robust_calc_Rt_C(phi=phi, rt=rt, S=S, cy2=cy2, chisq2=chisq2))
}

#' Evaluate cDCC(1, 1) and BIP-cDCC(1, 1) estimators on the estimation of the
#' weights for the minimum variance portfolio (MVP). The estimated weights are 
#' compared with the real ones using three performances criteria: 
#' Frobenius distance, Ratio of variance of the selected MVP and Percentage 
#' increase of the variance of the selected MVP in relation to the variance of 
#' the true MVP.
#' @param q_phi estimated parameter for the cDCC(1,1)
#' @param r_phi estimated parameter for the BIP-cDCC(1,1)
#' @param rt returns matrix
#' @param burn_rt burn-in returns (for simulation examples)
#' @param q_rt standardized returns by estimated GARCH(1,1) parameters
#' @param r_rt standardized returns by estimated BIP-GARCH(1,1) parameters
#' @param S unconditional covariance matrix
#' @param Dt matrix of conditional variance for each series 
#' @param q_Dt matrix of estimated conditional variance for each series by 
#' GARCH(1,1) Quasi-Maximum Likelihood estimator
#' @param r_Dt matrix of estimated conditional variance for each series by 
#' BIP-GARCH(1,1) robust estimator
#' @param q_S estimated unconditional covariance matrix for the 
#' Quasi-Maximum Likelihood estimator
#' @param r_S estimated unconditional covariance matrix for the 
#' robust estimator 
#' @param delta chi-squared quantile for the definition of the weighting scheme
#' 
#' @return MVP performance measures.

geral_calc_portfolio_variance = function(phi,
                                         q_phi,
                                         r_phi,
                                         rt,
                                         burn_rt,
                                         q_rt,
                                         r_rt,
                                         S,
                                         Dt,
                                         q_Dt,
                                         r_Dt,
                                         q_S,
                                         r_S,
                                         delta) {
  # Get dimensions and correction factors
  N = ncol(rt)
  
  if(N==2){
    robust_control = biv_robust_control(delta)  
    cy2 = robust_control[2]
  } else {
    cy2 = fc(delta, N)  
  }
  
  chisq2 = qchisq(delta, N)
  
  # Portfolio results
  results_lst = geral_calc_portfolio_variance_C(
    phi = phi,
    q_phi = q_phi,
    r_phi = r_phi,
    rt = rt,
    burn_rt = burn_rt,
    q_rt = q_rt,
    r_rt=r_rt,
    S = S,
    Dt = Dt,
    q_Dt = q_Dt,
    r_Dt = r_Dt,
    q_S = q_S,
    r_S = r_S,
    cy2 = cy2,
    chisq2 = chisq2
  )
  
  results_lst$portfolio_metrics = as.data.frame(results_lst$portfolio_metrics)
  results_lst$portfolio_metrics$metric = c(rep(c('evout', 'fro'), each=2), 
                                          rep('gmv', 3))
  results_lst$portfolio_metrics$model = c(rep(c('Q', 'R'), 2), 'real', 'Q', 'R')
  results_lst$portfolio_metrics = results_lst$portfolio_metrics %>%
    rename('value' = 'V1')
  return(results_lst)
}
