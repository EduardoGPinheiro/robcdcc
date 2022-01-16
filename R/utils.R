#' Calculate corretion factor
#'
#' @param delta Quantile for Chi-Squared 
#' @param N Dimension
#'
#' @return Correction factor for the wieghted variance estimator
#' @importFrom dplyr %>%
#' @importFrom stats rchisq, qchisq

fc = function(delta, N){
  set.seed(1)
  chisq = rchisq((10^5), N)
  w = as.list(chisq) %>% 
    lapply(function(x){x * min(1.0, qchisq(delta, N)/x)}) %>% 
    unlist
  
  res = mean(chisq) / mean(w)
  
  return(res)
}

#' Calculate corretion factor
#'
#' @param rt return matrix
#'
#' @return Fit results
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

calc_std_residuals = function(rt, eta_df){
  nobs = nrow(rt)
  std_rt = matrix(NA, ncol=ncol(rt), nrow=nobs)
  Dt = matrix(NA, ncol=ncol(rt), nrow=nobs+1)
  eta = eta_df %>% select(mu, omega, alpha, beta) %>% as.matrix
  
  for(i in 1:ncol(rt)){
    rt_i = rt[, i] - eta[i, 1]
    Dt[, i] = robcdcc::calc_ht_C(omega=eta[i, 2], alpha1=eta[i, 3], 
                                 beta1=eta[i, 4], rt=rt_i, nobs=nobs)
    std_rt[, i] = rt_i / sqrt(Dt[1:nobs, i])
  }
  return(list(std_rt=std_rt, Dt=Dt))
}

calc_ht = function(rt, eta) {
  nobs = length(rt)
  
  # parameters
  eta = as.numeric(eta)
  omega = eta[1]
  alpha = eta[2]
  beta = eta[3]
  
  return(calc_ht_C(omega, alpha, beta, rt, nobs))
}

calc_Dt = function(rt, eta){
  eta = eta %>% apply(FUN=list, MARGIN=1) %>% map(1) %>% lapply(unlist)
  rt = rt %>% apply(FUN=list, MARGIN=2) %>% map(1)
  
  Dt = mapply(
    calc_ht,
    rt,
    eta,
    SIMPLIFY = FALSE
  ) %>%
    rlist::list.cbind(.) 
  
  return(Dt)
}

# cDCC
estimateCDCC = function(rt){
  # GARCH estimation
  estimated_GARCH_result = estimateGARCH(rt)
  
  eta = estimated_GARCH_result$par_df
  epsilon = estimated_GARCH_result$std_residuals
  Dt = estimated_GARCH_result$Dt
  
  # cDCC estimation
  phi = optimCDCC(epsilon)
  
  return(list(eta=eta, phi=phi, epsilon=epsilon, Dt=Dt))
}

high_dimension_estimateCDCC = function(rt){
  phi = optimCDCC(rt)
  
  return(phi)
}

calc_Qs = function(rt, phi) {
  return(robust_calc_Qs_C(phi=phi, rt=rt, cy1=cy1, chisq1=chisq1))
}

calc_Rt = function(rt, phi) {
  nobs = nrow(rt)
  phi = as.numeric(phi)
  
  # parameters
  a = phi[1]
  b = phi[2]
  
  return(calc_Rt_C(a, b, rt, nobs))
}

# Robust GARCH
robust_estimateGARCH = function(rt, cy, chisq, k=30) {
  res = apply(rt, 2, robust_optimGARCH, cy, chisq, k)
  par_df = map(res, 1) %>% rlist::list.rbind(.)
  ri = seq(1, nrow(par_df))
  par_df$ri = ri
  std_residuals = map(res, 2) %>% rlist::list.cbind(.)
  Dt = map(res, 3) %>% rlist::list.cbind(.)
  
  return(list(par_df = par_df, std_residuals = std_residuals, Dt=Dt))
}

robust_calc_std_residuals = function(rt, eta_df, delta=.975){
  robust_control = case_when(
    delta == .99 ~ c(1.0185, 1.0101, .99), 
    delta == .975 ~ c(1.0465, 1.0257, .975), 
    delta == .95 ~ c(1.0953, 1.0526, .95), 
    delta == .90 ~ c(1.2030, 1.111, .90), 
    TRUE ~ c(fc(delta, 1), fc(delta, 2), delta))
  
  cy1 = robust_control[1]
  chisq1 = qchisq(delta, 1) # df = 1
  
  nobs = nrow(rt)
  std_rt = matrix(NA, ncol=ncol(rt), nrow=nobs)
  Dt = matrix(NA, ncol=ncol(rt), nrow=nobs+1)
  eta = eta_df %>% select(mu, omega, alpha, beta) %>% as.matrix
  
  for(i in 1:ncol(rt)){
    rt_i = rt[, i] - eta[i, 1]
    Dt[, i] = robcdcc::robust_calc_ht_C(omega=eta[i, 2], alpha1=eta[i, 3], 
                                        beta1=eta[i, 4], rt=rt_i, nobs=nobs, 
                                        cy=cy1, chisq=chisq1)
    std_rt[, i] = rt_i / sqrt(Dt[1:nobs, i])
  }
  return(list(std_rt=std_rt, Dt=Dt))
}

robust_calc_ht = function(rt, eta, rc_delta) {
  nobs = length(rt)
  
  # parameters
  eta = as.numeric(eta)
  omega = eta[1]
  alpha = eta[2]
  beta = eta[3]
  
  # robust parameters
  cy = rc_delta[1]
  chisq = qchisq(rc_delta[3], 1)
  
  return(robcdcc::robust_calc_ht_C(omega, alpha, beta, rt, nobs, cy, chisq))
}


robust_calc_Dt = function(rt, eta, rc_delta){
  eta = eta %>% apply(FUN=list, MARGIN=1) %>% map(1) %>% lapply(unlist)
  rt = rt %>% apply(FUN=list, MARGIN=2) %>% map(1)
  
  Dt = mapply(
    robust_calc_ht,
    rt,
    eta,
    MoreArgs = list(rc_delta = rc_delta),
    SIMPLIFY = FALSE
  ) %>%
    lapply(function(x){x[, 1]}) %>%
    rlist::list.cbind(.) 
  
  return(Dt)
}

# Robust cDCC
robust_estimateCDCC = function(rt, delta, k=30){
  robust_control = case_when(
    delta == .99 ~ c(1.0185, 1.0101, .99), 
    delta == .975 ~ c(1.0465, 1.0257, .975), 
    delta == .95 ~ c(1.0953, 1.0526, .95), 
    delta == .90 ~ c(1.2030, 1.111, .90), 
    TRUE ~ c(fc(delta, 1), fc(delta, 2), delta))
  
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
  phi = robust_optimCDCC(epsilon, cy1, chisq1, cy2, chisq2)
  
  return(list(eta=eta, phi=phi, epsilon=epsilon, Dt=Dt))
}

high_dimension_robust_estimateCDCC = function(rt, delta){
  robust_control = case_when(
    delta == .99 ~ c(1.0185, 1.0101, .99), 
    delta == .975 ~ c(1.0465, 1.0257, .975), 
    delta == .95 ~ c(1.0953, 1.0526, .95), 
    delta == .90 ~ c(1.2030, 1.111, .90), 
    TRUE ~ c(fc(delta, 1), fc(delta, 2), delta))

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
  robust_control = case_when(
    delta == .99 ~ c(1.0185, 1.0101, .99), 
    delta == .975 ~ c(1.0465, 1.0257, .975), 
    delta == .95 ~ c(1.0953, 1.0526, .95), 
    delta == .90 ~ c(1.2030, 1.111, .90),
    TRUE ~ c(fc(delta, 1), fc(delta, 2), delta))

  cy1 = robust_control[1]
  chisq1 = qchisq(delta, 1) 
  
  return(robust_calc_Qs_C(phi=phi, rt=rt, cy1=cy1, chisq1=chisq1))
}

robust_calc_Rt = function(rt, phi, S=S, delta=0.975) {
  N = ncol(rt)
  cy2 = fc(delta, N)
  chisq2 = qchisq(delta, N)
  
  return(robust_calc_Rt_C(phi=phi, rt=rt, S=S, cy2=cy2, chisq2=chisq2))
}


