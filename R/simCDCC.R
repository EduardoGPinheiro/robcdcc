#' Simulate GARCH(1,1)-cDCC(1,1)
#'
#' @param phi Parameters for cDCC(1,1)
#' @param eta Parameters for GARCH(1,1)
#' @param S Intercept matrix
#' @param nobs Number of observations
#' @param ndim Number of dimensions
#' @param seed seed for random sampling
#'
#' @return List with results
#' @importFrom dplyr %>%

simCDCC = function(phi, eta, S, nobs, ndim, seed) {
  alpha = phi[1]
  beta = phi[2]
  nobs = nobs + 500
  ht = matrix(0, ncol=ndim, nrow=nobs+1)
  ht[1,] = rep(1, ndim)
  R = list()
  
  rt = matrix(NA, ncol = ndim, nrow = nobs)
  Q = S
  Qs = diag(sqrt(diag(Q)))             # Qt_{*}^{1/2}
  IQs = diag(1/diag(Qs))               # Qt_{*}^{-1/2}
  R[[1]] = IQs %*% Q %*% IQs
  D = diag(ndim) 

  SS = S
  set.seed(seed)
  for (t in 1:nobs) {
    # iteration t:
    Q = (1 - alpha - beta) * S + alpha * Qs %*% SS %*% Qs + beta * Q
    Qs = diag(sqrt(diag(Q)))
    IQs = diag(1/diag(Qs)) 
    R[[t+1]] = IQs %*% Q %*% IQs

    rt[t, ] = mvnfast::rmvn(1, mu = rep(0, ndim), sigma = R[[t+1]])
    SS = rt[t, ] %o% rt[t, ]
    
    # iteration (t+1):
    rt[t, ] = rt[t,] %*% sqrt(D)
    ht[t,] = diag(D)
    diag(D) = eta[,1] + diag(eta[,2]) %*% (rt[t,]^2) + D %*% eta[,3]
  }
  
  ht[nobs+1, ] = diag(D)
  burnin_std_rt = rt[1:501, ] / sqrt(ht[1:501, ])
  std_rt = rt[501:nobs,] / sqrt(ht[501:nobs,])

  rt = rt[501:nobs,]
  ht = ht[501:(nobs+1),]
  R = R[501:nobs]
  
  return(list(rt=rt, ht=ht, Rt=R, burnin_std_rt=burnin_std_rt, std_rt=std_rt))
}


#' Simulate GARCH(1,1)-cDCC(1,1) with student t innovations
#'
#' @param phi Parameters for cDCC(1,1)
#' @param eta Parameters for GARCH(1,1)
#' @param S Intercept matrix
#' @param nobs Number of observations
#' @param ndim Number of dimensions
#' @param seed seed for random sampling
#'
#' @return List with results
#' @importFrom dplyr %>%

simCDCCt = function(phi, eta, S, nobs, ndim, seed) {
  alpha = phi[1]
  beta = phi[2]
  nobs = nobs + 500
  ht = matrix(0, ncol=ndim, nrow=nobs)
  ht[1,] = c(1,1)
  R = list()
  
  rt = matrix(NA, ncol = ndim, nrow = nobs)
  Q = S
  Qs = diag(sqrt(diag(Q)))             # Qt_{*}^{1/2}
  IQs = diag(1/diag(Qs))               # Qt_{*}^{-1/2}
  R[[1]] = IQs %*% Q %*% IQs
  D = diag(ndim) 
  
  SS = S
  set.seed(seed)
  for (t in 1:nobs) {
    # iteration t:
    Q = (1 - alpha - beta) * S + alpha * Qs %*% SS %*% Qs + beta * Q
    Qs = diag(sqrt(diag(Q)))
    IQs = diag(1/diag(Qs)) 
    R[[t+1]] = IQs %*% Q %*% IQs
    
    rt[t, ] = rt(ndim, df=5)
    rt[t, ] = sqrtm(R[[t+1]]) %*% rt[t, ]
    SS = rt[t, ] %o% rt[t, ]
    
    # iteration (t+1):
    rt[t, ] = rt[t,] %*% sqrt(D)
    diag(D) = eta[,1] + diag(eta[,2]) %*% (rt[t,]^2) + D %*% eta[,3]
    ht[t,] = diag(D)
  }
  
  rt = rt[501:nobs,]
  ht = ht[501:nobs,]
  R = R[501:nobs]
  
  return(list(rt=rt, ht=ht, Rt=R))
}



