loglikelihoodCDCC = function(param, rt){
  alpha <- param[1]
  beta <- param[2]
  
  ndim <- ncol(rt)
  nobs <- nrow(rt)
  
  res <- robcdcc::compositeCDCC_C(alpha, beta, rt, nobs, ndim)
  
  return(res)
}

