loglikelihoodGARCH = function(param, rt, h){
  alpha1 <- param[1]
  beta1 <- param[2]
  
  nobs <- length(rt)
  
  res <- robcdcc::loglikelihoodGARCH_C(alpha1, beta1, rt, nobs, h)
  
  return(res)
}
