robust_loglikelihoodCDCC = function(param, rt, cy1, chisq1, cy2, chisq2){
  alpha <- param[1]
  beta <- param[2]
  nobs <- nrow(rt)
  ndim <- ncol(rt)
  
  res <-robcdcc::robust_compositeCDCC_C(alpha, beta, rt, nobs, ndim, 
                                        cy1=cy1, chisq1=chisq1, cy2=cy2, 
                                        chisq2=chisq2)
  
  return(res)
}
