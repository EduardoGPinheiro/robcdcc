robust_gradientGARCH = function(param, rt, h, d=1e-5, cy, chisq){
  nobs = length(rt)
  n_param = length(param)
  Id = d*diag(n_param)
  param1 <- param + Id[,1]
  param2 <- param + Id[,2]
  
  lfc = robcdcc::robust_loglikelihoodGARCH_C(param[1], param[2], 
                                    rt, nobs, h, cy, chisq)
  lfc1 = robcdcc::robust_loglikelihoodGARCH_C(param1[1], param1[2], 
                                     rt, nobs, h, cy, chisq)
  lfc2 = robcdcc::robust_loglikelihoodGARCH_C(param2[1], param2[2],
                                     rt, nobs, h, cy, chisq)
  
  c(sum((lfc1 - lfc)/d), sum((lfc2 - lfc)/d) )
}
