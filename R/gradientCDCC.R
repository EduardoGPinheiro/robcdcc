gradientCDCC = function(param, rt, d=1e-5){
  ndim = ncol(rt)
  nobs = nrow(rt)
  n_param = length(param)
  Id = d*diag(n_param)
  param1 <- param + Id[,1]
  param2 <- param + Id[,2]
  
  lfc = robcdcc::compositeCDCC_C(param[1], param[2], rt, nobs, ndim)
  lfc1 = robcdcc::compositeCDCC_C(param1[1], param1[2], rt, nobs, ndim)
  lfc2 = robcdcc::compositeCDCC_C(param2[1], param2[2], rt, nobs, ndim)
  
  c(sum((lfc1 - lfc)/d), sum((lfc2 - lfc)/d))
}
