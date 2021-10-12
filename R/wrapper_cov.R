'Covariance matrix using wrapper.'
library(robustbase)
library(ellipse) 
library(svd) # for propack.svd
library(cellWise)
library(grid) # to plot images etc.
library(gridExtra) # has grid.arrange()
library(matrixStats)

wrapper_cov = function(X, estimate_loc_scale = FALSE) {
  cutoff <- sqrt(qchisq(0.975, ncol(X)))
  if (estimate_loc_scale == TRUE) {
    estX = estLocScale(X)
    Xw = wrap(X, estX$loc, estX$scale)$Xw
  }else{
    loc = rep(0, ncol(X))
    scale = rep(1, ncol(X))
    Xw = wrap(X, loc, scale)$Xw
  }
  
  Xw.cov = cov(Xw)
  return(list('cov' = Xw.cov, 'Xw' = Xw))
}

# robust_cor = function(data){
  wrapper_results = wrapper_cov(data)
  S = wrapper_results$cov
  rts = wrapper_results$Xw
  
  Ss = NA
  tryCatch({
    Ss = nlshrink_cov2(S=S, X=rts) %>% cov2cor
  }, error = function(e){cat("ERROR :",
                             conditionMessage(e))})
  return(Ss)
}

robust_cor2 = function(data){
  wrapper_results = wrapper_cov(data)
  S = wrapper_results$cov
  rts = wrapper_results$Xw
  
  Ss = NA
  tryCatch({
    Ss = nlshrink_cov(X=rts, method='nloptr') %>% cov2cor
  }, error = function(e){cat("ERROR :",
                             conditionMessage(e))})
  return(Ss)
}
