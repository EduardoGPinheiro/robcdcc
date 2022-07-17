optimGARCH = function(rt, ini_par=c(.05,.93)){
  mu_R = mean(rt)
  hat = var(rt)
  rt = rt - mu_R
  nobs = length(rt)
  
  ra <- matrix(c(1,0,0,1,-1,-1),ncol=2,byrow=TRUE)
  rb <- c(0.00001, 0.00001,-0.9999)
  
  opt = stats::constrOptim(theta=ini_par, 
                           f=loglikelihoodGARCH, 
                           grad=NULL, 
                           ui=ra, 
                           ci=rb, 
                           # mu=1e-5, 
                           outer.iterations=400, 
                           outer.eps=1e-07, 
                           rt=rt, 
                           h=hat)$par
  
  par = c(mu_R, hat*(1-opt[1]-opt[2]), opt[1], opt[2], hat)
  ht = calc_ht_C(par[2], par[3], par[4], rt, nobs)
  
  std_residuals = rt / sqrt(ht[1:nobs])
  par_df = data.frame(mu = par[1], omega = par[2], 
                      alpha = par[3], beta = par[4], h = par[5])
  return(list(par_df=par_df, std_residuals=std_residuals, ht=ht))
}
