robust_optimCDCC = function(rt, cy1, chisq1, cy2, chisq2, ini_par=c(.05, .93)){
  resta <- rbind(c(-1, -1), diag(2))
  restb <- c(-0.9999, 0.00001, 0.00001)

  opt = stats::constrOptim(theta=ini_par, 
                           f=robust_loglikelihoodCDCC, 
                           grad=NULL, 
                           ui=resta, 
                           ci=restb, 
                           # mu=1e-6, 
                           outer.iterations=400, 
                           outer.eps=1e-07, 
                           # control=list(maxit=100,reltol=1e-5), 
                           rt=rt, 
                           cy1=cy1, 
                           chisq1=chisq1,
                           cy2=cy2,
                           chisq2=chisq2)
  opt = opt$par
  
  return(opt)
}

