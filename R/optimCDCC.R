optimCDCC = function(rt, ini_par=c(.05, .93)){
  resta <- rbind(c(-1, -1), diag(2))
  restb <- c(-1, 0, 0)

  opt = stats::constrOptim(theta=ini_par,
                           f=loglikelihoodCDCC, 
                           grad=NULL, 
                           ui=resta, 
                           ci=restb,
                           # mu=1e-5, 
                           outer.iterations=400,
                           outer.eps=1e-07, 
                           # control=list(maxit=100,reltol=1e-5), 
                           rt=rt)
  opt = opt$par
  
  return(opt)
}
