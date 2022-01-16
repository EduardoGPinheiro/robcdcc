
# Simulation of GARCH(1,1)
simGARCH = function(omega, alpha1, beta1, nobs, seed){
  rt = vector(length = nobs)
  ht = omega/(1-alpha1-beta1)
  set.seed(seed)
  zt = rnorm(1)
  rt[1] =  sqrt(ht) * zt
  
  set.seed(seed)
  for(t in 1:(nobs-1)){
    ht = omega + alpha1 * rt[t]^2 + beta1 * ht
    zt = rnorm(1)
    rt[t+1] = sqrt(ht) * zt
  }
 
  return(rt[501:nobs])
}

# Simulation of GARCH(1,1) contaminated
simGARCH_outliers = function(omega, alpha1, beta1, nobs, nday, d, seed){
  rt = simGARCH(omega, alpha1, beta1, nobs, seed)
  
  nobs = nobs - 500
  times = seq(nobs / nday, nobs, nobs / nday)
  set.seed((seed * pi)^2)
  idx = sample(1:length(times), round(nday * .40))
  jumps = times[idx] 
  bin_jumps = times[-idx]
  
  set.seed((seed * pi)^3)
  cond = sample(c(0,1), length(bin_jumps), replace = TRUE)
  
  rt[jumps] = rt[jumps] + d * sign(rt[jumps])
  rt[bin_jumps] = rt[bin_jumps] + cond * d * sign(rt[bin_jumps])
  
  return(rt)
}
