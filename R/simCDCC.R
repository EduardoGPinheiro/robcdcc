#' Simulate cDCC(1,1)
#'
#' @param phi Parameters for cDCC(1,1)
#' @param eta Parameters for GARCH(1,1)
#' @param S Intercept matrix
#' @param nobs Number of observations
#' @param ndim Number of dimentions
#' @param seed seed for random sampling
#'
#' @return List with results
#' @importFrom dplyr %>%

simCDCC = function(phi, eta, S, nobs, ndim, seed) {
  'Simulation of cDCC(1,1)'
  alpha = phi[1]
  beta = phi[2]
  nobs = nobs + 500
  ht = matrix(0, ncol=ndim, nrow=nobs+1)
  ht[1,] = rep(1, ndim)
  R = list()
  
  rt = matrix(NA, ncol = ndim, nrow = nobs)
  Q = S
  Qs = diag(sqrt(diag(Q)))             # Qt_{*}^{1/2}
  IQs = diag(1/diag(Qs))               # Qt_{*}^{-1/2}
  R[[1]] = IQs %*% Q %*% IQs
  D = diag(ndim) 

  SS = S
  set.seed(seed)
  for (t in 1:nobs) {
    # iteraction t:
    Q = (1 - alpha - beta) * S + alpha * Qs %*% SS %*% Qs + beta * Q
    Qs = diag(sqrt(diag(Q)))
    IQs = diag(1/diag(Qs)) 
    R[[t+1]] = IQs %*% Q %*% IQs

    rt[t, ] = mvnfast::rmvn(1, mu = rep(0, ndim), sigma = R[[t+1]])
    SS = rt[t, ] %o% rt[t, ]
    
    # iteraction (t+1):
    rt[t, ] = rt[t,] %*% sqrt(D)
    ht[t,] = diag(D)
    diag(D) = eta[,1] + diag(eta[,2]) %*% (rt[t,]^2) + D %*% eta[,3]
  }
  
  ht[nobs+1, ] = diag(D)
  burnin_std_rt = rt[1:501, ] / sqrt(ht[1:501, ])
  std_rt = rt[501:nobs,] / sqrt(ht[501:nobs,])

  rt = rt[501:nobs,]
  ht = ht[501:(nobs+1),]
  R = R[501:nobs]
  
  return(list(rt=rt, ht=ht, Rt=R, burnin_std_rt=burnin_std_rt, std_rt=std_rt))
}

# Simulation cDCC wih t errors
simCDCCt = function(phi, eta, S, nobs, ndim, seed) {
  alpha = phi[1]
  beta = phi[2]
  nobs = nobs + 500
  ht = matrix(0, ncol=ndim, nrow=nobs)
  ht[1,] = c(1,1)
  R = list()
  
  rt = matrix(NA, ncol = ndim, nrow = nobs)
  Q = S
  Qs = diag(sqrt(diag(Q)))             # Qt_{*}^{1/2}
  IQs = diag(1/diag(Qs))               # Qt_{*}^{-1/2}
  R[[1]] = IQs %*% Q %*% IQs
  D = diag(ndim) 
  
  SS = S
  set.seed(seed)
  for (t in 1:nobs) {
    # iteraction t:
    Q = (1 - alpha - beta) * S + alpha * Qs %*% SS %*% Qs + beta * Q
    Qs = diag(sqrt(diag(Q)))
    IQs = diag(1/diag(Qs)) 
    R[[t+1]] = IQs %*% Q %*% IQs
    
    rt[t, ] = rt(ndim, df=5)
    rt[t, ] = sqrtm(R[[t+1]]) %*% rt[t, ]
    SS = rt[t, ] %o% rt[t, ]
    
    # iteraction (t+1):
    rt[t, ] = rt[t,] %*% sqrt(D)
    diag(D) = eta[,1] + diag(eta[,2]) %*% (rt[t,]^2) + D %*% eta[,3]
    ht[t,] = diag(D)
  }
  
  rt = rt[501:nobs,]
  ht = ht[501:nobs,]
  R = R[501:nobs]
  
  return(list(rt=rt, ht=ht, Rt=R))
}

# simulate cDCC with C++ function
high_dimension_simCDCC = function(phi, nobs, ndim, seed){
  nobs = nobs + 100
  # set.seed(sqrt(seed) + 1)
  set.seed(seed * 10)
  rho = rtruncnorm(
    ndim,
    a = .5 - 4 * .1,
    b = .5 + 4 * .1,
    mean = .5,
    sd = .1
  )
  S = rho %o% rho
  diag(S) = rep(1, ndim)
  
  rt = robcdcc::simCDCC_C(phi, S, nobs, seed)
  rt = rt[101:nobs, ]
  results = list(rt=rt, seed=seed)
  return(results)
}

burnin_high_dimension_simCDCC = function(phi, ndim, seed){
  set.seed(seed * 10)
  rho = rtruncnorm(
    ndim,
    a = .5 - 4 * .1,
    b = .5 + 4 * .1,
    mean = .5,
    sd = .1
  )
  S = rho %o% rho
  diag(S) = rep(1, ndim)
  
  rt = robcdcc::simCDCC_C(phi, S, 101, seed)
  results = list(rt=rt, seed=seed)
  return(results)
}

# contaminating series  
# type 1
contaminate1 = function(rt, ht, seed, epsilon, d, index=FALSE){
  nobs = nrow(rt)
  ndim = ncol(rt)
  nday = epsilon * nobs
  
  if(nday == 0){
    rt[,1] = rt[,1] + 0.05
    rt[,2] = rt[,2] - 0.05
    
    cond1 = jumps = cojumps = 0
  }else{
    rt[,1] = rt[,1] + 0.05
    rt[,2] = rt[,2] - 0.05
    
    aux1 = seq(nobs / nday, nobs, nobs / nday)
    aux2 = c(0, aux1[-length(aux1)])
    times = (aux1 + aux2) / 2 + 1
    
    set.seed(seed)
    idx = sample.int(length(times), round(nday * .40))
    
    cojumps = times[idx] %>% sort(.)
    jumps = times[-idx]
    
    cond1 = sample(c(0,1), length(jumps), replace = TRUE)
    cond2 = cond1-1
    
    rt[jumps, 1] = rt[jumps, 1] + cond1 * sqrt(ht[jumps, 1]) * d 
    rt[jumps, 2] = rt[jumps, 2] + cond2 * sqrt(ht[jumps, 2]) * d 
    
    rt[cojumps, 1] = rt[cojumps, 1] + sqrt(ht[cojumps, 1]) * d 
    rt[cojumps, 2] = rt[cojumps, 2] - sqrt(ht[cojumps, 2]) * d 
  }
  
  if(index == TRUE){
    return(c(cond1 * jumps, cojumps))
    break
  }
  
  return(rt)
}


# type 2
contaminate2 = function(rt, ht, seed, nday, d){
  nobs = nrow(rt)
  ndim = ncol(rt)
  sameday=TRUE
  perc = nday / nobs
  
  min_tempo = round(1/3 * nobs)
  max_tempo = round(2/3 * nobs)
  
  cond = 0
  if(sameday==TRUE){
    set.seed(seed)
    day = sample(min_tempo:max_tempo, nday, replace = FALSE)
    # selecting series
    while(cond!=1){
      idx = sample(1:(ndim-1), round(perc / 2 * ndim), replace = FALSE)
      idx = sort(idx)
      
      dif = idx - c( idx[-1], 0 )
      if(sum(as.numeric( abs( dif ) <= 2 )) == 0 ) cond = 1
    }
    
    # selecting series subsequent
    idx = sort(c( idx, idx +1 ))
    
    # generating outliers in data
    rt[day, idx] = rt[day, idx] + sign(rt[day, idx]) * d
    
    resultado = list(rt = rt, idx = idx, day = day)
  }else{
    set.seed(seed)
    day = sample(min_tempo:max_tempo, round(perc / 2 * ndim))
    day = rep(day, each = 2)
    while(cond!=1){
      idx = sample(1:(ndim-1), round(perc / 2 * ndim), replace = FALSE)
      idx = sort(idx)
      
      dif = idx - c( idx[-1], 0 )
      if(sum(as.numeric( abs( dif ) <= 2 )) == 0 ) cond = 1
    }
    
    # selecting series subsequent
    idx = sort( c( idx, idx +1 ) )
    coord = cbind(day, idx)
    
    #generating outliers in data
    rt[coord] = rt[coord] + sign(rt[coord]) * d 
  }
  
  return(rt)
}

# contamination for high dimension
contaminate_high1 = function(rt, epsilon, d, consecutive = FALSE, seed) {
  nseries = epsilon * ncol(rt)
  
  if (nseries == 0){
    return(rt)
  }
  
  # selecting series for contamination
  set.seed(seed)
  idx = sample(1:ncol(rt), nseries, replace = FALSE)
  
  # contaminating
  if (consecutive == FALSE) {
    rt[625, idx] = rt[625, idx] + d
    rt[1249, idx] = rt[1249, idx] + d
    
    return(rt)
  } else{
    rt[625, idx] = rt[625, idx] + d
    rt[626, idx] = rt[626, idx] + d

    rt[1249, idx] = rt[1249, idx] + d
    
    return(rt)
  }
}

contaminate_high2 = function(rt, epsilon, d, seed) {
  nseries = epsilon * ncol(rt)
  cont_t = seq(0, nrow(rt), 100)
  cont_t = cont_t[-1]
  
  if (nseries == 0){
    return(rt)
  }
  
  # selecting series for contamination
  set.seed(seed)
  idx = sample(1:ncol(rt), nseries, replace = FALSE)
  
  # contaminating
    rt[cont_t, idx] = rt[cont_t, idx] + d
    return(rt)
}