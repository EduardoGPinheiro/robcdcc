'Shrink eingvalues of a generic covariance matrix.'
library(nlshrink)

nlshrink_cov2 = function(S, X, k=0, method = 'nlminb', control = list()){
  n = nrow(X)
  p = ncol(X)
  
  if (k == 0) {
    X <- X - tcrossprod(rep(1,n), colMeans(X))
    k = 1
  }
  
  if (n > k) {
    effn <- n - k
  } else stop("k must be strictly less than nrow(X)")
  
  S_eigen = eigen(S)
  U = S_eigen$vectors
  lambda = S_eigen$values
  lambdasort = sort(lambda)
  lambdaorder = order(lambda)

  tau_est = tau_estimate2(S = S, X = X, k = k, method = method, control = control)
  nlshrink_tau = nlshrink_est(QuEST(tau_est, effn))
  return (U %*% (diag(nlshrink_tau[lambdaorder]) %*% t(U)) )
}

tau_estimate2 = function(S, X, k = 0, method = "nlminb", control = list()) {
  cat("Estimating population eigenvalues...")
  rho <- new.env()
  
  # initial set up
  rho$n = nrow(X)
  rho$p = ncol(X)
  rho$k = k
  rho$S = S
  
  if (rho$k == 0) {
    rho$X <- X - tcrossprod(rep(1,rho$n), colMeans(X))
    rho$k <- 1
  }  else rho$X <- X
  
  if (rho$n > rho$k) {
    rho$effn <- rho$n - rho$k
  } else stop("k must be strictly less than nrow(X)")
  
  rho$lambda <- sort(eigen(rho$S, only.values = TRUE)$values)
  rho$lambda[rho$lambda < 0] = 0
  if (rho$p > rho$effn)
    rho$lambda[1:(rho$p - rho$effn)] <- 0
  
  rho$lambda_ls <- linshrink(X = rho$X, k = rho$k)
  
  # Rest of function depends on method chosen
  if (method == "nloptr") {
    lambda_eval_f <- function(tau)
    {
      if (is.unsorted(tau)) {
        tausort <- sort(tau)
        tauorder <- order(tau)
      }
      else {
        tausort <- tau
        tauorder <- 1:rho$p
      }
      Q <- QuEST(tausort,rho$effn)
      lambda_est <- Q$lambda
      lambda_est_J <- get_lambda_J(Q)
      lambda_est_order <- lambda_est[tauorder]
      lambda_est_J_order <- lambda_est_J[tauorder,]
      return (list("objective" = mean((lambda_est_order - rho$lambda)^2) ,
                   "gradient" = 2/rho$p * crossprod((lambda_est_J_order), (lambda_est_order - rho$lambda) ) ) )
    }
    lambda_eq <- function(tau) sum(tau - rho$lambda)
    
    if (length(control) == 0) {
      control <- list(algorithm = "NLOPT_LD_AUGLAG",
                      xtol_rel=1e-8,
                      maxeval=2000,
                      local_opts=list(algorithm="NLOPT_LD_LBFGS", xtol_rel=1e-5))
    }
    nlopt_out <- nloptr(x0 = rho$lambda_ls,
                        eval_f = lambda_eval_f,
                        lb = rep(min(rho$lambda),rho$p),
                        ub = rep(max(rho$lambda), rho$p),
                        eval_g_eq = lambda_eq,
                        eval_jac_g_eq = function(x) t(rep(1,rho$p)),
                        opts=control)
    return (nlopt_out$solution)
  }
  else if (method == "nlminb") {
    lambda_obj <- function(tau)
    {
      if (is.unsorted(tau)) {
        tausort <- sort(tau)
        tauorder <- order(tau)
      }
      else {
        tausort <- tau
        tauorder <- 1:rho$p
      }
      Q <- QuEST(tausort,rho$effn)
      lambda_est <- Q$lambda
      lambda_est_order <- lambda_est[tauorder]
      return (mean((lambda_est_order - rho$lambda)^2) )
    }
    
    lambda_gr <- function(tau)
    {
      if (is.unsorted(tau)) {
        tausort <- sort(tau)
        tauorder <- order(tau)
      }
      else {
        tausort <- tau
        tauorder <- 1:rho$p
      }
      Q <- QuEST(tausort,rho$effn)
      lambda_est <- Q$lambda
      lambda_est_order <- lambda_est[tauorder]
      lambda_est_J <- get_lambda_J(Q)
      lambda_est_J_order <- lambda_est_J[tauorder,]
      return (as.numeric( 2/Q$p * crossprod((lambda_est_J_order), (lambda_est_order - rho$lambda) ) ) )
    }
    
    optim_out <- nlminb(start = rho$lambda_ls,
                        objective = lambda_obj,
                        gradient = lambda_gr,
                        lower = min(rho$lambda),
                        upper = max(rho$lambda),
                        control = control)
    return (optim_out$par)
  }
  else stop("method should be nloptr or nlminb")
}

QuEST <- function(tau, n)
{
  # first check inputs
  if (as.integer(n) <= 0) stop("n must be positive integer")
  
  # these parameters can be changed if needed
  left_tol <- 1e-8
  tol <- 1e-12
  
  if (is.unsorted(tau)) tau_ <- sort(tau) else tau_ <- tau
  tau_[tau_ < 0 & abs(tau_) < left_tol] = 0
  tau_[tau_ > 0 & abs(tau_) < tol] = 0
  if (any (tau_ < 0)) stop("all eigenvalues must be non-negative")
  if (all(tau_ == 0)) stop("atleast one eigenvalue must be non-zero")
  
  # set up environment
  Q <- new.env()
  Q$n <- as.integer(n)
  if (Q$n != n) print ("Note: n coerced to integer")
  Q$tau <- tau_
  Q$p <- length(Q$tau)
  Q$c <- Q$p / Q$n
  
  tau_nonzero <- Q$tau[Q$tau > 0]
  Q$t <- unique(tau_nonzero)
  Q$K <- length(Q$t)
  Q$pw <- vapply(Q$t, function(x) length(which(tau_nonzero == x)), integer(1))
  Q$pzw <- Q$p - sum(Q$pw)
  Q$pwt2 <- Q$pw * Q$t^2
  
  # main functions
  sup_fn <- function(Q)
  {
    #some short functions to be used in following code
    x_fn <- function(k)
      (Q$t[k]*Q$t[k+1])^(2/3) * ( (Q$pw[k]*Q$t[k+1])^(1/3) + (Q$pw[k+1]*Q$t[k])^(1/3) ) / ( Q$pwt2[k]^(1/3) + Q$pwt2[k+1]^(1/3) )
    
    p_theta_fn <- function(u, k) Q$pwt2[k]/(Q$t[k] - u)^2 + Q$pwt2[k+1]/(Q$t[k+1] - u)^2
    
    p_phi_L_fn <- function(u, k) if (k == 1) 0 else {
      sum(Q$pwt2[1:(k-1)]/(Q$t[1:(k-1)] - u)^2)
    }
    
    p_phi_R_fn <- function(u, k) if (k == (Q$K-1)) 0 else {
      sum(Q$pwt2[(k+2):Q$K]/(Q$t[(k+2):Q$K] - u)^2 )
    }
    
    #necessary (but not sufficient) condition for spectral separation between groups - returns boolean vector of length K-1, with element k indicating if there is spectral separation between eigenvectors t[k] and t[k+1]
    spec_sep <- function() {
      if (Q$K == 1) return (NULL) #if only one group nothing needed
      vapply(1:(Q$K-1), function(k)
        p_theta_fn(x_fn(k),k) + p_phi_L_fn(Q$t[k+1], k) + p_phi_R_fn(Q$t[k], k)  < Q$n, logical(1) )
    }
    
    p_phi_diff_fn <- function(u) 2*sum(Q$pwt2/(Q$t - u)^3)
    p_phi_fn <- function(u) sum(Q$pwt2/(Q$t - u)^2)
    
    #necessary and sufficient condition for spectral separation. returns vector with points of separation, with elements named by interval
    spec_sep2 <- function()
    {
      if (Q$K == 1) return (NULL) #if only one group nothing required
      
      spec_sep_out <- spec_sep()
      out <- c()
      for (k in which(spec_sep_out))
      {
        x_k <- x_fn(k)
        p_phi_diff_x_k <- p_phi_diff_fn(x_k)
        if (p_phi_diff_x_k == 0) out <- c(out, setNames(x_k,k) )
        else if (p_phi_diff_x_k > 0) {
          x_root <- optimize(f = function(x) (p_phi_diff_fn(x))^2, interval=c(Q$t[k], x_k), tol=1e-20)$minimum
          if (p_phi_fn(x_root) < Q$n) out <- c(out, setNames(x_root,k))
        }
        else {
          x_root <- optimize(f = function(x) (p_phi_diff_fn(x))^2, interval=c(x_k, Q$t[k+1]), tol=1e-20)$minimum
          if(p_phi_fn(x_root) < Q$n) out <- c(out, setNames(x_root,k))
        }
      }
      return ( out )
    }
    
    #find interior spectral boundaries, if any (and corresponding eigenvalue weights)
    spec_boundint <- function()
    {
      spec_sep_out <- spec_sep2()
      if (is.null(spec_sep_out)) return (NULL)
      K_subset <- as.integer(names(spec_sep_out))
      out <- do.call(rbind, lapply(K_subset, function(k)
      {
        x_mid <- spec_sep_out[as.character(k)]
        
        #bounds
        spec_k_lb <- optimize(f = function(x) (p_phi_fn(x) - Q$n)^2, interval=c(Q$t[k],x_mid), tol=1e-20)$minimum
        spec_k_ub <- optimize(f = function(x) (p_phi_fn(x) - Q$n)^2, interval=c(x_mid,Q$t[k+1]), tol=1e-20)$minimum
        
        #eigenvalue weights
        pw1 <- Q$pzw + sum(Q$pw[1:k])
        pw2 <- Q$p - pw1
        
        return(c(spec_k_lb, spec_k_ub, pw1, pw2))
      }) )
      rownames(out) <- K_subset
      return (out)
    }
    
    boundint <- spec_boundint()
    
    #all intervals of spectral support
    #returns matrix of intervals with rownames giving number of eigenvalues in each interval
    spec_supp <- function()
    {
      t1 <- min(Q$t)
      x_lb1 <- t1 - sqrt(1/Q$n * sum(Q$pwt2) ) - 1
      x_t1 <- optimize(f = function(x) (p_phi_fn(x) - Q$n)^2, interval=c(x_lb1, t1), tol=1e-20)$minimum
      
      tK <- max(Q$t)
      x_ubK <- tK + sqrt(1/Q$n * sum(Q$pwt2) ) + 1
      x_tK <- optimize(f = function(x) (p_phi_fn(x) - Q$n)^2, interval=c(tK, x_ubK), tol=1e-20)$minimum
      
      if (is.null(boundint)) {
        supp_intervals <- matrix(c(x_t1,x_tK), ncol=2)
        omega <- Q$p
      }
      else {
        supp_intervals <- matrix(c(x_t1, t(boundint[,1:2]), x_tK), byrow=TRUE, ncol=2)
        omega <- vapply(1:nrow(supp_intervals), function(i)
          sum(Q$pw[Q$t > supp_intervals[i,1] & Q$t < supp_intervals[i,2]] ), integer(1) )
      }
      
      rownames(supp_intervals) <- omega
      return(supp_intervals)
    }
    
    #support in u-space
    Q$u_Fbar <- spec_supp()
    Q$numint <- nrow(Q$u_Fbar)
    
    #modified Stieltjes transform of sample eigenvalues distribution
    m_LF_u_Fbar_fn <- function(u) 1/Q$p * sum(Q$pw*Q$t/(Q$t - u) )
    Q$m_LF_u_Fbar <- apply(Q$u_Fbar, c(1,2), m_LF_u_Fbar_fn)
    
    #support of sample eigenvalues distribution
    endpoint_fn <- function(u, m_LF_u) max(0, u - Q$c*u*m_LF_u)
    Q$endpoints <- matrix(mapply(endpoint_fn, Q$u_Fbar, Q$m_LF_u_Fbar), nrow=Q$numint)
    Q$endpoints[1,1] <- ifelse(Q$p==Q$n, 0, Q$endpoints[1,1] )
    
    #number of eigenvalues in each interval
    numeig_fn <- function()
    {
      if (Q$numint == 1) out <- Q$p - max(Q$pzw, Q$p - Q$n)
      else {
        out <- rep(0, Q$numint)
        out[1] <- boundint[1,3] - max(Q$pzw, Q$p - Q$n)
        if (Q$numint > 2) out[2:(Q$numint-1)] <- diff(boundint[1:(Q$numint - 1),3])
        out[Q$numint] <- boundint[(Q$numint-1),4]
      }
      return (as.integer(out) )
    }
    
    Q$numeig <- numeig_fn()
  }
  
  grid_fn <- function(Q)
  {
    minnxi = 100
    
    nxi0 = max(minnxi,min(sum(Q$pw),Q$n))
    Q$mult = ceiling(nxi0 / min(sum(Q$pw),Q$n) )
    Q$nxi = Q$mult*min(sum(Q$pw),Q$n)
    
    #grid points for each support (names give number of grid points in each interval - 2)
    Q$xi_list <- lapply(1:Q$numint, function(i) {
      lb <- Q$u_Fbar[i,1]; ub <- Q$u_Fbar[i,2]
      temp <- Q$numeig[i]*Q$mult
      return ( lb + (ub-lb)*(sin( pi/(2*(temp + 1)) * (1:temp) ) )^2  )
    })
  }
  
  sol_fn <- function(Q)
  {
    #function solves MP equation. y in [0, infty), xi is a grid point
    Gamma_fn <- function(y, xi)
      sum(Q$tau^2/( (Q$tau - xi)^2 + y^2 ) ) - Q$n
    
    #solve for roots at each grid point
    Q$zxi_list <- lapply(1:length(Q$xi_list), function(i) {
      vapply(1:(length(Q$xi_list[[i]]) ), function(j) {
        xi <- Q$xi_list[[i]][j]
        delta <- min( (Q$t - xi)^2 )
        y_ub <- sqrt(max(1/Q$n * sum(Q$pwt2) - delta, 0) ) + 1
        root <- optimize(f = function(y) (Gamma_fn(y, xi))^2, interval=c(0,y_ub), tol=1e-20)$minimum
        return (complex(real = xi, imaginary = root))}, complex(1))
    })
  }
  
  den_fn <- function(Q)
  {
    Q$a = 4 #exponent of power transform of abscissa
    
    m_LF_fn <- Vectorize ( function(z) 1/Q$p * sum(Q$tau/(Q$tau - z) ), vectorize.args = "z" )
    Q$m_LF_zxi_list <- lapply(1:length(Q$zxi_list), function(i) m_LF_fn(Q$zxi_list[[i]]) )
    
    Q$x_list <- lapply(1:length(Q$zxi_list), function(i) Re (Q$zxi_list[[i]] * (1 - Q$c*Q$m_LF_zxi_list[[i]]) ) )
    
    Q$f_list <- lapply(1:length(Q$zxi_list), function(i) 1/(pi*Q$c) * Im(Q$zxi_list[[i]]) / (abs(Q$zxi_list[[i]])^2) )
    
    Q$zeta_list <- lapply(1:length(Q$x_list), function(i) Q$x_list[[i]]^(1/Q$a) )
    
    Q$g_list <- lapply(1:length(Q$zeta_list), function(i) Q$a*Q$zeta_list[[i]]^(Q$a-1) * Q$f_list[[i]])
  }
  
  dis_fn <- function(Q)
  {
    Q$dis_x_list <- lapply(1:Q$numint, function(i) c(Q$endpoints[i,1], Q$x_list[[i]], Q$endpoints[i,2]))
    Q$dis_zeta_list <- lapply(1:Q$numint, function(i) c(Q$endpoints[i,1]^(1/Q$a), Q$zeta_list[[i]], Q$endpoints[i,2]^(1/Q$a) ) )
    Q$dis_m_LF_list <- lapply(1:Q$numint, function(i) c(Q$m_LF_u_Fbar[i,1], Q$m_LF_zxi_list[[i]], Q$m_LF_u_Fbar[i,2]))
    Q$dis_g_list <- lapply(1:Q$numint, function(i) c(0, Q$g_list[[i]], 0) )
    Q$F0 <- 1/Q$p * max(Q$pzw, (Q$p - Q$n) )
    Q$Fstart <- Q$F0 + c(0, cumsum(Q$numeig[-Q$numint]) ) / Q$p
    Q$Fend <- Q$F0 + cumsum(Q$numeig) / Q$p
    Q$intlength <- sapply(1:Q$numint, function(i) length(Q$dis_x_list[[i]]))
    Q$dis_G_raw <- lapply(1:Q$numint, function(i)
      c(0, cumsum(diff(Q$dis_zeta_list[[i]]) * (Q$dis_g_list[[i]][-1] + Q$dis_g_list[[i]][-Q$intlength[i]]) * 0.5 ) ) )
    Q$dis_G_list <- lapply(1:Q$numint, function(i) Q$Fstart[i] +(Q$Fend[i] - Q$Fstart[i]) * Q$dis_G_raw[[i]] / Q$dis_G_raw[[i]][Q$intlength[i]] )
  }
  
  lambda_fn <- function(Q)
  {
    Q$F <- lapply(1:Q$numint, function(i) unique(Q$dis_G_list[[i]]))
    Q$F_idx <- lapply(1:Q$numint, function(i) which(!duplicated(Q$dis_G_list[[i]] ) ) )
    Q$nidx <- sapply(1:Q$numint, function(i) length(Q$F_idx[[i]]))
    Q$x_F <- lapply(1:Q$numint, function(i) (Q$dis_zeta_list[[i]][Q$F_idx[[i]]])^(Q$a) )
    Q$x_F_mean <- lapply(1:Q$numint, function(i) 0.5*(Q$x_F[[i]][-1] + Q$x_F[[i]][-length(Q$x_F[[i]])] ) )
    Q$x_F_diff <- lapply(1:Q$numint, function(i) diff(Q$x_F[[i]]))
    Q$F_diff <- lapply(1:Q$numint, function(i) diff(Q$F[[i]]))
    Q$nquant <- Q$numeig + 1
    Q$quant <- lapply(1:Q$numint, function(i) seq(Q$F[[i]][1], Q$F[[i]][Q$nidx[i]],length.out = Q$nquant[i]))
    Q$bins <- lapply(1:Q$numint, function(i) findInterval(Q$quant[[i]], Q$F[[i]], rightmost.closed=TRUE))
    if (any(sapply(1:length(Q$bins), function(i) c(Q$bins[[i]][1], Q$bins[[i]][Q$nquant[i]] ) ) != rbind(rep(1,Q$numint), Q$nidx-1) ) ) stop("Unexpected bin allocation")
    Q$integral_indic <- lapply(1:Q$numint, function(i)
      1*( tcrossprod(rep(1,Q$nquant[[i]]), Q$F_idx[[i]][-1]) <= tcrossprod(Q$bins[[i]], rep(1, (Q$nidx[i]-1) ) ) ) )
    
    Q$lambda <- rep(0, Q$p)
    X_integral <- lapply(1:Q$numint, function(i) Q$F_diff[[i]] * Q$x_F_mean[[i]] )
    integral_j_kappa <- lapply(1:Q$numint, function(i) (Q$quant[[i]] - Q$F[[i]][Q$bins[[i]]])* ( Q$x_F[[i]][Q$bins[[i]]] +
                                                                                                   0.5* (Q$quant[[i]] - Q$F[[i]][Q$bins[[i]]]) * Q$x_F_diff[[i]][Q$bins[[i]]] / Q$F_diff[[i]][Q$bins[[i]]]) )
    X_kappa_integral <- lapply(1:Q$numint, function(i)
      rowSums ( tcrossprod(rep(1,Q$nquant[[i]]), X_integral[[i]]) * Q$integral_indic[[i]] ) + integral_j_kappa[[i]] )
    for (i in 1:Q$numint)
      Q$lambda[round(Q$F[[i]][1] * Q$p + 1):round(Q$F[[i]][Q$nidx[i]] * Q$p)] <- diff(X_kappa_integral[[i]]) * Q$p
  }
  
  sup_fn(Q)
  grid_fn(Q)
  sol_fn(Q)
  den_fn(Q)
  dis_fn(Q)
  lambda_fn(Q)
  
  return (Q)
}

#compute gradient of QuEST function
get_lambda_J <- function(Q)
{
  rep1p <- rep(1,Q$p)
  
  #Jacobian of support endpoints
  sup_J_fn <- function(Q)
  {
    u_vec <- c(t(Q$u_Fbar))
    Dr <- apply(Q$u_Fbar, c(1,2), function(u) sum(Q$tau^2 / (Q$tau-u)^3 ) )
    temp <- (array(1,dim=dim(Q$u_Fbar)) %o% Q$tau) - (Q$u_Fbar %o% rep1p)
    Q$sup_J <- (Q$u_Fbar %o% Q$tau) / (temp^3 * (Dr %o% rep1p) )
    
    Q$endpoints_J <- 1/Q$n * (Q$u_Fbar^2 %o% rep(1,Q$p)) / temp^2
    if (Q$p == Q$n) Q$endpoints_J[1,1,] <- 0
  }
  
  sup_J_fn(Q)
  
  #Jacobian of grid points in u-space
  Q$xi_J_list <- lapply(1:Q$numint, function(i) {
    length_temp <- length(Q$xi_list[[i]])
    temp <- (sin( pi/(2*(length_temp + 1)) * (1:length_temp) ) )^2
    return( tcrossprod((1 - temp), Q$sup_J[i,1,]) + tcrossprod(temp, Q$sup_J[i,2,]) )
  })
  
  den_J_fn <- function(Q)
  {
    #some objects to be used to in following functions
    TAU_list <- lapply(1:Q$numint, function(i) tcrossprod(rep(1,length(Q$xi_list[[i]])), Q$tau) )
    TAU2_list <- lapply(1:Q$numint, function(i) tcrossprod(rep(1,length(Q$xi_list[[i]])), Q$tau^2) )
    ZXI_list <- lapply(1:Q$numint, function(i) tcrossprod(Q$zxi_list[[i]], rep1p ) )
    XI_list <- lapply(1:Q$numint, function(i) tcrossprod(Re(Q$zxi_list[[i]]), rep1p ) )
    YXI_list <- lapply(1:Q$numint, function(i) tcrossprod(Im(Q$zxi_list[[i]]), rep1p ) )
    YXI2_list <- lapply(1:Q$numint, function(i) tcrossprod(Im(Q$zxi_list[[i]])^2, rep1p ) )
    TAU_XI_list <- lapply(1:Q$numint, function(i) TAU_list[[i]] - XI_list[[i]])
    
    #Jacobian of yxi := imag(zxi)
    Q$yxi_J_list <- lapply(1:Q$numint, function(k) {
      NORM <- TAU_XI_list[[k]]^2 + YXI2_list[[k]]; NORM2 <- NORM*NORM
      Nr <- rowSums(TAU2_list[[k]] * TAU_XI_list[[k]] / NORM2)
      Dr <- rowSums(TAU2_list[[k]] * YXI_list[[k]] / NORM2)
      return ( (TAU_list[[k]] / NORM - TAU2_list[[k]] * TAU_XI_list[[k]] / NORM2 + tcrossprod(Nr, rep1p) * Q$xi_J_list[[k]])  / tcrossprod(Dr, rep1p) )
    })
    
    Q$zxi_J_list <- lapply(1:Q$numint, function(k) Q$xi_J_list[[k]] + 1i*Q$yxi_J_list[[k]])
    
    Q$f_J_list <- lapply(1:Q$numint, function(k) 
      1/(pi*Q$c) * Im(Q$zxi_J_list[[k]] * tcrossprod(1/(Q$zxi_list[[k]]*Q$zxi_list[[k]]), rep1p) ) )
    
    Q$m_LF_J_list <- lapply(1:Q$numint, function(k) {
      TAU_ZXI2 <- (TAU_list[[k]] - ZXI_list[[k]])^2
      1/Q$p * (-ZXI_list[[k]] / TAU_ZXI2 + Q$zxi_J_list[[k]] * tcrossprod(rowSums(TAU_list[[k]] / TAU_ZXI2 ),rep1p) )
    })
    
    Q$x_J_list <- lapply(1:Q$numint, function(k) 
      Re ( Q$zxi_J_list[[k]] * (1 - Q$c*tcrossprod(Q$m_LF_zxi_list[[k]], rep1p) ) - Q$c * ZXI_list[[k]] * Q$m_LF_J_list[[k]] ) )
    
    Q$zeta_J_list <- lapply(1:Q$numint, function(k) 
      1/Q$a * Q$x_J_list[[k]] * tcrossprod(Q$x_list[[k]]^(1/Q$a - 1), rep1p) )
    
    Q$g_J_list <- lapply(1:Q$numint, function(k) 
      Q$a * ( (Q$a - 1) * Q$zeta_J_list[[k]] * tcrossprod(Q$f_list[[k]]*Q$zeta_list[[k]]^(Q$a - 2), rep1p) + 
                tcrossprod(Q$zeta_list[[k]]^(Q$a - 1), rep1p) * Q$f_J_list[[k]] ) )
  }
  
  den_J_fn(Q)
  
  dis_J_fn <- function(Q)
  {
    Q$dis_g_J_list <- lapply(1:Q$numint, function(k) rbind(rep(0,Q$p), Q$g_J_list[[k]], rep(0,Q$p)))
    
    Q$dis_zeta_J_list <- lapply(1:Q$numint, function(k) 
      rbind(1/Q$a * Q$endpoints[k,1]^(1/Q$a - 1) * Q$endpoints_J[k,1,],
            Q$zeta_J_list[[k]],
            1/Q$a * Q$endpoints[k,2]^(1/Q$a - 1) * Q$endpoints_J[k,2,]) )
    if (Q$p == Q$n) Q$dis_zeta_J_list[[1]][1,] <- 0
    
    dis_zeta_J_diff <- lapply(1:Q$numint, function(k) Q$dis_zeta_J_list[[k]][-1,] - Q$dis_zeta_J_list[[k]][-Q$intlength[k],])
    dis_zeta_diff <- lapply(1:Q$numint, function(k) diff(Q$dis_zeta_list[[k]]) )
    dis_g_mean <- lapply(1:Q$numint, function(k) 0.5*(Q$dis_g_list[[k]][-1] + Q$dis_g_list[[k]][-Q$intlength[k]]) )
    dis_g_J_mean <- lapply(1:Q$numint, function(k) 0.5*(Q$dis_g_J_list[[k]][-1,] + Q$dis_g_J_list[[k]][-Q$intlength[k],] ) )
    
    G_J_raw_list <- lapply(1:Q$numint, function(k) rbind(rep(0,Q$p), 
                                                         apply(dis_zeta_J_diff[[k]] * tcrossprod(dis_g_mean[[k]], rep1p), 2, cumsum) + 
                                                           apply(tcrossprod(dis_zeta_diff[[k]], rep1p) * dis_g_J_mean[[k]], 2, cumsum) ) )
    
    Fends_diff <- Q$Fend - Q$Fstart
    Q$dis_G_J <- lapply(1:Q$numint, function(k)
      Fends_diff[k] / Q$dis_G_raw[[k]][Q$intlength[k]] * (G_J_raw_list[[k]] - 
                                                            tcrossprod(Q$dis_G_raw[[k]], rep1p) * tcrossprod(rep(1,Q$intlength[k]), G_J_raw_list[[k]][Q$intlength[k],]) ) )
  }
  
  dis_J_fn(Q)
  
  lambda_J_fn <- function(Q)
  {
    Q$lambda_J <- matrix(0, nrow=Q$p, ncol = Q$p)
    if ((Q$p - Q$n) <= Q$pzw & Q$pzw > 0) 
      Q$lambda_J[max(1,Q$p - Q$n + 1):Q$pzw, 1:Q$pzw] <- 1 - Q$c
    Q$F_J <- lapply(1:Q$numint, function(k) Q$dis_G_J[[k]][Q$F_idx[[k]],])
    Q$x_J_F <- lapply(1:Q$numint, function(k) 
      tcrossprod((Q$a*Q$dis_zeta_list[[k]][Q$F_idx[[k]]]^(Q$a-1)), rep1p) * Q$dis_zeta_J_list[[k]][Q$F_idx[[k]],] )
    
    F_J_diff <- lapply(1:Q$numint, function(k) Q$F_J[[k]][-1,] - Q$F_J[[k]][-Q$nidx[k],])
    x_J_F_mean <- lapply(1:Q$numint, function(k) 0.5*(Q$x_J_F[[k]][-1,] + Q$x_J_F[[k]][-Q$nidx[k],]) )
    x_J_F_diff <- lapply(1:Q$numint, function(k) Q$x_J_F[[k]][-1,] - Q$x_J_F[[k]][-Q$nidx[k],])
    X_J_integral <- lapply(1:Q$numint, function(k) 
      F_J_diff[[k]] * tcrossprod(Q$x_F_mean[[k]], rep1p) + tcrossprod(Q$F_diff[[k]], rep1p) * x_J_F_mean[[k]] )
    
    integral_j_kappa <- lapply(1:Q$numint, function(k) {
      q <- Q$quant[[k]]
      b <- Q$bins[[k]]
      tcrossprod(q,rep1p) * Q$x_J_F[[k]][b,] - Q$x_J_F[[k]][b,] * tcrossprod(Q$F[[k]][b], rep1p) -
        Q$F_J[[k]][b,] * tcrossprod(Q$x_F[[k]][b], rep1p) - 
        Q$F_J[[k]][b,] * tcrossprod(((q - Q$F[[k]][b])*Q$x_F_diff[[k]][b] / Q$F_diff[[k]][b]), rep1p) +
        x_J_F_diff[[k]][b,] * tcrossprod( (0.5 * (q - Q$F[[k]][b])^2 / Q$F_diff[[k]][b]), rep1p) -
        F_J_diff[[k]][b,] * tcrossprod( (0.5 * (q - Q$F[[k]][b])^2 * Q$x_F_diff[[k]][b] / Q$F_diff[[k]][b]^2), rep1p) })
    
    X_J_kappa_integral <- lapply(1:Q$numint, function(i) Q$integral_indic[[i]] %*% X_J_integral[[i]] + integral_j_kappa[[i]])
    for (i in 1:Q$numint)
      Q$lambda_J[round(Q$F[[i]][1] * Q$p + 1):round(Q$F[[i]][Q$nidx[i]] * Q$p), ] <- apply(X_J_kappa_integral[[i]], 2, diff) * Q$p
  }
  
  lambda_J_fn(Q)
  return (Q$lambda_J)
}

nlshrink_est <- function(Q) {
  delta <- rep(0,Q$p)
  if ((Q$p - Q$n) > Q$pzw) {
    lb <- min(Q$tau) - 1/ Q$n * sum(Q$tau) - 1
    ub <- 1/Q$p * (Q$p - Q$pw[1])*Q$t[1]
    optim_fn <- function(u) (sum(Q$pw*Q$t / (Q$t - u)) - Q$n)^2
    u_Fbar0 <- optimize(f = optim_fn, interval = c(lb,ub))$minimum
    null_adj <- u_Fbar0 / (1 - Q$c)
    delta[1:(Q$p - Q$n)] <- null_adj
  }
  
  lim0 <- ifelse(Q$p == Q$n & all(Q$tau > 0), Q$n / sum(Q$pw / Q$t), 0)
  
  out <- lapply(1:Q$numint, function(i) {
    y <- rep(0, Q$nidx[i])
    Dr <- abs(1 - Q$c*Q$dis_m_LF_list[[i]][Q$F_idx[[i]]])^2
    if (Q$p == Q$n) {
      j = (Q$x_F[[i]] == 0 & Dr == 0)
      y[j] <- lim0
      y[!j] <- Q$x_F[[i]][!j] / Dr[!j]
    } else {
      y <- Q$x_F[[i]] / Dr
    }
    F_temp <- Q$F[[i]]
    b <- Q$bins[[i]]
    F_temp_b <- F_temp[b]
    integral_ydF <- diff(F_temp) * 0.5 * (y[-1] + y [-Q$nidx[i]])
    integral_ydF2 <- (Q$quant[[i]] - F_temp_b) * (y[b] +
                                                    0.5 * (Q$quant[[i]] - F_temp_b) * (y[b + 1] - y[b]) / (F_temp[b+1] - F_temp_b) )
    integral_ydF3 <- rowSums(tcrossprod(rep(1,Q$nquant[i]), integral_ydF) * Q$integral_indic[[i]]) + integral_ydF2
    return (diff(integral_ydF3) * Q$p)
  })
  for (i in 1:Q$numint)
    delta[round(Q$F[[i]][1]*Q$p + 1):round(Q$F[[i]][Q$nidx[i]] * Q$p)] <- out[[i]]
  if (Q$p == Q$n)
    delta[delta < lim0] <- lim0
  
  return (delta)
}
