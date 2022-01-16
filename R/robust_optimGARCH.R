medianB = function(r, k=30) {
  n = length(r)
  med = c()
  mad = c()
  w = matrix(NA, ncol = 2, nrow = n)
  for (i in 1:n) {
    med[i] = ifelse(i < k / 2, median(r[1:(k + 1)]),
                    ifelse(i > n - k / 2, median(r[(n - k):n]),
                           median(r[(i - k / 2):(i + k / 2)])))
    mad[i] = ifelse(i < k / 2, median(abs(r[1:(k + 1)] - med[i])),
                    ifelse(i > n - k / 2, median(abs(r[(n - k):n] - med[i])),
                           median(abs(r[(i - k / 2):(i + k / 2)] - med[i]))))
    w[i, 1] = ifelse(i <= k / 2, 1, ifelse(i > n - k / 2, (n - k), (i -
                                                                      k / 2)))
    w[i, 2] = ifelse(i <= k / 2, k + 1, ifelse(i > n - k / 2, n, (i + k /
                                                                    2)))
  }
  return(list(med, mad, w))
}

robust_optimGARCH = function(rt, cy, chisq, k=30) {
  AUX = medianB(rt, k=k)
  Med = AUX[[1]]
  MAD = AUX[[2]]
  I = (rt - Med) ^ 2 / (1.486 * MAD) ^ 2 <= 3.841459
  mu_R = sum(rt * I) / sum(I)
  J = (rt - mu_R) ^ 2 / (1.486 * MAD) ^ 2 <= 3.841459
  hat = 1.318 * sum((rt - mu_R) ^ 2 * J) / sum(J)
  
  rt = rt - mu_R
  ra <- matrix(c(1, 0, 0, 1, -1, -1), ncol = 2, byrow = TRUE)
  rb <- c(0.00001, 0.00001, -0.9999)
  nobs = length(rt)
  
  opt = stats::constrOptim(
    theta = c(.05, .93),
    f = robust_loglikelihoodGARCH,
    grad = robust_gradientGARCH,
    ui = ra,
    ci = rb,
    mu = 1e-5,
    outer.iterations = 400,
    outer.eps = 1e-07,
    rt = rt,
    h = hat,
    cy = cy,
    chisq = chisq
  )$par
  
  par = c(mu_R, hat * (1 - opt[1] - opt[2]), opt[1], opt[2], hat)
  ht = robcdcc::robust_calc_ht_C(par[2], par[3], par[4], rt, nobs, cy, chisq)
  std_residuals = rt / sqrt(ht[1:nobs])
  par_df = data.frame(
    mu = par[1],
    omega = par[2],
    alpha = par[3],
    beta = par[4],
    h = par[5]
  )
  return(list(par_df = par_df, std_residuals = std_residuals, ht=ht))
}
