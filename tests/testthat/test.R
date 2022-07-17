# Parameters
phi = c(.1, .8)

# Estimation
rt = simCDCC_C(phi=phi, S=diag(2), nobs=1000, seed=1)
q_results = estimateCDCC(rt)
r_results = robust_estimateCDCC(rt, delta=.975)

test_that("phi estimation", {
  q_phi_ape = abs(q_results$phi - phi) / phi
  r_phi_ape = abs(r_results$phi - phi) / phi
  
  expect_true(max(q_phi_ape) < .5)
  expect_true(max(r_phi_ape) < .5)
})

test_that("ht estimation", {
  Dt = calc_Dt(rt=rt, eta=q_results$eta)
  rDt_lst = robust_calc_Dt(rt=rt, eta=r_results$eta, delta = .975) 
  rDt = rDt_lst$Dt
  wmin = rDt_lst$w %>% apply(min, MARGIN=2)
  
  expect_true((q_results$Dt == Dt) %>% apply(all, MARGIN = 2) %>% all)
  expect_true((r_results$Dt == rDt) %>% apply(all, MARGIN = 2) %>% all)
  
  expect_true((nrow(Dt) - nrow(rt)) == 1)
  expect_true((nrow(rDt) - nrow(rt)) == 1)
})

test_that("Rt estimation", {
  Rt = calc_Rt(rt=q_results$epsilon, phi=q_results$phi, S=diag(2))
  
  robust_Rt = robust_calc_Rt(
    rt = r_results$epsilon,
    phi = r_results$phi,
    S = diag(2),
    delta = .975
  )
  
  expect_equal(diag(Rt), c(1, 1))
  expect_equal(diag(robust_Rt$Rt), c(1, 1))
})

test_that("portfolio results", {
 results_lst = geral_calc_portfolio_variance(
    phi = c(.1, .8),
    q_phi = q_results$phi,
    r_phi = r_results$phi,
    rt = rt,
    burn_rt = rt[1:2, ],
    r_rt = r_results$epsilon,
    q_rt = q_results$epsilon,
    S = diag(2),
    Dt = diag(2),
    q_Dt = q_results$Dt %>% tail(1) %>% as.vector %>% diag,
    r_Dt = r_results$Dt %>% tail(1) %>% as.vector %>% diag,
    q_S = diag(2),
    r_S = diag(2),
    delta = .975
  )
 
 q_Dt = q_results$Dt %>% tail(1) %>% as.vector %>% diag %>% sqrt
 r_Dt = r_results$Dt %>% tail(1) %>% as.vector %>% diag %>% sqrt
 
 # qmvn
 qRt = diag(1/diag(q_Dt)) %*% results_lst$q_Ht %*% diag(1/diag(q_Dt))
 validate_qRt = calc_Rt(rt=q_results$epsilon, phi=q_results$phi, S=diag(2))
 
 # robust
 rRt = diag(1/diag(r_Dt)) %*% results_lst$r_Ht %*% diag(1/diag(r_Dt))
 validate_rRt = robust_calc_Rt(rt=r_results$epsilon, 
                               phi=r_results$phi, S=diag(2))$Rt
 
 expect_equal(qRt, validate_qRt)
 expect_equal(rRt, validate_rRt)
})

test_that("devolatilized returns", {
  std_rt = calc_devolatilized_returns(rt=rt, eta_df=q_results$eta)
  robust_std_rt = calc_robust_devolatilized_returns(
    rt=rt, eta_df=r_results$eta)

  expect_equal(std_rt, q_results$epsilon)
  expect_equal(robust_std_rt, r_results$epsilon)
})

test_that("real portfolio variance", {
  results_lst = geral_calc_portfolio_variance(
    phi = c(.1, .8),
    q_phi = q_results$phi,
    r_phi = r_results$phi,
    rt = rt,
    burn_rt = rt[1:2, ],
    r_rt = r_results$epsilon,
    q_rt = q_results$epsilon,
    S = diag(2),
    Dt = diag(2),
    q_Dt = q_results$Dt %>% tail(1) %>% as.vector %>% diag,
    r_Dt = r_results$Dt %>% tail(1) %>% as.vector %>% diag,
    q_S = diag(2),
    r_S = diag(2),
    delta = .975
  )
  
  real_var = calc_real_portfolio_variance_C(
    phi = c(.1, .8),
    rt = rt,
    burn_rt = rt[1:2, ],
    S = diag(2),
    Dt = diag(2)
  )
  
  Ht = results_lst$Ht
  validate_real_var = 1 / sum(1/2 * diag(solve(Ht)))
  expect_equal(validate_real_var, real_var)
})
