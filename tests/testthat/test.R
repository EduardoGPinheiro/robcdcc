rt = simCDCC_C(phi=c(.1, .8), S=diag(2), nobs=1000, seed=1)
phi_results = estimateCDCC(rt)
robust_phi_results = robust_estimateCDCC(rt, delta=.975)

test_that("parameters estimation", {
  expect_equal(phi_results$phi, c(0.09859258, 0.74636124))
  expect_equal(phi_results$eta$alpha, c(0.0343791365913, 0.0000100205095))
  expect_equal(phi_results$eta$beta, c(0.740706442, 0.463080309))
  
  expect_equal(robust_phi_results$phi, c(0.07297929, 0.79259606))
  expect_equal(robust_phi_results$eta$alpha, c(0.001300310921, 0.000010032098))
  expect_equal(robust_phi_results$eta$beta, c(0.961287836, 0.957906203))
})

test_that("ht estimation", {
  Dt = calc_Dt(rt=rt, eta=phi_results$eta)
  rDt_lst = robust_calc_Dt(rt=rt, eta=robust_phi_results$eta, delta = .975) 
  rDt = rDt_lst$Dt
  wmin = rDt_lst$w %>% apply(min, MARGIN=2)
  
  expect_true((phi_results$Dt == Dt) %>% apply(all, MARGIN = 2) %>% all)
  expect_true((robust_phi_results$Dt == rDt) %>% apply(all, MARGIN = 2) %>% all)
  expect_equal(wmin, c(0.38883892, 0.44211281))
  
  expect_true((nrow(Dt) - nrow(rt)) == 1)
  expect_true((nrow(rDt) - nrow(rt)) == 1)
})

test_that("Rt estimation", {
  Rt = calc_Rt(rt=phi_results$epsilon, phi=phi_results$phi, S=diag(2))
  
  robust_Rt = robust_calc_Rt(
    rt = robust_phi_results$epsilon,
    phi = robust_phi_results$phi,
    S = diag(2),
    delta = .975
  )
  
  expect_equal(Rt[1,2], -0.099229167)
  expect_equal(diag(Rt), c(1, 1))
  
  expect_equal(min(robust_Rt$w), 0.54527522)
  expect_equal(robust_Rt$Rt[1,2], -0.082190407)
  expect_equal(diag(robust_Rt$Rt), c(1, 1))
})

test_that("portfolio results", {
 results_lst = geral_calc_portfolio_variance(
    phi = c(.1, .8),
    q_phi = phi_results$phi,
    r_phi = robust_phi_results$phi,
    rt = rt,
    burn_rt = rt[1:2, ],
    r_rt = robust_phi_results$epsilon,
    q_rt = phi_results$epsilon,
    S = diag(2),
    Dt = diag(2),
    q_Dt = phi_results$Dt %>% tail(1) %>% as.vector %>% diag,
    r_Dt = robust_phi_results$Dt %>% tail(1) %>% as.vector %>% diag,
    q_S = diag(2),
    r_S = diag(2),
    delta = .975
  )
 
 q_Dt = phi_results$Dt %>% tail(1) %>% as.vector %>% diag %>% sqrt
 r_Dt = robust_phi_results$Dt %>% tail(1) %>% as.vector %>% diag %>% sqrt
 
 # qmvn
 qRt = diag(1/diag(q_Dt)) %*% results_lst$q_Ht %*% diag(1/diag(q_Dt))
 validate_qRt = calc_Rt(rt=phi_results$epsilon, phi=phi_results$phi, S=diag(2))
 
 # robust
 rRt = diag(1/diag(r_Dt)) %*% results_lst$r_Ht %*% diag(1/diag(r_Dt))
 validate_rRt = robust_calc_Rt(rt=robust_phi_results$epsilon, 
                               phi=robust_phi_results$phi, S=diag(2))$Rt
 
 expect_equal(qRt, validate_qRt)
 expect_equal(rRt, validate_rRt)
})

test_that("devolatilized returns", {
  std_rt = calc_devolatilized_returns(rt=rt, eta_df=phi_results$eta)
  robust_std_rt = calc_robust_devolatilized_returns(
    rt=rt, eta_df=robust_phi_results$eta)

  expect_equal(std_rt, phi_results$epsilon)
  expect_equal(robust_std_rt, robust_phi_results$epsilon)
})
