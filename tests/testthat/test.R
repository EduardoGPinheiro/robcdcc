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
  
  expect_equal(min(robust_Rt$w), 0.54528087)
  expect_equal(robust_Rt$Rt[1,2], -0.082167733)
  expect_equal(diag(robust_Rt$Rt), c(1, 1))
})