test_that("'ebnm' calls give same results as 'ebnm_point_normal' or 'ebnm_point_laplace' calls",{
  n = 100
  set.seed(1)
  s = rgamma(n,1,1)
  x = rnorm(n,0,1+s)
  
  ### test point_normal
  # fix nothing
  ebnm.res1 = ebnm(x, s, "point_normal", g = NULL, fix_mu = F)
  ebnm.pn.res1 = ebnm_point_normal(x, s)
  expect_equal(ebnm.res1, ebnm.pn.res1)
  # fix pi0
  ebnm.res2 = ebnm(x, s, "point_normal", g = list(pi0 = .2), fix_pi0 = T, fix_mu = F)
  ebnm.pn.res2 = ebnm_point_normal(x, s, g = list(pi0 = .2), fix_pi0 = T)
  expect_equal(ebnm.res2, ebnm.pn.res2)
  # fix mu
  ebnm.res3 = ebnm(x, s, "point_normal")
  ebnm.pn.res3 = ebnm_point_normal(x, s, g = list(mu = 0), fix_mu = T)
  expect_equal(ebnm.res3, ebnm.pn.res3)
  # fix both
  ebnm.res4 = ebnm(x, s, "point_normal", g = list(pi0 = .2, mu = 0), fix_pi0 = T, fix_mu = T)
  ebnm.pn.res4 = ebnm_point_normal(x, s, g = list(pi0 = .2, mu = 0), fix_pi0 = T, fix_mu = T)
  expect_equal(ebnm.res4, ebnm.pn.res4)
  
  ### test point_laplace
  ebnm.res5 = ebnm(x, s, "point_laplace")
  ebnm.pl.res5 = ebnm_point_laplace(x, s)
  ebnm.pl.res5$fitted_g$mu = 0 # have to manually add this for test. NOTE: point_laplace will likely change, and so too will this check
  expect_equal(ebnm.res5, ebnm.pl.res5)
})
