test_that("postmean gives same result as ashr",{
  n = 100
  set.seed(1)
  s = rgamma(n,1,1)
  x = rnorm(n,0,1+s)
  ebnm.res = ebnm_point_normal(x,s)
  pi0 = ebnm.res$fitted_g$pi0
  a = ebnm.res$fitted_g$a
  ash.res = ashr::ash(x,s,outputlevel = 5,method = "shrink",
                      g = ashr::normalmix(c(pi0,1-pi0),c(0,0),c(0,sqrt(1/a))),
                      fixg = TRUE)
  expect_equal(ebnm.res$result$PosteriorMean,
               ash.res$flash_data$postmean,tol = 1e-6)
  expect_equal(ebnm.res$result$PosteriorMean2,
               ash.res$flash_data$postmean2, tol = 1e-6)
  expect_equal(ebnm.res$loglik,ash.res$flash_data$penloglik,tol = 1e-6)

  ebnm.res2 = ebnm_point_normal(x,s,norm=1)
  expect_equal(ebnm.res2$loglik,ebnm.res$loglik)
  expect_equal(ebnm.res2$result,ebnm.res$result)
  expect_equal(ebnm.res2$fitted_g$a,ebnm.res$fitted_g$a)

  # check that fixing g works as intended
  g = list(pi0 = 0, a = 0.5)
  ebnm.res3 = ebnm_point_normal(x, s, g, fixg=T)
  expect_identical(ebnm.res3$fitted_g, g)

  # check fixing pi0
  g = list(pi0 = 0.2, a = 0.5)
  ebnm.res4 = ebnm_point_normal(x, s, g, fix_pi0 = TRUE)
  expect_identical(ebnm.res4$fitted_g$pi0, g$pi0)
  expect_false(ebnm.res4$fitted_g$a == g$a)
  g = list(pi0 = 0)
  ebnm.res5 = ebnm_point_normal(x, s, g, fix_pi0 = TRUE)
  ebnm.res6 = ebnm_normal(x, s)
  expect_identical(ebnm.res5, ebnm.res6)

  # test control parameter
  ebnm.res7 = ebnm_point_normal(x, s, control=list(factr=1000))

  # test infinite and zero SEs
  x = c(rep(0, 5), rep(1, 5))
  s = rep(1, 10)
  s[10] = Inf
  expect_error(ebnm_point_normal(x, s))
  s[5] = s[10] = 0
  ebnm.res8 = ebnm_point_normal(x, s)
  expect_equal(ebnm.res8$result$PosteriorMean[c(5, 10)], x[c(5, 10)])
  expect_equal(ebnm.res8$result$PosteriorMean2[c(5, 10)], x[c(5, 10)]^2)
})
