test_that("postmean gives same result as ashr",{
  n=100
  set.seed(1)
  s = rgamma(n,1,1)
  x=rnorm(n,0,1+s)
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
})
