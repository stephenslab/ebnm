test_that("sampler gives reasonable results", {
  set.seed(1)
  mu = c(rep(0,50), rexp(50))
  s = rgamma(50,1,1)
  x = mu + rnorm(50,0,s)

  res = ebnm_point_normal(x,s, output=c("result", "fitted_g", "post_sampler"))

  samp <- res$post_sampler(1000)
  expect_equal(colMeans(samp), res$result$PosteriorMean, tolerance = 0.1)
  expect_equal(colMeans(samp^2), res$result$PosteriorMean2, tolerance = 0.1)
})
