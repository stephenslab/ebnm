context("Point normal sampler")

set.seed(1)
mu <- c(rep(0, 25), rexp(25))
s <- rgamma(50, 1, 1)
x <- mu + rnorm(50, 0, s)

res <- ebnm_point_normal(x, s,
                         output = c("result", "fitted_g", "post_sampler"))
samp <- res$post_sampler(1000)

test_that("point-normal sampler gives reasonable results", {
  expect_equal(colMeans(samp), res$result$posterior_mean, tol = 0.1)
  expect_equal(colMeans(samp^2), res$result$posterior_mean2, tol = 0.1)
})
