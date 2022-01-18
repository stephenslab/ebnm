context("NPMLE")

n <- 100
set.seed(1)
s <- rnorm(n, 1, 0.1)
x <- runif(n, -10, 10) + rnorm(n, sd = s)

cdf_grid <- seq(-10, 10, 0.1)

true_g <- ashr::unimix(1, -10, 10)
true_cdf <- ashr::comp_cdf(true_g, cdf_grid)

test_that("Basic functionality works", {
  npmle.res <- ebnm(x, s, prior_family = "npmle")
  npmle.res2 <- ebnm_npmle(x, s)
  npmle.res$call <- npmle.res2$call <- NULL
  expect_identical(npmle.res, npmle.res2)

  est_cdf <- drop(npmle.res$fitted_g$pi %*% ashr::comp_cdf(npmle.res$fitted_g, cdf_grid))
  expect_equal(true_cdf, est_cdf, tolerance = 0.2)
})

test_that("Fixing the scale works", {
  npmle.res <- ebnm_npmle(x, s, scale = 1)
  g_scale <- npmle.res$fitted_g$mean[2] - npmle.res$fitted_g$mean[1]
  expect_equal(1, g_scale, tolerance = 0.1)
})

test_that("Fixing g works", {
  g_init = normalmix(rep(0.2, 5), seq(-10, 10, by = 5), 0)
  npmle.res <- ebnm_npmle(x, s, g_init = g_init, fix_g = TRUE)
  expect_identical(npmle.res[[g_ret_str()]], g_init)
})

test_that("Gaussian grid is selected when the range of x is large", {
  x <- 100 * rcauchy(n)
  npmle.res <- ebnm_npmle(x, s)
  expect_true(npmle.res$fitted_g$sd[1] > 0)
})
