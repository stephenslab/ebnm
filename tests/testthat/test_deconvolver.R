context("deconvolveR")

n <- 1000
set.seed(1)
x <- runif(n, -10, 10) + rnorm(n)
s <- 1

cdf_grid <- seq(-10, 10, 0.1)

true_g <- ashr::unimix(1, -10, 10)
true_cdf <- ashr::comp_cdf(true_g, cdf_grid)

test_that("Basic functionality works", {
  deconv.res <- ebnm(x, s, prior_family = "deconvolver")
  deconv.res2 <- ebnm_deconvolver(x, s)
  deconv.res$call <- deconv.res2$call <- NULL
  expect_identical(deconv.res, deconv.res2)

  est_cdf <- drop(deconv.res$fitted_g$pi %*% ashr::comp_cdf(deconv.res$fitted_g, cdf_grid))
  expect_equal(true_cdf, est_cdf, tolerance = 0.1)
})

test_that("Fixing the scale works", {
  deconv.res <- ebnm_deconvolver(x, s, scale = 1)
  g_scale <- deconv.res$fitted_g$mean[2] - deconv.res$fitted_g$mean[1]
  expect_equal(1, g_scale, tolerance = 0.1)
})

test_that("Fixing g works", {
  g_init = normalmix(rep(0.2, 5), seq(-10, 10, by = 5), 0)
  deconv.res <- ebnm_deconvolver(x, s, g_init = g_init, fix_g = TRUE)
  expect_identical(deconv.res[[g_ret_str()]], g_init)
})

test_that("deconv and nlm parameters get passed in", {
  deconv.res <- ebnm_deconvolver(x, s)
  deconv.res2 <- ebnm_deconvolver(x, s, c0 = 2, pDegree = 4, control = list(steptol = 1e-4))
  expect_false(identical(deconv.res$fitted_g$pi, deconv.res2$fitted_g$pi))
})
