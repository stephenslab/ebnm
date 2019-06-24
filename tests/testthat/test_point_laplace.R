context("Point Laplace")

n <- 1000
set.seed(1)
s <- rnorm(n, 1, 0.1)
x <- c(rexp(n / 2, rate = 0.1), rep(0, n / 2)) + rnorm(n, sd = s)

true_pi0 <- 0.5
true_scale <- 10
true_g <- laplacemix(pi = c(true_pi0, 1 - true_pi0),
                     mean = rep(0, 2),
                     scale = c(0, true_scale))

pl.res <- ebnm(x, s, prior_type = "point_laplace")

test_that("Basic functionality works", {
  pl.res2 <- ebnm_point_laplace(x, s)
  expect_identical(pl.res, pl.res2)
  expect_equal(pl.res$fitted_g, true_g, tolerance = 0.1)
})

test_that("Fixing the scale works", {
  pl.res2 <- ebnm_point_laplace(x, s, scale = true_scale)
  expect_equal(pl.res2$fitted_g, true_g, tolerance = 0.1)
  expect_identical(pl.res2$fitted_g$scale[2], true_scale)
})

test_that("Fixing g works", {
  pl.res2 <- ebnm_point_laplace(x, s, g_init = pl.res$fitted_g, fix_g = TRUE)
  expect_identical(pl.res$fitted_g, pl.res2$fitted_g)
  expect_equal(pl.res$loglik, pl.res2$loglik)
})

test_that("Initializing g works", {
  pl.res2 <- ebnm_point_laplace(x, s, g_init = true_g)
  expect_equal(pl.res$loglik, pl.res2$loglik)
})

test_that("Output parameter works", {
  pl.res <- ebnm_point_laplace(x, s, output = "fitted_g")
  expect_identical(names(pl.res), "fitted_g")
})
