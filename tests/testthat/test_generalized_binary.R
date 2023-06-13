context("generalized binary")

set.seed(1)
theta <- c(pmax(rnorm(100, 2, 0.1), 1), rep(0, 400))
s <- rep(0.5, length(theta))
x <- rnorm(n=length(theta), mean=theta, sd=s)

test_that("Basic functionality works", {
  gb.res <- ebnm(x, s, prior_family = "generalized_binary")
  gb.res2 <- ebnm_generalized_binary(x, s)
  gb.res$call <- gb.res2$call <- NULL
  expect_identical(gb.res, gb.res2)
  expect_equal(gb.res$fitted_g$pi[2], 0.2, tolerance = 0.1)
  expect_equal(gb.res$fitted_g$mean[2], 2, tolerance = 0.1)
})

test_that("Changing the scale works", {
  gb.res3 <- ebnm_generalized_binary(x, s, scale = 0.01)
  g <- gb.res3$fitted_g
  expect_equal(g$sd[2] / g$mean[2], 0.01)
})

gb.res3 <- ebnm_generalized_binary(x, s, mode = 3)

test_that("Setting the mode works", {
  expect_identical(gb.res3$fitted_g$mean[2], 3)
})

test_that("Fixing g works", {
  gb.res4 <- ebnm_generalized_binary(x, s, g_init = gb.res3$fitted_g, fix_g = FALSE)
  expect_equal(gb.res4$fitted_g$mean[2], 2, tolerance = 0.1)
  gb.res5 <- ebnm_generalized_binary(x, s, g_init = gb.res3$fitted_g, fix_g = TRUE)
  expect_identical(gb.res5$fitted_g$mean[2], 3)
})

test_that("Additional parameters can be passed in", {
  gb.res4 <- ebnm_generalized_binary(x, s, mu_range = c(2.8, 3.2))
  expect_true(gb.res4$fitted_g$mean[2] >= 2.8)
  expect_true(gb.res4$fitted_g$mean[2] <= 3.2)
})
