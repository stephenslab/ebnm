context("point mass")

x <- 1:10
s <- 1

test_that("Basic functionality works", {
  pm.res <- ebnm(x, s, prior_family = "point_mass", mode = 0)
  pm.res2 <- ebnm_point_mass(x, s, mode = 0)
  pm.res$call <- pm.res2$call <- NULL
  expect_identical(pm.res, pm.res2)
  expect_identical(pm.res$posterior$mean, rep(0, length(x)))
  expect_identical(pm.res$posterior$sd, rep(0, length(x)))
})

test_that("Parameters that do nothing are ignored", {
  expect_warning(pm.res <- ebnm(x, s, prior_family = "point_mass", scale = 1))
})

test_that("Fixing g works", {
  g <- ashr::normalmix(1, 0, 0)
  pm.res <- ebnm_point_mass(x, s, g_init = g, fix_g = TRUE)
  expect_true(all(pm.res$posterior$mean == 0))
})

test_that("Setting the mode works", {
  pm.res <- ebnm_point_mass(x, s, mode = 5)
  expect_equal(pm.res$fitted_g$mean, 5)
})

test_that("Mode estimation works", {
  pm.res <- ebnm_point_mass(x, s, mode = "estimate")
  expect_equal(pm.res$fitted_g$mean, mean(x))
})
