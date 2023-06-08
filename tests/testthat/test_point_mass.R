context("point mass")

x <- 1:10
s <- 1

test_that("Basic functionality works", {
  pm.res <- ebnm(x, s, prior_family = "point_mass")
  pm.res2 <- ebnm_point_mass(x, s)
  pm.res$call <- pm.res2$call <- NULL
  expect_identical(pm.res, pm.res2)
  expect_identical(pm.res$posterior$mean, rep(0, length(x)))
  expect_identical(pm.res$posterior$sd, rep(0, length(x)))
})

test_that("Parameters that do nothing are ignored", {
  expect_warning(pm.res <- ebnm(x, s, prior_family = "point_mass", scale = 1))
  expect_warning(pm.res <- ebnm(x, s, prior_family = "point_mass",
                                g_init = ashr::normalmix(1, 1, 0)))
})

test_that("Setting the mode works", {
  pm.res <- ebnm_point_mass(x, s, mode = 5)
  expect_equal(pm.res$fitted_g$mean, 5)
})

test_that("Mode estimation works", {
  pm.res <- ebnm_point_mass(x, s, mode = "estimate")
  expect_equal(pm.res$fitted_g$mean, mean(x))
})