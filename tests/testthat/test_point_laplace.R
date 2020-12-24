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

pl.res <- ebnm(x, s, prior_family = "point_laplace")

test_that("Basic functionality works", {
  pl.res2 <- ebnm_point_laplace(x, s)
  expect_identical(pl.res, pl.res2)
  expect_equal(pl.res[[g_ret_str()]], true_g, tolerance = 0.1)
})

test_that("Mode estimation works", {
  pl.res2 <- ebnm_point_laplace(x, s, mode = "est")
  expect_equal(pl.res2[[g_ret_str()]], true_g, tolerance = 0.5)
  expect_false(identical(pl.res2[[g_ret_str()]]$mean[1], true_mean))
})

test_that("Fixing the scale works", {
  pl.res2 <- ebnm_point_laplace(x, s, scale = true_scale)
  expect_equal(pl.res2[[g_ret_str()]], true_g, tolerance = 0.1)
  expect_identical(pl.res2[[g_ret_str()]]$scale[2], true_scale)
})

test_that("Fixing g works", {
  pl.res2 <- ebnm_point_laplace(x, s, g_init = pl.res[[g_ret_str()]], fix_g = TRUE)
  expect_identical(pl.res[[g_ret_str()]], pl.res2[[g_ret_str()]])
  expect_equal(pl.res[[llik_ret_str()]], pl.res2[[llik_ret_str()]])
})

test_that("Initializing g works", {
  pl.res2 <- ebnm_point_laplace(x, s, g_init = true_g)
  expect_equal(pl.res[[llik_ret_str()]], pl.res2[[llik_ret_str()]])
})

test_that("Output parameter works", {
  pl.res <- ebnm_point_laplace(x, s, output = c("fitted_g"))
  expect_identical(names(pl.res), "fitted_g")
})

test_that("Infinite and zero SEs give expected results", {
  x <- c(rep(0, 5), rep(5, 5))
  s <- rep(1, 10)
  # s[6] <- 0
  s[10] <- Inf

  pn.res <- ebnm_point_laplace(x, s, output = output_all())

  # expect_equal(pn.res[[df_ret_str()]][[pm_ret_str()]][6], x[6])
  # expect_equal(pn.res[[df_ret_str()]][[pm_ret_str()]]2[6], x[6]^2)
  expect_equal(pn.res[[df_ret_str()]][[pm_ret_str()]][10], 0)
  expect_equal(pn.res[[df_ret_str()]][[pm2_ret_str()]][10],
               2 * pn.res[[g_ret_str()]]$scale[2]^2 * pn.res[[g_ret_str()]]$pi[2])
  expect_equal(pn.res[[df_ret_str()]][[lfsr_ret_str()]][10],
               pn.res[[g_ret_str()]]$pi[1] + pn.res[[g_ret_str()]]$pi[2] / 2)
})
