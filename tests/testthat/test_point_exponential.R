context("Point Exponential")

n <- 1000
set.seed(1)
s <- rnorm(n, 1, 0.1)
x <- c(rexp(n / 2, rate = 0.1), rep(0, n / 2)) + rnorm(n, sd = s)

true_pi0 <- 0.5
true_scale <- 10
true_mode <- 0

true_g <- exponentialmix(pi = c(true_pi0, 1 - true_pi0),
                         mode = rep(true_mode, 2),
                         scale = c(0, true_scale))

pe.res <- ebnm(x, s, prior_family = "point_exponential")

test_that("Basic functionality works", {
  pe.res2 <- ebnm_point_exponential(x, s)
  expect_identical(pe.res, pe.res2)
  expect_equal(pe.res[[g_ret_str()]], true_g, tolerance = 0.1)
})

test_that("Mode estimation works", {
  pe.res2 <- ebnm_point_exponential(x, s, mode = "est")
  expect_equal(pe.res2[[g_ret_str()]], true_g, tolerance = 0.5)
  expect_false(identical(pe.res2[[g_ret_str()]]$mean[1], true_mean))
})

test_that("Fixing the scale works", {
  pe.res2 <- ebnm_point_exponential(x, s, scale = true_scale)
  expect_equal(pe.res2[[g_ret_str()]], true_g, tolerance = 0.1)
  expect_equal(pe.res2[[g_ret_str()]]$scale[2], true_scale)
})

test_that("Fixing g works", {
  pe.res2 <- ebnm_point_exponential(x, s, g_init = pe.res[[g_ret_str()]], fix_g = TRUE)
  expect_identical(pe.res[[g_ret_str()]], pe.res2[[g_ret_str()]])
  expect_equal(pe.res[[llik_ret_str()]], pe.res2[[llik_ret_str()]])
})

test_that("Initializing g works", {
  pe.res2 <- ebnm_point_exponential(x, s, g_init = true_g)
  expect_equal(pe.res[[llik_ret_str()]], pe.res2[[llik_ret_str()]])
})

test_that("Output parameter works", {
  pe.res <- ebnm_point_exponential(x, s, output = c("fitted_g"))
  expect_identical(names(pe.res), "fitted_g")
})
