context("Point Exponential")

n <- 1000
set.seed(1)
s <- rnorm(n, 1, 0.1)
x <- c(rexp(n / 2, rate = 0.1), rep(0, n / 2)) + rnorm(n, sd = s)

true_pi0 <- 0.5
true_scale <- 10
true_mode <- 0

true_g <- gammamix(pi = c(true_pi0, 1 - true_pi0),
                   shape = c(1, 1),
                   scale = c(0, true_scale),
                   shift = rep(true_mode, 2))

pe.res <- ebnm(x, s, prior_family = "point_exponential")

test_that("Basic functionality works", {
  pe.res2 <- ebnm_point_exponential(x, s)
  pe.res$call <- pe.res2$call <- NULL
  expect_identical(pe.res, pe.res2)
  expect_equal(pe.res[[g_ret_str()]], true_g, tolerance = 0.1)
})

test_that("Mode estimation works", {
  pe.res2 <- ebnm_point_exponential(x, s, mode = "est")
  expect_equal(pe.res2[[g_ret_str()]], true_g, tolerance = 0.5)
  expect_false(identical(pe.res2[[g_ret_str()]]$mean[1], true_mode))
})

test_that("Fixing the scale works", {
  pe.res2 <- ebnm_point_exponential(x, s, scale = true_scale)
  expect_equal(pe.res2[[g_ret_str()]], true_g, tolerance = 0.1)
  expect_equal(pe.res2[[g_ret_str()]]$scale[2], true_scale)
})

test_that("Fixing g works", {
  pe.res2 <- ebnm_point_exponential(x, s, g_init = pe.res[[g_ret_str()]], fix_g = TRUE)
  expect_identical(pe.res[[g_ret_str()]], pe.res2[[g_ret_str()]])
  expect_equal(as.numeric(pe.res[[llik_ret_str()]]), as.numeric(pe.res2[[llik_ret_str()]]))
})

test_that("Initializing g works", {
  pe.res2 <- ebnm_point_exponential(x, s, g_init = true_g)
  expect_equal(pe.res[[llik_ret_str()]], pe.res2[[llik_ret_str()]])
})

test_that("Output parameter works", {
  pe.res <- ebnm_point_exponential(x, s, output = c("fitted_g"))
  pe.res$call <- NULL
  expect_identical(names(pe.res), "fitted_g")
})

test_that("Can fix g with one component", {
  g_init <- gammamix(pi = 1,
                     shape = 1,
                     scale = true_scale,
                     shift = true_mode)
  pe.res <- ebnm_point_exponential(x, s, g_init = g_init, fix_g = TRUE)

  g_init2 <- gammamix(pi = c(0, 1),
                      shape = c(1, 1),
                      scale = c(0, true_scale),
                      shift = rep(true_mode, 2))
  pe.res2 <- ebnm_point_exponential(x, s, g_init = g_init2, fix_g = TRUE)

  expect_equal(pe.res$log_likelihood, pe.res2$log_likelihood)
})

test_that("Null case estimates pi0 = 1", {
  x <- rnorm(n)
  pe.res <- ebnm_point_exponential(x, s = 1, scale = 1)
  expect_equal(pe.res[[g_ret_str()]]$pi[1], 1)
})

# test_that("predict method works as expected", {
#   pe.res <- ebnm_point_exponential(x, s)
#   pe.res2 <- predict(pe.res, list(x = 1:10, s = 1))
#   expect_equal(pe.res$fitted_g, pe.res2$fitted_g)
# })
