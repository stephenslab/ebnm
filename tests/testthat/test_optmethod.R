# test_n <- function(n, seed = 666, mb_times = 100L) {
#   set.seed(seed)
#   # true g is 0.5 delta_0 + 0.5 N(0, 2^2)
#   theta <- c(rep(0, n), rnorm(n, sd = 2))
#   s <- sqrt(rexp(2 * n))
#   x <- theta + rnorm(2 * n, sd = s)
#
#   test_res <- microbenchmark::microbenchmark(
#     ebnm_point_normal(x, s, optmethod = "nlm"),
#     ebnm_point_normal(x, s, optmethod = "lbfgsb"),
#     ebnm_point_normal(x, s, optmethod = "trust"),
#     ebnm_point_normal(x, s, optmethod = "nograd_nlm"),
#     ebnm_point_normal(x, s, optmethod = "nograd_lbfgsb"),
#     ebnm_point_normal(x, s, optmethod = "nohess_nlm"),
#     times = mb_times)
#
#   return(test_res)
# }

context("Optimization Methods")

n <- 500
set.seed(666)
theta <- c(rep(0, 3 * n), rnorm(n, sd = 2))
s <- sqrt(rexp(4 * n))
x <- theta + rnorm(4 * n, sd = s)
ebnm.res <- ebnm_point_normal(x, s, optmethod = "nlm")

test_that("Different optimization methods work", {
  ebnm.res2 <- ebnm_point_normal(x, s, optmethod = "lbfgsb")
  expect_equal(ebnm.res$fitted_g$pi, ebnm.res2$fitted_g$pi, tol = 1e-4)

  ebnm.res2 <- ebnm_point_normal(x, s, optmethod = "trust")
  expect_equal(ebnm.res$fitted_g$pi, ebnm.res2$fitted_g$pi, tol = 1e-4)

  ebnm.res2 <- ebnm_point_normal(x, s, optmethod = "nograd_nlm")
  expect_equal(ebnm.res$fitted_g$pi, ebnm.res2$fitted_g$pi, tol = 1e-4)

  ebnm.res2 <- ebnm_point_normal(x, s, optmethod = "nograd_lbfgsb")
  expect_equal(ebnm.res$fitted_g$pi, ebnm.res2$fitted_g$pi, tol = 1e-4)

  ebnm.res2 <- ebnm_point_normal(x, s, optmethod = "nohess_nlm")
  expect_equal(ebnm.res$fitted_g$pi, ebnm.res2$fitted_g$pi, tol = 1e-4)

  ebnm.res3 <- ebnm_point_normal(x, s, scale = 2)
  ebnm.res4 <- ebnm_point_normal(x, s, scale = 2, optmethod = "optimize",
                                 control = list(interval = c(0, 1)))
  expect_equal(ebnm.res$fitted_g$pi, ebnm.res2$fitted_g$pi, tol = 1e-4)
})
