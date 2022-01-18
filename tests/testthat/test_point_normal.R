context("Point normal")

n <- 1000
set.seed(1)
s <- rnorm(n, 1, 0.1)
x <- c(rnorm(n / 2, 0, 10 + s), rnorm(n / 2, 0, s))

true_pi0 <- 0.5
true_mean <- 0
true_sd <- 10
true_g <- ashr::normalmix(pi = c(true_pi0, 1 - true_pi0),
                         mean = rep(true_mean, 2),
                         sd = c(0, true_sd))

test_that("Basic functionality works", {
  pn.res <- ebnm(x, s, prior_family = "point_normal")
  pn.res2 <- ebnm_point_normal(x, s)
  pn.res$call <- pn.res2$call <- NULL
  expect_identical(pn.res, pn.res2)
  expect_equal(pn.res[[g_ret_str()]], true_g, tolerance = 0.1)
})

test_that("Mode estimation works", {
  pn.res <- ebnm_point_normal(x, s, mode = "est")
  expect_equal(pn.res[[g_ret_str()]], true_g, tolerance = 0.1)
  expect_false(identical(pn.res[[g_ret_str()]]$mean[1], true_mean))
})

test_that("Fixing the sd works", {
  pn.res <- ebnm_point_normal(x, s, scale = true_sd)
  expect_equal(pn.res[[g_ret_str()]], true_g, tolerance = 0.1)
  expect_equal(pn.res[[g_ret_str()]]$sd[2], true_sd, tolerance = 1e-8)
})

test_that("Fixing g works", {
  pn.res <- ebnm_point_normal(x, s, g_init = true_g, fix_g = TRUE)
  expect_equal(pn.res[[g_ret_str()]], true_g)
})

test_that("Output parameter works", {
  pn.res <- ebnm_point_normal(x, s, output = samp_arg_str())
  pn.res$call <- NULL
  expect_identical(names(pn.res), samp_ret_str())
})

test_that("compute_summary_results gives same results as ashr", {
  pn.res <- ebnm_point_normal(x, s, output = output_all())
  ash.res <- ebnm_ash(x, s, g_init = pn.res[[g_ret_str()]], fix_g = TRUE,
                      output = output_all(), method = "shrink")

  expect_equal(pn.res[[df_ret_str()]][[pm_ret_str()]],
               ash.res[[df_ret_str()]][[pm_ret_str()]], tol = 1e-6)
  expect_equal(pn.res[[df_ret_str()]][[psd_ret_str()]],
               ash.res[[df_ret_str()]][[psd_ret_str()]], tol = 1e-6)
  expect_equal(pn.res[[df_ret_str()]]$lfsr,
               ash.res[[df_ret_str()]]$lfsr, tol = 1e-6)
  expect_equal(pn.res[[llik_ret_str()]],
               ash.res[[llik_ret_str()]], tol = 1e-6)
})

# test_that("Infinite and zero SEs give expected results", {
#   x <- c(rep(0, 5), rep(10, 5))
#   s <- rep(1, 10)
#   # s[6] <- 0
#   s[10] <- Inf
#
#   pn.res <- ebnm_point_normal(x, s, output = output_all())
#
#   # expect_equal(pn.res[[df_ret_str()]][[pm_ret_str()]][6], x[6])
#   # expect_equal(pn.res[[df_ret_str()]][[pm2_ret_str()]][6], x[6]^2)
#   expect_equal(pn.res[[df_ret_str()]][[pm_ret_str()]][10], 0)
#   expect_equal(pn.res[[df_ret_str()]][[pm2_ret_str()]][10],
#                pn.res[[g_ret_str()]]$pi[2] * pn.res[[g_ret_str()]]$sd[2]^2)
#
#   # Zero SEs should throw an error if mu is not fixed.
#   # expect_error(ebnm_point_normal(x, s, mode = "est"))
# })

test_that("Null case estimates pi0 = 1", {
  x <- rnorm(n, s = 0.5)
  pn.res <- ebnm_point_normal(x, s = 1)
  expect_equal(pn.res[[g_ret_str()]]$pi[1], 1)
})

test_that("Very large observations give reasonable results", {
  scl <- 1e8
  pn.res <- ebnm_point_normal(x, s, mode = "estimate")
  pn.res.lg <- ebnm_point_normal(scl*x, scl*s, mode = "estimate")

  expect_equal(pn.res[[g_ret_str()]]$pi[1], pn.res.lg[[g_ret_str()]]$pi[1])
  expect_equal(scl * pn.res[[g_ret_str()]]$sd[2], pn.res.lg[[g_ret_str()]]$sd[2])
  expect_equal(scl * pn.res[[g_ret_str()]]$mean[1], pn.res.lg[[g_ret_str()]]$mean[1])
})

test_that("Very small observations give reasonable results", {
  scl <- 1e-8
  pn.res <- ebnm_point_normal(x, s, mode = "estimate")
  pn.res.sm <- ebnm_point_normal(scl*x, scl*s, mode = "estimate")

  expect_equal(pn.res[[g_ret_str()]]$pi[1], pn.res.sm[[g_ret_str()]]$pi[1])
  expect_equal(scl * pn.res[[g_ret_str()]]$sd[2], pn.res.sm[[g_ret_str()]]$sd[2])
  expect_equal(scl * pn.res[[g_ret_str()]]$mean[1], pn.res.sm[[g_ret_str()]]$mean[1])
})

test_that("g_init with pi0 = 0 or pi0 = 1 isn't a dealbreaker", {
  pn.res <- ebnm_point_normal(x, s)

  bad_g <- normalmix(c(1, 0), c(0, 0), c(0, true_sd))
  pn.res2 <- ebnm_point_normal(x, s, g_init = bad_g, fix_g = FALSE)
  expect_equal(pn.res[[llik_ret_str()]], pn.res2[[llik_ret_str()]])

  bad_g <- normalmix(c(0, 1), c(0, 0), c(0, true_sd))
  pn.res3 <- ebnm_point_normal(x, s, g_init = bad_g, fix_g = FALSE)
  expect_equal(pn.res[[llik_ret_str()]], pn.res3[[llik_ret_str()]])

})
