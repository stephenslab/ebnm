context("Normal")

n <- 1000
set.seed(1)
s <- rnorm(n, 1, 0.1)
x <- rnorm(n, 0, 10 + s)

true_mean <- 0
true_sd <- 10
true_g <- ashr::normalmix(pi = 1, mean = true_mean, sd = true_sd)

test_that("Basic functionality works", {
  norm.res <- ebnm(x, s, prior_family = "normal")
  norm.res2 <- ebnm_normal(x, s)
  norm.res$call <- norm.res2$call <- NULL
  expect_identical(norm.res, norm.res2)
  expect_equal(norm.res[[g_ret_str()]], true_g, tolerance = 0.2)
})

test_that("Mode estimation works", {
  norm.res <- ebnm_normal(x, s, mode = "est")
  expect_equal(norm.res[[g_ret_str()]], true_g, tolerance = 0.2)
  expect_false(identical(norm.res[[g_ret_str()]]$mean, true_mean))
})

test_that("Fixing the sd works", {
  norm.res <- ebnm_normal(x, s, scale = true_sd)
  expect_equal(norm.res[[g_ret_str()]], true_g, tolerance = 0.1)
  expect_equal(norm.res[[g_ret_str()]]$sd, true_sd)
})

test_that("Fixing g works", {
  norm.res <- ebnm_normal(x, s, g_init = true_g, fix_g = TRUE)
  expect_equal(norm.res[[g_ret_str()]], true_g)
})

test_that("Output parameter works", {
  norm.res <- ebnm_normal(x, s, output = samp_arg_str())
  norm.res$call <- NULL
  expect_identical(names(norm.res), samp_ret_str())
})

test_that("compute_summary_results gives same results as ashr", {
  norm.res <- ebnm_normal(x, s, output = ebnm_output_all())
  ash.res <- ebnm_ash(x, s, g_init = norm.res[[g_ret_str()]], fix_g = TRUE,
                      output = ebnm_output_all(), method = "shrink")

  expect_equal(norm.res[[df_ret_str()]][[pm_ret_str()]],
               ash.res[[df_ret_str()]][[pm_ret_str()]], tol = 1e-6)
  expect_equal(norm.res[[df_ret_str()]][[psd_ret_str()]],
               ash.res[[df_ret_str()]][[psd_ret_str()]], tol = 1e-6)
  expect_equal(norm.res[[df_ret_str()]][[lfsr_ret_str()]],
               ash.res[[df_ret_str()]][[lfsr_ret_str()]], tol = 1e-6)
  expect_equal(as.numeric(norm.res[[llik_ret_str()]]),
               as.numeric(ash.res[[llik_ret_str()]]), tol = 1e-6)
})

# test_that("Infinite and zero SEs give expected results", {
#   x <- c(rep(0, 5), rep(10, 5))
#   s <- rep(1, 10)
#   # s[6] <- 0
#   s[10] <- Inf
#
#   norm.res <- ebnm_normal(x, s)
#
#   # expect_equal(norm.res[[df_ret_str()]][[pm_ret_str()]][6], x[6])
#   # expect_equal(norm.res[[df_ret_str()]][[psd_ret_str()]][6], x[6]^2)
#   expect_equal(norm.res[[df_ret_str()]][[pm_ret_str()]][10], 0)
#   expect_equal(norm.res[[df_ret_str()]][[psd_ret_str()]][10], norm.res[[g_ret_str()]]$sd)
# })
