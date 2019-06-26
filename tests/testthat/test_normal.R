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
  expect_identical(norm.res, norm.res2)
  expect_equal(norm.res$fitted_g, true_g, tolerance = 0.2)
})

test_that("Mode estimation works", {
  norm.res <- ebnm_normal(x, s, mode = "est")
  expect_equal(norm.res$fitted_g, true_g, tolerance = 0.2)
  expect_false(identical(norm.res$fitted_g$mean, true_mean))
})

test_that("Fixing the sd works", {
  norm.res <- ebnm_normal(x, s, scale = true_sd)
  expect_equal(norm.res$fitted_g, true_g, tolerance = 0.1)
  expect_identical(norm.res$fitted_g$sd, true_sd)
})

test_that("Fixing g works", {
  norm.res <- ebnm_normal(x, s, g_init = true_g, fix_g = TRUE)
  expect_identical(norm.res$fitted_g, true_g)
})

test_that("Output parameter works", {
  norm.res <- ebnm_normal(x, s, output = "post_sampler")
  expect_identical(names(norm.res), "post_sampler")
})

test_that("compute_summary_results gives same results as ashr", {
  output <- c("result", "lfsr", "fitted_g")
  norm.res <- ebnm_normal(x, s, output = output)
  ash.res <- ebnm_ash(x, s, g_init = norm.res$fitted_g, fix_g = TRUE,
                      output = output, method = "shrink")

  expect_equal(norm.res$result$posterior_mean, ash.res$result$posterior_mean,
               tol = 1e-6)
  expect_equal(norm.res$result$posterior_mean2, ash.res$result$posterior_mean2,
               tol = 1e-6)
  expect_equal(norm.res$result$lfsr, ash.res$result$lfsr,
               tol = 1e-6)
  expect_equal(norm.res$loglik, ash.res$loglik,
               tol = 1e-6)
})

test_that("Infinite and zero SEs give expected results", {
  x <- c(rep(0, 5), rep(1, 5))
  s <- rep(1, 10)
  # s[6] <- 0
  s[10] <- Inf

  norm.res <- ebnm_normal(x, s)

  # expect_equal(norm.res$result$posterior_mean[6], x[6])
  # expect_equal(norm.res$result$posterior_mean2[6], x[6]^2)
  expect_equal(norm.res$result$posterior_mean[10], 0)
  expect_equal(norm.res$result$posterior_mean2[10],
               norm.res$fitted_g$sd^2)
})
