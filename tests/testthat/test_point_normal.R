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
  pn.res <- ebnm(x, s, prior_type = "point_normal")
  pn.res2 <- ebnm_point_normal(x, s)
  expect_identical(pn.res, pn.res2)
  expect_equal(pn.res$fitted_g, true_g, tolerance = 0.1)
})

test_that("Mode estimation works", {
  pn.res <- ebnm_point_normal(x, s, mode = "est")
  expect_equal(pn.res$fitted_g, true_g, tolerance = 0.1)
  expect_false(identical(pn.res$fitted_g$mean[1], true_mean))
})

test_that("Fixing the sd works", {
  pn.res <- ebnm_point_normal(x, s, sd = true_sd)
  expect_equal(pn.res$fitted_g, true_g, tolerance = 0.1)
  expect_identical(pn.res$fitted_g$sd[2], true_sd)
})

test_that("Fixing g works", {
  pn.res <- ebnm_point_normal(x, s, g_init = true_g, fix_g = TRUE)
  expect_identical(pn.res$fitted_g, true_g)
})

test_that("Output parameter works", {
  pn.res <- ebnm_point_normal(x, s, output = "post_sampler")
  expect_identical(names(pn.res), "post_sampler")
})

test_that("compute_summary_results gives same results as ashr", {
  output <- c("result", "lfsr", "fitted_g")
  pn.res <- ebnm_point_normal(x, s, output = output)
  ash.res <- ebnm_ash(x, s, g_init = pn.res$fitted_g, fix_g = TRUE,
                      output = output, method = "shrink")

  expect_equal(pn.res$result$PosteriorMean, ash.res$result$PosteriorMean,
               tol = 1e-6)
  expect_equal(pn.res$result$PosteriorMean2, ash.res$result$PosteriorMean2,
               tol = 1e-6)
  expect_equal(pn.res$result$lfsr, ash.res$result$lfsr,
               tol = 1e-6)
  expect_equal(pn.res$loglik, ash.res$loglik,
               tol = 1e-6)
})

test_that("Infinite and zero SEs give expected results", {
  x <- c(rep(0, 5), rep(1, 5))
  s <- rep(1, 10)
  s[6] <- 0
  s[10] <- Inf

  pn.res <- ebnm_point_normal(x, s)

  expect_equal(pn.res$result$PosteriorMean[6], x[6])
  expect_equal(pn.res$result$PosteriorMean2[6], x[6]^2)
  expect_equal(pn.res$result$PosteriorMean[10], 0)
  expect_equal(pn.res$result$PosteriorMean2[10],
               pn.res$fitted_g$sd[2]^2 * pn.res$fitted_g$pi[2])

  # Should get an error if mu is not fixed.
  expect_error(ebnm_point_normal(x, s, mode = "est"))
})
