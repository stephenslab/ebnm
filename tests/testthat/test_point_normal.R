context("Point normal")

n <- 100
set.seed(1)
s <- rgamma(n, 1, 1)
x <- rnorm(n, 0, s + 1)

test_that("compute_summary_results gives same results as ashr", {
  ebnm.res <- ebnm_point_normal(x, s, g = list(mu = 0), fix_mu = TRUE,
                                output = c("result", "lfsr", "fitted_g", "loglik"))
  pi0 <- ebnm.res$fitted_g$pi0
  a <- ebnm.res$fitted_g$a
  ash.res <- ashr::ash(x, s, method = "shrink",
                      g = ashr::normalmix(c(pi0, 1 - pi0),
                                          c(0, 0),
                                          c(0, sqrt(1 / a))),
                      fixg = TRUE)

  expect_equal(ebnm.res$result$PosteriorMean,
               ash.res$result$PosteriorMean,
               tol = 1e-6)
  expect_equal(ebnm.res$result$PosteriorMean2,
               ash.res$result$PosteriorMean^2 + ash.res$result$PosteriorSD^2,
               tol = 1e-6)
  expect_equal(ebnm.res$result$lfsr,
               ash.res$result$lfsr,
               tol = 1e-6)
  expect_equal(ebnm.res$loglik,
               ash.res$loglik,
               tol = 1e-6)
})

test_that("fixing g works as intended", {
  g <- list(pi0 = 0, a = 0.5, mu = 0)
  ebnm.res <- ebnm_point_normal(x, s, g, fixg = TRUE)

  expect_identical(ebnm.res$fitted_g, g)
})

g <- list(pi0 = 0.2, a = 0.5, mu = 0)

test_that("fixing pi0 works as intended", {
  ebnm.res <- ebnm_point_normal(x, s, g, fix_pi0 = TRUE, fix_mu = FALSE)

  expect_identical(ebnm.res$fitted_g$pi0, g$pi0)
  expect_false(ebnm.res$fitted_g$a == g$a)
  expect_false(ebnm.res$fitted_g$mu == g$mu)
})

test_that("fixing a works as intended", {
  ebnm.res <- ebnm_point_normal(x, s, g, fix_a = TRUE, fix_mu = FALSE)

  expect_identical(ebnm.res$fitted_g$a, g$a)
  expect_false(ebnm.res$fitted_g$pi0 == g$pi0)
  expect_false(ebnm.res$fitted_g$mu == g$mu)
})

test_that("fixing mu works as intended", {
  ebnm.res <- ebnm_point_normal(x, s, g, fix_mu = TRUE)

  expect_identical(ebnm.res$fitted_g$mu, g$mu)
  expect_false(ebnm.res$fitted_g$a == g$a)
  expect_false(ebnm.res$fitted_g$pi0 == g$pi0)
})

test_that("fixing pi0 and mu together works as intended", {
  ebnm.res <- ebnm_point_normal(x, s, g, fix_pi0 = TRUE, fix_mu = TRUE)

  expect_identical(ebnm.res$fitted_g$pi0, g$pi0)
  expect_identical(ebnm.res$fitted_g$mu, g$mu)
  expect_false(ebnm.res$fitted_g$a == g$a)
})

test_that("infinite and zero SEs give expected results", {
  x <- c(rep(0, 5), rep(1, 5))
  s <- rep(1, 10)
  s[6] <- 0
  s[10] <- Inf

  # first, fix mu = 0
  ebnm.res <- ebnm_point_normal(x, s, g = list(mu = 0), fix_mu = TRUE)

  expect_equal(ebnm.res$result$PosteriorMean[6], x[6])
  expect_equal(ebnm.res$result$PosteriorMean2[6], x[6]^2)
  expect_equal(ebnm.res$result$PosteriorMean[10], 0)
  expect_equal(ebnm.res$result$PosteriorMean2[10],
               1 / ebnm.res$fitted_g$a * (1 - ebnm.res$fitted_g$pi0))

  # now, don't fix mu
  ebnm.res = ebnm_point_normal(x, s)
  expect_equal(ebnm.res$result$PosteriorMean[6], x[6])
  expect_equal(ebnm.res$result$PosteriorMean2[6], x[6]^2)
  expect_equal(ebnm.res$result$PosteriorMean[10], ebnm.res$fitted_g$mu)
  expect_equal(ebnm.res$result$PosteriorMean2[10],
               (1 / ebnm.res$fitted_g$a * (1 - ebnm.res$fitted_g$pi0)) + ebnm.res$fitted_g$mu^2)
})

test_that("removing null component gives reasonable estimates for mu and a", {
  n = 1000
  mu = 3
  a = 1/25
  theta = rnorm(n, mu, 1 / sqrt(a))
  s = rgamma(n, 1, 1)
  x = rnorm(n, theta, s)
  ebnm.res = ebnm_point_normal(x, s, g = list(pi0 = 0), fix_pi0 = TRUE,
                                 fix_mu = FALSE)
  expect_equal(ebnm.res$fitted_g$mu, mu, tol = .25)
  expect_equal(ebnm.res$fitted_g$a, a, tol = .02)
})
