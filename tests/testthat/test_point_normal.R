context("Point normal")

n <- 100
set.seed(1)
s <- rgamma(n, 1, 1)
x <- rnorm(n, 0, s + 1)

ebnm.res <- ebnm_point_normal(x, s, g = list(mu = 0), fix_mu = TRUE)

test_that("compute_summary_results gives same results as ashr", {
  pi0 <- ebnm.res$fitted_g$pi0
  a <- ebnm.res$fitted_g$a
  ash.res <- ashr::ash(x, s, outputlevel = 5, method = "shrink",
                      g = ashr::normalmix(c(pi0, 1 - pi0),
                                          c(0, 0),
                                          c(0, sqrt(1 / a))),
                      fixg = TRUE)

  expect_equal(ebnm.res$result$PosteriorMean,
               ash.res$flash_data$postmean, tol = 1e-6)
  expect_equal(ebnm.res$result$PosteriorMean2,
               ash.res$flash_data$postmean2, tol = 1e-6)
  expect_equal(ebnm.res$loglik,ash.res$flash_data$penloglik, tol = 1e-6)
})

test_that("setting norm parameter does not change results", {
  ebnm.res2 <- ebnm_point_normal(x, s)
  ebnm.res3 <- ebnm_point_normal(x, s, norm = 1)

  expect_equal(ebnm.res2$loglik, ebnm.res3$loglik, tol = 1e-4)
  expect_equal(ebnm.res2$result, ebnm.res3$result, tol = 1e-4)
  expect_equal(ebnm.res2$fitted_g$a, ebnm.res3$fitted_g$a, tol = 1e-4)
  expect_equal(ebnm.res2$fitted_g$pi0, ebnm.res3$fitted_g$pi0, tol = 1e-4)
  expect_equal(ebnm.res2$fitted_g$mu, ebnm.res3$fitted_g$mu, tol = 1e-4)
})

test_that("fixing g works as intended", {
  g <- list(pi0 = 0, a = 0.5, mu = 0)
  ebnm.res4 <- ebnm_point_normal(x, s, g, fixg = TRUE)

  expect_identical(ebnm.res4$fitted_g, g)
})

test_that("fixing pi0 works as intended", {
  g <- list(pi0 = 0.2, a = 0.5, mu = 0)
  ebnm.res5 <- ebnm_point_normal(x, s, g, fix_pi0 = TRUE, fix_mu = FALSE)

  expect_identical(ebnm.res5$fitted_g$pi0, g$pi0)
  expect_false(ebnm.res5$fitted_g$a == g$a)
  expect_false(ebnm.res5$fitted_g$mu == g$mu)
})

test_that("fixing mu works as intended", {
  g <- list(pi0 = 0.2, a = 0.5, mu = 0)
  ebnm.res6 <- ebnm_point_normal(x, s, g, fix_mu = TRUE)

  expect_identical(ebnm.res6$fitted_g$mu, g$mu)
  expect_false(ebnm.res6$fitted_g$a == g$a)
  expect_false(ebnm.res6$fitted_g$pi0 == g$pi0)
})

test_that("fixing pi0 and mu together works as intended", {
  g <- list(pi0 = 0.2, a = 0.5, mu = 0)
  ebnm.res7 <- ebnm_point_normal(x, s, g, fix_pi0 = TRUE, fix_mu = TRUE)

  expect_identical(ebnm.res7$fitted_g$pi0, g$pi0)
  expect_identical(ebnm.res7$fitted_g$mu, g$mu)
  expect_false(ebnm.res7$fitted_g$a == g$a)
})

test_that("infinite and zero SEs give expected results", {
  x <- c(rep(0, 5), rep(1, 5))
  s <- rep(1, 10)
  s[6] <- 0
  s[10] <- Inf
  # first, fix mu = 0
  ebnm.res8 <- ebnm_point_normal(x, s, g = list(mu = 0), fix_mu = TRUE)

  expect_equal(ebnm.res8$result$PosteriorMean[6], x[6])
  expect_equal(ebnm.res8$result$PosteriorMean2[6], x[6]^2)
  expect_equal(ebnm.res8$result$PosteriorMean[10], 0)
  expect_equal(ebnm.res8$result$PosteriorMean2[10],
               1 / ebnm.res8$fitted_g$a * (1 - ebnm.res8$fitted_g$pi0))

  # now, don't fix mu
  ebnm.res9 = ebnm_point_normal(x, s)
  expect_equal(ebnm.res9$result$PosteriorMean[6], x[6])
  expect_equal(ebnm.res9$result$PosteriorMean2[6], x[6]^2)
  expect_equal(ebnm.res9$result$PosteriorMean[10], ebnm.res9$fitted_g$mu)
  expect_equal(ebnm.res9$result$PosteriorMean2[10],
               (1 / ebnm.res9$fitted_g$a * (1 - ebnm.res9$fitted_g$pi0)) + ebnm.res9$fitted_g$mu^2)
})

test_that("removing null component gives reasonable estimates for mu and a", {
  n = 1000
  mu = 3
  a = 1/25
  theta = rnorm(n, mu, 1 / sqrt(a))
  s = rgamma(n, 1, 1)
  x = rnorm(n, theta, s)
  ebnm.res10 = ebnm_point_normal(x, s, g = list(pi0 = 0), fix_pi0 = TRUE,
                                 fix_mu = FALSE)
  expect_equal(ebnm.res10$fitted_g$mu, mu, tol = .25)
  expect_equal(ebnm.res10$fitted_g$a, a, tol = .02)
})
