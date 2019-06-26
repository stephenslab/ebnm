context("ebnm_normal_scale_mixture")

n <- 1000
set.seed(1)
s <- rnorm(n, 1, 0.1)
x <- c(rnorm(n / 2, 0, 10 + s), rnorm(n / 2, 0, s))

test_that("Basic ebnm_normal_scale_mixture functionality works", {
  ebnm.res <- ebnm(x, s, prior_family = "normal_scale_mixture")
  ebnm.res2 <- ebnm_normal_scale_mixture(x, s)
  expect_identical(ebnm.res, ebnm.res2)
})

test_that("ebnm_normal_scale_mixture returns the same results as ashr", {
  ebnm.res <- ebnm_normal_scale_mixture(x, s, output = output_all())
  ash.res  <- ebnm_ash(x, s, mixcompdist = "normal", prior = "uniform",
                       output = output_all())

  # The ash grid is the same:
  expect_equal(ebnm.res$fitted_g$sd, ash.res$fitted_g$sd)

  # The estimated mixture probabilities are the same:
  expect_equal(ebnm.res$fitted_g$pi, ash.res$fitted_g$pi)

  # The posterior quantities are the same:
  expect_equal(ebnm.res$result, ash.res$result)

  # Likelihoods are the same:
  expect_equal(ebnm.res$loglik, ash.res$loglik)

  # Now repeat with scalar s.
  s <- 1
  ebnm.res <- ebnm_normal_scale_mixture(x, s, output = output_all())
  ash.res  <- ebnm_ash(x, s, mixcompdist = "normal", prior = "uniform",
                       output = output_all())
  expect_equal(ebnm.res$fitted_g$sd, ash.res$fitted_g$sd)
  expect_equal(ebnm.res$fitted_g$pi, ash.res$fitted_g$pi, tol = 1e-6)
  expect_equal(ebnm.res$result, ash.res$result)
  expect_equal(ebnm.res$loglik, ash.res$loglik)
})

test_that("Mode estimation works", {
  x <- x + 3
  ebnm.res <- ebnm_normal_scale_mixture(x, s, mode = "est")
  expect_equal(ebnm.res$fitted_g$mean[1], 3, tol = .01)
})
