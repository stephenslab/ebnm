context("ebnm_normal_scale_mixture")

n <- 1000
set.seed(1)
s <- rnorm(n, 1, 0.1)
x <- c(rnorm(n / 2, 0, 10 + s), rnorm(n / 2, 0, s))

test_that("Basic ebnm_normal_scale_mixture functionality works", {
  ash.res <- ebnm(x, s, prior_type = "normal_scale_mixture")
  ash.res2 <- ebnm_normal_scale_mixture(x, s)
  expect_identical(ash.res, ash.res2)
})

test_that("ebnm_normal_scale_mixture returns the same results as ashr", {
  ash.res  <- ebnm_normal_scale_mixture(x, s, output = output_all())
  ash.res2 <- ebnm_ash(x, s, mixcompdist = "normal", prior = "uniform",
                       output = output_all())

  # The ash grid is the same:
  expect_equal(ash.res$fitted_g$sd, ash.res2$fitted_g$sd)

  # The estimated mixture probabilities are the same:
  expect_equal(ash.res$fitted_g$pi, ash.res2$fitted_g$pi)

  # The posterior quantities are the same:
  expect_equal(ash.res$result, ash.res2$result)

  # Likelihoods are the same:
  expect_equal(ash.res$loglik, ash.res2$loglik)

  # Now repeat with scalar s.
  s <- 1
  ash.res  <- ebnm_normal_scale_mixture(x, s, output = output_all())
  ash.res2 <- ebnm_ash(x, s, mixcompdist = "normal", prior = "uniform",
                       output = output_all())
  expect_equal(ash.res$fitted_g$sd, ash.res2$fitted_g$sd)
  expect_equal(ash.res$fitted_g$pi, ash.res2$fitted_g$pi, tol = 1e-6)
  expect_equal(ash.res$result, ash.res2$result)
  expect_equal(ash.res$loglik, ash.res2$loglik)
})
