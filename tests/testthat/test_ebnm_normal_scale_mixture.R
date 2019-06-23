context("ebnm_normal_scale_mixture")

n <- 1000
set.seed(1)
s <- rnorm(n, 1, 0.1)
x <- c(rnorm(n / 2, 0, 10 + s), rnorm(n / 2, 0, s))

test_that("Basic ebnm_normal_scale_mixture functionality works", {
  ash.res <- ebnm(x, s, prior_type = "normal_scale_mixture")
  ash.res2 <- ebnm_normal_scale_mixture(x, s)
  expect_identical(ash.res, ash.res2)
  expect_true(inherits(ash.res2$fitted_g, "normalmix"))
})
