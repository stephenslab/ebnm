context("Prior family inference")

n <- 1000
set.seed(1)
x <- runif(n, -10, 10) + rnorm(n)
s <- 1

test_family <- function(family) {
  ebnm.res <- ebnm(x, s, prior_family = family)
  expect_identical(infer_prior_family(ebnm.res$fitted_g), family)
}

test_that("all families are inferred correctly", {
  test_family("point_normal")
  test_family("point_laplace")
  test_family("point_exponential")
  test_family("normal")
  test_family("horseshoe")
  test_family("normal_scale_mixture")
  test_family("unimodal")
  test_family("unimodal_symmetric")
  test_family("unimodal_nonnegative")
  test_family("unimodal_nonpositive")
  test_family("generalized_binary")
  test_family("npmle")
  test_family("flat")
  test_family("point_mass")
})
