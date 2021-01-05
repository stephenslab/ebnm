context("NPMLE")

n <- 1000
set.seed(1)
s <- rnorm(n, 1, 0.1)
x <- runif(n, -10, 10)

cdf_grid <- seq(-10, 10, 0.1)

true_g <- ashr::unimix(1, -10, 10)
true_cdf <- ashr::comp_cdf(true_g, cdf_grid)

test_that("Basic functionality works", {
  npmle.res <- ebnm(x, s, prior_family = "npmle")
  npmle.res2 <- ebnm_npmle(x, s)
  expect_identical(npmle.res, npmle.res2)

  est_cdf <- drop(npmle.res$fitted_g$pi %*% ashr::comp_cdf(npmle.res$fitted_g, cdf_grid))
  expect_equal(true_cdf, est_cdf, tolerance = 0.1)
})

test_that("Fixing the scale works", {
  npmle.res <- ebnm_npmle(x, s, scale = 1)
  g_scale <- npmle.res$fitted_g$a[2] - npmle.res$fitted_g$a[1]
  expect_equal(1, g_scale, tolerance = 0.1)
})

test_that("Fixing g works", {
  npmle.res <- ebnm_npmle(x, s, g_init = true_g, fix_g = TRUE)
  expect_identical(npmle.res[[g_ret_str()]], true_g)
})
