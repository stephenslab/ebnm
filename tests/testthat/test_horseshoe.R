context("Horseshoe")

n <- 100
set.seed(1)
s <- 2
tau <- 0.3
lambda <- abs(rcauchy(n))
x <- rnorm(n, sd = lambda * s * tau) + rnorm(n, sd = s)

true_g <- horseshoe(s * tau)

hs.res <- ebnm(x, s, prior_family = "horseshoe")

test_that("Basic functionality works", {
  hs.res2 <- ebnm_horseshoe(x, s)
  hs.res$call <- hs.res2$call <- NULL
  expect_identical(hs.res, hs.res2)
  expect_equal(hs.res[[g_ret_str()]], true_g, tolerance = 0.5)
})

hs.res2 <- ebnm_horseshoe(x, s, scale = 0.5)

test_that("Fixing the scale works", {
  expect_equal(hs.res2[[g_ret_str()]], horseshoe(scale = 0.5))
})

test_that("Fixing g works", {
  hs.res3 <- ebnm_horseshoe(x, s, g_init = hs.res2[[g_ret_str()]], fix_g = TRUE)
  expect_identical(hs.res2[[g_ret_str()]], hs.res3[[g_ret_str()]])
  expect_equal(hs.res2[[llik_ret_str()]], hs.res3[[llik_ret_str()]])
})

test_that("Initializing g works", {
  hs.res3 <- ebnm_horseshoe(x, s, g_init = horseshoe(scale = 1))
  expect_equal(hs.res[[llik_ret_str()]], hs.res3[[llik_ret_str()]])
})

test_that("Any mode other than zero results in an error", {
  expect_error(ebnm(x, s, mode = 1, prior_family = "horseshoe"))
  expect_error(ebnm(x, s, mode = "estimate", prior_family = "horseshoe"))
})

test_that("Heteroskedastic standard errors do not work", {
  s_hetero <- rexp(n)
  expect_error(ebnm_horseshoe(x, s_hetero))
})

test_that("Likelihood sanity checks", {
  pn.res <- ebnm_point_normal(x, s)
  expect_true(hs.res[[llik_ret_str()]] > pn.res[[llik_ret_str()]])

  normal.x <- rnorm(n, s = 10)
  pn.res2 <- ebnm_point_normal(normal.x, s = 1)
  hs.res3 <- ebnm_horseshoe(normal.x, s = 1)
  expect_true(hs.res3[[llik_ret_str()]] < pn.res2[[llik_ret_str()]])
})

test_that("Posterior sampler works", {
  hs.res <- ebnm_horseshoe(x, s, output = c("posterior_sampler", "posterior_mean"))
  sampler <- hs.res$posterior_sampler
  zz <- capture_output(sampres <- sampler(1000))
  # The analytic posterior means are not very accurate:
  expect_equal(hs.res$posterior$mean, colMeans(sampres), tol = 1)
})

test_that("df is correct for returned logLik", {
  hs.res <- ebnm_horseshoe(x, s)
  expect_equal(attr(logLik(hs.res), "df"), 1)
  hs.res2 <- ebnm_horseshoe(x, s, scale = 1)
  expect_equal(attr(logLik(hs.res2), "df"), 0)
})

# test_that("predict method works as expected", {
#   hs.res <- ebnm_horseshoe(x, s)
#   hs.res2 <- predict(hs.res, list(x = 1:10, s = 1))
#   expect_equal(hs.res$fitted_g, hs.res2$fitted_g)
# })
