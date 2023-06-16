context("Posterior samplers")

set.seed(1)
mu <- c(rep(0, 25), rexp(25))
s <- rgamma(50, 1, 1)
x <- mu + rnorm(50, 0, s)

test_that("point-normal sampler gives reasonable results", {
  res <- ebnm_point_normal(x, s, output = ebnm_output_all())
  samp <- res$posterior_sampler(1000)

  expect_equal(colMeans(samp), res[[df_ret_str()]][[pm_ret_str()]], tol = 0.1)
  expect_equal(colMeans(samp^2), res[[df_ret_str()]][[pm2_ret_str()]], tol = 0.1)
})

test_that("point-Laplace sampler gives reasonable results", {
  res <- ebnm_point_laplace(x, s, output = ebnm_output_all())
  samp <- res$posterior_sampler(1000)

  expect_equal(colMeans(samp), res[[df_ret_str()]][[pm_ret_str()]], tol = 0.1)
  expect_equal(colMeans(samp^2), res[[df_ret_str()]][[pm2_ret_str()]], tol = 0.1)
})

test_that("point-exponential sampler gives reasonable results", {
  res <- ebnm_point_exponential(x, s, output = ebnm_output_all())
  samp <- res$posterior_sampler(1000)

  expect_equal(colMeans(samp), res[[df_ret_str()]][[pm_ret_str()]], tol = 0.1)
  expect_equal(colMeans(samp^2), res[[df_ret_str()]][[pm2_ret_str()]], tol = 0.1)
})

test_that("normal-mixture sampler gives reasonable results", {
  res <- ebnm_normal_scale_mixture(x, s, output = ebnm_output_all())
  samp <- res$posterior_sampler(1000)

  expect_equal(colMeans(samp), res[[df_ret_str()]][[pm_ret_str()]], tol = 0.1)
  expect_equal(colMeans(samp^2), res[[df_ret_str()]][[pm2_ret_str()]], tol = 0.1)
})
