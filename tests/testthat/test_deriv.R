test_that("derivatives of point normal loglik are right",{
  library(numDeriv)
  n=100
  set.seed(1)
  s = rgamma(n,1,1)
  x=rnorm(n,0,1+s)
  a = 0.5
  w = 0.2
  d = ebnm:::grad_negloglik_point_normal(x,s,w,a)
  expect_equal(d[1],grad(function(w) -ebnm:::loglik_point_normal(x,s,w,a),w),
               tol = 1e-4)
  expect_equal(d[2],grad(function(a) -ebnm:::loglik_point_normal(x,s,w,a),a),
               tol = 1e-4)
})

test_that("derivatives of laplace loglik are right",{
  library(numDeriv)
  n=100
  set.seed(1)
  s = rgamma(n,1,1)
  x=rnorm(n,0,1+s)
  a = 0.5
  w = 0.2
  d = ebnm:::grad_negloglik_laplace(x,s,w,a)
  expect_equal(d[1],grad(function(w) -ebnm:::loglik_laplace(x,s,w,a),w),
               tol = 1e-4)
  expect_equal(d[2],grad(function(a) -ebnm:::loglik_laplace(x,s,w,a),a),
               tol = 1e-4)
})
