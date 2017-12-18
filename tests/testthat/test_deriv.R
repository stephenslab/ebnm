test_that("derivatives of normal loglik are right",{
  n=100
  set.seed(1)
  s = rgamma(n,1,1)
  x=rnorm(n,0,1+s)
  a = 0.5
  w = 0.2
  d = grad_negloglik_normal(x,s,w,a)
  expect_equal(d[1],numDeriv::grad(function(w){-loglik_normal(x,s,w,a)},w),tol=1e-4)
  expect_equal(d[2],numDeriv::grad(function(a){-loglik_normal(x,s,w,a)},a),tol=1e-4)
})

test_that("derivatives of laplace loglik are right",{
  n=100
  set.seed(1)
  s = rgamma(n,1,1)
  x=rnorm(n,0,1+s)
  a = 0.5
  w = 0.2
  d = grad_negloglik_laplace(x,s,w,a)
  expect_equal(d[1],numDeriv::grad(function(w){-loglik_laplace(x,s,w,a)},w),tol=1e-4)
  expect_equal(d[2],numDeriv::grad(function(a){-loglik_laplace(x,s,w,a)},a),tol=1e-4)
})
