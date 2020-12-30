context("Point Laplace optimization functions")

n = 100
set.seed(1)
s = rnorm(n, 1, 0.1)
x = c(rexp(n / 4, rate = 0.1), rep(0, 3 * n / 4)) + s + 1

true_pi0 = 0.25
true_a = 0.1
true_mu = 1

par = c(log(1 / true_pi0 - 1), log(true_a), true_mu)

optval = pl_nllik(par, x, s, par_init = NULL, fix_par = c(FALSE, FALSE, FALSE),
                  calc_grad = TRUE, calc_hess = TRUE)

test_that("pn_nlm_fn value agrees with loglik_point_normal value", {
  true_llik = loglik_point_laplace(x, s,
                                   w = 1 - true_pi0,
                                   a = true_a,
                                   mu = true_mu)
  expect_equivalent(-true_llik, optval)
})

true_llik_fn <- function(par) {
  w = 1 / (1 + exp(-par[1]))
  a = exp(par[2])
  mu = par[3]
  return(-loglik_point_laplace(x, s, w = w, a = a, mu = mu))
}

true_grad <- numDeriv::grad(true_llik_fn, x = par)
test_that("pl_nlm_fn returns correct gradient", {
  expect_equal(true_grad, attr(optval, "gradient"))
})

true_hess <- numDeriv::hessian(true_llik_fn, x = par)
test_that("pl_nlm_fn returns correct Hessian", {
  expect_equal(true_hess, attr(optval, "hessian"))
})
