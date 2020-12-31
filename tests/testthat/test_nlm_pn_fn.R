context("Point normal optimization functions")

n = 100
set.seed(1)
s = rnorm(n, 1, 0.1)
x = 2 + c(rnorm(n / 4, 0, 10 + s), rnorm(3 * n / 4, 0, s))

true_pi0 = 0.25
true_a = 0.01
true_mu = 2

z = (x - true_mu)^2 / s^2
sum_z = sum(z)

par = c(-log(1 / true_pi0 - 1), -log(true_a), true_mu)

optval = pn_nllik(par, x, s, par_init = NULL,
                  fix_par = c(FALSE, FALSE, FALSE),
                  n0 = 0, n1 = 0, sum1 = 0, n2 = n,
                  s2 = s^2, z = z, sum_z = sum_z,
                  calc_grad = TRUE, calc_hess = TRUE)

test_that("pn_nlm_fn value agrees with loglik_point_normal value", {
  llik_from_optval = pn_llik_from_optval(optval, n1 = 0, n2 = n, s2 = s^2)
  true_llik = loglik_point_normal(x, s,
                                  w = 1 - true_pi0,
                                  a = true_a,
                                  mu = true_mu)
  expect_equivalent(true_llik, llik_from_optval)
})

true_llik_fn <- function(par) {
  pi0 = 1 / (1 + exp(-par[1]))
  a = exp(-par[2])
  mu = par[3]
  return(-loglik_point_normal(x, s, w = 1 - pi0, a = a, mu = mu))
}

true_grad <- numDeriv::grad(true_llik_fn, x = par)
test_that("pn_nlm_fn returns correct gradient", {
  expect_equal(true_grad, attr(optval, "gradient"))
})

true_hess <- numDeriv::hessian(true_llik_fn, x = par)
test_that("pn_nlm_fn returns correct Hessian", {
  expect_equal(true_hess, attr(optval, "hessian"))
})

s[1:2] <- 0
x[(n - 1):n] <- true_mu
s[(n - 1):n] <- 0
sum1 <- sum((x[1:2] - true_mu)^2)
xsub <- x[-which(s == 0)]
ssub <- s[-which(s == 0)]
z <- z[-which(s == 0)]
sum_z <- sum(z)

par <- par[1:2]
optval = pn_nllik(par, xsub, ssub, par_init = NULL, fix_par = c(FALSE, FALSE, TRUE),
                  n0 = 2, n1 = 2, sum1 = sum1, n2 = n - 4,
                  s2 = ssub^2, z = z, sum_z = sum_z,
                  calc_grad = TRUE, calc_hess = TRUE)

test_that("pn_nlm_fn and loglik_point_normal agree when some SEs are zero", {
  llik_from_optval = pn_llik_from_optval(optval, n1 = 2, n2 = n - 4,
                                         s2 = ssub^2)
  true_llik = loglik_point_normal(x, s,
                                  w = 1 - true_pi0,
                                  a = true_a,
                                  mu = true_mu)
  expect_equivalent(true_llik, llik_from_optval)
})

true_llik_fn <- function(par) {
  pi0 = 1 / (1 + exp(-par[1]))
  a = exp(-par[2])
  return(-loglik_point_normal(x, s, w = 1 - pi0, a = a, mu = true_mu))
}

true_grad <- numDeriv::grad(true_llik_fn, x = par)
test_that("pn_nlm_fn returns correct gradient when some SEs are zero", {
  expect_equal(true_grad, attr(optval, "gradient"))
})

true_hess <- numDeriv::hessian(true_llik_fn, x = par)
test_that("pn_nlm_fn returns correct Hessian when some SEs are zero", {
  expect_equal(true_hess, attr(optval, "hessian"))
})
