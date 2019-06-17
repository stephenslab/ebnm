context("Point normal optimization functions")

n = 100
set.seed(1)
s = rnorm(n, 1, 0.1)
x = 2 + c(rnorm(n / 4, 0, 10 + s), rnorm(3 * n / 4, 0, s))

true_pi0 = 0.25
true_a = 0.01
true_mu = 2

z = (x - true_mu)^2 / s^2
sum.z = sum(z)

test_that("pn_fn value agrees with loglik_point_normal value", {
  llik = loglik_point_normal(x, s, w = 1 - true_pi0, a = true_a, mu = true_mu)
  optval = pn_fn(NULL, fix_pi0 = TRUE, fix_a = TRUE, fix_mu = TRUE,
                 alpha = -log(1 / true_pi0 - 1), beta = -log(true_a), mu = true_mu,
                 n0 = 0, n1 = 0, sum1 = 0, n2 = n, x = x, s2 = s^2, z = z,
                 sum.z = sum.z)
  llik_from_optval = pn_llik_from_optval(optval, n1 = 0, n2 = n, s2 = s^2)
  expect_equal(llik, llik_from_optval)
})

test_that("pn_gr returns correct gradient", {
  # alpha is defined as -log(1 / pi0 - 1).
  llik_alpha = function(alpha) {
    pi0 = 1 / (1 + exp(-alpha))
    return(loglik_point_normal(x, s, w = 1 - pi0, a = true_a, mu = true_mu))
  }
  grad_alpha = numDeriv::grad(llik_alpha, x = -log(1 / true_pi0 - 1))
  grad_alpha_from_optim = pn_gr(-log(1 / true_pi0 - 1),
                                fix_pi0 = FALSE, fix_a = TRUE, fix_mu = TRUE,
                                alpha = NULL, beta = -log(true_a), mu = true_mu,
                                n0 = 0, n1 = 0, sum1 = 0, n2 = n,
                                x = x, s2 = s^2, z = z, sum.z = sum.z)

  # beta is defined as -log(a).
  llik_beta = function(beta) {
    a = exp(-beta)
    return(loglik_point_normal(x, s, w = 1 - true_pi0, a = a, mu = true_mu))
  }
  grad_beta = numDeriv::grad(llik_beta, x = -log(true_a))
  grad_beta_from_optim = pn_gr(-log(true_a),
                               fix_pi0 = TRUE, fix_a = FALSE, fix_mu = TRUE,
                               alpha = -log(1 / true_pi0 - 1),
                               beta = NULL, mu = true_mu,
                               n0 = 0, n1 = 0, sum1 = 0, n2 = n,
                               x = x, s2 = s^2, z = z, sum.z = sum.z)

  llik_mu = function(mu) {
    return(loglik_point_normal(x, s, w = 1 - true_pi0, a = true_a, mu = mu))
  }
  grad_mu = numDeriv::grad(llik_mu, x = true_mu)
  grad_mu_from_optim = pn_gr(true_mu,
                             fix_pi0 = TRUE, fix_a = TRUE, fix_mu = FALSE,
                             alpha = -log(1 / true_pi0 - 1),
                             beta = -log(true_a), mu = NULL,
                             n0 = 0, n1 = 0, sum1 = 0, n2 = n,
                             x = x, s2 = s^2, z = z, sum.z = sum.z)

  all_grad_from_optim = pn_gr(c(-log(1 / true_pi0 - 1), -log(true_a), true_mu),
                              fix_pi = FALSE, fix_a = FALSE, fix_mu = FALSE,
                              alpha = NULL, beta = NULL, mu = NULL,
                              n0 = 0, n1 = 0, sum1 = 0, n2 = n,
                              x = x, s2 = s^2, z = z, sum.z = sum.z)

  # optim works on the negative log likelihood.
  expect_equal(grad_alpha, -grad_alpha_from_optim)
  expect_equal(grad_beta, -grad_beta_from_optim)
  expect_equal(grad_mu, -grad_mu_from_optim)
  expect_equal(all_grad_from_optim,
               -c(grad_alpha, grad_beta, grad_mu))
})

s[1:2] <- 0
x[(n - 1):n] <- true_mu
s[(n - 1):n] <- 0
sum1 <- sum((x[1:2] - true_mu)^2)
xsub <- x[-which(s == 0)]
ssub <- s[-which(s == 0)]
z <- z[-which(s == 0)]
sum.z <- sum(z)

test_that("pn_fn and loglik_point_normal agree when some SEs are zero", {
  llik = loglik_point_normal(x, s, w = 1 - true_pi0, a = true_a, mu = true_mu)

  optval = pn_fn(NULL, fix_pi0 = TRUE, fix_a = TRUE, fix_mu = TRUE,
                 alpha = -log(1 / true_pi0 - 1), beta = -log(true_a), mu = true_mu,
                 n0 = 2, n1 = 2, sum1 = sum1, n2 = n - 4, x = xsub, s2 = ssub^2,
                 z = z, sum.z = sum.z)
  llik_from_optval = pn_llik_from_optval(optval, n1 = 2, n2 = n - 4, s2 = ssub^2)
  expect_equal(llik, llik_from_optval)
})

test_that("pn_gr returns correct gradient when some SEs are zero", {
  # alpha is defined as -log(1 / pi0 - 1).
  llik_alpha = function(alpha) {
    pi0 = 1 / (1 + exp(-alpha))
    return(loglik_point_normal(x, s, w = 1 - pi0, a = true_a, mu = true_mu))
  }
  grad_alpha = numDeriv::grad(llik_alpha, x = -log(1 / true_pi0 - 1))
  grad_alpha_from_optim = pn_gr(-log(1 / true_pi0 - 1),
                                fix_pi0 = FALSE, fix_a = TRUE, fix_mu = TRUE,
                                alpha = NULL, beta = -log(true_a), mu = true_mu,
                                n0 = 2, n1 = 2, sum1 = sum1, n2 = n - 4,
                                x = xsub, s2 = ssub^2, z = z, sum.z = sum.z)

  # beta is defined as -log(a).
  llik_beta = function(beta) {
    a = exp(-beta)
    return(loglik_point_normal(x, s, w = 1 - true_pi0, a = a, mu = true_mu))
  }
  grad_beta = numDeriv::grad(llik_beta, x = -log(true_a))
  grad_beta_from_optim = pn_gr(-log(true_a),
                               fix_pi0 = TRUE, fix_a = FALSE, fix_mu = TRUE,
                               alpha = -log(1 / true_pi0 - 1),
                               beta = NULL, mu = true_mu,
                               n0 = 2, n1 = 2, sum1 = sum1, n2 = n - 4,
                               x = xsub, s2 = ssub^2, z = z, sum.z = sum.z)

  expect_equal(grad_alpha, -grad_alpha_from_optim)
  expect_equal(grad_beta, -grad_beta_from_optim)
})
