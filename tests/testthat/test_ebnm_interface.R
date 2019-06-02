context("'ebnm' interface function")

n = 1000
set.seed(1)
s = rnorm(n, 1, 0.1)
x = c(rnorm(n / 2, 0, 10 + s), rnorm(n / 2, 0, s))

true_pi0 = 0.5
true_a = 0.01
true_mu = 0

test_that("point_normal works with nothing fixed", {
  ebnm.res = ebnm(x, s, "point_normal", fix_mu = FALSE)
  ebnm.pn.res = ebnm_point_normal(x, s, fix_mu = FALSE)
  expect_equal(ebnm.res, ebnm.pn.res)
  expect_equal(ebnm.res$fitted_g, list(pi0 = true_pi0, a = true_a, mu = true_mu),
               tolerance = 0.1)
})

test_that("point_normal works with pi0 fixed", {
  ebnm.res = ebnm(x, s, "point_normal", g = list(pi0 = true_pi0),
                  fix_pi0 = TRUE, fix_mu = FALSE)
  ebnm.pn.res = ebnm_point_normal(x, s, g = list(pi0 = true_pi0),
                                  fix_pi0 = TRUE, fix_mu = FALSE)
  expect_equal(ebnm.res, ebnm.pn.res)
  expect_equal(ebnm.res$fitted_g, list(pi0 = true_pi0, a = true_a, mu = true_mu),
               tolerance = 0.1)
})

test_that("point_normal works with pi0 fixed", {
  ebnm.res = ebnm(x, s, "point_normal", g = list(a = true_a),
                  fix_a = TRUE, fix_mu = FALSE)
  ebnm.pn.res = ebnm_point_normal(x, s, g = list(a = true_a),
                                  fix_a = TRUE, fix_mu = FALSE)
  expect_equal(ebnm.res, ebnm.pn.res)
  expect_equal(ebnm.res$fitted_g, list(pi0 = true_pi0, a = true_a, mu = true_mu),
               tolerance = 0.1)
})

test_that("point_normal works with mu fixed", {
  ebnm.res = ebnm(x, s, "point_normal")
  ebnm.pn.res = ebnm_point_normal(x, s)
  expect_equal(ebnm.res, ebnm.pn.res)
  expect_equal(ebnm.res$fitted_g, list(pi0 = true_pi0, a = true_a, mu = true_mu),
               tolerance = 0.1)
})

test_that("point_normal works with pi0 and mu fixed", {
  ebnm.res = ebnm(x, s, "point_normal", g = list(pi0 = true_pi0, mu = true_mu),
                  fix_pi0 = TRUE, fix_mu = TRUE)
  ebnm.pn.res = ebnm_point_normal(x, s, g = list(pi0 = true_pi0, mu = true_mu),
                                  fix_pi0 = TRUE, fix_mu = TRUE)
  expect_equal(ebnm.res, ebnm.pn.res)
  expect_equal(ebnm.res$fitted_g, list(pi0 = true_pi0, a = true_a, mu = true_mu),
               tolerance = 0.1)
})

test_that("point_normal works with a and mu fixed", {
  ebnm.res = ebnm(x, s, "point_normal", g = list(a = true_a, mu = true_mu),
                  fix_a = TRUE, fix_mu = TRUE)
  ebnm.pn.res = ebnm_point_normal(x, s, g = list(a = true_a, mu = true_mu),
                                  fix_a = TRUE, fix_mu = TRUE)
  expect_equal(ebnm.res, ebnm.pn.res)
  expect_equal(ebnm.res$fitted_g, list(pi0 = true_pi0, a = true_a, mu = true_mu),
               tolerance = 0.1)
})

test_that("point_normal works with a and pi0 fixed", {
  ebnm.res = ebnm(x, s, "point_normal", g = list(a = true_a, pi0 = true_pi0),
                  fix_a = TRUE, fix_pi0 = TRUE, fix_mu = FALSE)
  ebnm.pn.res = ebnm_point_normal(x, s, g = list(a = true_a, pi0 = true_pi0),
                                  fix_a = TRUE, fix_pi0 = TRUE, fix_mu = FALSE)
  expect_equal(ebnm.res, ebnm.pn.res)
  expect_equal(ebnm.res$fitted_g, list(pi0 = true_pi0, a = true_a, mu = true_mu),
               tolerance = 0.1)
})

test_that("control parameter works", {
  ebnm.res = ebnm(x, s)
  ebnm.res.ctrl = ebnm(x, s, control = list(ndigit = 12))
  expect_equal(ebnm.res, ebnm.res.ctrl)
})

test_that("point_laplace works", {
  ebnm.res = ebnm(x, s, "point_laplace")
  ebnm.pl.res = ebnm_point_laplace(x, s)
  ebnm.pl.res$fitted_g$mu = 0 # have to manually add this for test.
  expect_equal(ebnm.res, ebnm.pl.res)
})
