context("Optimization Methods")

set.seed(666)
x <- rnorm(30, sd = 2)

test_that("Different optimization methods work", {
  n_res <- ebnm_normal(x, optmethod = "nlm", control = list(ndigit = 8))
  n_res2 <- ebnm(x, prior_family = "normal", optmethod = "nograd_lbfgsb", control = list(ndeps = 1e-2))
  expect_equal(n_res$fitted_g$sd, n_res2$fitted_g$sd, 1e-4)

  pn_res <- ebnm_point_normal(x, optmethod = "nograd_nlm", control = list(gradtol = 1e-4))
  pn_res2 <- ebnm(x, prior_family = "point_normal", optmethod = "trust", control = list(fterm = 1e-6))
  expect_equal(pn_res$fitted_g$sd[2], pn_res2$fitted_g$sd[2], 1e-4)

  pl_res <- ebnm_point_laplace(x, optmethod = "trust", control = list(iterlim = 50))
  pl_res2 <- ebnm(x, prior_family = "point_laplace", optmethod = "nohess_nlm", control = list(steptol = 1e-4))
  expect_equal(pl_res$fitted_g$scale[2], pl_res2$fitted_g$scale[2], 1e-4)

  pe_res <- ebnm_point_exponential(x, optmethod = "nlm", control = list(check.analyticals = TRUE))
  pe_res2 <- ebnm_point_exponential(x, optmethod = "lbfgsb", control = list(maxit = 10))
  expect_equal(pe_res$fitted_g$scale[2], pe_res2$fitted_g$scale[2], 1e-4)

  # smn_res <- ebnm_normal_scale_mixture(x, optmethod = "mixIP", control = list(rtol = 1e-4))
  smn_res <- ebnm_normal_scale_mixture(x, optmethod = "mixSQP", control = list(numiter.em = 10))
  smn_res2 <- ebnm(x, prior_family = "normal_scale_mixture", optmethod = "mixEM", control = list(tol = 1e-6))
  expect_equal(smn_res$log_likelihood, smn_res2$log_likelihood, 1e-1)

  uni_res <- ebnm_unimodal(x, optmethod = "mixVBEM", control = list(tol = 1e-5))
  uni_res2 <- ebnm(x, prior_family = "unimodal", optmethod = "mixSQP", control = list(maxiter.activeset = 10))
  expect_equal(uni_res$log_likelihood, uni_res2$log_likelihood, 1e-2)

  symmuni_res <- ebnm_unimodal_symmetric(x, optmethod = "w_mixEM", control = list(tol = 1e-5))
  symmuni_res2 <- ebnm(x, prior_family = "unimodal_symmetric", optmethod = "mixSQP", control = list(maxiter.sqp = 5))
  expect_equal(symmuni_res$log_likelihood, symmuni_res2$log_likelihood, 1e-2)

  # Remove cxxMixSquarem test until fixed within ashr:
  # nnuni_res <- ebnm_unimodal_nonnegative(x, optmethod = "cxxMixSquarem", control = list(tol = 1e-4))
  # nnuni_res2 <- ebnm(x, prior_family = "unimodal_nonnegative", optmethod = "mixEM", control = list(maxiter = 100))
  # expect_equal(nnuni_res$log_likelihood, nnuni_res2$log_likelihood, 1e-2)

  npuni_res <- ebnm_unimodal_nonpositive(x, optmethod = "mixSQP", control = list(numiter.em = 2))
  npuni_res2 <- ebnm(x, prior_family = "unimodal_nonpositive", optmethod = "mixVBEM", control = list(step.max0 = 3))
  expect_equal(npuni_res$log_likelihood, npuni_res2$log_likelihood, 1e-2)

  npmle_res <- ebnm_npmle(x, optmethod = "mixSQP", control = list(numiter.em = 5))
  npmle_res2 <- ebnm_npmle(x, optmethod = "mixEM", control = list(maxiter = 500))
  # npmle_res2 <- ebnm(x, prior_family = "npmle", optmethod = "REBayes", control = list(rtol = 1e-8))
  expect_equal(npmle_res$log_likelihood, npmle_res2$log_likelihood, 1e-1)

  deconv_res <- ebnm_deconvolver(x, control = list(ndigit = 6))
  deconv_res2 <- ebnm(x, prior_family = "deconvolver")
  expect_equal(deconv_res$log_likelihood, deconv_res2$log_likelihood)

  hs_res <- ebnm_horseshoe(x, control = list(tol = 1e-8))
  hs_res2 <- ebnm(x, prior_family = "horseshoe")
  expect_equal(hs_res$fitted_g$scale, hs_res2$fitted_g$scale, 1e-4)

  gb_res <- ebnm_generalized_binary(x, control = list(tol = 1e-8))
  gb_res2 <- ebnm(x, prior_family = "generalized_binary")
  expect_equal(gb_res$log_likelihood, gb_res2$log_likelihood, 1e-4)

  pm_res <- ebnm_point_mass(x, mode = "estimate", control = list(tol = 1e-8))
  pm_res2 <- ebnm(x, prior_family = "point_mass", mode = "estimate", control = list(interval = c(-1, 1)))
  expect_equal(pm_res$fitted_g$mean, pm_res2$fitted_g$mean, 1e-4)
})

test_that("Warnings get thrown when they should be", {
  expect_warning(ebnm(x, prior_family = "horseshoe", optmethod = "nlm"))
  expect_error(ebnm_generalized_binary(x, optmethod = "nlm"))
  expect_warning(ebnm(x, prior_family = "generalized_binary", optmethod = "nlm"))
  expect_warning(ebnm(x, prior_family = "flat", optmethod = "nlm"))
  expect_warning(ebnm(x, prior_family = "point_mass", optmethod = "nlm"))
})
