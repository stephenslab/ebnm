context("ebnm_normal_scale_mixture")

n <- 1000
set.seed(1)
s <- rnorm(n, 1, 0.1)
x <- c(rnorm(n / 2, 0, 10 + s), rnorm(n / 2, 0, s))

# The control parameters are different depending on whether or not mixsqp
#   is called from ashr or not, so just pass them in here.
control = list(eps = 1e-6, numiter.em = 20)

test_that("Basic ebnm_normal_scale_mixture functionality works", {
  ebnm.res <- ebnm(x, s, prior_family = "normal_scale_mixture", control = control)
  ebnm.res2 <- ebnm_normal_scale_mixture(x, s, control = control)
  ebnm.res$call <- ebnm.res2$call <- NULL
  expect_identical(ebnm.res, ebnm.res2)
})

test_that("ebnm_normal_scale_mixture returns the same results as ashr", {
  ebnm.res <- ebnm_normal_scale_mixture(x, s, output = ebnm_output_all(),
                                        scale = ebnm_scale_normalmix(x, s, mode = 0),
                                        control = control)
  ash.res  <- ebnm_ash(x, s, mixcompdist = "normal", prior = "uniform",
                       scale = ebnm_scale_normalmix(x, s, mode = 0)[-1],
                       output = ebnm_output_all(), control = control)

  # The ash grid is the same:
  expect_equal(ebnm.res[[g_ret_str()]]$sd, ash.res[[g_ret_str()]]$sd)

  # Control defaults are different depending on whether mixsqp is called from
  #   ashr or called directly, so the following quantities can be slightly
  #   different.

  # The estimated mixture probabilities are the same:
  # expect_equal(ebnm.res[[g_ret_str()]]$pi, ash.res[[g_ret_str()]]$pi, tol = 1e-6)

  # The posterior quantities are the same:
  # expect_equal(ebnm.res[[df_ret_str()]], ash.res[[df_ret_str()]], tol = 1e-6)

  # Likelihoods are the same:
  expect_equal(ebnm.res[[llik_ret_str()]], ash.res[[llik_ret_str()]],
               tol = 1e-2, scale = 1)

  # Now repeat with scalar s.
  s <- 1
  ebnm.res <- ebnm_normal_scale_mixture(x, s, output = ebnm_output_all(),
                                        scale = ebnm_scale_normalmix(x, s, mode = 0),
                                        control = control)
  ash.res  <- ebnm_ash(x, s, mixcompdist = "normal", prior = "uniform",
                       scale = ebnm_scale_normalmix(x, s, mode = 0)[-1],
                       output = ebnm_output_all(), control = control)
  expect_equal(ebnm.res[[g_ret_str()]]$sd, ash.res[[g_ret_str()]]$sd)
  # expect_equal(ebnm.res[[g_ret_str()]]$pi, ash.res[[g_ret_str()]]$pi, tol = 1e-6)
  # expect_equal(ebnm.res[[df_ret_str()]], ash.res[[df_ret_str()]], tol = 1e-6)
  expect_equal(ebnm.res[[llik_ret_str()]], ash.res[[llik_ret_str()]],
               tol = 1e-2, scale = 1)
})

test_that("Mode estimation works", {
  x <- x + 3
  ebnm.res <- ebnm_normal_scale_mixture(x, s, mode = "est")
  expect_equal(ebnm.res[[g_ret_str()]]$mean[1], 3, tol = .01)
})
