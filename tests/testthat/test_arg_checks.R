context("Argument checks")

n <- 1000
set.seed(1)
s <- rnorm(n, 1, 0.1)
x <- c(rnorm(n / 2, 0, 10 + s), rnorm(n / 2, 0, s))


test_that("g_init is valid", {
  g_init <- laplacemix(c(0.5, 0.5), c(0, 0), c(0, 1))
  expect_error(ebnm_point_normal(x, s, g_init = g_init))

  g_init <- normalmix(c(1/3, 1/3, 1/3), c(0, 0, 0), c(0, 1, 2))
  expect_error(ebnm_point_normal(x, s, g_init = g_init))

  g_init <- normalmix(c(0.5, 0.5), c(0, 0), c(1, 2))
  expect_error(ebnm_point_normal(x, s, g_init = g_init))

  g_init <- normalmix(c(0.5, 0.5), c(1, 1), c(0, 10))
  ebnm.res <- ebnm_point_normal(x, s, g_init = g_init)
  expect_error(ebnm_point_normal(x, s, g_init = g_init, mode = 2))
  expect_error(ebnm_point_normal(x, s, g_init = g_init, scale = 1))
})
