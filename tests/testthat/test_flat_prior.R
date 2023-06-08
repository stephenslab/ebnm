context("flat")

x <- 1:10
s <- 1

test_that("Basic functionality works", {
  flat.res <- ebnm(x, s, prior_family = "flat")
  flat.res2 <- ebnm_flat(x, s)
  flat.res$call <- flat.res2$call <- NULL
  expect_identical(flat.res, flat.res2)
  expect_identical(x, flat.res$posterior$mean)
})

test_that("Parameters that do nothing are ignored", {
  expect_warning(flat.res <- ebnm(x, s, prior_family = "flat", mode = 1))
  expect_warning(flat.res <- ebnm(x, s, prior_family = "flat", scale = 1))
  ginit <- ashr::normalmix(1, 0, 1)
  expect_warning(flat.res <- ebnm(x, s, prior_family = "flat", g_init = ginit))
})
