context("ebnm_unimodal")

n <- 1000
set.seed(1)
s <- rnorm(n, 1, 0.1)
x <- c(rnorm(n / 2, 0, 10 + s), rnorm(n / 2, 0, s))

test_that("Basic ebnm_unimodal functionality works", {
  ash.res <- ebnm(x, s, prior_family = "unimodal")
  ash.res2 <- ebnm_unimodal(x, s)
  ash.res$call <- ash.res2$call <- NULL
  expect_identical(ash.res, ash.res2)
  expect_true(inherits(ash.res2[[g_ret_str()]], "unimix"))
  expect_false(identical(ash.res2[[g_ret_str()]]$a, -ash.res2[[g_ret_str()]]$b))
})

test_that("Basic ebnm_unimodal_symmetric functionality works", {
  ash.res <- ebnm(x, s, prior_family = "unimodal_symmetric")
  ash.res2 <- ebnm_unimodal_symmetric(x, s)
  ash.res$call <- ash.res2$call <- NULL
  expect_identical(ash.res, ash.res2)
  expect_true(inherits(ash.res2[[g_ret_str()]], "unimix"))
  expect_identical(ash.res2[[g_ret_str()]]$a, -ash.res2[[g_ret_str()]]$b)
})

test_that("Basic ebnm_unimodal_nonnegative functionality works", {
  ash.res <- ebnm(x, s, prior_family = "unimodal_nonnegative")
  ash.res2 <- ebnm_unimodal_nonnegative(x, s)
  ash.res$call <- ash.res2$call <- NULL
  expect_identical(ash.res, ash.res2)
  expect_true(inherits(ash.res2[[g_ret_str()]], "unimix"))
  expect_identical(ash.res2[[g_ret_str()]]$a, rep(0, length(ash.res2[[g_ret_str()]]$a)))
})

test_that("Basic ebnm_unimodal_nonpositive functionality works", {
  ash.res <- ebnm(x, s, prior_family = "unimodal_nonpositive")
  ash.res2 <- ebnm_unimodal_nonpositive(x, s)
  ash.res$call <- ash.res2$call <- NULL
  expect_identical(ash.res, ash.res2)
  expect_true(inherits(ash.res2[[g_ret_str()]], "unimix"))
  expect_identical(ash.res2[[g_ret_str()]]$b, rep(0, length(ash.res2[[g_ret_str()]]$b)))
})
