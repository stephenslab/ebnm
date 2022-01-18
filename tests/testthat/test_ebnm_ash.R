context("ebnm_ash")

n <- 1000
set.seed(1)
s <- rnorm(n, 1, 0.1)
x <- c(rnorm(n / 2, 0, 10 + s), rnorm(n / 2, 0, s))

true_pi0 <- 0.5
true_mean <- 0
true_sd <- 10
true_g <- ashr::normalmix(pi = c(true_pi0, 1 - true_pi0),
                          mean = rep(true_mean, 2),
                          sd = c(0, true_sd))

test_that("Basic functionality works", {
  ash.res <- ebnm(x, s, prior_family = "ash")
  ash.res2 <- ebnm_ash(x, s)
  ash.res$call <- ash.res2$call <- NULL
  expect_identical(ash.res, ash.res2)
})

test_that("Mode estimation works", {
  ash.res <- ebnm_ash(x, s, mode = "est", mixcompdist = "normal")
  expect_false(identical(ash.res[[g_ret_str()]]$mean[1], true_mean))
})

test_that("Fixing the sd works", {
  ash.res <- ebnm_ash(x, s, scale = true_sd, mixcompdist = "normal")
  expect_equal(ash.res[[g_ret_str()]], true_g, tolerance = 0.1)
  expect_identical(ash.res[[g_ret_str()]]$sd[2], true_sd)
})

test_that("Fixing g works", {
  ash.res <- ebnm_ash(x, s, g_init = true_g, fix_g = TRUE)
  expect_identical(ash.res[[g_ret_str()]], true_g)
  ash.res2 <- ebnm(x, s, g_init = true_g, fix_g = TRUE, prior_family = "ash")
  ash.res$call <- ash.res2$call <- NULL
  expect_identical(ash.res, ash.res2)
})

test_that("Output parameter works", {
  ash.res <- ebnm_ash(x, s, output = samp_arg_str())
  ash.res$call <- NULL
  expect_identical(names(ash.res), samp_ret_str())
})
