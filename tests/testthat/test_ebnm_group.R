context("ebnm_group")

n <- 1000
grps <- c("A", "B", "C")
sds <- c(3, 10, 100)
pi0s <- c(0, 0.5, 0.9)

set.seed(1)

theta <- rep(0, n)
s <- rnorm(n, 1, 0.1)
group <- sample(grps, n, replace = TRUE)

for (i in 1:length(grps)) {
  n_grp <- sum(group == grps[i])
  x_grp <- rnorm(n_grp, sd = sds[i])
  x_grp <- x_grp * rbinom(n_grp, 1, prob = 1 - pi0s[i])
  theta[group == grps[i]] <- x_grp
}
x <- theta + rnorm(n)
ebnm.res <- ebnm_group(x, s, group)

test_that("Basic functionality works", {
  rmse.ebnm <- sqrt(mean((ebnm.res$posterior$mean - theta)^2))
  expect_true(rmse.ebnm < 1)
})

test_that("Works with only one group", {
  ebnm.res <- ebnm_group(x, s, rep("A", length(x)))
  expect_identical(names(ebnm.res$fitted_g), "A")
})

test_that("Works when group has a single observation", {
  group1 <- group
  group1[1] <- "Z"
  ebnm.res <- ebnm_group(x, s, group1)
  expect_true("Z" %in% names(ebnm.res$fitted_g))
})

test_that("Named list arguments work", {
  pf_list <- list(
    A = "point_normal",
    B = "point_laplace",
    C = "normal"
  )
  mode_list <- list(
    A = "estimate",
    B = 0,
    C = 1
  )
  scale_list <- list(
    A = 3,
    B = "estimate",
    C = "estimate"
  )
  ebnm.res <- ebnm_group(x, s, group,
                         prior_family = pf_list,
                         mode = mode_list,
                         scale = scale_list)
  expect_identical(class(ebnm.res$fitted_g[["A"]]), "normalmix")
  expect_identical(class(ebnm.res$fitted_g[["B"]]), "laplacemix")
  expect_failure(expect_equal(ebnm.res$fitted_g[["A"]]$mean[1], 0))
  expect_equal(ebnm.res$fitted_g[["B"]]$mean[1], 0)
  expect_equal(ebnm.res$fitted_g[["A"]]$sd[2], 3)
  expect_failure(expect_equal(ebnm.res$fitted_g[["A"]]$sd[2], 10))
})

test_that("Argument g_init works", {
  g_init <- ebnm.res$fitted_g
  ebnm.res <- ebnm_group(x, s, group, g_init = g_init)
  rmse.ebnm <- sqrt(mean((ebnm.res$posterior$mean - theta)^2))
  expect_true(rmse.ebnm < 1)
})
