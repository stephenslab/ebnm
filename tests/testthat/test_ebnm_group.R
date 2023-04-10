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

test_that("Basic functionality works", {
  ebnm.res <- ebnm_group(x, s, group)
  rmse.ebnm <- sqrt(mean((ebnm.res$posterior$mean - theta)^2))
  expect_true(rmse.ebnm < 1)
})

# TODO: One group; group with single observation; mode, scale, sampler, init_g, fixg, etc.
