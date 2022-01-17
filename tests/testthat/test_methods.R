context("methods")

n <- 1000
set.seed(1)
s <- 1
x <- c(rnorm(n / 2, 0, 10 + s), rnorm(n / 2, 0, s))

test_that("Plotting method returns a ggplot object", {
  pn.res <- ebnm(x, s, prior_family = "point_normal")
  plt <- plot(pn.res)
  expect_s3_class(plt, "ggplot")
})
