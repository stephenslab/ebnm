context("Speed Comparison")

n <- 1000
set.seed(1)
s <- rgamma(n, 1, 1)
x <- rnorm(n, 0, s + 1)

bench <- microbenchmark::microbenchmark(ebnm(x, s, "point_normal"), ashr::ash(x, s, "normal"), ebnm(x, s, "point_laplace"))

test_that("Speed of 'point_normal' relative to ash didn't change much", {
  t_mean_pn <- mean(bench$time[bench$expr == 'ebnm(x, s, "point_normal")']) # avg ebnm time
  t_mean_ash <- mean(bench$time[bench$expr == 'ashr::ash(x, s, "normal")']) # avg ash time
  
  t_ratio_pn <- t_mean_ash / t_mean_pn # calculated ~7.5 on my desktop
  
  expect_gte(t_ratio_pn, 7)
})

test_that("Speed of 'point_laplace' relative to ash didn't change much", {
  t_mean_pl <- mean(bench$time[bench$expr == 'ebnm(x, s, "point_laplace")']) # avg ebnm time
  t_mean_ash <- mean(bench$time[bench$expr == 'ashr::ash(x, s, "normal")']) # avg ash time
  
  t_ratio_pl <- t_mean_ash / t_mean_pl # calculated ~2.5 on my desktop
  
  expect_gte(t_ratio_pl, 2)
})
