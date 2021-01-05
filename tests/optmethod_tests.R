library(microbenchmark)

test_n <- function(n, seed = 666, mb_times = 100L) {
  set.seed(seed)
  # true g is 0.5 delta_0 + 0.5 N(0, 2^2)
  theta <- c(rep(0, n), rnorm(n, sd = 2))
  s <- sqrt(rexp(2 * n))
  x <- theta + rnorm(2 * n, sd = s)

  test_res <- microbenchmark(
    ebnm_point_normal(x, s, optmethod = "nlm"),
    ebnm_point_normal(x, s, optmethod = "lbfgsb"),
    ebnm_point_normal(x, s, optmethod = "trust"),
    ebnm_point_normal(x, s, optmethod = "nograd_nlm"),
    ebnm_point_normal(x, s, optmethod = "nograd_lbfgsb"),
    ebnm_point_normal(x, s, optmethod = "nohess_nlm"),
    times = mb_times)

  return(test_res)
}

test_n(1000)
test_n(10000)
test_n(100000, mb_times = 10L)
