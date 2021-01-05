init_g_for_npmle <- function(x, s, scale) {
  min_K <- 10
  max_K <- 300
  KLdiv_target <- 0.5

  xrange <- diff(range(x))

  # See ebnm paper for rationale.
  if (xrange / max_K < 2.3 * min(s)) {
    # Use point-mass mixture.
    d <- min(s) * (64 * KLdiv_target)^0.25
    ncomp <- min(max(min_K, ceiling(xrange / d)), max_K) + 1
    grid <- seq(min(x), max(x), length.out = ncomp)

    g_init <- ashr::unimix(pi = rep(1 / ncomp, ncomp),
                           a = grid,
                           b = grid)
  } else {
    # Use Gaussian mixture.
    d <- 2 * min(s) * sqrt(exp(2 * KLdiv_target) - 1)
    ncomp <- min(max(min_K, ceiling(xrange / d)), max_K) + 1
    grid <- seq(min(x), max(x), length.out = ncomp)

    g_init <- ashr::normalmix(pi = rep(1 / ncomp, ncomp),
                              mean = grid,
                              sd = d / 2)
  }

  return(g_init)
}
