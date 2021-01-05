init_g_for_npmle <- function(x, scale) {
  if (identical(scale, "estimate")) {
    ncomp <- min(max(10, ceiling(sqrt(length(x)))), 300)
  } else {
    ncomp <- scale
  }

  grid <- quantile(x, probs = c(0, 1:ncomp / ncomp))

  g_init <- ashr::unimix(pi = rep(1 / ncomp, ncomp),
                         a = grid[1:ncomp],
                         b = grid[2:(ncomp + 1)])

  return(g_init)
}
