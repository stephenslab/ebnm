init_g_for_npmle <- function(x, scale) {
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)

  if (identical(scale, "estimate")) {
    ncomp <- max(10, ceiling(sqrt(length(x))))
  } else {
    ncomp <- max(1, ceiling((x_max - x_min) / scale))
  }

  comp_len <- (x_max - x_min) / ncomp

  g_init <- ashr::unimix(pi = rep(1 / ncomp, ncomp),
                         a = x_min + (0:(ncomp - 1)) * comp_len,
                         b = x_min + (1:ncomp) * comp_len)

  return(g_init)
}
