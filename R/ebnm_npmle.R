ebnm_npmle <- function(x,
                       s = 1,
                       g_init = NULL,
                       fix_g = FALSE,
                       output = output_default(),
                       ...) {
  # Create ash grid.
  if (is.null(g_init)) {
    if (any(is.infinite(s))) {
      grid.x <- x[is.finite(s)]
      grid.s <- s[is.finite(s)]
    } else {
      grid.x <- x
      grid.s <- s
    }

    n <- length(grid.x)
    z <- qnorm(0.5 / n, lower.tail = FALSE)
    grid.min    <- min(grid.x - z * grid.s)
    grid.max    <- max(grid.x + z * grid.s)
    grid.nseg   <- max(8, ceiling((grid.max - grid.min)
                                  / (min(grid.s[grid.s > 0]) / 2)))
    grid.seglen <- (grid.max - grid.min) / grid.nseg
    g_init <- ashr::unimix(pi = rep(0, grid.nseg),
                           a = grid.min + (0:(grid.nseg - 1)) * grid.seglen,
                           b = grid.min + (1:grid.nseg) * grid.seglen)
  }

  return(ebnm_ash_workhorse(x = x,
                            s = s,
                            mode = NULL,
                            scale = NULL,
                            g_init = g_init,
                            fix_g = fix_g,
                            output = output,
                            call = match.call(),
                            ...))
}
