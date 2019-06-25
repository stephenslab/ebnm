ebnm_npmle <- function(x,
                       s = 1,
                       g_init = NULL,
                       fix_g = FALSE,
                       output = output_default(),
                       ...) {
  max_ncomp <- 20

  # Create ash grid.
  if (is.null(g_init)) {
    pts <- sort(x)
    if (length(x) > max_ncomp) {
      idx <- c(1, ceiling(length(x) / (max_ncomp - 1) * 1:(max_ncomp - 1)))
      idx[max_ncomp] <- length(x)
      pts <- pts[idx]
    }
    ncomp <- length(pts)
    mids <- (pts[1:(ncomp - 1)] + pts[2:ncomp]) / 2
    mids <- c(pts[1] - (mids[1] - pts[1]),
              mids,
              pts[ncomp] + (pts[ncomp] - mids[ncomp - 1]))
    g_init <- ashr::unimix(pi = rep(0, ncomp),
                           a = mids[1:ncomp],
                           b = mids[2:(ncomp + 1)])
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
