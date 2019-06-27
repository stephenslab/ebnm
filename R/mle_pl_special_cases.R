# The scale of the slab is fixed, so only pi0 needs to be estimated.
#
#' @importFrom stats pnorm optimize
#'
mle_point_laplace_fixa <- function(x, s, g, control) {
  a <- g$a

  lf <- calc_lf(x, s)

  # These calculations are from nlm_fn_point_laplace.
  xleft <- x / s + s * a
  lpnormleft <- pnorm(xleft, log.p = TRUE, lower.tail = FALSE)
  lgleft <- a * x + lpnormleft
  xright <- x / s - s * a
  lpnormright <- pnorm(xright, log.p = TRUE)
  lgright <- -a * x + lpnormright
  lg <- log(a / 2) + s^2 * a^2 / 2 + logscale_add(lgleft, lgright)

  optargs <- list(f = pl_fixa_llik, lf = lf, lg = lg, interval = c(0, 1),
                  maximum = TRUE)
  optres <- do.call(optimize, c(optargs, control))

  g$pi0 <- 1 - optres$maximum
  g$val <- optres$objective

  return(g)
}

pl_fixa_llik <- function(w, lf, lg) {
  if (w == 0) {
    return(sum(lf))
  } else if (w == 1) {
    return(sum(lg))
  } else {
    return(sum(logscale_add(log(1 - w) + lf, log(w) + lg)))
  }
}
