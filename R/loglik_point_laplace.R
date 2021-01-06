# Functions to compute the log likelihood under the point-laplace prior.
#   Currently used only for testing/reference.

loglik_point_laplace = function(x, s, w, a, mu) {
  return(sum(vloglik_point_laplace(x, s, w, a, mu)))
}

# Return log((1 - w)f + wg) as a vector (deal with cases w = 1 and w = 0
#   separately for stability).
#
#' @importFrom stats dnorm
#'
vloglik_point_laplace = function(x, s, w, a, mu) {
  if (w <= 0) {
    return(dnorm(x - mu, sd = s, log = TRUE))
  }

  lg <- logg_laplace(x - mu, s, a)
  if (w == 1) {
    return(lg)
  }

  lf <- dnorm(x - mu, sd = s, log = TRUE)
  lfac <- pmax(lg, lf)
  return(lfac + log((1 - w) * exp(lf - lfac) + w * exp(lg - lfac)))
}
