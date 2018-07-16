# Functions to compute the log likelihood under the normal prior.


# The log density of normal(0, 1/a) convolved with normal(0, s).
#
#' @importFrom stats dnorm
#'
logg_normal = function(x, s, a) {
  dnorm(x, 0, sqrt(s^2 + 1/a), log=TRUE)
}


# Return log((1-w)f + wg) as a vector.
# Deal with cases w = 1 and w = 0 separately for stability.
#
#' @importFrom stats dnorm
#'
vloglik_normal = function(x, s, w, a) {
  if (w <= 0) {
    return(dnorm(x, 0, s, log = TRUE))
  }
  lg = logg_normal(x, s, a)
  if (w >= 1) {
    return(lg)
  }

  lf = dnorm(x, 0, s, log = TRUE)
  lfac = pmax(lg, lf)
  result = lfac + log((1 - w) * exp(lf - lfac) + w * exp(lg - lfac))
  # Handle the case where s = 0, so that the density is infinite:
  result[is.infinite(lf)] = log(1 - w)

  return(result)
}


loglik_normal = function(x, s, w, a) {
  sum(vloglik_normal(x, s, w, a))
}
